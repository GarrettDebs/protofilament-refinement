try:
    nope
    import cupy as np
except:
    import numpy as np

from pf_refinement import InfoFile
from gui.base import gdGui
import os
import pdb_util
from EMImage import EMImage
from util3d import spherical_cosmask
from tqdm import tqdm
from scipy.ndimage.interpolation import affine_transform


class WedgeMasks(InfoFile):
    def __init__(self, info='info.txt', shift=True):
        super(WedgeMasks, self).__init__(info)
        #self.subunits_to_use=subunitsToUse(self.rise_per_subunit, self.num_pfs, \
        #                       self.num_starts, shift)
        
    def getValues(self):
        vals=self.startingValues()
        temp=gdGui('Generate Wedge Masks', **vals)
        self.vals=temp.sendValues()
        #self.vals={'input_volume':'test.mrc', 'mask':'test_mask.mrc', 'pdb': 'fit.pdb'}
        #parser.add_argument('-w','--patch_width',default=3,type=int)
        #parser.add_argument('-l','--patch_length',default=3,type=int)
        #parser.add_argument('-s','--shift',action='store_true')
        #parser.add_argument('-f','--frealign_input_file',default='info.txt')
        
    def startingValues(self):
        vals={'input_volume':'XXXMRCFILEXXX', 'mask':'XXXMRCFILEXXX', \
              'pdb': 'XXXFILEXXX'}#, 'output_volume':'XXXMRCFILEXXX'}
        return vals
    
    def readFiles(self):
        self.pdb=pdb_util.PDB()
        self.pdb.read_file(self.vals['pdb'])
        self.com=self.pdb.calc_center()
        self.vol=EMImage(self.vals['input_volume'])
        self.vol_dim=np.asarray(self.vol.data.shape)
        mask=EMImage(self.vals['mask'])
        self.vol.mult(mask)
        
    def initCoords(self):
        size=self.vol_dim[0]/2
        x, y=np.meshgrid(np.arange(-size+1,size+1),np.arange(-size+1,size+1))
        #temp=np.arctan2(y, x)
        #temp=np.reshape(temp,(size*2,size*2,1))
        #self.radmatrix=np.repeat(temp, size*2, axis=2)
        self.radmatrix=np.remainder(np.arctan2(x, y)+2*np.pi,2*np.pi)-2*np.pi
        self.zline=np.arange(-size+1,size+1)
        
    def createWedge(self, subunit_num=0):
        ###Initialize the location of the protofilament to remove
        theta0=np.arctan2(self.com[0], self.com[1])-\
        np.deg2rad(subunit_num*self.twist_per_subunit)
            
        z0=(self.com[2]+self.rise_per_subunit*subunit_num)/self.pixel_size
        
        ###Define the length along the protofilament in terms of subunits 
        zsubunits=(self.zline-z0)/(self.rise_per_subunit*self.num_pfs/\
                                   (self.num_starts*self.pixel_size))
        
        ###Define the angle of the center of the protofilament along the length of the segment
        theta=np.deg2rad((self.twist_per_subunit*self.num_pfs-360)*zsubunits)-theta0
        
        ###Initialize the wedge mask
        wedge=np.zeros(self.vol_dim.tolist())
        
        ###Define the size of the wedgemask
        #fudge=np.deg2rad(360.0/(self.num_pfs*2))
        fudge=np.deg2rad(16)
        
        ###Generate the wedge mask
        for i in range(len(theta)):
            temp1=np.remainder(theta[i]-fudge+2*np.pi,2*np.pi)-2*np.pi
            temp2=np.remainder(theta[i]+fudge+2*np.pi,2*np.pi)-2*np.pi
            angles=[temp1, temp2]
            if max(angles)-min(angles)>2*fudge+.2:
                above=max(angles)
                below=min(angles)
                inds=np.logical_or(self.radmatrix>above,self.radmatrix<below)
            else:
                above=min(angles)
                below=max(angles)
                inds=np.logical_and(self.radmatrix>above,self.radmatrix<below)
                
            wedge[i,:,:][inds]=1
            
        print above, below
        print np.rad2deg(above), np.rad2deg(below)
        return wedge
    
    def cosmask_filter(self):
        edge_resolution=20
        edge_width = self.pixel_size * np.ceil(edge_resolution/(2*self.pixel_size))
        cosmask_filter = np.fft.fftshift(spherical_cosmask(self.vol_dim, 0, edge_width / self.pixel_size))
        self.cosmask_filter = np.fft.fftn(cosmask_filter) / np.sum(cosmask_filter)
    
    def rotshift3D_spline(self, v, eulers, shifts=np.array([0,0,0]), mode='wrap'):
    # With a nod to:
    #  http://stackoverflow.com/questions/20161175/how-can-i-use-scipy-ndimage-interpolation-affine-transform-to-rotate-an-image-ab
        print shifts
        print eulers
        rot_origin = 0.5*np.array(v.shape)
        rot_rad = -np.deg2rad(eulers)
        phi = np.array([[np.cos(rot_rad[0]), np.sin(rot_rad[0]), 0],
                               [-np.sin(rot_rad[0]),np.cos(rot_rad[0]), 0],
                               [0                , 0              , 1]])
        theta=np.array([[np.cos(rot_rad[1]), 0,-np.sin(rot_rad[1])],[0,1,0],\
                            [np.sin(rot_rad[1]), 0, np.cos(rot_rad[1])]])
        psi=np.array([[np.cos(rot_rad[2]), np.sin(rot_rad[2]), 0],\
                          [-np.sin(rot_rad[2]),np.cos(rot_rad[2]), 0],\
                          [0,0,1]])
        rot_matrix=np.dot(np.dot(phi,theta),psi)
        offset = -(rot_origin-rot_origin.dot(rot_matrix)).dot(np.linalg.inv(rot_matrix))
        offset = offset - shifts
        transformed_v = affine_transform(v,rot_matrix,offset=offset,mode=mode)
        return transformed_v
    
    def makeMasks(self):
        self.cosmask_filter()
        #hard_wedge=self.createWedge()
        rot=self.twist_per_subunit*np.arange(self.num_pfs)
        
        for pf in tqdm(range(self.num_pfs)):
            if not os.path.isdir('pf%g'%pf):
                os.mkdir('pf%g'%pf)

            rot_wedge=self.createWedge(pf)
            soft_m = np.real(np.fft.ifftn(self.cosmask_filter * np.fft.fftn(rot_wedge)))
            soft_m[soft_m<0]=0
            wedge=EMImage(soft_m)            
            wedge.write_mrc('pf%g/pf_wedge.mrc'%pf)
            wedge.add(-1)
            wedge.mult(-1)
            wedge.mult(self.vol)
            wedge.write_mrc('pf%g/pf_masked.mrc'%pf)
        
    def run(self):
        self.getValues()
        self.readFiles()
        self.initCoords()
        self.makeMasks()