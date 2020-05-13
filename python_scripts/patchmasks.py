try:
    nope
    import cupy as np
except:
    import numpy as np

from info_file import InfoFile
from gui.base import gdGui
import os
import pdb_util
from helix.basics import helical_wedge_mask, subunitsToUse
from EMImage import EMImage
from util3d import spherical_cosmask
from cupyx.scipy.ndimage.interpolation import affine_transform

class WedgeMasks(InfoFile):
    def __init__(self, info='info.txt', shift=True):
        super(WedgeMasks, self).__init__(info)
        self.subunits_to_use=subunitsToUse(self.rise_per_subunit, self.num_pfs, \
                               self.num_starts, shift)
        
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
        self.radmatrix=np.arctan2(x, y)
        self.zline=np.arange(-size+1,size+1)
        
    def createWedge(self, subunit_num=0):
        ###Initialize the location of the protofilament to remove
        theta0=np.arctan2(self.com[0], self.com[1])+\
        np.deg2rad(subunit_num*self.twist_per_subunit)
            
        z0=(self.com[2]+self.rise_per_subunit*subunit_num)/self.pixel_size
        
        ###Define the length along the protofilament in terms of subunits 
        zsubunits=(self.zline-z0)/(self.rise_per_subunit*self.num_pfs/\
                                   (self.num_starts*self.pixel_size))
        
        ###Define the angle of the center of the protofilament along the length of the segment
        theta=np.deg2rad((360-self.twist_per_subunit*self.num_pfs)*zsubunits)-theta0
        
        ###Initialize the wedge mask
        wedge=np.zeros(self.vol_dim.tolist())
        
        ###Define the size of the wedgemask
        fudge=np.deg2rad(360.0/(self.num_pfs*2))
        
        ###Generate the wedge mask
        for i in range(len(theta)):
            above=theta[i]-fudge
            below=theta[i]+fudge
            inds=np.logical_and(self.radmatrix>above,self.radmatrix<below)
            wedge[i,:,:][inds]=1
            
        ###Make a soft mask
        edge_resolution=20
        edge_width = self.pixel_size * np.ceil(edge_resolution/(2*self.pixel_size))
        cosmask_filter = np.fft.fftshift(spherical_cosmask(self.vol_dim, 0, edge_width / self.pixel_size))
        cosmask_filter = np.fft.fftn(cosmask_filter) / np.sum(cosmask_filter)
        soft_m = np.real(np.fft.ifftn( cosmask_filter * np.fft.fftn(wedge)))
        soft_m[soft_m<0]=0
        
        
        return EMImage(soft_m)
    
    def makeMasks(self):
        for pf in range(self.num_pfs)[:1]:
            if not os.path.isdir('pf%g'%pf):
                os.mkdir('pf%g'%pf)
                
            wedge=self.createWedge(pf)
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