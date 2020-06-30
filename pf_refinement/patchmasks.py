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


class WedgeMasks(InfoFile):
    def __init__(self, info='info.txt', shift=True):
        super(WedgeMasks, self).__init__(info)
        #self.subunits_to_use=subunitsToUse(self.rise_per_subunit, self.num_pfs, \
        #                       self.num_starts, shift)
        
    def getValues(self):
        vals=self.startingValues()
        temp=gdGui('Generate Wedge Masks', **vals)
        self.vals=temp.sendValues()
        
    def startingValues(self):
        vals={
            'microtubule_volume':'XXXMRCFILEXXX', 
            'microtubule_mask':'XXXMRCFILEXXX', 
            'fit_tubulin_pdb': 'XXXFILEXXX'
            }
        return vals
    
    def readFiles(self):
        self.pdb=pdb_util.PDB()
        self.pdb.read_file(self.vals['fit_tubulin_pdb'])
        self.com=self.pdb.calc_center()
        self.vol=EMImage(self.vals['microtubule_volume'])
        self.vol_dim=np.asarray(self.vol.data.shape)
        mask=EMImage(self.vals['microtubule_mask'])
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
        theta0=np.arctan2(self.com[0], self.com[1])+\
        np.deg2rad(subunit_num*self.twist_per_subunit)
            
        z0=(self.com[2]+self.rise_per_subunit*subunit_num)/self.pixel_size
        
        ###Define the length along the protofilament in terms of subunits 
        zsubunits=(self.zline.copy()-z0)*self.pixel_size/self.dimer_repeat_dist
        
        ###Define the angle of the center of the protofilament along the length of the segment
        theta=np.deg2rad((-self.helical_twist)*zsubunits)+theta0
        
        ###Initialize the wedge mask
        wedge=np.zeros(self.vol_dim.tolist())
        
        ###Define the size of the wedgemask
        fudge=np.deg2rad(360.0/(self.num_pfs*2))
        
        ###Generate the wedge mask
        for i in range(len(theta)):
            ### Find the range of angles for the  wedge mask to fall between
            temp1=np.remainder(theta[i]-fudge+2*np.pi,2*np.pi)-2*np.pi
            temp2=np.remainder(theta[i]+fudge+2*np.pi,2*np.pi)-2*np.pi
            angles=[temp1, temp2]
            ### Account for the 360 degree wrapping
            if max(angles)-min(angles)>2*fudge+.2:
                above=max(angles)
                below=min(angles)
                inds=np.logical_or(self.radmatrix>above,self.radmatrix<below)
            else:
                above=min(angles)
                below=max(angles)
                inds=np.logical_and(self.radmatrix>above,self.radmatrix<below)
                
            wedge[i,:,:][inds]=1
            
        return wedge
    
    def cosmask_filter(self):
        ### Add soft edge to the mask
        edge_resolution=20
        edge_width = self.pixel_size * np.ceil(edge_resolution/(2*self.pixel_size))
        cosmask_filter = np.fft.fftshift(spherical_cosmask(self.vol_dim, 0, edge_width / self.pixel_size))
        self.cosmask_filter = np.fft.fftn(cosmask_filter) / np.sum(cosmask_filter)
    
    
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