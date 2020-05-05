from info_file import InfoFile
from gui.base import gdGui
import os
import numpy as np
#import cupy as cp
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
        #self.startingValues()
        #temp=gdGui()
        #self.vals=temp.sendValues()
        self.vals={'input_volume':'test.mrc', 'mask':'test_mask.mrc', 'pdb': 'fit.pdb'}
        #parser.add_argument('-w','--patch_width',default=3,type=int)
        #parser.add_argument('-l','--patch_length',default=3,type=int)
        #parser.add_argument('-s','--shift',action='store_true')
        #parser.add_argument('-f','--frealign_input_file',default='info.txt')
        
    def makeMasks(self):
        self.readFiles()
        
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
        x, y=np.meshgrid(np.arange(-size,size),np.arange(-size,size))
        #temp=np.arctan2(y, x)
        #temp=np.reshape(temp,(size*2,size*2,1))
        #self.radmatrix=np.repeat(temp, size*2, axis=2)
        self.radmatrix=np.arctan2(y, x)
        self.zline=np.arange(-size,size)
        
    def something(self):
        theta0=np.arctan2(self.com[1], self.com[0])    
        z0=self.com[2]/self.pixel_size
        theta=self.rise_per_subunit/self.twist_per_subunit*(self.zline-z0)-theta0
        wedge=np.zeros(self.vol_dim)
        
        fudge=np.deg2rad(20)
        for i in range(len(theta)):
            above=theta[i]-fudge
            below=theta[i]+fudge
            inds=np.logical_and(self.radmatrix>above,self.radmatrix<below)
            wedge[:,:,i][inds]=1
            
        test=EMImage(wedge)
        test.write_mrc('puhleeze.mrc')
        
    def createWedge(self, vol_dim, num_subunits, width, ref_pdb, pf_num=6, edge_resolution=20):
        ###Create a wedge mask centered at pf6 to be applied to reference volume
        wedge=np.zeros((vol_dim[0],vol_dim[1],vol_dim[2]))
        subunit_num=self.subunits_to_use[pf_num,0]*self.num_pfs+pf_num
        mask=wedge
        
        edge_width = self.pixel_size * np.ceil(edge_resolution/(2*self.pixel_size))
        cosmask_filter = np.fft.fftshift(spherical_cosmask(vol_dim, 0, edge_width / self.pixel_size))
        cosmask_filter = np.fft.fftn(cosmask_filter) / np.sum(cosmask_filter)
        soft_m = np.real(np.fft.ifftn( cosmask_filter * np.fft.fftn(mask)))
        
        soft_m[soft_m<0]=0
        ##Need to transpose to be on the same axes as before
        self.soft_wedge=EMImage(soft_m.transpose())
        
    def run(self):
        self.getValues()
        self.makeMasks()