from starfile import StarFile
from helix import subunitsToUse
import numpy as np
import pandas as pd

class StarOp(StarFile):
    def __init__(self, file, info='info.txt'):
        super(StarOp, self).__init__(file)
        self.readInfo(info)
        
    def mtGroup(self):
        ###Create a means of easily indexing each microtubule
        self.mt=self.df.groupby(['MicrographName', 'HelicalTubeID'])
        
    def delPhiXs(self):
        ###Make sure that the microtubule groups have been made
        if not hasattr(self, 'mt'):
            self.mtGroup()
            
        self.dphixs={}
        for key in self.mt.groups.keys():   
            ##Calculate delta phi for for all protofilaments in the microtubule
            phis=self.mt.get_group(key).loc[:,'AngleRot'].to_numpy(dtype=float)
            dphi=phis[1:]-phis[:-1]
            dphi=np.append(dphi,0)
            dphi[self.num_pfs-1::self.num_pfs]=phis[::self.num_pfs]-\
            phis[self.num_pfs-1::self.num_pfs]
            dphi[dphi>100]-=360
            dphi[dphi<-100]+=360
            temp=np.zeros((len(dphi)/self.num_pfs,self.num_pfs))
            for i in range(self.num_pfs):
                temp[:,i]=dphi[i::self.num_pfs]
                
            self.dphixs[key]=temp
        
    def symStar(self):
        try:
            self.mt.keys[0]
        except Exception:
            self.mtGroup()
            
        use=subunitsToUse(self.rise_per_subunit, num_pfs=self.num_pfs, 
                          num_starts=self.num_starts)
        
        num_particles=len(self.df)
        
        rise=np.tile(use[:,2]/self.pixel_size,num_particles)
        pfnums=np.tile(np.arange(self.num_pfs),num_particles)

        temp=self.df.iloc[np.repeat(np.arange(num_particles),self.num_pfs)].copy()
        temp.reset_index(inplace=True)
        temp.drop(columns='index',inplace=True)
        
        phi=np.deg2rad(temp.AngleRot.values.astype(float))
        theta=np.deg2rad(temp.AngleTilt.values.astype(float))
        psi=np.deg2rad(temp.AnglePsi.values.astype(float))

        temp.OriginX=(temp.OriginX.values.astype(float) + 
                                  rise*np.sin(theta)*np.cos(psi)).astype(str)           
        temp.OriginY=(temp.OriginY.values.astype(float) - 
                                  rise*np.sin(theta)*np.sin(psi)).astype(str)
        temp.AngleRot=(temp.AngleRot.values.astype(float) +
                                   pfnums*self.twist_per_subunit).astype(str)
                                   
        self.df=temp
        self.writeStar(self.name.split('.')[0]+'_sub.star')
        
        self.df.rename(columns={'ImageName': 'ImageOriginalName'}, inplace=True)
        
        newnames=np.empty(len(self.df), dtype='S100')
        for i in range(1, num_particles+1):
            for j in range(self.num_pfs):
                name='%g@pf%g/proto_particles.mrcs'%(i, j)
                newnames[(i-1)*self.num_pfs+j]=name
        
        self.df=self.df.assign(ImageName=newnames)
        self.writeStar('proto_particle_stack.star')
            