from starfile import StarFile
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
        