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
            # Calculate delta phi for for all protofilaments in the microtubule
            
            ### Get the phi values
            phis=self.mt.get_group(key).loc[:,'AngleRot'].to_numpy(dtype=float)
            
            ### Caclulate dphi for all adjacent phi values
            dphi=phis[1:]-phis[:-1]
            
            ### Because the seam dphi is not correctly calculated, we add
            ### an extra column, and specifically calculate the seam dphis
            dphi=np.append(dphi,0)
            dphi[self.num_pfs-1::self.num_pfs]=phis[::self.num_pfs]-\
            phis[self.num_pfs-1::self.num_pfs]
            
            ### Check for any 'wrapping errors' (i.e. 359-1 should be 2 degrees
            ### not 358 degrees. 
            dphi[dphi>100]-=360
            dphi[dphi<-100]+=360
            
            ### Next reshape it so that we can easily access all dphis for a
            ### given particle number and/or protofilament number
            temp=np.zeros((len(dphi)/self.num_pfs,self.num_pfs))
            for i in range(self.num_pfs):
                temp[:,i]=dphi[i::self.num_pfs]
                
            ### Save the numpy array in a dictionary.
            self.dphixs[key]=temp
            
        ### Flatten all data for all microtubules into a single array.
        ### All dphi values associated with protofilament i are accessible
        ### with by indexing as [i::num_pfs].
        temp=[]
        blah=[temp.extend(self.dphixs[key].flatten().tolist()) 
              for key in self.dphixs.keys()]

        self.dphi_flat=np.array(temp)
        
    def symStar(self):
        ### Defines the shift for each protofilament to ensure centering
        ### (i.e. if the shift is greater than the helical rise, shift down)
        use=subunitsToUse(self.rise_per_subunit, num_pfs=self.num_pfs, 
                          num_starts=self.num_starts)
        
        num_particles=len(self.df)
        
        ### Calculate the rise and rotation for each particle
        rise=np.tile(use[:,2]/self.pixel_size,num_particles)
        pfnums=np.tile(np.arange(self.num_pfs),num_particles)
        shifted_rot=np.tile(use[:,0]*-self.helical_twist, num_particles)

        ### Duplicate the DataFrame at the correct size
        temp=self.df.iloc[np.repeat(np.arange(num_particles),self.num_pfs)].copy()
        temp.reset_index(inplace=True)
        temp.drop(columns='index',inplace=True)
        
        ### Get the Euler Angles
        phi=np.deg2rad(temp.AngleRot.values.astype(float))
        theta=np.deg2rad(temp.AngleTilt.values.astype(float))
        psi=np.deg2rad(temp.AnglePsi.values.astype(float))

        ### Add what we need to
        temp.OriginX=(temp.OriginX.values.astype(float) + 
                                  rise*np.sin(theta)*np.cos(psi)).astype(str)           
        temp.OriginY=(temp.OriginY.values.astype(float) - 
                                  rise*np.sin(theta)*np.sin(psi)).astype(str)
        temp.AngleRot=(temp.AngleRot.values.astype(float) +
                                   pfnums*self.twist_per_subunit +
                                   shifted_rot).astype(str)
                                   
        ### Save file
        self.df=temp
        self.writeStar('sym_coords.star')
        
        ### Now need to save with the new name accounting for the subtraction
        self.df.rename(columns={'ImageName': 'ImageOriginalName'}, inplace=True)
        
        newnames=np.empty(len(self.df), dtype='S100')
        for i in range(1, num_particles+1):
            for j in range(self.num_pfs):
                name='%g@pf%g/proto_particles.mrcs'%(i, j)
                newnames[(i-1)*self.num_pfs+j]=name
        
        self.df=self.df.assign(ImageName=newnames)
        self.writeStar('proto_particle_stack.star')
        
    def truncStar(self, size, output):
        try:
            self.mt.keys[0]
        except Exception:
            self.mtGroup()
        
        i=1
        keys=self.mt.groups.keys()
        
        newdf=self.mt.get_group(keys[0])
        
        while len(newdf)<size and i<len(keys):
            newdf=pd.concat([newdf, self.mt.get_group(keys[i])])
            i+=1
        
        newdf=newdf.sort_index()
        self.writeStar(output, newdf)
            