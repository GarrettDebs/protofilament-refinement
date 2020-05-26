import numpy as np

class Basics():
    def __init__(self,info_file):
        self.readInfo(info_file)
    
    def readInfo(self,info_file='info.txt'):
        f=open(info_file)
        info=f.readlines()
        f.close()
        for i in range(len(info)):
            temp=info[i].split()
            if temp[0].startswith('voxel'):
                self.pixel_size=float(temp[1])
                continue
            elif temp[0].startswith('bin'):
                self.bin_factor=float(temp[1])
                continue
            elif temp[0].startswith('num_pfs'):
                self.num_pfs=int(temp[1])
                continue
            elif temp[0].startswith('helical_rise'):
                self.dimer_repeat_dist=float(temp[1])
                continue
            elif temp[0].startswith('num_starts'):
                self.num_starts=float(temp[1])/2
                continue
            elif temp[0].startswith('helical_twist'):
                self.helical_twist=float(temp[1])
            elif temp[0].startswith('micrograph_pixels_per'):
                self.dimer_repeat_dist=float(temp[1])*self.pixel_size/self.bin_factor
        
        self.rise_per_subunit=(self.dimer_repeat_dist*self.num_starts)/self.num_pfs
        self.twist_per_subunit=(360-self.helical_twist/self.num_starts)/self.num_pfs
        print 'voxel: %g bin: %d num_pfs: %d dimer_repeat: %g twist: %g'\
        %(self.pixel_size, self.bin_factor, self.num_pfs, self.dimer_repeat_dist, self.helical_twist)
    
    def subunitsToUse(self, shift=True):
        index=0
        out=np.zeros([len(range(int(-2*self.num_starts),int(2*self.num_starts+1)))*self.num_pfs,3])
        
        for i in range(int(-2*self.num_starts),int(2*self.num_starts+1)):
            for j in range(0,self.num_pfs):
                z = i*self.num_pfs*self.rise_per_subunit/self.num_starts + (j)*self.rise_per_subunit;
                out[index,0]=i
                out[index,1]=j
                out[index,2]=z
                index += 1
                
        out=out[out[:,2].argsort()]
        if shift:
            lowerb=np.nonzero(out[:,2]>=-(self.rise_per_subunit*self.num_pfs/(self.num_starts*2)))
            upperb=np.nonzero(out[:,2]<(self.rise_per_subunit*self.num_pfs/(self.num_starts*2)))
        else:
            lowerb=np.nonzero(out[:,2]>=0)
            upperb=np.nonzero(out[:,2]<(self.rise_per_subunit*self.num_pfs/(self.num_starts)))
                              
        inter=np.intersect1d(lowerb,upperb)
        use=out[inter]
        use=use[use[:,1].argsort()]
        self.subunits_to_use=use
        
    def rotateMT(self,mt,pfshift,halfshift):
        phiInd=self.key.index('AngleRot')
        psiInd=self.key.index('AnglePsi')
        thetaInd=self.key.index('AngleTilt')
        xInd=self.key.index('OriginX')
        yInd=self.key.index('OriginY')
        
        if not hasattr(self,'subunits_to_use'):
            self.subunitsToUse(True)
            self.halfdist=self.dimer_repeat_dist*.5/self.pixel_size
                        
        use=self.subunits_to_use[:,2]/self.pixel_size
                
        data=np.array(self.mt[mt])
        phi=data[:,phiInd].astype(float)
        psi=np.deg2rad(data[:,psiInd].astype(float))
        theta=np.deg2rad(data[:,thetaInd].astype(float))
        x=data[:,xInd].astype(float)
        y=data[:,yInd].astype(float)
        
        x+=use[pfshift]*np.cos(psi)*np.sin(theta)-\
        self.halfdist*halfshift*np.cos(psi)*np.sin(theta)
        #print self.halfdist*halfshift*np.cos(psi)*np.sin(theta), self.halfdist*halfshift*np.sin(psi)*np.sin(theta)  
        
        y-=use[pfshift]*np.sin(psi)*np.sin(theta)+\
        self.halfdist*halfshift*np.sin(psi)*np.sin(theta)
        
        phi+=self.twist_per_subunit*pfshift
        
        data[:,phiInd]=phi.astype(str)
        data[:,xInd]=x.astype(str)
        data[:,yInd]=y.astype(str)
        
        self.mt[mt]=data.tolist()
        