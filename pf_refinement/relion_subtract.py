from pf_refinement import InfoFile, StarFile
from gui import gdGui
import numpy as np
import os
from distutils.log import info

class SubStar(InfoFile):
    def __init__(self, info='info.txt'):
        super(SubStar, self).__init__(info)
        
    def getVals(self, protosub=False):
        if protosub:
            vals={
                'input_star': 'XXXSTARFILEXXX',
                'input_volume': 'XXXMRCFILEXXX',
                'protofilament_mask': 'XXXMRCFILE',
                'cores': 1
                }
        else:
            vals={
                'input_star': 'XXXSTARFILEXXX',
                'cores': 1
                }
        temp=gdGui('Generate Protofilament Particles', **vals)
        self.vals=temp.sendValues()
        
    def __call__(self, protosub=False):
        self.getVals(protosub)
        self.star_file=self.vals['input_star']
        self.cores=int(self.vals['cores'])
        if protosub:
            self.makeProto()
            self.volume=self.vals['input_volume']
            self.mask=self.vals['protofilament_mask']
            self.splitPfs()
            
    def makeProto(self):
        from pf_refinement import EMImage
        
        proto=EMImage(self.vals['input_volume'])
        mask=EMImage(self.vals['protofilament_mask'])
        
        proto.mult(mask)
        proto.write_mrc('protofilament_for_sub.mrc')
    
    def splitPfs(self):        
        star=StarFile(self.star_file)
        star.df.ImageName=star.df.ImageOriginalName
        for i in range(self.num_pfs):
            data=star.df.iloc[i::self.num_pfs]
            star.writeStar('pf%g_for_sub.star'%i, data)
    
class Sub(InfoFile):
    def __init__(self, star_file, info='info.txt'):
        super(Sub, self).__init__(info)
        self.star_file = star_file
        
    def __call__(self, pfnum, rank=0):
        command='relion_project --i pf%g/pf_masked.mrc --o '\
        'pf%g/proto_particles --ctf --angpix %g --ang %s --subtract_exp'%\
        (pfnum, pfnum, self.pixel_size, self.star_file)
        code=os.system(command)
        if code!=0:
            raise RuntimeError('RELION did not run properly. Try running the '
                            'following command to troubleshoot \n\n'
                            '%s \n'%command)
        
    def parallel(self):
        from mpi4py import MPI
        
        pfs=np.arange(self.num_pfs)
        comm=MPI.COMM_WORLD
        rank=comm.Get_rank()
        size=comm.Get_size()
        nums=pfs[pfs%size==rank]
        
        for i in nums:
            self(i, rank=rank)
            
class ProtoSub(InfoFile):
    def __init__(self, info='info.txt'):
        super(ProtoSub, self).__init__(info)
        
        from mpi4py import MPI

        self.comm=MPI.COMM_WORLD
        self.rank=self.comm.Get_rank()
        self.size=self.comm.Get_size()
        
    def __call__(self):
        pfs=np.arange(self.num_pfs)
        nums=pfs[pfs%self.size==self.rank]
        
        for start_pf in nums:
            self.neighbors(start_pf)
            self.subProtoPfs(start_pf)
            self.finalSub(start_pf,True)
            
        if rank==0:
            data=1
            for i in range(1, self.size):
                req=self.comm.Isend(data, dest=i)
                req.wait()
            
            self.writeProtoStar()
                
        else:
            data=0
            req=self.comm.Irecv(data,source=0)
            req.wait()


    def neighbors(self,pf_num):
        self.n1=pf_num-1
        self.n2=pf_num+1
        if pf_num==self.num_pfs-1:
            self.n2=0
        elif pf_num==0:
            self.n1=self.num_pfs-1
        
    def subProtoPfs(self, start_pf):
        dirname='pf%g_protosubbed'%start_pf
        if not os.path.isdir(dirname):
            os.mkdir(dirname)
        
        start=start_pf+1
        if start==self.num_pfs:
            start=0
        
        command='relion_project --i protofilament_for_sub.mrc '\
        '--o %s/pf%g_subbed --ctf ' \
        '--angpix %g --ang pf%g_for_sub.star --subtract_exp'%\
        (dirname, start, self.pixel_size, start)
        #os.system(command)
        for i in range(start+1, self.num_pfs):
            if i==self.num_pfs-1 and start_pf+1==self.num_pfs:
                break
            self.subPf(i, dirname)
            if i>start+1 and i!=self.n1:
                #os.system('rm %s/pf%g_subbed.mrcs'%(dirname, i-1))
                
        if start_pf!=self.num_pfs-1:
            for i in range(start_pf):
                self.subPf(i, dirname)
                if i>2 and i!=self.n1:
                    #os.system('rm %s/pf%g_subbed.mrcs'%(dirname, i-1))

            
    def subPf(self, pf_num, dirname, vol='protofilament_for_sub.mrc'):
        if pf_num==0:
            sub=self.num_pfs-1
        else:
            sub=pf_num-1

        def changename(x):
            return '%g@%s/pf%g_subbed.mrcs'%(x.name+1, dirname, sub)
        
        star=StarFile('pf%g_for_sub.star'%pf_num) 
          
        star.df.ImageName=star.df.apply(changename, axis=1)
        
        star.writeStar('%s/pf%g_for_sub.star'%(dirname,pf_num))
        
        command='relion_project --i %s --o %s/pf%g_subbed --ctf ' \
        '--angpix %g --ang %s/pf%g_for_sub.star --subtract_exp'\
        %(vol, dirname, pf_num, self.pixel_size, dirname, pf_num)
        #os.system(command)
    
    def finalSub(self,start_pf,clean=False):
        dirname='pf%g_protosubbed'%start_pf
        
        def changename(x):
            return '%g@%s/pf%g_subbed.mrcs'%(x.name+1, dirname, self.n1)
        
        star=StarFile('pf%g_for_sub.star'%start_pf)  
        star.df.ImageName=star.df.apply(changename, axis=1)
        
        star.writeStar('%s/pf%g_finalsub.star'%(dirname,start_pf))
        if clean:
            os.system('rm %s/*subbed.star'%(dirname))
            os.system('rm %s/*for_sub.star'%(dirname))
            os.system('mv %s/pf%g_subbed.mrcs %s/safe'%\
                      (dirname, self.n1, dirname))
            os.system('rm %s/*subbed.mrcs'%(dirname))
            os.system('mv %s/safe %s/pf%g_subbed.mrcs'%\
                      (dirname, dirname, self.n1))
            
    def writeProtoStar(self):
        dfs=[]
        star=StarFile('pf0_protosubbed/pf0_finalsub.star')
        star.readInfo('info.txt')
        
        num_pfs=self.num_pfs
        for i in range(num_pfs):
            star=StarFile('pf%g_protosubbed/pf%g_finalsub.star'%(i,i))
            index=pd.Index(range(i, len(star.df)*num_pfs, num_pfs))
            dfs.append(star.df.set_index(index))
        
        
        df=pd.concat(dfs).sort_index()
        star.writeStar('protosubbed.star',df)
    
        