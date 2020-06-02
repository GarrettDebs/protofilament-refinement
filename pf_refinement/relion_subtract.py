from pf_refinement import InfoFile
from gui import gdGui
import numpy as np
import os

class SubStar(InfoFile):
    def __init__(self, info='info.txt'):
        super(SubStar, self).__init__(info)
        
    def getVals(self):
        vals={'input_star': 'XXXSTARFILEXXX',
              'cores': 1}
        temp=gdGui('Generate Protofilament Particles', **vals)
        self.vals=temp.sendValues()
        
    def __call__(self):
        self.getVals()
        self.star_file=self.vals['input_star']
        self.cores=int(self.vals['cores'])
    
class Sub(InfoFile):
    def __init__(self, star_file, cores=1, info='info.txt'):
        super(Sub, self).__init__(info)
        self.star_file = star_file
        self.cores=cores
        
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
        
        b=Sub(self.star_file)
        for i in nums:
            b(i, rank=rank)    