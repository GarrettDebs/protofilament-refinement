import gui.base as bs
import os


class InfoFile(object):
    def __init__(self, file='info.txt'):
        #if os.path.isfile('info.txt'):
        #    print 'what'
        try:
            self.readInfo(file)
        except:
            raise NameError('No info file, please initialize the project.')
            
    def readInfo(self, file='info.txt'):
        f=open(file)
        lines=f.readlines()
        for line in lines:
            temp=line.split()
            if temp[0].startswith('pix'):
                self.pixel_size=float(temp[1])
            elif temp[0].startswith('sym'):
                self.num_pfs=int(temp[1])
            elif temp[0].startswith('num'):
                self.num_starts=float(temp[1])
            elif temp[0].startswith('dimer'):
                self.dimer_repeat_dist=float(temp[1])
            elif temp[0].startswith('hel'):
                self.helical_twist=float(temp[1])
        
        f.close()
        
        self.rise_per_subunit=(self.dimer_repeat_dist*self.num_starts)/self.num_pfs
        self.twist_per_subunit=(360-self.helical_twist)/self.num_pfs
                
        
        