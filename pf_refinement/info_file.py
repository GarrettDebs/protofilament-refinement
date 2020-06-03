import gui.base as bs
import logging
import os


class InfoFile(object):
    def __init__(self, file='info.txt'):
        #if os.path.isfile('info.txt'):
        #    print 'what'
        try:
            self.readInfo(file)
        except:
            raise NameError('No info file, please initialize the project.')
            
    def readInfo(self, file):
        f=open(file)
        lines=f.readlines()
        for line in lines:
            temp=line.split()
            if temp[0].startswith('micro'):
                self.micrograph_pixel_size=float(temp[1])
            elif temp[0].startswith('sym'):
                self.num_pfs=int(temp[1])
            elif temp[0].startswith('num'):
                self.num_starts=float(temp[1])
            elif temp[0].startswith('dimer'):
                self.dimer_repeat_dist=float(temp[1])
            elif temp[0].startswith('hel'):
                self.helical_twist=float(temp[1])
        
        f.close()
        
        self.pixel_size=self.micrograph_pixel_size
        
        self.rise_per_subunit=(self.dimer_repeat_dist*self.num_starts)/self.num_pfs
        self.twist_per_subunit=(360-self.helical_twist)/self.num_pfs
        print 'micrograph_pixel_size: %g num_pfs: %d dimer_repeat: %g twist: %g'\
        %(self.micrograph_pixel_size, self.num_pfs, self.dimer_repeat_dist, self.helical_twist)
        
        
        
class NewInfo():
    def __init__(self):
        ###Start by getting the values to be written out
        self.getValues()
        ###Write out the info file
        self.writeInfo()
        
    def getValues(self):
        ###Initialize the starting parameters with the defaults
        vals=self.startingValues()
        ###Call the GUI to get user input values
        temp=bs.gdGui('Initialize Project', **vals)
        ###Retrieve the values that have been input from the user
        self.vals=temp.sendValues()
        
    def startingValues(self):
        ###Default starting values      
        vals={}
        vals['micrograph_pixel_size']=1.3
        vals['symmetry_type']=13
        vals['num_starts']=1.5
        vals['dimer_repeat_distance']=82
        vals['helical_twist']=0
        return vals
    
    def writeInfo(self):
        f=open('info.txt','w')
        ###Iterate through the variables and write them out
        for key, item in self.vals.items():
            ###Make sure there are no missing items
            if not item:
                logging.error('Missing input parameter, rerun initialization')
                exit()
                
            f.write(' '.join([key, item+ '\n']))
            
        f.close()
                
        
        