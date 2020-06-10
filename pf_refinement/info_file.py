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
        
        
        
class NewInfo():
    def __call__(self):
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
        vals['pixel_size']=1.3
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
                raise NameError('Missing input parameter, rerun initialization')
                
            f.write(' '.join([key, item+ '\n']))
            
        f.close()
                
        
        