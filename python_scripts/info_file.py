import gui.base as bs
import logging
import os


class infoFile():
    def __init__(self):
        #if os.path.isfile('info.txt'):
        #    print 'what'
        try:
            self.readInfo()
        except:
            logging.debug('No info file, initialize the project.')
            
    def readInfo(self):
        f=open('info.txt')
        f.readlines()
        f.close()
        
class newInfo():
    def __init__(self):
        self.getValues()
        self.writeInfo()
        
    def getValues(self):
        vals=self.startingValues()
        bs.gdGui('Initialize Project', **vals)
        
    def startingValues(self):
        vals={}
        vals['dimer_repeat_distance']=82
        vals['helical_twist']=0
        vals['num_starts']=3
        vals['pixel_size']=1.3
        return vals
    
    def writeInfo(self):
        return
                
        
        