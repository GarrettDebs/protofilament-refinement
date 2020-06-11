from gui import gdGui
from collections import OrderedDict
from pf_refinement import StarFile

class Preprocess(StarFile):
    def __call__(self):
        ###Start by getting the values to be written out
        self.getValues()
        ###Write out the info file
        self.writeInfo()
        if self.vals['preprocess']:
            self.preprocess()
        
    def getValues(self):
        ###Initialize the starting parameters with the defaults
        vals=self.startingValues()
        ###Call the GUI to get user input values
        temp=gdGui('Initialize Project', **vals)
        ###Retrieve the values that have been input from the user
        self.vals=temp.sendValues()
        
    def startingValues(self):
        ###Default starting values      
        vals={
            'pixel_size': 1,
            'symmetry_type': 13,
            'num_starts': 1.5,
            'dimer_repeat_distance': 82,
            'helical_twist': 0,
            'preprocess': False,
            }
            
        return vals
    
    def writeInfo(self):
        f=open('info.txt','w')
        keys=['pixel_size', 'symmetry_type', 'num_starts', 
              'dimer_repeat_distance', 'helical_twist']
        ###Iterate through the variables and write them out
        for key in keys:
            item=self.vals[key]
            ###Make sure there are no missing items
            if not item:
                raise NameError('Missing input parameter, rerun initialization')
                
            f.write(' '.join([key, item+ '\n']))
        f.close()
    
    def preprocess(self):
        vals={
        'input_star_file': 'XXXSTARFILEXXX',
            'output_star_file': 'OUTPUT.star',
            'invert_polarity': False,
            'remove_symmetry': False
            }
        temp=gdGui('Initialize Project', **vals)
        ###Retrieve the values that have been input from the user
        self.vals=temp.sendValues()
        
        self.readStar(self.vals['input_star_file'])
        if self.vals['remove_symmetry']:
            self.df=self.df.iloc[::self.num_pfs]
            
        if self.vals['invert_polarity']:
            self.df.AnglePsi=(self.df.AnglePsi.astype(float)+180).astype(str)
        
        self.writeStar(self.vals['output_star_file'])
