from pf_refinement import InfoFile, StarFile
from gui import gdGui
import numpy as np
import os

### Class for gathering the user input, and setting up your working directory
class SmoothGui(InfoFile):
    def __init__(self, info='info.txt'):
        super(SmoothGui, self).__init__()
    
    def __call__(self):
        ### Get Variables
        self.getVals()
        self.cores=int(self.vals['cores'])
        
        ### Set up working directory which will also be passed on
        self.vals['workdir'] = self.get_next_workdir('CFSmoothCoords');
    
    def getVals(self):
        ### If you'd like to add more user input options, add them here.
        ### Note you'll also need to add them to the argument parser under
        ### the exact same name in pf_smooth as well.
        vals={
            'star_file': 'XXXStarFileXXX',
            'micrograph_pixel_size': 1,
            'subunits_per_repeat': 2,
            'fit_order': 2,
            'keep_filament_id': False,
            'only_interpolate': False,
            'cores': 1
            }
        temp=gdGui('Smooth Coordinates', **vals)
        self.vals=temp.sendValues()
        
    def get_next_workdir(self, task_name, create=True):
        if not create:
            return "."
        if not os.path.exists(task_name):
            os.mkdir(task_name)
        for i in range(1,1000):
            if not os.path.exists(os.path.join(task_name, 'job%03d' % i)):
                try:
                    os.mkdir(os.path.join(task_name, 'job%03d' % i))
                except OSError as e:
                    if e.errno != 17: #  File exists
                        raise e
                break
        return os.path.join(task_name, 'job%03d' % i)
        

class SmoothCoords(StarFile):
    """ Smooth Coords instance. For command line processing
    """
    def __init__(self, **kwargs):
        """ Constructor. No description
        """
        #### Reads the star file and info file
        super(SmoothCoords, self).__init__(kwargs['star_file'])
        ### This will read variables like rise_per_subunit, num_pfs etc...
        self.readInfo()
        
        self._pdf = None
        self.project = None
        self.filament_offset = 100
        self.is_2d = True
        ### This is where all the input variables (as defined in getVals)
        ### will be stored. Just a heads up, all values will be strings
        self.args = kwargs
        
        ### Finally, initialize the mpi stuff
        try:
            from mpi4py import MPI
            
            self.comm=MPI.COMM_WORLD
            self.rank=self.comm.Get_rank()
            self.size=self.comm.Get_size()
        except Exception:
            self.rank=0
            self.size=1
            
        

    def __call__(self):
        ### Group microtubules
        try:
            self.mt=self.df.groupby(['MicrographName', 'HelicalTubeID'])
        except Exception:
            raise NameError('No metadata for HelicalTubeID, please add column '
                            'so we can properly group filaments.')
        
        ### Get the keys for each microtubule    
        keys=np.array(self.mt.groups.keys())
        fil_ids=np.arange(len(keys))
        ### Get the keys for the microtubule that this rank will run
        mt_keys=keys[fil_ids%self.size==self.rank]
        
        for mt_key in mt_keys:
            self.SmoothFunction(mt_key)
            
        if self.size>1:
            ### This is just a way to wait for all filaments to finish
            ### processing before moving onto the next step.
            if self.rank==0:
                ### Arbitrary array, just need to have something to send
                data=np.arange(2)
                for i in range(1, self.size):
                    req=self.comm.Isend(data, dest=i)
                    req.wait()
                
                ### This is where you would gather all the intermediate files
                self.combineIntermediateFiles()
                    
            else:
                data=np.empty(2)
                req=self.comm.Irecv(data,source=0)
                req.wait()
                
        else:
            ### If you're not running in parallel no need to wait for anything
            self.combineIntermediateFiles()
            
    def SmoothFunction(self, mt_key):
        '''
        Smooth one filaments
        '''
        
        metadata=self.mt.get_group(tuple(mt_key))
        data=self.get_values(metadata)
        
        ### Do the actual calculations with your stuff here.
        
    def get_values(self, metadata):
        """ Get coord values from m_in, as ndarray
        """
        scale=self.pixel_size/float(self.args['micrograph_pixel_size'])
        
        coord_x = (metadata.CoordinateX.astype(float) - 
                   metadata.OriginX.astype(float) * scale).values
                   
        coord_y = (metadata.CoordinateY.astype(float) - 
                   metadata.OriginY.astype(float) * scale).values
    
        try:
            coord_z = (metadata.CoordinateZ.astype(float) - 
                   metadata.OriginZ.astype(float) * scale).values
        except:
            coord_z = np.zeros(coord_y.shape)
        phi = metadata.AngleRot.astype(float).values
        theta = metadata.AngleTilt.astype(float).values
        psi = metadata.AnglePsi.astype(float).values
        return np.vstack((coord_x, coord_y, coord_z, phi, theta, psi)).T
    
    def combineIntermediateFiles(self):
        testfile=self.args['workdir']+'/test.txt'
        f=open(testfile,'w')
        f.write('hey it works!')
        f.close()



