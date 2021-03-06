###Make sure the python path is clear
import os

try: 
    os.environ['PYTHONPATH']
    print 'PYTHONPATH is set which can lead to some issues. Please run the '\
    'following command and try again.\n'
    print 'unset PYTHONPATH'
    exit()
except Exception:
    pass


from info_file import *
from starfile import *

from starmethods import *
from .helix.basics import (expand_helical_parameters, get_helical_transform, 
                           get_helix_subunit_list, helical_wedge_mask, 
                           symmetrize_filament_by_subunit, symmetrize_filament_by_repeat, symmetrize_filament)
from .helix.mt import mt_lattice_params
from .util3d import (vol_intersect, spherical_cosmask, rotshift3D_spline,read_pdb,
                     find_binary_threshold, soft_mask)

from patchmasks import *
from relion_subtract import *
from masked_classification import *