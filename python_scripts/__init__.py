from .helix.basics import (expand_helical_parameters, get_helical_transform, 
                           get_helix_subunit_list, helical_wedge_mask, 
                           symmetrize_filament_by_subunit, symmetrize_filament_by_repeat, symmetrize_filament)
from .helix.mt import mt_lattice_params
from .util3d import (vol_intersect, spherical_cosmask, rotshift3D_spline,read_pdb,
                     find_binary_threshold, soft_mask)

import logging

logger = logging.getLogger(__name__)

def warn(msg):
    print "Warning:", msg

def info(msg):
    print "Info:", msg

def report_error(e):
    pass

# to Fix TclError: no display name and no $DISPLAY environment variable,
# with no gui interface needed
import __main__ as tmp_main
import matplotlib
if hasattr(tmp_main, '__file__'):
    matplotlib.use('Agg')
else:
    matplotlib.use('TkAgg')

