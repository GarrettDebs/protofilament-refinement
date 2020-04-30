import cupy as np
import cupyx.scipy.ndimage

import cupyx.scipy.fftpack

from tqdm import tqdm

def vol_intersect(v1_dim, v2_dim, v1_origin, v2_origin):
    """ This function is designed to help with summing, windowing and other operations with
    two volumes where the volumes are of different sizes and have different origins.

    v1_origin are the absolute coordinates of the center of v1 (i.e. v1_dim/2+1 in the v1 coordinate system)
    v2_origin are the absolute coordinates of the center of v2
    
    It finds the intersection of the two input volumes, and then returns:
      v_dim_intersect, the dimensions of the intersection volume
      v1_topleft, the coordinates of the top left corner of the intersection in v1 (where 0,0,0 is the v1 top left corner)
      v2_topleft, the coordinates of the top left corner of the intersection in v2 (where 0,0,0 is the v2 top left corner)
    """

    v_dim_intersect = np.zeros(3)
    v1_topleft = np.zeros(3)
    v2_topleft = np.zeros(3)
    v1_botright = np.zeros(3)
    v2_botright = np.zeros(3)

    for dim in np.arange(3):
#################
# The idea here is to redefine the boundaries of v1 as the boundaries for v2 and visa versa. 
# Then, truncate these as necessary to fit within the actual boundaries of v1 and v2
##################
#                       = (Vector from v2 ctr -> v1 ctf ) - v2 half width           + shift origin from v1 center to corner
        v1_topleft[dim] = v2_origin[dim] - v1_origin[dim] - np.floor(v2_dim[dim]/2) + np.floor(v1_dim[dim]/2)
#                       = (Vector from v1 ctr -> v2 ctf ) - v1 half width           + shift origin from v2 center to corner
        v2_topleft[dim] = v1_origin[dim] - v2_origin[dim] - np.floor(v1_dim[dim]/2) + np.floor(v2_dim[dim]/2)

        v1_botright[dim] = v1_topleft[dim] + v2_dim[dim]
        v2_botright[dim] = v2_topleft[dim] + v1_dim[dim]

        if v1_topleft[dim] < 0:
            v1_topleft[dim] = 0
        if v2_topleft[dim] < 0:
            v2_topleft[dim] = 0

        if v1_botright[dim] > v1_dim[dim]:
            v1_botright[dim] = v1_dim[dim]
        if v2_botright[dim] > v2_dim[dim]:
            v2_botright[dim] = v2_dim[dim]

        v_dim_intersect[dim] = v1_botright[dim] - v1_topleft[dim]
        
    return v_dim_intersect.astype(int), v1_topleft.astype(int), v2_topleft.astype(int)

def spherical_cosmask(n,mask_radius, edge_width, origin=None):
    """mask = spherical_cosmask(n, mask_radius, edge_width, origin)
    """

    if type(n) is int:
        n = np.array([n])

    sz = np.array([1, 1, 1])
    sz[0:np.size(n)] = n[:]

    szl = -np.floor(sz/2)
    szh = szl + sz

    x,y,z = np.meshgrid( np.arange(szl[0],szh[0]), 
                         np.arange(szl[1],szh[1]), 
                         np.arange(szl[2],szh[2]), indexing='ij', sparse=True)

    r = np.sqrt(x*x + y*y + z*z)

    m = np.zeros(sz)

#    edgezone = np.where( (x*x + y*y + z*z >= mask_radius) & (x*x + y*y + z*z <= np.square(mask_radius + edge_width)))

    edgezone = np.all( [ (x*x + y*y + z*z >= mask_radius), (x*x + y*y + z*z <= np.square(mask_radius + edge_width))], axis=0)
    m[edgezone] = 0.5 + 0.5*np.cos( 2*np.pi*(r[edgezone] - mask_radius) / (2*edge_width))
    m[ np.all( [ (x*x + y*y + z*z <= mask_radius*mask_radius) ], axis=0 ) ] = 1

#    m[ np.where(x*x + y*y + z*z <= mask_radius*mask_radius)] = 1

    return m

def rotshift3D_spline(v, phi=0, shifts=np.array([0,0,0]), mode='wrap'):
# With a nod to:
#  http://stackoverflow.com/questions/20161175/how-can-i-use-scipy-ndimage-interpolation-affine-transform-to-rotate-an-image-ab
    rot_origin = 0.5*np.array(v.shape)
    rot_rad = -phi*np.pi/180.0
    rot_matrix = np.array([[np.cos(rot_rad), np.sin(rot_rad), 0],
                           [-np.sin(rot_rad),np.cos(rot_rad), 0], 
                           [0                , 0              , 1]])
    offset = -(rot_origin-rot_origin.dot(rot_matrix)).dot(np.linalg.inv(rot_matrix))
    offset = offset - shifts

    transformed_v = scipy.ndimage.interpolation.affine_transform(v,rot_matrix,offset=offset,mode=mode)

    return transformed_v

# NOTE: annoyingly, we can't used ndimage.rotate because the origin is off by 0.5 pixel units...
#  the following works but is slow and fuzzy:
#                vol_rotshift[:,:,:] = scipy.ndimage.rotate(v_small, -rot, reshape=False, mode='wrap')
#                v_small[:,:,:] = scipy.ndimage.shift(vol_rotshift, [0, 0, frac_adjust], mode='wrap')
#                vol_rotshift[:,:,:] = v_small[:,:,:]

def read_pdb(infilename, chains=None):

    """chains option: String list of chain identifiers to include, eg 'ABEFG'
    """

#    atoms = read_pdb('kin_alf8A_fit.pdb')
#    atoms = e2.PointArray()
#    atoms.read_from_pdb('kin_alf8A_fit.pdb')

#    print len(atoms.get_points())

    try : infile=open(infilename,"r")
    except : raise IOError("%s is an invalid file name" %infilename)

    atomdefs={'H':(1.0,1.00794),'C':(6.0,12.0107),'A':(7.0,14.00674),'N':(7.0,14.00674),'O':(8.0,15.9994),'P':(15.0,30.973761),
              'S':(16.0,32.066),'W':(18.0,1.00794*2.0+15.9994),'AU':(79.0,196.96655) }

    aavg=[0,0,0]	# calculate atomic center
    amin=[1.0e20,1.0e20,1.0e20]		# coords of lower left front corner of bounding box
    amax=[-1.0e20,-1.0e20,-1.0e20]	# coords
    natm=0
    atoms=[]		# we store a list of atoms to process to avoid multiple passes
    nelec=0
    mass=0

	# parse the pdb file and pull out relevant atoms
    for line in infile:
        if (line[:4]=='ATOM'): # or (line[:6]=='HETATM' and options.het)) :
            if chains and not (line[21] in chains) : continue
            
            try:
                a=line[12:14].strip()
                aseq=int(line[6:11].strip())
                resol=int(line[22:26].strip())
                
                x=float(line[30:38])
                y=float(line[38:46])
                z=float(line[46:54])
            except:
                print "PDB Parse error:\n%s\n'%s','%s','%s'  '%s','%s','%s'\n"%(
                    line,line[12:14],line[6:11],line[22:26],line[30:38],line[38:46],line[46:54])
                print a,aseq,resol,x,y,z

            atoms.append((a,x,y,z))
                
            aavg[0]+=x
            aavg[1]+=y
            aavg[2]+=z
            natm+=1
				
            amin[0]=min(x,amin[0])
            amin[1]=min(y,amin[1])
            amin[2]=min(z,amin[2])
            amax[0]=max(x,amax[0])
            amax[1]=max(y,amax[1])
            amax[2]=max(z,amax[2])
	
            try:
                nelec+=atomdefs[a.upper()][0]
                mass+=atomdefs[a.upper()][1]
            except:
                print("Unknown atom %s ignored at %d"%(a,aseq))
							
    infile.close()
		
    print "%d atoms used with a total charge of %d e- and a mass of %d kDa"%(natm,nelec,mass/1000)
    print "atomic center at %1.1f,%1.1f,%1.1f (center of volume at 0,0,0)"%(aavg[0]/natm,aavg[1]/natm,aavg[2]/natm)
    print "Bounding box: x: %7.2f - %7.2f"%(amin[0],amax[0])
    print "              y: %7.2f - %7.2f"%(amin[1],amax[1])
    print "              z: %7.2f - %7.2f"%(amin[2],amax[2])

    return atoms

def get_center_of_pdb(pdb_file):
    from pdb_util import PDB
    pdb=PDB()
    pdb.read_file(pdb_file)
    com = pdb.center
    return com

# volume Tools

def write_vol(vol_file, vol_data):
    '''
        write numpy volume, or eman2 volume.
        Note the numpy volume is the transpose of EMAn2 volume. 
        Do not do the transpose of numpy volume before write
        
        Note: Numpy has indexing method called view. 
        arr.transpose() only return new x y z axis but the memory layout does not change.
        When numpy2em is called, if the array data type is float32(f), the memory is copied without changing the axis.
        But if the array data type is not float32, (np.zeros() return float64(d) array)
        the data is copied to a temporary array, 
        which reorders memory layout according to the view of the original array,
        in this case the coordinate is transposed.
        
        
    '''
    from EMAN2 import EMNumPy, IMAGE_SINGLE_SPIDER, EMData
    args = [vol_file]
    if vol_file.endswith(".spi"): # spider
        args.append(0)
        args.append(IMAGE_SINGLE_SPIDER)
    
    if type(vol_data) is np.ndarray:
        # the transpose is only affect the array that is not float32 type
        EMNumPy.numpy2em(vol_data.transpose()).write_image(*args)
    elif type(vol_data) is EMData:
        vol_data.write_image(*args)
    
def find_binary_threshold(v, vol_frac):

    frac_ind = np.int(np.ceil(v.size * (1-vol_frac)))
    return np.partition(v.ravel(), frac_ind)[frac_ind]

def soft_mask(v, voxel_size, num_subunit_residues,
              helical_repeat_distance=None, repeats_to_include=0,
              filter_resolution=20, expansion_factor=1.2, 
              expansion_radius=0, print_progress=True, return_mask=False):

    full_expansion_radius = expansion_radius + filter_resolution/2

# avg AA mol wt. in g/mol, density in g/cm3
    avg_aa_molwt = 110
    protein_density = 1.4

# 2
    print helical_repeat_distance
    v_thresh = np.zeros(v.shape)
    v_thresh[:] = v[:]

    sz = np.array(v.shape).astype(int)

    total_molwt = num_subunit_residues*avg_aa_molwt/6.023e23
    if helical_repeat_distance != None:
        total_molwt = total_molwt * sz[2]*voxel_size / helical_repeat_distance
    total_vol = np.prod(sz) * voxel_size**3                  # vol in A3
    mol_vol = total_molwt/protein_density / (1.0e-24)        # vol in A3
    mol_vol_frac = mol_vol/total_vol
    target_vol_frac = mol_vol_frac*expansion_factor

    thresh = find_binary_threshold(v_thresh, target_vol_frac)
    true_frac = (0.0 + np.sum(v_thresh >= thresh)) / v_thresh.size

    if repeats_to_include != 0:
        zdim = np.round(repeats_to_include * helical_repeat_distance/voxel_size)
    else:
        zdim = sz[2]

    if zdim > sz[2] - 4*np.ceil(filter_resolution/voxel_size):
        zdim = sz[2] - 4*np.ceil(filter_resolution/voxel_size)

    zdim = zdim.astype(int)

    v_thresh[:,:,0:np.floor(sz[2]/2).astype(int) - np.floor(zdim/2).astype(int)] = 0
    v_thresh[:,:,np.floor(sz[2]/2).astype(int) - np.floor(zdim/2).astype(int) + 1 + zdim - 1:] = 0
    v_thresh[v_thresh < thresh] = 0

    if print_progress:
        print 'Target volume fraction: {}'.format(target_vol_frac)
        print 'Achieved volume fraction: {}'.format(true_frac)
        print 'Designated threshold: {}'.format(thresh)

    progress_bar = tqdm(total=5)

    v_thresh = scipy.fftpack.fftn(v_thresh)
    progress_bar.update(1)
# 3
    cosmask_filter = np.fft.fftshift(spherical_cosmask(sz, 0, np.ceil(filter_resolution/voxel_size)))
    cosmask_filter = scipy.fftpack.fftn(cosmask_filter) / np.sum(cosmask_filter)
    progress_bar.update(1)

    v_thresh = v_thresh * cosmask_filter
    v_thresh = scipy.fftpack.ifftn(v_thresh)
    progress_bar.update(1)
    v_thresh = np.real(v_thresh)
    v_thresh[np.abs(v_thresh) < 10*np.finfo(type(v_thresh.ravel()[0])).eps] = 0

    v_thresh[v_thresh != 0] = 1

# The extent of blurring is equal to the diameter of the cosmask sphere; 
#  if we want this to equal the expected falloff for filter_resolution, 
#  we therefore need to divide filter_res by 4 to get the 
#  desired radius for spherical_cosmask.

    v_thresh = scipy.fftpack.fftn(v_thresh)
    progress_bar.update(1)

    v_thresh = v_thresh * cosmask_filter
    v_thresh = scipy.fftpack.ifftn(v_thresh)
    progress_bar.update(1)
    v_thresh = np.real(v_thresh)
    v_thresh[np.abs(v_thresh) < 10*np.finfo(type(v_thresh.ravel()[0])).eps] = 0

    if return_mask:
        v[:,:,:] = v_thresh
    else:
        v *= v_thresh

    return v_thresh

def synthetic_volume_helix(pdb_file, symm):
    '''
        synthetic volume from pdb, using symmetry
        @param symm: The symmetry matrix to apply to the pdb. 
            Symm should have N entries to symmetrize n subunits. 
            [[1,0,0], [0,1,0], [0,0,1]] is the unit 
        @param pdb: The pdb file used to synthesize the volume
    '''
    from cryofilia.pdb_util import PDB
    ps = []
    p = PDB.readFile(pdb_file) 
    for sym in symm:
        ps.append(p.symmetrize(sym))
    p = sum(ps)
    mrc = pdb2mrc(pdb)
    return mrc

