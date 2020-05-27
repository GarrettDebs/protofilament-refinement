import numpy as np
import cupy as cp

from tqdm import tqdm

import pf_refinement.util3d as cf

def subunitsToUse(rise_per_subunit, num_pfs=13, num_starts=1.5, shift=True):
    index=0
    out=np.zeros([len(range(int(-2*num_starts),int(2*num_starts+1)))*num_pfs,3])
    
    for i in range(int(-2*num_starts),int(2*num_starts+1)):
        for j in range(0,num_pfs):
            z = i*num_pfs*rise_per_subunit/num_starts + (j)*rise_per_subunit;
            out[index,0]=i
            out[index,1]=j
            out[index,2]=z
            index += 1
            
    out=out[out[:,2].argsort()]
    if shift:
        lowerb=np.nonzero(out[:,2]>=-(rise_per_subunit*num_pfs/(num_starts*2)))
        upperb=np.nonzero(out[:,2]<(rise_per_subunit*num_pfs/(num_starts*2)))
    else:
        lowerb=np.nonzero(out[:,2]>=0)
        upperb=np.nonzero(out[:,2]<(rise_per_subunit*num_pfs/(num_starts)))
                          
    inter=np.intersect1d(lowerb,upperb)
    use=out[inter]
    use=use[use[:,1].argsort()]
    return use

def expand_helical_parameters(rise_per_subunit, twist_per_subunit, num_pfs=None, num_starts=1):
    """axial_repeat_dist, twist_per_repeat, supertwist, num_pfs, num_starts = expand_helical_parameters(...)

    This function returns extra parameters related to a specified helical symmetry. Specifically:

      - The helix is divided into "repeats", where each repeat has a specified number of protofilaments, 
        the repeats are separated by axial_repeat_dist, and rotate by twist_per_repeat.

    A minimal description of the helical lattice can be as simple as the rise and the twist.
    However, by introducing an additional parameter, num_starts, we can generalize the description of helices
    to include two additional features that are sometimes present but that cannot be captured simply
    by the rise and twist:
      (1) n-fold symmetry
      (2) the presence of a (single) seam

    It is important to note that if num_starts is not equal to 1, the resulting lattice geometry 
    (including, notably, the unit cell shape) will strongly depend on the
    number of protofilaments, which is in fact arbitrary!

    When dealing with the latter case (num_starts not equal to 1), it is therefore 
    important to explicitly specify the number of protofilaments when calling
    this function. If not specified, the number of protofilaments will be guessed by using 
    the specified twist to give the straightest possible protofilament (i.e. supertwist
    closest to zero).  However, this may or may not agree with the whims of nature...
    """

    if num_pfs is None:
        supertwist = None
        n_tmp = np.ceil(360 / abs(twist_per_subunit)) + 1
        for i in np.arange(1,n_tmp+1):
            dif = twist_per_subunit * i
            dif = dif - np.round(dif/360)*360
            if supertwist is None or abs(dif) < abs(supertwist):
                supertwist = dif
                num_pfs = i

    supertwist = twist_per_subunit * num_pfs
    supertwist = supertwist - round(float(supertwist)/360)*360

    axial_repeat_dist = rise_per_subunit * num_pfs / num_starts
    twist_per_repeat = supertwist / num_starts

    return axial_repeat_dist, twist_per_repeat, supertwist, num_pfs, num_starts

def get_helical_transform(n, rise_per_subunit, twist_per_subunit, num_pfs=None, num_starts=1):
    """ Given a particular helical symmetry, defined by the rise, twist, number of protofilaments and number of 
    starts, compute the shift and rotation corresponding to the n'th subunit.
    """
    axial_repeat_dist, twist_per_repeat, supertwist = \
        expand_helical_parameters(rise_per_subunit, twist_per_subunit, num_pfs, num_starts)[0:3]

    repeat_num = np.floor(n/num_pfs)
    pf_num = np.mod(n,num_pfs)

    z_shift = repeat_num * axial_repeat_dist + pf_num * rise_per_subunit
    rot =   repeat_num * twist_per_repeat  + pf_num * twist_per_subunit 

    return z_shift, rot, repeat_num, pf_num

def get_helix_subunit_list(v_dim, voxel_size, rise_per_subunit, twist_per_subunit, num_pfs, num_starts, ref_com=None, edge_width=0):
    """ Compute the list of subunits that lie within the z-range of a specified
    output volume.
    edge_width (in Angstroms) specifies extra padding to be included above and below the boundary subunits.
    Return values: 
       subunits: list of signed integers
       full_subunit_height: positive real number (units of voxels)

    """

    if ref_com is None:
        ref_com_z = 0
    else:
        ref_com_z = ref_com[2]

    axial_repeat_dist, twist_per_repeat, supertwist = \
        expand_helical_parameters(rise_per_subunit, twist_per_subunit, num_pfs, num_starts)[0:3]

###############
# Below, full_subunit_height captures the complete vertical dimension of a subunit wedge,
#  also accounting for its skewed shape (i.e. parallelogram rather than square when one
#  unrolls the cylindrical lattice onto a flat surface)
# Extra padding is also added to account for "edge_width", which specifies a soft-edged cosine mask
#  which broadens the full extend of the wedge at both the top and bottom.
###############
    full_subunit_height = axial_repeat_dist + (twist_per_subunit + twist_per_repeat)/twist_per_subunit * rise_per_subunit + 2*edge_width

    repeat_num = 0
    min_shift = 0
    max_shift = 0
    found_new_bounds = True

    while found_new_bounds:
        found_new_bounds = False
        for dir in np.array([-1, 1]):
            for pf_num in np.arange( 0,num_pfs):
                subunit_num = pf_num + dir * repeat_num * num_pfs
                shift = get_helical_transform(subunit_num, rise_per_subunit, twist_per_subunit, num_pfs, num_starts)[0]

                if ( (ref_com_z + shift - full_subunit_height/2) / voxel_size <= -np.floor(v_dim[2]/2) + v_dim[2]) and \
                        ( (ref_com_z + shift + full_subunit_height/2) / voxel_size >= -np.floor(v_dim[2]/2)):
                    if shift > max_shift:
                        found_new_bounds = True
                        repeat_num_maxz = dir*repeat_num
                        max_shift = shift
                    if shift < min_shift:
                        found_new_bounds = True
                        repeat_num_minz = dir*repeat_num
                        min_shift = shift

        repeat_num = repeat_num + 1

    if(repeat_num_minz <= repeat_num_maxz):
        first_subunit = repeat_num_minz*num_pfs
        last_subunit = (repeat_num_maxz + 1)*num_pfs
    else:
        first_subunit = repeat_num_maxz*num_pfs
        last_subunit = (repeat_num_minz + 1)*num_pfs

    generous_subunits = np.arange(first_subunit, last_subunit)

    tot_subunit_count = 0
    subunits = []
    for subunit_num in generous_subunits:
        shift = get_helical_transform(subunit_num, rise_per_subunit, twist_per_subunit, num_pfs, num_starts)[0]

        dest_shift, dest_rot = get_helical_transform(subunit_num, rise_per_subunit, twist_per_subunit, num_pfs, num_starts)[0:2]

        dest_center_z_pix = (ref_com_z + dest_shift)/voxel_size

        if ( dest_center_z_pix - full_subunit_height/2/voxel_size <= -np.floor(v_dim[2]/2) + v_dim[2]) and \
                (dest_center_z_pix + full_subunit_height/2/voxel_size >= -np.floor(v_dim[2]/2)):
            subunits.append(subunit_num)

    return subunits, full_subunit_height

def helical_wedge_mask(v_dim,voxel_size,rise_per_subunit,twist_per_subunit,
                       subunit_num=0,num_pfs=1,num_starts=1,
                       edge_resolution=20, outer_radius=None, ref_com=None, 
                       shifts=None, eulers=None):
    """Generates a wedge-shaped mask surrounding the i'th subunit of a helical lattice.
    The origin of the helical lattice is defined by ref_com; if not specified, the origin
    will lie on the z = 0 xy-plane, falling on the y-axis (Note: x-axis may make more sense?).
    subunit_num = 0 gives the wedge centered around ref_com.
    """

    fudge = 0.000001

    if type(v_dim) is int:
        v_dim = np.array([v_dim])

    if v_dim.size < 2:
        v_dim = np.array([v_dim[0], v_dim[0], v_dim[0]])
    if v_dim.size < 3:
        v_dim = np.array([v_dim[0], v_dim[1], v_dim[1]])

    if outer_radius is None:
        outer_radius = v_dim[0]/2 * voxel_size - edge_resolution

    if ref_com is None:
        ref_com = np.array([0, v_dim[1]/4, 0]) * voxel_size

    axial_repeat_dist, twist_per_repeat, supertwist = \
        expand_helical_parameters(rise_per_subunit, twist_per_subunit, num_pfs, num_starts)[0:3]

    shift, rot, repeat_num, pf_num = get_helical_transform(subunit_num, rise_per_subunit, twist_per_subunit, num_pfs, num_starts)

    wedge_angle_rad = 2*np.pi / num_pfs
    supertwist_rad = supertwist*np.pi/180;

    subunit_height_pix = axial_repeat_dist / voxel_size
    z_begin = shift / voxel_size

    m = np.zeros(v_dim)

    szl = -np.floor(v_dim/2)
    szh = szl + v_dim

#    X,Y,Z = np.mgrid[ szl[0]:szh[0], szl[1]:szh[1], szl[2]:szh[2]]
# NOTE: I have not yet figured out why 'xy' indexing does not work here!
#  The octave 'gold standard' (test_wedge2.m) confirms that 'ij' indexing gives the correct result
    X,Y,Z = np.meshgrid(np.arange(szl[0],szh[0]), 
                        np.arange(szl[1],szh[1]), 
                        np.arange(szl[2],szh[2]), indexing='ij')
#                        np.arange(szl[2],szh[2]), indexing='ij',sparse=True)

    if shifts is not None and (shifts != np.array([0,0,0])).any():
        X = X + shifts[0]/voxel_size
        Y = Y + shifts[1]/voxel_size
        Z = Z + shifts[2]/voxel_size

    if eulers is not None and (eulers != np.array([0,0,0])).any():
        rot = e2.Transform()
        rot.set_rotation({"type":"spider","phi":-eulers[2],"theta":-eulers[1],"psi":-eulers[0]})
        rot = np.array(rot.get_matrix()).reshape([3,4])[0:3, 0:3]

        xyz = np.reshape([X.flatten(), Y.flatten(), Z.flatten()], [3, X.size]).transpose()
        xyz = np.dot(xyz, rot)
        x = xyz[:,0].reshape(X.shape)
        y = xyz[:,1].reshape(Y.shape)
        z = xyz[:,2].reshape(Z.shape)
    else:
        x = X
        y = Y
        z = Z

    theta = np.arctan2(y,x)
    theta_begin = -np.arctan2(ref_com[1], ref_com[0]) + pf_num*wedge_angle_rad
    theta = np.mod(theta + theta_begin + np.pi, 2*np.pi) - np.pi;
    theta = theta + (z - ref_com[2]/voxel_size) / subunit_height_pix * supertwist_rad;
    theta = theta + fudge;            # This keeps away a weird 'fencepost' problem
                                      #  where single rows of pixels escape every wedge.
    m[ np.abs(theta) <= wedge_angle_rad/2] = 1;

    z = z - ref_com[2]/voxel_size + rise_per_subunit/voxel_size * theta/wedge_angle_rad

    m[ np.where(x*x + y*y > np.square(outer_radius/voxel_size)) ] = 0
    m[ np.where(np.abs(z - z_begin) > subunit_height_pix/2) ] = 0

    return m

def symmetrize_filament(vol, voxel_size, rise_per_subunit, twist_per_subunit, 
                        output_dim=None,
                        num_pfs=None, num_starts=1, 
                        source_subunit_num=None,
                        ref_com=None,subunit_pdb_model=None, outer_radius=None,
                        edge_resolution=20, interpolation_method='spline'):
    """ Uses wedge masks to rebuild a filament based on the subunit located
    at the 'origin' of the helical lattice, which is by default positioned
    in the xy plane at x=0, y=some positive value.

    ref_com can be used to redefine where the helix origin should start .
    
    In the current convention, subunit_num = 0 corresponds to the left side
    of any seam that is present, if twist_per_subunit is positive.

    source_subunit_num defines where, within the helical lattice, a wedge
    shaped copy should be taken, for replication in the symmetrized output.

    If a seam is present, you probably don't want to specify source_subunit_num = 0
    or source_subunit_num = num_pfs-2, because those particular subunits lie next to
    the seam and are thus likely to differ more strongly from the remaining num_pfs-2 subunits
    going around the helix.

    If num_pfs is specified, source_subunit_num by default is set to floor(num_pfs/2)
    because this is in some respects the 'safest' choice: far from the seam.  This 
    choice also coincides with where FREALIGN places the one 'good' protofilament
    when applying helical symmetry to filaments that have a seam.
    """

    v_dim = np.array(vol.shape)

    if output_dim is None:
        output_dim = v_dim
    else:
        if type(output_dim) is int:
            output_dim = np.array([output_dim]).astype(int)
        output_dim = output_dim.astype(int)

        if output_dim.size < 2:
            output_dim = np.array([output_dim[0], output_dim[0], output_dim[0]]).astype(int)
        if output_dim.size < 3:
            output_dim = np.array([output_dim[0], output_dim[1], output_dim[1]]).astype(int)

    if outer_radius is None:
        outer_radius = np.floor(output_dim[0]/2) * voxel_size

    edge_width = voxel_size * np.ceil(edge_resolution/(2*voxel_size))

    axial_repeat_dist, twist_per_repeat, supertwist = \
        expand_helical_parameters(rise_per_subunit, twist_per_subunit, num_pfs, num_starts)[0:3]

    diameter_voxels = 2 * (np.ceil(outer_radius / voxel_size))
    v_repeat_dim = np.array([diameter_voxels, diameter_voxels, 4 + 2*np.ceil(axial_repeat_dist/voxel_size/2 + 1)]).astype(int) + \
        2*(np.ceil(edge_width/voxel_size) + 10)

    v_repeat = symmetrize_filament_by_subunit(vol, voxel_size, rise_per_subunit, twist_per_subunit, 
                                              output_dim=v_repeat_dim,
                                              num_pfs=num_pfs, num_starts=num_starts, 
                                              source_subunit_num=source_subunit_num,
                                              ref_com=ref_com,subunit_pdb_model=subunit_pdb_model, outer_radius=outer_radius,
                                              edge_resolution=edge_resolution, interpolation_method=interpolation_method)[0]

    v_symm = symmetrize_filament_by_repeat(v_repeat, voxel_size, rise_per_subunit, twist_per_subunit, 
                                           output_dim=output_dim,
                                           num_pfs=num_pfs, num_starts=num_starts, 
                                           interpolation_method=interpolation_method)

    return v_symm

def symmetrize_filament_by_repeat(vol, voxel_size, rise_per_subunit, twist_per_subunit, 
                                  output_dim=None,
                                  num_pfs=None, num_starts=1, 
                                  interpolation_method='spline'):

    v_dim = np.array(vol.shape)

    if output_dim is None:
        output_dim = v_dim

    axial_repeat_dist, twist_per_repeat, supertwist = \
        expand_helical_parameters(rise_per_subunit, twist_per_subunit, num_pfs, num_starts)[0:3]

    v_repeat_dim = np.array([v_dim[0], v_dim[1], 2*10 + 4 + 2*np.ceil(axial_repeat_dist/voxel_size/2 + 1)]).astype(int)
    tight_repeat_dim = np.array([v_dim[0], v_dim[1], 4 + 2*np.ceil(axial_repeat_dist/voxel_size/2 + 1)]).astype(int)

    v_repeat = np.zeros(v_repeat_dim)
    v_temp = np.zeros(tight_repeat_dim)
    vol_rotshift = np.zeros(v_repeat_dim)

    source_dim_trim, source_topleft, source_small_topleft = \
        cf.vol_intersect(v_dim, v_repeat_dim, np.array([0,0,0]), np.array([0,0,0]))

    v_repeat[source_small_topleft[0]:source_small_topleft[0] + source_dim_trim[0],
             source_small_topleft[1]:source_small_topleft[1] + source_dim_trim[1],
             source_small_topleft[2]:source_small_topleft[2] + source_dim_trim[2]] = \
             vol[source_topleft[0] : source_topleft[0] + source_dim_trim[0], 
                 source_topleft[1] : source_topleft[1] + source_dim_trim[1], 
                 source_topleft[2] : source_topleft[2] + source_dim_trim[2]]
#    v_repeat_em = e2.EMNumPy.numpy2em(v_repeat.transpose())
#
    v_symm = np.zeros(output_dim)

    half_n_repeats = np.ceil(output_dim[2]*voxel_size / axial_repeat_dist / 2) + 1

    repeats = np.arange(-half_n_repeats, half_n_repeats + 1)

    with tqdm(total=len(repeats)) as progress_bar:
     for repeat in repeats:
         progress_bar.update(1)

         shift = repeat * axial_repeat_dist / voxel_size
         rot =   repeat * twist_per_repeat

         shift_int = round(shift)
         shift_frac = shift - shift_int

         dest_dim_trim, dest_topleft, dest_repeat_topleft = \
             cf.vol_intersect(tight_repeat_dim, v_repeat_dim, np.array([0,0,0]), np.array([0,0,0]))
         dest2_dim_trim, dest2_topleft, dest2_repeat_topleft = \
             cf.vol_intersect(output_dim, tight_repeat_dim, np.array([0,0,0]), np.array([0,0,shift_int]))

         if (dest2_dim_trim > 0).all():

             if interpolation_method is 'gridding':
                 vol_rotshift_em = sx.rot_shift3D_grid(v_repeat_em, phi=rot, sz=shift_frac)
                 vol_rotshift[:,:,:] = e2.EMNumPy.em2numpy(vol_rotshift_em).transpose()
             elif interpolation_method is 'nearest_neighbor':
                 vol_rotshift_em = sx.rot_shift3D(v_repeat_em, phi=rot, sz=shift_frac)
                 vol_rotshift[:,:,:] = e2.EMNumPy.em2numpy(vol_rotshift_em).transpose()
             elif interpolation_method is 'spline':
                 vol_rotshift[:,:,:] = cf.rotshift3D_spline(v_repeat, phi=rot, shifts=np.array([0, 0, shift_frac]), mode='reflect')

             v_temp[dest_topleft[0] : dest_topleft[0] + dest_dim_trim[0], 
                    dest_topleft[1] : dest_topleft[1] + dest_dim_trim[1], 
                    dest_topleft[2] : dest_topleft[2] + dest_dim_trim[2]] = \
                    vol_rotshift[dest_repeat_topleft[0] : dest_repeat_topleft[0] + dest_dim_trim[0], 
                                 dest_repeat_topleft[1] : dest_repeat_topleft[1] + dest_dim_trim[1], 
                                 dest_repeat_topleft[2] : dest_repeat_topleft[2] + dest_dim_trim[2]]
             v_symm[dest2_topleft[0] : dest2_topleft[0] + dest2_dim_trim[0], 
                    dest2_topleft[1] : dest2_topleft[1] + dest2_dim_trim[1], 
                    dest2_topleft[2] : dest2_topleft[2] + dest2_dim_trim[2]] = \
                    v_temp[dest2_repeat_topleft[0] : dest2_repeat_topleft[0] + dest2_dim_trim[0], 
                           dest2_repeat_topleft[1] : dest2_repeat_topleft[1] + dest2_dim_trim[1], 
                           dest2_repeat_topleft[2] : dest2_repeat_topleft[2] + dest2_dim_trim[2]]

    return v_symm

def symmetrize_filament_by_subunit(vol, voxel_size, rise_per_subunit, twist_per_subunit, 
                                   output_dim=None,
                                   num_pfs=None, num_starts=1, 
                                   source_subunit_num=None,
                                   ref_com=None,subunit_pdb_model=None, outer_radius=None,
                                   edge_resolution=20, interpolation_method='spline'):
    """ Uses wedge masks to rebuild a filament based on the subunit located
    at the 'origin' of the helical lattice, which is by default positioned
    in the xy plane at x=0, y=some positive value.

    ref_com can be used to redefine where the helix origin should start .
    
    In the current convention, subunit_num = 0 corresponds to the left side
    of any seam that is present, if twist_per_subunit is positive.

    source_subunit_num defines where, within the helical lattice, a wedge
    shaped copy should be taken, for replication in the symmetrized output.

    If a seam is present, you probably don't want to specify source_subunit_num = 0
    or source_subunit_num = num_pfs-2, because those particular subunits lie next to
    the seam and are thus likely to differ more strongly from the remaining num_pfs-2 subunits
    going around the helix.

    If num_pfs is specified, source_subunit_num by default is set to floor(num_pfs/2)
    because this is in some respects the 'safest' choice: far from the seam.  This 
    choice also coincides with where FREALIGN places the one 'good' protofilament
    when applying helical symmetry to filaments that have a seam.
    """

    v_dim = np.array(vol.shape).astype(int)

    if output_dim is None:
        output_dim = v_dim

    output_dim = output_dim.astype(int)

    if source_subunit_num is None:
        if num_pfs is None:
            source_subunit_num = 0
        else:
            source_subunit_num = np.floor(num_pfs/2)
    
    if ref_com is None:
        if (subunit_pdb_model is None):
            ref_com = np.array([0, output_dim[1]/4, 0]) * voxel_size
        else:
            ref_com = np.array([0, output_dim[1]/4, 0]) * voxel_size

    ref_com_z = ref_com[2]

    edge_width = voxel_size * np.ceil(edge_resolution/(2*voxel_size))

    if outer_radius is None:
        outer_radius = np.floor(output_dim[0]/2) * voxel_size

    axial_repeat_dist, twist_per_repeat, supertwist = \
        expand_helical_parameters(rise_per_subunit, twist_per_subunit, num_pfs, num_starts)[0:3]

################
# Compute the list of subunits that lie within the z-range of a specified output volume, 
#  and also the full z extent of a subunit wedge:
################

    subunits, full_subunit_height = \
        get_helix_subunit_list(output_dim, voxel_size, rise_per_subunit, twist_per_subunit, num_pfs, num_starts, ref_com, edge_width)

    print subunits
    if interpolation_method is 'gridding':
#  For gridding: need cubic volume
        v_dim_small = np.array([2*np.ceil( (outer_radius + 2*edge_width)/voxel_size), 
                                2*np.ceil( (outer_radius + 2*edge_width)/voxel_size),
                                2*np.ceil( (outer_radius + 2*edge_width)/voxel_size)])
    else:
# Non-gridding: can make smaller volume
        v_dim_small = np.array([2*np.ceil( (outer_radius + 2*edge_width)/voxel_size), 
                                2*np.ceil( (outer_radius + 2*edge_width)/voxel_size),
                                2*np.ceil( (full_subunit_height + 4*edge_width) / (2*voxel_size)) ])
    v_dim_small = v_dim_small.astype(int)

    v_symm = np.zeros(output_dim)
    v_temp = np.zeros(output_dim)
    vol_rotshift = np.zeros(v_dim_small)
    wedge = np.zeros(v_dim_small)
    wedge_sum = np.zeros(output_dim)

    cosmask_filter = np.fft.fftshift(cf.spherical_cosmask(v_dim_small, 0, edge_width / voxel_size))
    cosmask_filter = np.fft.fftn(cosmask_filter) / np.sum(cosmask_filter)

    source_shift, source_rot = get_helical_transform(source_subunit_num, rise_per_subunit, twist_per_subunit, num_pfs, num_starts)[0:2]
    source_z_pix = (ref_com_z + source_shift)/voxel_size
    source_z_pix_int = round(source_z_pix)

    source_dim_trim, source_topleft, source_small_topleft = \
        cf.vol_intersect(v_dim, v_dim_small, np.array([0,0,0]), np.array([0,0,source_z_pix_int]))
            
    v_small = np.zeros(v_dim_small)
    v_small[source_small_topleft[0]:source_small_topleft[0] + source_dim_trim[0],
            source_small_topleft[1]:source_small_topleft[1] + source_dim_trim[1],
            source_small_topleft[2]:source_small_topleft[2] + source_dim_trim[2]] = \
            vol[source_topleft[0] : source_topleft[0] + source_dim_trim[0], 
                source_topleft[1] : source_topleft[1] + source_dim_trim[1], 
                source_topleft[2] : source_topleft[2] + source_dim_trim[2]]
#    v_small_em = e2.EMNumPy.numpy2em(v_small.transpose())

    with tqdm(total=len(subunits)) as progress_bar:
     for subunit_num in subunits:

        dest_shift, dest_rot = get_helical_transform(subunit_num, rise_per_subunit, twist_per_subunit, num_pfs, num_starts)[0:2]
        dest_center_z_pix = (ref_com_z + dest_shift)/voxel_size
        dest_center_z_pix_int = round(dest_center_z_pix)

        progress_bar.update(1)

        net_shift = dest_shift - source_shift
        net_rot = dest_rot - source_rot
        frac_adjust = net_shift/voxel_size - (dest_center_z_pix_int - source_z_pix_int)
        
        dest_dim_trim, dest_topleft, dest_small_topleft = \
            cf.vol_intersect(output_dim, v_dim_small, np.array([0,0,0]), np.array([0,0,dest_center_z_pix_int]))

        if (dest_dim_trim > 0).all():

            if interpolation_method is 'gridding':
                vol_rotshift_em = sx.rot_shift3D_grid(v_small_em, phi=net_rot, sz=frac_adjust)
                vol_rotshift[:,:,:] = e2.EMNumPy.em2numpy(vol_rotshift_em).transpose()
            elif interpolation_method is 'spline':
                vol_rotshift[:,:,:] = cf.rotshift3D_spline(v_small, phi=net_rot, shifts=np.array([0, 0, frac_adjust]))
            elif interpolation_method is 'nearest_neighbor':
                vol_rotshift_em = sx.rot_shift3D(v_small_em, phi=net_rot, sz=frac_adjust)
                vol_rotshift[:,:,:] = e2.EMNumPy.em2numpy(vol_rotshift_em).transpose()
               # 
            wedge[:,:,:] = helical_wedge_mask(v_dim_small, 
                                              ref_com=ref_com,
                                              voxel_size=voxel_size,
                                              subunit_num=subunit_num,
                                              num_pfs=num_pfs,
                                              num_starts=num_starts, 
                                              rise_per_subunit=rise_per_subunit,
                                              twist_per_subunit=twist_per_subunit,
                                              outer_radius = outer_radius, 
                                              eulers = np.array([0, 0, 0]),
                                              shifts = np.array([0, 0, dest_center_z_pix_int*voxel_size]))

            wedge[:,:,:] = np.real(np.fft.ifftn( cosmask_filter * np.fft.fftn(wedge)))

            v_temp[:,:,:] = 0
            v_temp[dest_topleft[0] : dest_topleft[0] + dest_dim_trim[0], 
                   dest_topleft[1] : dest_topleft[1] + dest_dim_trim[1], 
                   dest_topleft[2] : dest_topleft[2] + dest_dim_trim[2]] = \
                   wedge[dest_small_topleft[0] : dest_small_topleft[0] + dest_dim_trim[0], 
                         dest_small_topleft[1] : dest_small_topleft[1] + dest_dim_trim[1], 
                         dest_small_topleft[2] : dest_small_topleft[2] + dest_dim_trim[2]] * \
                         vol_rotshift[dest_small_topleft[0] : dest_small_topleft[0] + dest_dim_trim[0], 
                                      dest_small_topleft[1] : dest_small_topleft[1] + dest_dim_trim[1], 
                                      dest_small_topleft[2] : dest_small_topleft[2] + dest_dim_trim[2]]

            wedge_sum[dest_topleft[0] : dest_topleft[0] + dest_dim_trim[0], 
                      dest_topleft[1] : dest_topleft[1] + dest_dim_trim[1], 
                      dest_topleft[2] : dest_topleft[2] + dest_dim_trim[2]] += \
                      wedge[dest_small_topleft[0] : dest_small_topleft[0] + dest_dim_trim[0], 
                            dest_small_topleft[1] : dest_small_topleft[1] + dest_dim_trim[1], 
                            dest_small_topleft[2] : dest_small_topleft[2] + dest_dim_trim[2]]
            v_symm += v_temp

    return v_symm, wedge_sum
