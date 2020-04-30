import cupy as np

def mt_lattice_params(num_pfs, num_starts, dimer_repeat_dist, radius13pf):

    """ Returns canonical microtubule parameters according to chretien and wade:

    twist_per_subunit, rise_per_subunit, mt_radius, axial_repeat_dist, twist_per_repeat = 
      mt_lattice_params(...)

    PLEASE NOTE:

    Angles are returned in degrees.

    The number of starts MUST be input in terms of 80A repeats, so 
    canonical microtubules will have 1.5 starts in this notation!

    Another point of confusion: twist_per_subunit and twist_per_repeat are returned
    according to a ***right-handed*** coordinate convention, which contrasts with the
    left-handed convention used by SPIDER, etc. for Euler angles. 
    """

    monomer_starts = num_starts*2.

    subunit_angle13 = np.arctan2(3./2.*dimer_repeat_dist, 2.*np.pi*radius13pf)
    delta_x = 2.*np.pi*radius13pf/(13.*np.cos(subunit_angle13))
    delta_n = dimer_repeat_dist * (3./2. * num_pfs / 13. - num_starts)
    theta = np.arctan2(delta_n,num_pfs*delta_x)

    pitch = num_starts*dimer_repeat_dist*np.cos(theta)

    rise_per_subunit = pitch/num_pfs
    subunit_angle = subunit_angle13 - theta

# mt_radius = num_pfs * delta_x * cos(subunit_angle13) / (2*pi);
#           = num_pfs * 2*pi*radius13pf/(13*cos(subunit_angle13))/(2*pi) / cos(subunit_angle13)
    mt_radius = num_pfs * radius13pf/13.
    twist_per_subunit = delta_x * np.cos(subunit_angle) / mt_radius

    twist_per_repeat = (twist_per_subunit * num_pfs - 2.*np.pi) / num_starts
    axial_repeat_dist = np.cos(theta) * dimer_repeat_dist

    return 180/np.pi * twist_per_subunit, rise_per_subunit, mt_radius, axial_repeat_dist, \
        180/np.pi * twist_per_repeat

