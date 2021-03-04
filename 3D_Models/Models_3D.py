import numpy
import moments
import time


#Simultaneous Split Models


#1. Simultaneous Split No Migration
def sim_split_no_mig(params, ns):
    """
    Model with simultaneous split between pop 1, 2 and 3, gene flow does not occur.
    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    nu3: Size of population 3 after split.
	T1: The scaled time between the split and the present (in units of 2*Na generations).
	"""
    #4 parameters
    nu1,nu2,nu3, T1 = params
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1] + ns[2])
    fs = moments.Spectrum(sts)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1] + ns[2])
    fs = moments.Manips.split_2D_to_3D_2(fs, ns[1], ns[2])
    nomig = numpy.array([[0, 0], [0, 0]])
    fs.integrate([nu1, nu2, nu3], T1, m=nomig, dt_fac=0.01)
    return fs

def sim_split_sym_mig_adjacent(params, ns):
    """
    Model with simultaneous split between pop 1, 2 and 3, symmetric gene flow occurs between adjacent populations.
    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    nu3: Size of population 3 after split.
    m12: Migration rate between populations 1 and 2 (2*Na*m)
    m23: Migration rate between populations 2 and 3
	T1: The scaled time between the split and the present (in units of 2*Na generations).
	"""
    #6 parameters
    nu1,nu2,nu3,m12,m23,T1 = params
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1] + ns[2])
    fs = moments.Spectrum(sts)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1] + ns[2])
    fs = moments.Manips.split_2D_to_3D_2(fs, ns[1], ns[2])
    sym_mig = numpy.array([[0, m12, 0], [m12, 0, m23], [0, m23, 0]])
    fs.integrate([nu1, nu2, nu3], T1, m=sym_mig, dt_fac=0.01)
    return fs

def sim_split_asym_mig_adjacent(params, ns):
    """
    Model with simultaneous split between pop 1, 2 and 3, asymmetric gene flow occurs between adjacent populations.
    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    nu3: Size of population 3 after split.
    m12: Migration rate from population 2 to 1
    m21: Migration rate from population 1 to 2
    m23: Migration rate from population 3 to 2
    m32: Migration rate from population 2 to 3
    T1: The scaled time between the split and the present (in units of 2*Na generations).
    """
    # 8 parameters
    nu1, nu2, nu3, m12, m21, m23, m32, T1 = params
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1] + ns[2])
    fs = moments.Spectrum(sts)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1] + ns[2])
    fs = moments.Manips.split_2D_to_3D_2(fs, ns[1], ns[2])
    asym_mig_adj = numpy.array([[0, m12, 0],
                              [m21, 0, m23],
                              [0, m32, 0]])
    fs.integrate([nu1, nu2, nu3], T1, m=asym_mig_adj, dt_fac=0.01)
    return fs

def sim_split_sym_mig_all(params, ns):
    """
    Model with simultaneous split between pop 1, 2 and 3, symmetric gene flow occurs between all populations.
    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    nu3: Size of population 3 after split.
    m12: Migration rate from population 2 to 1
    m21: Migration rate from population 1 to 2
    m13: Migration rate from population 3 to 1
    m31: Migration rate from population 1 to 3
    m23: Migration rate from population 3 to 2
    m32: Migration rate from population 2 to 3
	T1: The scaled time between the split and the present (in units of 2*Na generations).
	"""
    #7 parameters
    nu1, nu2, nu3, m12, m13, m23, T1 = params
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1] + ns[2])
    fs = moments.Spectrum(sts)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1] + ns[2])
    fs = moments.Manips.split_2D_to_3D_2(fs, ns[1], ns[2])
    sym_mig_all = numpy.array([[0, m12, m13],
                             [m12, 0, m23],
                             [m13, m23, 0]])
    fs.integrate([nu1, nu2, nu3], T1, m=sym_mig_all, dt_fac=0.01)
    return fs

def sim_split_asym_mig_all(params, ns):
    """
    Model with simultaneous split between pop 1, 2 and 3, asymmetric gene flow occurs between all populations.
    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    nu3: Size of population 3 after split.
    m12: Migration rate from population 2 to 1
    m21: Migration rate from population 1 to 2
    m13: Migration rate from population 3 to 1
    m31: Migration rate from population 1 to 3
    m23: Migration rate from population 3 to 2
    m32: Migration rate from population 2 to 3
	T1: The scaled time between the split and the present (in units of 2*Na generations).
	"""
    #10 parameters
    nu1, nu2, nu3, m12, m21, m13, m31, m23, m32, T1 = params

    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1] + ns[2])
    fs = moments.Spectrum(sts)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1] + ns[2])
    fs = moments.Manips.split_2D_to_3D_2(fs, ns[1], ns[2])
    asym_mig_all = numpy.array([[0, m12, m13],
                               [m21, 0, m23],
                               [m31, m32, 0]])
    fs.integrate([nu1, nu2, nu3], T1, m=asym_mig_all, dt_fac=0.01)
    return fs


#Split models

#1. Split with no migration
def split_nomig(params, ns):
    """
    Model with split between pop 1 and (2,3), then split between 2 and 3.
    Migration does not occur between any population pair.
    nu1: Size of population 1 after split.
    nuA: Size of population (2,3) after split from 1.
    nu2: Size of population 2 after split.
    nu3: Size of population 3 after split.
    mA: Migration between ancestral populations
    T1: The scaled time between the split of pops 1 vs 2 and 3 (in units of 2*Na generations).
    T2: The scaled time between the split of pops 2 and 3 (in units of 2*Na generations).
    """
    # 7 parameters
    nu1, nuA, nu2, nu3, mA, T1, T2 = params

    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1] + ns[2])
    fs = moments.Spectrum(sts)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1] + ns[2])

    sym_mig = numpy.array([[0, mA], [mA, 0]])

    fs.integrate([nu1, nuA], T1, m=sym_mig)

    fs = moments.Manips.split_2D_to_3D_2(fs, ns[1], ns[2])
    
    nomig = numpy.array([[0,0,0],
                         [0,0,0],
                         [0,0,0]])

    fs.integrate([nu1, nu2, nu3], T2, m=nomig)

    return fs


#2. Split with migration the entire time
def split_sym_mig_all(params, ns):
    """
     Model with split between pop 1 and (2,3), then split between 2 and 3.
     Migration is symmetrical between all population pairs (ie 1<->2, 2<->3, and 1<->3).
     nu1: Size of population 1 after split.
     nuA: Size of population (2,3) after split from 1.
     nu2: Size of population 2 after split.
     nu3: Size of population 3 after split.
     mA: Migration rate between population 1 and population (2,3)
     m1: Migration rate between populations 1 and 2 (2*Na*m)
     m2: Migration rate between populations 2 and 3
     m3: Migration rate between populations 1 and 3
     T1: The scaled time between the split of pops 1 vs 2 and 3 (in units of 2*Na generations).
     T2: The scaled time between the split of pops 2 and 3 (in units of 2*Na generations).
     """

    # 10 parameters
    nu1, nuA, nu2, nu3, mA, m12, m23, m13, T1, T2 = params

    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1] + ns[2])
    fs = moments.Spectrum(sts)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1] + ns[2])

    sym_mig_1 = numpy.array([[0, mA], [mA, 0]])

    fs.integrate([nu1, nuA], T1, m=sym_mig_1)

    fs = moments.Manips.split_2D_to_3D_2(fs, ns[1], ns[2])

    sym_mig_2 = numpy.array([[0, m12, m13],
                             [m12, 0, m23],
                             [m13, m23, 0]])

    fs.integrate([nu1, nu2, nu3], T2, m=sym_mig_2)

    return fs

#3. Split with asymmetric migration
def split_asym_mig_all(params, ns):
    """
     Model with split between pop 1 and (2,3), then split between 2 and 3.
     Migration is symmetrical between all population pairs (ie 1<->2, 2<->3, and 1<->3).
     nu1: Size of population 1 after split.
     nuA: Size of population (2,3) after split from 1.
     nu2: Size of population 2 after split.
     nu3: Size of population 3 after split.
     mA: Migration rate between population 1 and population (2,3)
     m1: Migration rate between populations 1 and 2 (2*Na*m)
     m2: Migration rate between populations 2 and 3
     m3: Migration rate between populations 1 and 3
     T1: The scaled time between the split of pops 1 vs 2 and 3 (in units of 2*Na generations).
     T2: The scaled time between the split of pops 2 and 3 (in units of 2*Na generations).
     """

    # 13 parameters
    nu1, nuA, nu2, nu3, mA, m12, m13, m21, m23, m31, m32, T1, T2 = params

    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1] + ns[2])
    fs = moments.Spectrum(sts)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1] + ns[2])

    sym_mig_1 = numpy.array([[0, mA], [mA, 0]])

    fs.integrate([nu1, nuA], T1, m=sym_mig_1, dt_fac=0.01)

    fs = moments.Manips.split_2D_to_3D_2(fs, ns[1], ns[2])

    sym_mig_2 = numpy.array([[0, m12, m13],
                             [m21, 0, m23],
                             [m31, m32, 0]])

    fs.integrate([nu1, nu2, nu3], T2, m=sym_mig_2, dt_fac=0.01)

    return fs

#4. Split with symmetric migration between adjacent populations
def split_symmig_adjacent(params, ns):
    """
        Model with split between pop 1 and (2,3), then split between 2 and 3. Assume 2 occurs
        in between populations 1 and 3, which do not come in to contact with one another.
        Migration is symmetrical between 'adjacent' population pairs (ie 1<->2, 2<->3, but not 1<->3).
        nu1: Size of population 1 after split.
        nuA: Size of population (2,3) after split from 1.
        nu2: Size of population 2 after split.
        nu3: Size of population 3 after split.
        mA: Migration rate between population 1 and population (2,3)
        m12: Migration rate between populations 1 and 2 (2*Na*m)
    	m23: Migration rate between populations 2 and 3
        T1: The scaled time between the split of pops 1 vs 2 and 3 (in units of 2*Na generations).
        T2: The scaled time between the split of pops 2 and 3 (in units of 2*Na generations).
       """
    # 8 parameters
    nu1, nuA, nu2, nu3, mA, m23, T1, T2 = params
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1] + ns[2])
    fs = moments.Spectrum(sts)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1] + ns[2])

    sym_mig_1 = numpy.array([[0, mA], [mA, 0]])

    fs.integrate([nu1, nuA], T1, m=sym_mig_1, dt_fac=0.01)

    fs = moments.Manips.split_2D_to_3D_2(fs, ns[1], ns[2])

    sym_mig_2 = numpy.array([[0, 0, 0], [0, 0, m23], [0, m23, 0]])

    fs.integrate([nu1, nu2, nu3], T2, m=sym_mig_2, dt_fac=0.01)

    return fs

#5. Split with asymmetric migration between adjacent populations
def split_asymmig_adjacent(params, ns):
    """
        Model with split between pop 1 and (2,3), then split between 2 and 3. Assume 2 occurs
        in between populations 1 and 3, which do not come in to contact with one another.
        Migration is symmetrical between 'adjacent' population pairs (ie 1<->2, 2<->3, but not 1<->3).
       """
    # 9 parameters
    nu1, nuA, nu2, nu3, mAB, m23, m32, T1, T2 = params
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1] + ns[2])
    fs = moments.Spectrum(sts)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1] + ns[2])

    sym_mig_1 = numpy.array([[0, mAB], [mAB, 0]])

    fs.integrate([nu1, nuA], T1, m=sym_mig_1, dt_fac=0.01)

    fs = moments.Manips.split_2D_to_3D_2(fs, ns[1], ns[2])

    asym_mig_2 = numpy.array([[0, 0, 0],
                              [0, 0, m23],
                              [0, m32, 0]])

    fs.integrate([nu1, nu2, nu3], T2, m=asym_mig_2, dt_fac=0.01)

    return fs

