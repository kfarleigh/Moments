# Authors: Keaka Farleigh and Nicolas Finger
# email: keakafarleigh@gmail.com
# If you use these models please cite the appropriate papers (see https://github.com/kfarleigh/Moments)
import numpy
import moments
import time


#Simultaneous Split Models

def sim_split_no_mig(params, ns):
    """
    Model with simultaneous split between pop 1, 2 and 3, gene flow does not occur.
	"""
    #4 parameters
    nu1,nu2,nu3, T1 = params
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1] + ns[2])
    fs = moments.Spectrum(sts)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1] + ns[2])
    fs = moments.Manips.split_2D_to_3D_2(fs, ns[1], ns[2])
    nomig = numpy.array([[0, 0], [0, 0]])
    fs.integrate([nu1, nu2, nu3], T1, m=nomig)
    return fs

def sim_split_sym_mig_adjacent(params, ns):
    """
    Model with simultaneous split between pop 1, 2 and 3, symmetric gene flow occurs between adjacent populations.
	"""
    #6 parameters
    nu1,nu2,nu3,m12,m23,T1 = params
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1] + ns[2])
    fs = moments.Spectrum(sts)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1] + ns[2])
    fs = moments.Manips.split_2D_to_3D_2(fs, ns[1], ns[2])
    sym_mig = numpy.array([[0, m12, 0], [m12, 0, m23], [0, m23, 0]])
    fs.integrate([nu1, nu2, nu3], T1, m=sym_mig)
    return fs

def sim_split_asym_mig_adjacent(params, ns):
    """
    Model with simultaneous split between pop 1, 2 and 3, asymmetric gene flow occurs between adjacent populations.
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
    fs.integrate([nu1, nu2, nu3], T1, m=asym_mig_adj)
    return fs

def sim_split_sym_mig_all(params, ns):
    """
    Model with simultaneous split between pop 1, 2 and 3, symmetric gene flow occurs between all populations.
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
    fs.integrate([nu1, nu2, nu3], T1, m=sym_mig_all)
    return fs

def sim_split_asym_mig_all(params, ns):
    """
    Model with simultaneous split between pop 1, 2 and 3, asymmetric gene flow occurs between all populations.
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
    fs.integrate([nu1, nu2, nu3], T1, m=asym_mig_all)
    return fs


# Split models

def split_nomig(params, ns):
    """
    Model with split between pop 1 and (2,3), then split between 2 and 3.
    """
    # 7 parameters
    nu1, nuA, nu2, nu3, mA, T1, T2 = params

    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1] + ns[2])
    fs = moments.Spectrum(sts)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1] + ns[2])

    sym_mig_1 = numpy.array([[0, mA], [mA, 0]])

    fs.integrate([nu1, nuA], T1, m=sym_mig_1)

    fs = moments.Manips.split_2D_to_3D_2(fs, ns[1], ns[2])
    
    nomig = numpy.array([[0, 0], [0, 0]])

    fs.integrate([nu1, nu2, nu3], T2, m=nomig)

    return fs


def split_sym_mig_all(params, ns):
    """
     Model with split between pop 1 and (2,3), then split between 2 and 3.
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

def split_asym_mig_all(params, ns):
    """
     Model with split between pop 1 and (2,3), then split between 2 and 3.
     """

    # 13 parameters
    nu1, nuA, nu2, nu3, mA, m12, m13, m21, m23, m31, m32, T1, T2 = params

    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1] + ns[2])
    fs = moments.Spectrum(sts)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1] + ns[2])

    sym_mig_1 = numpy.array([[0, mA], [mA, 0]])

    fs.integrate([nu1, nuA], T1, m=sym_mig_1)

    fs = moments.Manips.split_2D_to_3D_2(fs, ns[1], ns[2])

    sym_mig_2 = numpy.array([[0, m12, m13],
                             [m21, 0, m23],
                             [m31, m32, 0]])

    fs.integrate([nu1, nu2, nu3], T2, m=sym_mig_2)

    return fs

def split_symmig_adjacent(params, ns):
    """
        Model with split between pop 1 and (2,3), then split between 2 and 3. Assume 2 occurs
        in between populations 1 and 3, which do not come in to contact with one another.
        Migration is symmetrical between 'adjacent' population pairs (ie 1<->2, 2<->3, but not 1<->3).
             """
    # 8 parameters
    nu1, nuA, nu2, nu3, mA, m23, T1, T2 = params
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1] + ns[2])
    fs = moments.Spectrum(sts)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1] + ns[2])

    sym_mig_1 = numpy.array([[0, mA], [mA, 0]])

    fs.integrate([nu1, nuA], T1, m=sym_mig_1)

    fs = moments.Manips.split_2D_to_3D_2(fs, ns[1], ns[2])

    sym_mig_2 = numpy.array([[0, 0, 0], [0, 0, m23], [0, m23, 0]])

    fs.integrate([nu1, nu2, nu3], T2, m=sym_mig_2)

    return fs


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

    fs.integrate([nu1, nuA], T1, m=sym_mig_1)

    fs = moments.Manips.split_2D_to_3D_2(fs, ns[1], ns[2])

    asym_mig_2 = numpy.array([[0, 0, 0],
                              [0, 0, m23],
                              [0, m32, 0]])

    fs.integrate([nu1, nu2, nu3], T2, m=asym_mig_2)

    return fs

