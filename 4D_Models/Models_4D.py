import moments
import numpy

### Simultaneous Split Models ###

def sim_split_nomig_4D(params, ns):
    """
    Model with simultaneous split between populations, gene flow occurs btw all population.
    Split between pops 3 and 4 and gene flow occurs only between pops 1 and 2, 1 and 3, and 2 and 4
    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    nuA: Size of population (3,4) after split.
    nu3: Size of population 3 after split.
    nu4: Size of population 4 after split
    T1: The scaled time between the split of pops in units of 2*Na generations).
    """
    #5 parameters

    nu1, nu2, nu3, nu4 , T1 = params

    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1] + ns[2] +ns[3])
    fs = moments.Spectrum(sts)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1] + ns[2] + ns[3])
    fs = moments.Manips.split_2D_to_3D_2(fs, ns[1], ns[2] + ns[3])
    fs = moments.Manips.split_3D_to_4D_3(fs, ns[2], ns[3])

    no_mig = numpy.array([[0, 0, 0, 0],
                             [0, 0, 0, 0],
                             [0, 0, 0, 0],
                             [0, 0, 0, 0]])

    fs.integrate([nu1, nu2, nu3, nu4], T1, m=no_mig)

    return fs

def sim_split_all_sym_mig_4D(params, ns):
    """
    Model with simultaneous split between populations, gene flow occurs btw all population.
    Split between pops 3 and 4 and gene flow occurs only between pops 1 and 2, 1 and 3, and 2 and 4
    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    nuA: Size of population (3,4) after split.
    nu3: Size of population 3 after split.
    nu4: Size of population 4 after split
    m12: Migration rate between populations 1 and 2 during T2
    m13: Migration rate between populations 1 and 3 during T2
    m24: Migration rate between populations 2 and 4 during T2
    T1: The scaled time between the split of pops (in units of 2*Na generations).
    """
    # 11 parameters

    nu1, nu2, nu3, nu4, m12, m13, m14, m23, m24, m34, T1 = params

    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1] + ns[2] +ns[3])
    fs = moments.Spectrum(sts)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1] + ns[2] + ns[3])
    fs = moments.Manips.split_2D_to_3D_2(fs, ns[1], ns[2] + ns[3])
    fs = moments.Manips.split_3D_to_4D_3(fs, ns[2], ns[3])

    sym_mig = numpy.array([[0, m12, m13, m14],
                             [m12, 0, m23, m24],
                             [m13, m23, 0, m34],
                             [m14, m24, m34, 0  ]])

    fs.integrate([nu1, nu2, nu3, nu4], T1, m=sym_mig)

    return fs

def sim_split_all_asym_mig_4D(params, ns):
    """
    Model with simultaneous split between populations, gene flow occurs btw all populations.
    Split between pops 3 and 4 and gene flow occurs only between pops 1 and 2, 1 and 3, and 2 and 4
    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    nuA: Size of population (3,4) after split.
    nu3: Size of population 3 after split.
    nu4: Size of population 4 after split
    mAB: Migration rate between populations 1 and 2 during T1
    mAC: Migration rate between populations 1 and 3 during T1
    mBC: Migration rate between populations 2 and 3 during T1
    m12: Migration rate between populations 1 and 2 during T2
    m13: Migration rate between populations 1 and 3 during T2
    m24: Migration rate between populations 2 and 4 during T2
    T1: The scaled time between the split of pops (in units of 2*Na generations).
    """
    # 17 parameters

    nu1, nu2, nu3, nu4, m12, m13, m14, m21, m23, m24, m31, m32, m34, m41, m42, m43, T1 = params

    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1] + ns[2] +ns[3])
    fs = moments.Spectrum(sts)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1] + ns[2] + ns[3])
    fs = moments.Manips.split_2D_to_3D_2(fs, ns[1], ns[2] + ns[3])
    fs = moments.Manips.split_3D_to_4D_3(fs, ns[2], ns[3])

    asym_mig = numpy.array([[0, m12, m13, m14],
                             [m21, 0, m23, m24],
                             [m31, m32, 0, m34],
                             [m41, m42, m43, 0  ]])

    fs.integrate([nu1, nu2, nu3, nu4], T1, m=asym_mig)

    return fs

def sim_split_sym_mig_barrier_4D(params, ns):
    """
    Simaltaneous split, symmetric migration only occurs between Great Basin and Southern Populations 
    """
    # 7 parameters

    nu1, nu2, nu3, nu4, m12, m34, T1 = params

    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1] + ns[2] +ns[3])
    fs = moments.Spectrum(sts)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1] + ns[2] + ns[3])
    fs = moments.Manips.split_2D_to_3D_2(fs, ns[1], ns[2] + ns[3])
    fs = moments.Manips.split_3D_to_4D_3(fs, ns[2], ns[3])
    sym_mig_2 = numpy.array([[0, m12, 0,0],
			     [m12, 0, 0, 0],
                             [0, 0, m34,0],
                             [0, 0, 0,m34]])

    fs.integrate([nu1, nu2, nu3, nu4], T1, m=sym_mig_2)

    return fs

def sim_split_asym_mig_barrier_4D(params, ns):
    """
    Simaltaneous split, asymmetric migration only occurs between Great Basin and Southern Populations 
    """
    # 9 parameters

    nu1, nu2, nu3, nu4, m12, m21, m34, m43, T1 = params

    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1] + ns[2] +ns[3])
    fs = moments.Spectrum(sts)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1] + ns[2] + ns[3])
    fs = moments.Manips.split_2D_to_3D_2(fs, ns[1], ns[2] + ns[3])
    fs = moments.Manips.split_3D_to_4D_3(fs, ns[2], ns[3])
    asym_mig_2 = numpy.array([[0, m12, 0, 0],
                             [m21, 0, 0, 0],
                             [0, 0, m34, 0  ],
                             [0, 0, 0, m43  ]])

    fs.integrate([nu1, nu2, nu3, nu4], T1, m=asym_mig_2)

    return fs




### Split Models  ###

def split_sym_mig_4D(params, ns):
    """
    Symmetrical Migration between all modern populations 
    """

    # 19 Parameters 
    nu1, nu2, nuA, nu3, nuB, nu4, nuC, mAB, mAC, mBC, m12, m13, m14, m23, m24, m34, T1, T2, T3 = params

    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1] + ns[2] + ns[3])
    fs = moments.Spectrum(sts)
    # split S1 from the other populations, will result in [A,B], A will eventually become S1 
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1] + ns[2] + ns[3])
 
    # Symmetric migration between ancestral populations 
    sym_mig_1 = numpy.array([[0, mAB], [mAB, 0]])
    
    fs.integrate([nuA, nuB], T1, m=sym_mig_1)
    
    # the current order is [A,B], B = S1
    # Now we will split S2 from population A
    fs = moments.Manips.split_2D_to_3D_2(fs, ns[1], ns[2] + ns[3])
    sym_mig_2 = numpy.array([[0, mAB, mAC],
                             [mAB, 0, mBC],
                             [mAC, mBC, 0]])
    fs.integrate([nuA, nuB, nuC], T2, m=sym_mig_2) 
   
    # Split the E and West Pops 
    fs = moments.Manips.split_3D_to_4D_3(fs, ns[2], ns[3])

 
    sym_mig_3 = numpy.array([[0, m12, m13, m14],
                             [m12, 0, m23, m24],
                             [m13, m23, 0, m34  ],
                             [m14, m24, m34, 0  ]])

    # and now the final time epoch is integrated:
    fs.integrate([nu1, nu2, nu3, nu4], T3, m=sym_mig_3)
    
    return fs

def split_asym_mig_all_4D(params, ns):
    """
    Asymmetrical migration between all modern populations 
    """

        #25 Parameters 
    nu1, nu2, nuA, nu3, nuB, nu4, nuC, mAB, mAC, mBC, m12, m13, m14, m21, m23, m24, m31, m32, m34, m41, m42, m43, T1, T2, T3 = params

    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1] + ns[2] + ns[3])
    fs = moments.Spectrum(sts)
    # split S1 from the other populations, will result in [A,B], A will eventually become S1 
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1] + ns[2] + ns[3])
 
    # Symmetric migration between ancestral populations 
    sym_mig_1 = numpy.array([[0, mAB], [mAB, 0]])
    
    fs.integrate([nuA, nuB], T1, m=sym_mig_1)
    
    # the current order is [A,B], B = S1
    # Now we will split S2 from population A
    fs = moments.Manips.split_2D_to_3D_2(fs, ns[1], ns[2] + ns[3])
    sym_mig_2 = numpy.array([[0, mAB, mAC],
                             [mAB, 0, mBC],
                             [mAC, mBC, 0]])
    fs.integrate([nuA, nuB, nuC], T2, m=sym_mig_2) 
   
    # Split the E and West Pops 
    fs = moments.Manips.split_3D_to_4D_3(fs, ns[2], ns[3])

 
    asym_mig = numpy.array([[0, m12, m13, m14],
                             [m21, 0, m23, m24],
                             [m31, m32, 0, m34  ],
                             [m41, m42, m43, 0  ]])

    # and now the final time epoch is integrated:
    fs.integrate([nu1, nu2, nu3, nu4], T3, m=asym_mig)
    
    return fs

def split_nomig_4D(params, ns):
    """
    No migration between modern populations 
    """
        # 13 Parameters 
    nu1, nu2, nuA, nu3, nuB, nu4, nuC, mAB, mAC, mBC, T1, T2, T3 = params

    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1] + ns[2] + ns[3])
    fs = moments.Spectrum(sts)
    # split S1 from the other populations, will result in [A,B], A will eventually become S1 
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1] + ns[2] + ns[3])
 
    # Symmetric migration between ancestral populations 
    sym_mig_1 = numpy.array([[0, mAB], [mAB, 0]])
    
    fs.integrate([nuA, nuB], T1, m=sym_mig_1)
    
    # the current order is [A,B], B = S1
    # Now we will split S2 from population A
    fs = moments.Manips.split_2D_to_3D_2(fs, ns[1], ns[2] + ns[3])
    sym_mig_2 = numpy.array([[0, mAB, mAC],
                             [mAB, 0, mBC],
                             [mAC, mBC, 0]])
    fs.integrate([nuA, nuB, nuC], T2, m=sym_mig_2) 
   
    # Split the E and West Pops 
    fs = moments.Manips.split_3D_to_4D_3(fs, ns[2], ns[3])

 
    no_mig = numpy.array([[0, 0, 0, 0],
                             [0, 0, 0, 0],
                             [0, 0, 0, 0  ],
                             [0, 0, 0, 0  ]])

    # and now the final time epoch is integrated:
    fs.integrate([nu1, nu2, nu3, nu4], T3, m=no_mig)
    
    return fs

def split_symmig_barrier_4D(params, ns):
    """
    Symmetric migration between pops 3 & 4 and between pops 1 & 2 
    """
        # 15 Parameters 
    nu1, nu2, nuA, nu3, nuB, nu4, nuC, mAB, mAC, mBC, m12, m34, T1, T2, T3 = params

    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1] + ns[2] + ns[3])
    fs = moments.Spectrum(sts)
    # split S1 from the other populations, will result in [A,B], A will eventually become S1 
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1] + ns[2] + ns[3])
 
    # Symmetric migration between ancestral populations 
    sym_mig_1 = numpy.array([[0, mAB], [mAB, 0]])
    
    fs.integrate([nuA, nuB], T1, m=sym_mig_1)
    
    # the current order is [A,B], B = S1
    # Now we will split S2 from population A
    fs = moments.Manips.split_2D_to_3D_2(fs, ns[1], ns[2] + ns[3])
    sym_mig_2 = numpy.array([[0, mAB, mAC],
                             [mAB, 0, mBC],
                             [mAC, mBC, 0]])
    fs.integrate([nuA, nuB, nuC], T2, m=sym_mig_2) 
   
    # Split the E and West Pops 
    fs = moments.Manips.split_3D_to_4D_3(fs, ns[2], ns[3])

 
    sym_migbarrier = numpy.array([[0, m12, 0, 0],
                             [m12, 0, 0, 0],
                             [0, 0, 0, m34  ],
                             [0, 0, m34, 0  ]])

    # and now the final time epoch is integrated:
    fs.integrate([nu1, nu2, nu3, nu4], T3, m=sym_migbarrier)
    
    return fs

def split_asymmig_barrier_4D(params, ns):
    """
    Asymmetric migration between pops 3 & 4 and between pops 1 & 2 
    """
        # 17 Parameters 
    nu1, nu2, nuA, nu3, nuB, nu4, nuC, mAB, mAC, mBC, m12, m21, m34, m43, T1, T2, T3 = params

    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1] + ns[2] + ns[3])
    fs = moments.Spectrum(sts)
    # split S1 from the other populations, will result in [A,B], A will eventually become S1 
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1] + ns[2] + ns[3])
 
    # Symmetric migration between ancestral populations 
    sym_mig_1 = numpy.array([[0, mAB], [mAB, 0]])
    
    fs.integrate([nuA, nuB], T1, m=sym_mig_1)
    
    # the current order is [A,B], B = S1
    # Now we will split S2 from population A
    fs = moments.Manips.split_2D_to_3D_2(fs, ns[1], ns[2] + ns[3])
    sym_mig_2 = numpy.array([[0, mAB, mAC],
                             [mAB, 0, mBC],
                             [mAC, mBC, 0]])
    fs.integrate([nuA, nuB, nuC], T2, m=sym_mig_2) 
   
    # Split the E and West Pops 
    fs = moments.Manips.split_3D_to_4D_3(fs, ns[2], ns[3])

 
    asym_migbarrier = numpy.array([[0, m12, 0, 0],
                             [m21, 0, 0, 0],
                             [0, 0, 0, m34  ],
                             [0, 0, m43, 0  ]])

    # and now the final time epoch is integrated:
    fs.integrate([nu1, nu2, nu3, nu4], T3, m=asym_migbarrier)
    
    return fs




