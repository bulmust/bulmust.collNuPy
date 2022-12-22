"""
-----
ID      : initial.py
Author  : Taygun Bulmus
E-Mail  : bulmust@gmail.com
-----
"""
# ============================
import numpy as np
from math import gamma
from readFiles import technicalParameters_read, physicalParameters_read
from unitConv import erg_s_TO_MeV_km, muB_TO_Gauss__km
# ============================

# ============================
# Colors For Printing
class tLog:
    INFO  = '\033[94m'+'[INFO]    '+'\033[0m'
    OK      = '\033[92m'+'[OK]        '+'\033[0m'
    WARNING = '\033[93m'+'[WARNING]   '+'\033[0m'
    ERROR   = '\033[91m'+'[ERROR]     '+'\033[0m'
    EXIT    = '\033[91m'+'[EXIT]      initial.py'+'\033[0m'
# ============================

# ============================
# Read Parameters
technicalParametersDic= technicalParameters_read()
physicalParametersDic= physicalParameters_read()
# ============================

# ============================
# New Variables For Shorter Notation
hamOsc= physicalParametersDic['hamiltonian_oscillation']
hamMat= physicalParametersDic['hamiltonian_matter']
hamSelfSA= physicalParametersDic['hamiltonian_selfInteraction_singleAngle']
hamEM= physicalParametersDic['hamiltonian_electromagnetic']

# ============================
flavNum=physicalParametersDic['flavor_number']; 
# flavNum might be different from total number of flavors in the model 
# (e.g. 2+1 sterile model)
if flavNum in [1,2,3,4]:
    totFlav= int(flavNum)
# SterileBool 
if flavNum == 4:
    sterileBool = True
    print(tLog.WARNING+'Sterile neutrino is considered. '\
        'Flavors are "e", "mu", "tau", "sterile".')
else:
    sterileBool = False
# ============================

# ============================
Emod= physicalParametersDic['numberOf_energyMode']

deltaCP_Rad= physicalParametersDic['CP_Phase']

nuDistParam= physicalParametersDic['neutrino_distributionParameter']

#! ============================
if not sterileBool and nuDistParam == 5:
    # ERROR Sterile neutrinos are not considered, but neutrino_distributionParameter=5
    print(tLog.INFO+'Sterile neutrinos are not considered'\
        ', but neutrino_distributionParameter=5.')
    exit(tLog.EXIT)
elif sterileBool and nuDistParam != 5:
    # WARNING 
    print(tLog.WARNING+'Sterile neutrino is considered but neutrino distribution'
              ' is not "Fermi-Dirac for active neutrinos, Zero for sterile neutrino".')
#! ============================

# ============================
# Neutrinosphere Radius (also initial distance) [km]
ri_km= physicalParametersDic['distance_initial']; 
# Final Distance [km]
rf_km= physicalParametersDic['distance_final']
# All Distance
distAll_km = np.append(np.arange(\
    ri_km, rf_km, technicalParametersDic['holdData_Every']), rf_km)  # [km]
# To show intermediate distance
distTemp = 0;distTemp= distAll_km[0]
# ============================

# ============================
if hamEM:
    # Neutrino Magnetic Moment
    nuMagMom__Gauss_km=\
        physicalParametersDic['neutrino_MagneticMoment']* muB_TO_Gauss__km
    # (2*totFlav x 2*totFlav) BOOL
    dim_rho_2totFlav_Bool= True
else:
    # (2*totFlav x 2*totFlav) BOOL
    dim_rho_2totFlav_Bool= False
# ============================

# ============================
# Zero Arrays
E_MeV= np.zeros([int(Emod)])
# Luminosity and Temperature Arrays
luminosityArray_MeV_km= np.zeros(2*totFlav)
tempArray_MeV= np.zeros(2*totFlav)
# Distribution Matrices
distributionInit_nu_flav = np.zeros((Emod, totFlav, totFlav), dtype=complex)
distributionInit_nub_flav= np.zeros((Emod, totFlav, totFlav), dtype=complex)
rho= np.zeros((Emod, totFlav, totFlav), dtype=complex)
rhob= np.zeros((Emod, totFlav, totFlav), dtype=complex)
# TraceTerm [MeV^-1 * km^-3]
# Trace terms that should multiply to self interaction Hamiltonian
#   \sum_alpha[ (L_alpha/<E_alpha>)*f_alpha*(1/2*np.pi*SN_R**2) ]
TraceTerm  = np.zeros(Emod, dtype= np.complex)
TraceTermb = np.zeros(Emod, dtype= np.complex)
# ============================

# ============================
if dim_rho_2totFlav_Bool:
    # dimRhoMat = (1 x 2*(flavor) x 2*(flavor))
    dimRhoMat= (Emod, 2*totFlav, 2*totFlav)
    
    # If hamEM and hamSelfSA is considered, We will define gFactorSelf
    # See [Abbar:2020ggq] eq (10)
    if hamSelfSA:
        #gFactorSelfDagger = np.conjugate(np.transpose(gFactorSelf))    
        gFactorSelf= np.block([[ np.eye(totFlav, totFlav)\
            , np.zeros((totFlav, totFlav))], [np.zeros((totFlav, totFlav))\
            , -np.eye(totFlav, totFlav)]])    
else:
    # dimRhoMat = (neutrino/antineutrino) x (flavor) x (flavor))
    dimRhoMat= (2, Emod, totFlav, totFlav)
rhoInit_flav= np.zeros(dimRhoMat, dtype=complex)
# ============================

# ============================
# Energy Array [MeV]
E_MeV[0]= float(physicalParametersDic['energy_initial'])
for i1 in range(1, int(Emod)):
    # interval :(float(SNParam[2]) -float(SNParam[1]))/(energyMode-1)
    E_MeV[i1]= E_MeV[0]\
        + (((float(physicalParametersDic['energy_final']))- E_MeV[0])/ (Emod - 1)) * i1
# ============================

# ============================
if flavNum == 2:
    # Vacuum Mixing Angle
    theta_Rad= physicalParametersDic['theta12']

    # Mixing matrix
    U_mix = np.array([[np.cos(theta_Rad), np.sin(theta_Rad)],
                      [-np.sin(theta_Rad), np.cos(theta_Rad)]])
    
    # Delta M Square
    dMSqrHi_MeV2= physicalParametersDic['neutrino_hierarchy']\
        * physicalParametersDic['deltaM_square21']
    
    # Luminosity And Temperatures
    if nuDistParam == 1 or nuDistParam == 2:
        luminosityArray_MeV_km[0]= physicalParametersDic['luminosity_e']\
            * erg_s_TO_MeV_km
        luminosityArray_MeV_km[1]= physicalParametersDic['luminosity_mu']\
            * erg_s_TO_MeV_km
        luminosityArray_MeV_km[2]= physicalParametersDic['luminosity_eb']\
            * erg_s_TO_MeV_km
        luminosityArray_MeV_km[3]= physicalParametersDic['luminosity_mub']\
            * erg_s_TO_MeV_km
        tempArray_MeV[0]= physicalParametersDic['temperature_e']
        tempArray_MeV[1]= physicalParametersDic['temperature_mu']
        tempArray_MeV[2]= physicalParametersDic['temperature_eb']
        tempArray_MeV[3]= physicalParametersDic['temperature_mub']    
elif flavNum == 3:
    # Vacuum Mixing Angle
    theta12_Rad= physicalParametersDic['theta12']
    theta23_Rad= physicalParametersDic['theta23']
    theta13_Rad= physicalParametersDic['theta13']

    # Mixing matrix
    U_mix = (np.array([[1, 0, 0],
                       [0,  np.cos(theta23_Rad), np.sin(theta23_Rad)],
                       [0, -np.sin(theta23_Rad), np.cos(theta23_Rad)]])) @\
            (np.array([[np.cos(theta13_Rad), 0, np.sin(theta13_Rad)\
                * np.exp(-1j * deltaCP_Rad)],
                       [0, 1, 0],
                       [-np.sin(theta13_Rad)\
                * np.exp(1j * deltaCP_Rad), 0, np.cos(theta13_Rad)]])) @ \
            (np.array([[ np.cos(theta12_Rad), np.sin(theta12_Rad), 0],
                       [-np.sin(theta12_Rad), np.cos(theta12_Rad), 0],
                       [0, 0, 1]]))
    
    # Delta M Square
    dMSqr21_MeV2= physicalParametersDic['deltaM_square21']
    dMSqr32Hi_MeV2= physicalParametersDic['neutrino_hierarchy']\
        * physicalParametersDic['deltaM_square32']
    
    # Luminosity And Temperatures
    if nuDistParam == 1 or nuDistParam == 2:
        luminosityArray_MeV_km[0]= physicalParametersDic['luminosity_e']\
            * erg_s_TO_MeV_km
        luminosityArray_MeV_km[1]= physicalParametersDic['luminosity_mu']\
            * erg_s_TO_MeV_km
        luminosityArray_MeV_km[2]= physicalParametersDic['luminosity_tau']\
            * erg_s_TO_MeV_km
        luminosityArray_MeV_km[3]= physicalParametersDic['luminosity_eb']\
            * erg_s_TO_MeV_km
        luminosityArray_MeV_km[4]= physicalParametersDic['luminosity_mub']\
            * erg_s_TO_MeV_km
        luminosityArray_MeV_km[5]= physicalParametersDic['luminosity_taub']\
            * erg_s_TO_MeV_km
        tempArray_MeV[0]= physicalParametersDic['temperature_e']
        tempArray_MeV[1]= physicalParametersDic['temperature_mu']
        tempArray_MeV[2]= physicalParametersDic['temperature_tau']
        tempArray_MeV[3]= physicalParametersDic['temperature_eb']
        tempArray_MeV[4]= physicalParametersDic['temperature_mub']
        tempArray_MeV[5]= physicalParametersDic['temperature_taub']    
elif flavNum == 4:
    # Vacuum Mixing Angle
    theta12_Rad= physicalParametersDic['theta12']
    theta23_Rad= physicalParametersDic['theta23']
    theta13_Rad= physicalParametersDic['theta13']
    theta14_Rad= physicalParametersDic['theta_14']
    theta24_Rad= physicalParametersDic['theta_24']
    theta34_Rad= physicalParametersDic['theta_34']

    # Delta M Square
    dMSqr21_MeV2= physicalParametersDic['deltaM_square21']
    dMSqr32Hi_MeV2= physicalParametersDic['neutrino_hierarchy']\
        * physicalParametersDic['deltaM_square32']
    dMSqr41_MeV2= physicalParametersDic['deltaM_square41']
    
    # Luminosity And Temperatures
    if nuDistParam == 1 or nuDistParam == 2 or nuDistParam == 5:
        luminosityArray_MeV_km[0]= physicalParametersDic['luminosity_e']\
            * erg_s_TO_MeV_km
        luminosityArray_MeV_km[1]= physicalParametersDic['luminosity_mu']\
            * erg_s_TO_MeV_km
        luminosityArray_MeV_km[2]= physicalParametersDic['luminosity_tau']\
            * erg_s_TO_MeV_km
        luminosityArray_MeV_km[3]= physicalParametersDic['luminosity_sterile']\
            * erg_s_TO_MeV_km
        luminosityArray_MeV_km[4]= physicalParametersDic['luminosity_eb']\
            * erg_s_TO_MeV_km
        luminosityArray_MeV_km[5]= physicalParametersDic['luminosity_mub']\
            * erg_s_TO_MeV_km
        luminosityArray_MeV_km[6]= physicalParametersDic['luminosity_taub']\
            * erg_s_TO_MeV_km
        luminosityArray_MeV_km[7]= physicalParametersDic['luminosity_sterileb']\
            * erg_s_TO_MeV_km

        tempArray_MeV[0]= physicalParametersDic['temperature_e']
        tempArray_MeV[1]= physicalParametersDic['temperature_mu']
        tempArray_MeV[2]= physicalParametersDic['temperature_tau']
        tempArray_MeV[3]= physicalParametersDic['temperature_sterile']
        tempArray_MeV[4]= physicalParametersDic['temperature_eb']
        tempArray_MeV[5]= physicalParametersDic['temperature_mub']
        tempArray_MeV[6]= physicalParametersDic['temperature_taub']
        tempArray_MeV[7]= physicalParametersDic['temperature_sterileb']
    
    # Mixing matrix
    # U_mix = R34 @ R24 @ R14 @ R23 @ R13 @ R12 => [Barry:2011wb] eqn (2)
    U_mix = (np.array([[1, 0, 0, 0],\
                       [0, 1, 0, 0],\
                       [0, 0,  np.cos(theta34_Rad), np.sin(theta34_Rad)],\
                       [0, 0, -np.sin(theta34_Rad), np.cos(theta34_Rad)]])) @\
            (np.array([[1, 0, 0, 0],\
                       [0, np.cos(theta24_Rad), 0, np.sin(theta24_Rad)],\
                       [0, 0, 1, 0],\
                       [0, -np.sin(theta24_Rad), 0, np.cos(theta24_Rad)]])) @\
            (np.array([[np.cos(theta14_Rad), 0, 0, np.sin(theta14_Rad)],\
                       [0, 1, 0, 0],\
                       [0, 0, 1, 0],\
                       [-np.sin(theta14_Rad), 0, 0, np.cos(theta14_Rad)]])) @\
            (np.array([[1, 0, 0, 0],\
                       [0,  np.cos(theta23_Rad), np.sin(theta23_Rad), 0],\
                       [0, -np.sin(theta23_Rad), np.cos(theta23_Rad), 0],\
                       [0, 0, 0, 1]])) @\
            (np.array([[np.cos(theta13_Rad), 0, np.sin(theta13_Rad)\
                * np.exp(-1j * deltaCP_Rad), 0],\
                       [0, 1, 0, 0],\
                       [-np.sin(theta13_Rad)\
                * np.exp(1j * deltaCP_Rad), 0, np.cos(theta13_Rad), 0],\
                       [0, 0, 0, 1]])) @\
            (np.array([[ np.cos(theta12_Rad), np.sin(theta12_Rad), 0, 0],\
                       [-np.sin(theta12_Rad), np.cos(theta12_Rad), 0, 0],\
                       [0, 0, 1, 0],\
                       [0, 0, 0, 1]]))
# ============================

#! ============================
# ERROR Temperature Errors
if 0 in tempArray_MeV and nuDistParam != 5:
    if nuDistParam == 1:
        print(tLog.INFO+'Fermi Dirac Distribution'\
        ': (E^2)/(F_2(0) T^3 (exp(E/T) +1)) <== [Duan:2006an] eqn (12)')
        print(tLog.ERROR+'Temperature of one of the neutrinos is zero. '\
            'Any of temperatures can not be zero '\
            'for Fermi Dirac Distribution.')
        exit(tLog.EXIT)
    elif nuDistParam == 2:
        print(tLog.INFO+'Pinched Dirac Distribution'\
            ': (((alpha+1)/<E>)^alpha)*(E^alpha/Gamma(alpha+1))*exp(-(alpha+1)*E/<E>) '\
            '[Keil:2002in]')
        print(tLog.ERROR+'Temperature of one of the neutrinos is zero. '\
            'Any of temperatures can not be zero '
            'for Pinched Dirac Distribution.')
        exit(tLog.EXIT)
if nuDistParam == 5:
    if physicalParametersDic['temperature_e'] == 0\
        or physicalParametersDic['temperature_mu'] == 0\
        or physicalParametersDic['temperature_tau'] == 0\
        or physicalParametersDic['temperature_eb'] == 0\
        or physicalParametersDic['temperature_mub'] == 0\
        or physicalParametersDic['temperature_taub'] == 0:
        print(tLog.INFO+'Fermi Dirac Distribution'\
            ': (E^2)/(F_2(0) T^3 (exp(E/T) +1)) <== [Duan:2006an] eqn (12)')
        print(tLog.ERROR+'Temperature of one of the active neutrinos is zero. '\
            'Any of temperatures can not be zero '
            'for Fermi Dirac Distribution.')
        exit(tLog.EXIT)
    if physicalParametersDic['temperature_sterile'] != 0\
        or physicalParametersDic['temperature_sterileb'] != 0:
        print(tLog.ERROR+'One of the sterile neutrinos temperature is not zero.')
        exit(tLog.EXIT)
#! ============================

# ============================
# Distribution, Density Matrix and its trace
# \rho = (Distribution x (L/E) x (1/(4\pi^2 ri_km^{2})))/ (Diagonal_Term)
# Diagonal term is for normalization.
for i2 in range(Emod):
    # Distributions
    for i3 in range(totFlav):
        if physicalParametersDic['neutrino_distributionParameter'] == 1:
            # Fermi Dirac Distribution [1/MeV] <= [Duan:2006an] eqn (12),
            # Fermi Dirac Integral F2=1.803;
            distributionInit_nu_flav [i2, i3, i3]=\
                (1/ (1.803* (tempArray_MeV[i3]**3)))\
                * (E_MeV[i2]**2) / (1 + np.exp(E_MeV[i2]/ tempArray_MeV[i3]))
            distributionInit_nub_flav[i2, i3, i3]=\
                (1/ (1.803* (tempArray_MeV[i3+totFlav]**3)))\
                * (E_MeV[i2]**2) / (1 + np.exp(E_MeV[i2]/ tempArray_MeV[i3+totFlav]))

            rho[i2, i3, i3] =\
                (1 / (4 * (ri_km**2) * np.pi**2))* distributionInit_nu_flav[i2,i3, i3]\
                * luminosityArray_MeV_km[i3]/ (3.1514 * tempArray_MeV[i3])
            rhob[i2, i3, i3]=\
                (1/ (4 * (ri_km**2) * np.pi**2))* distributionInit_nub_flav[i2, i3, i3]\
                *luminosityArray_MeV_km[i3+totFlav]/ (3.1514* tempArray_MeV[i3+totFlav])
        elif physicalParametersDic['neutrino_distributionParameter'] == 2:
            # Pinched Distribution [1/MeV]
            # Mean Energy MeanE=<E>=F(3)/F(2) T = 3.1514 T
            # Gamma Function Gamma(alpha+1)= alpha!
            # alphaRho = (<E^2> - 2* <E>^2)/(<E>^2 - <E^2>)
            #          = (F4 - 2*(F3^2)/F2)/((F3^2)/F2 - F4)
            # Fermi Dirac Integrals:F2= 1.803, F3= 5.682, F4= 23.331
            # alpha = 2.302
            alphaRho = 3
            Gamma = gamma(alphaRho + 1)
            # Distribution w.r.t E
            distributionInit_nu_flav[i2, i3, i3]=\
                 (((alphaRho + 1)/ (3.1514 * tempArray_MeV[i3]))**(alphaRho + 1))\
                * ((E_MeV[i2]**(alphaRho)) / Gamma)* (np.exp(- (alphaRho + 1)\
                * E_MeV[i2]/ (3.1514 * tempArray_MeV[i3])))
            distributionInit_nub_flav[i2, i3, i3]=\
                 (((alphaRho + 1)/(3.1514* tempArray_MeV[i3+totFlav]))**(alphaRho + 1))\
                * ((E_MeV[i2]**(alphaRho)) / Gamma)* (np.exp(- (alphaRho + 1)\
                * E_MeV[i2]/ (3.1514 * tempArray_MeV[i3+totFlav])))

            rho[i2, i3, i3] =\
                (1 / (4 * (ri_km**2) * np.pi**2))* distributionInit_nu_flav[i2,i3, i3]\
                * luminosityArray_MeV_km[i3] / (3.1514 * tempArray_MeV[i3])
            rhob[i2, i3, i3]=\
                (1/ (4* (ri_km**2) * np.pi**2))* distributionInit_nub_flav[i2, i3, i3]\
                *luminosityArray_MeV_km[i3+totFlav]/ (3.1514* tempArray_MeV[i3+totFlav])
        elif physicalParametersDic['neutrino_distributionParameter'] == 3:
            # Electron Box Distribution [Dimensionless]
            distributionInit_nu_flav[i2, 0, 0] = 1
            distributionInit_nub_flav[i2, 0, 0] = 1
            # Self interaction needs coefficient
            if hamSelfSA:
                rho[i2, 0, 0] =\
                    (1/ (4 * (ri_km**2) * np.pi**2))* distributionInit_nu_flav[i2,0, 0]\
                    * luminosityArray_MeV_km[i3] / E_MeV[i2]
                rhob[i2, 0, 0]=\
                    (1/ (4 * (ri_km**2) * np.pi**2))* distributionInit_nu_flav[i2,0, 0]\
                    * luminosityArray_MeV_km[i3+totFlav] / E_MeV[i2]
            else:
                rho[i2, 0, 0]  = 1
                rhob[i2, 0, 0] = 1
        elif physicalParametersDic['neutrino_distributionParameter'] == 4:
            # Muon Box Distribution [Dimensionless]
            distributionInit_nu_flav[i2, 1, 1] = 1
            distributionInit_nub_flav[i2, 1, 1] = 1
            # Self interaction needs coefficient
            if hamSelfSA:
                rho[i2, 1, 1]=\
                    (1/ (4 * (ri_km**2) * np.pi**2))* distributionInit_nu_flav[i2,0, 0]\
                    * luminosityArray_MeV_km[i3] / E_MeV[i2]
                rhob[i2, 1, 1]=\
                    (1/ (4 * (ri_km**2) * np.pi**2))* distributionInit_nu_flav[i2,0, 0]\
                    * luminosityArray_MeV_km[i3+totFlav] / E_MeV[i2]
            else:
                rho[i2, 1, 1]  = 1
                rhob[i2, 1, 1] = 1
        elif physicalParametersDic['neutrino_distributionParameter'] == 5:
            # Sterile neutrinos values are set to zero
            # Sterile neutrinos are always assumed to be at the end of the array
            if i3 != totFlav-1:
                # Active Fermi-Dirac [1/MeV], Sterile Zero Distribution
                distributionInit_nu_flav[i2, i3, i3] = (1 / (1.803\
                    * (tempArray_MeV[i3]**3)))* (E_MeV[i2]**2)\
                    / (1 + np.exp(E_MeV[i2] / tempArray_MeV[i3]))
                distributionInit_nub_flav[i2, i3, i3] = (1 / (1.803\
                    * (tempArray_MeV[i3+totFlav]**3)))* (E_MeV[i2]**2)\
                    / (1 + np.exp(E_MeV[i2] / tempArray_MeV[i3+totFlav]))

                rho[i2, i3, i3]  = (1 / (4 * (ri_km**2) * np.pi**2))\
                    * distributionInit_nu_flav[i2,i3, i3]\
                    * luminosityArray_MeV_km[i3]/ (3.1514 * tempArray_MeV[i3])
                rhob[i2, i3, i3] = (1 / (4 * (ri_km**2) * np.pi**2))\
                    * distributionInit_nub_flav[i2, i3, i3]\
                    * luminosityArray_MeV_km[i3+totFlav]\
                    / (3.1514 * tempArray_MeV[i3+totFlav])
        elif physicalParametersDic['neutrino_distributionParameter'] == 6:
            # Only Electron Box Distribution (Not ebar) [Dimensionless]
            distributionInit_nu_flav[i2, 0, 0] = 1
            if hamSelfSA:
                rho[i2, 0, 0]= (1 / (4 * (ri_km**2) * np.pi**2))\
                    * distributionInit_nu_flav[i2,0, 0]\
                    * luminosityArray_MeV_km[i3] / E_MeV[i2]
            else:
                rho[i2, 0, 0]= 1
        elif physicalParametersDic['neutrino_distributionParameter'] == 7:
            # Only Electron Box Distribution (Not ebar) [Dimensionless]
            distributionInit_nu_flav[i2, 0, 0] = 1
            if hamSelfSA:
                rhob[i2, 0, 0] = (1 / (4 * (ri_km**2) * np.pi**2))\
                    * distributionInit_nu_flav[i2,0, 0]\
                    * luminosityArray_MeV_km[i3]/ E_MeV[i2]
            else:
                rhob[i2, 0, 0] = 1
    
    # TraceTerm
    # Trace terms that should multiply to self interaction Hamiltonian
    # \sum_i[(L_i/<E_i>)*f_i] *(1/(4*(SN_R^2)*mt.pi^2))
    if dim_rho_2totFlav_Bool:
        # If the rho matrix is totFlav*2 x totFlav*2, normalization must include
            # both neutrinos and anti-neutrinos together
            # I keep it two different variable of traceTerms which are equal Due
            # to adequate run for the rest of code
        TraceTerm[i2]  = np.trace(rho[i2] + rhob[i2])
        TraceTermb[i2] = np.trace(rho[i2] + rhob[i2])

        # Normalize density matrix
        rho[i2]  = rho[i2]  / TraceTerm[i2]
        rhob[i2] = rhob[i2] / TraceTermb[i2]
    
        # Set Initial Density Matrix
        rhoInit_flav[:, 0:totFlav, 0:totFlav]= rho[i2]
        rhoInit_flav[:, totFlav:totFlav* 2, totFlav:totFlav* 2]= rhob[i2]
    else:
        TraceTerm[i2]  = np.trace(rho[i2])
        TraceTermb[i2] = np.trace(rhob[i2])

        # Normalize density matrix
        rho[i2]  = rho[i2]  / TraceTerm[i2]
        rhob[i2] = rhob[i2] / TraceTermb[i2]
        
        # Set Initial Density Matrix
        rhoInit_flav[0,i2]= rho[i2]
        rhoInit_flav[1,i2]= rhob[i2]
# ============================

# ============================
# Warnings
if hamSelfSA:
    if Emod < 100:
        print(tLog.WARNING+'Self Interaction is considered but number of energy'
              ' mode is '+str(Emod)+'. Due to the integration over'
              ' energy in self interaction Hamiltonian, energyMode should'
              ' be increased regarding precision.')
# ============================
