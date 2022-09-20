"""
-----
ID      : hamilt.py
Author  : Taygun Bulmus
E-Mail  : bulmust@gmail.com
-----
"""
import math as mt
import numpy as np
import initial as init
from unitConv import MeV_TO_km_1, avogadroNum_g_1, G_f_km2
import os as osCommand
from scipy.interpolate import interp1d
from scipy.integrate import simps

# ============================
# Colors For Printing
class tLog:
    STATUS  = '\033[94m'+'[STATUS]    '+'\033[0m'
    OK      = '\033[92m'+'[OK]        '+'\033[0m'
    WARNING = '\033[93m'+'[WARNING]   '+'\033[0m'
    ERROR   = '\033[91m'+'[ERROR]     '+'\033[0m'
    EXIT    = '\033[91m'+'[EXIT]      hamilt.py'+'\033[0m'
# ============================

# ============================
# Oscillation Hamiltonian
 # It is constant, so it should be outside of totalHamCalc class due to
 # speed consideration.
if init.hamOsc:
    # Define Initial Oscillation Hamiltonian [1/km]
    HAM_OSC_FLAV_KM_1 = np.zeros([init.Emod, init.totFlav, init.totFlav], dtype = np.complex128)
    # ============================
    # Two Flavor
    if init.flavNum == 2:
        for i1 in range(init.Emod):
            # Oscillation Hamiltonian In Flavor Base
            HAM_OSC_FLAV_KM_1[i1] = (MeV_TO_km_1 / (4.0 * init.E_MeV[i1]))* (init.U_mix @\
                np.array([[-init.dMSqrHi_MeV2, 0], [0, init.dMSqrHi_MeV2]]) @\
                np.transpose(init.U_mix))
    # ============================
    # Three Flavor
    elif init.flavNum == 3:
        for i2 in range(init.Emod):
            # Oscillation Hamiltonian In Flavor Base
            dMSqr31_MeV2= init.dMSqr21_MeV2+ init.dMSqr32Hi_MeV2
            HAM_OSC_FLAV_KM_1[i2] = (MeV_TO_km_1 / (6.0 * init.E_MeV[i2]))* (init.U_mix @\
                np.array([[-init.dMSqr21_MeV2 - dMSqr31_MeV2, 0, 0]\
                        , [0, init.dMSqr21_MeV2 - init.dMSqr32Hi_MeV2, 0]\
                        , [0, 0, dMSqr31_MeV2 + init.dMSqr32Hi_MeV2] ], dtype=complex)@\
                np.conjugate(np.transpose(init.U_mix)))
    # ============================
    # Four Flavor
    elif init.flavNum == 4:
        # Normal | Inverted
        # =4     | =4
        # =3     | =2
        # =2     | =1
        # =1     | =3
        for i2 in range(init.Emod):
            # Oscillation Hamiltonian In Flavor Base
            dMSqr31_MeV2= init.dMSqr21_MeV2+ init.dMSqr32Hi_MeV2
            dMSqr42_MeV2= init.dMSqr41_MeV2- init.dMSqr21_MeV2
            dMSqr43_MeV2= init.dMSqr41_MeV2- dMSqr31_MeV2
            
            HAM_OSC_FLAV_KM_1[i2] = (MeV_TO_km_1 / (8.0 * init.E_MeV[i2]))* (init.U_mix @\
                np.array([[-init.dMSqr21_MeV2- dMSqr31_MeV2- init.dMSqr41_MeV2, 0, 0, 0 ]\
                        , [0, init.dMSqr21_MeV2- init.dMSqr32Hi_MeV2- dMSqr42_MeV2, 0, 0]\
                        , [0, 0, dMSqr31_MeV2 + init.dMSqr32Hi_MeV2- dMSqr43_MeV2, 0    ]\
                        , [0, 0, 0, init.dMSqr41_MeV2+ dMSqr42_MeV2+ dMSqr43_MeV2       ] ], dtype=complex)@\
                np.conjugate(np.transpose(init.U_mix)))
    # ============================
    # If Majorana Interaction is used
    # HAM_OSC_FLAV_KM_1[Emod,totFlav*2,totFlav*2]\
    # = np.block([[ HAM_OSC_FLAV_KM_1   ,np.zeros((totFlav, totFlav))], [np.zeros((totFlav, totFlav)),HAM_OSC_FLAV_KM_1] ])
    if init.dim_rho_2totFlav_Bool:
        tmp = np.zeros([init.Emod, init.totFlav*2\
         , init.totFlav*2], dtype=np.complex128)
        tmp[:, 0:init.totFlav, 0:init.totFlav] = HAM_OSC_FLAV_KM_1
        tmp[:, init.totFlav:init.totFlav*2, init.totFlav:init.totFlav*2] = HAM_OSC_FLAV_KM_1
        HAM_OSC_FLAV_KM_1= tmp; del tmp
# ============================

# ============================
# Matter Interaction Hamiltonian
# It is outside due to speed consideration.
if init.hamMat:
    # ============================
    # Default Values
    if init.physicalParametersDic['use_defaultMatterProfile']:
        # ============================
        # Constant Profile
        if   init.physicalParametersDic['matterDensity_profile'] == 0:
            nBaryon_1_km3_LAMBDA =lambda dist : init.physicalParametersDic['matterDensity_initial']*avogadroNum_g_1/ ((1e-5)**3)
        # ============================
        # Exponential Profile
        elif init.physicalParametersDic['matterDensity_profile'] == 1.0:
            nBaryon_1_km3_LAMBDA =lambda dist : (init.physicalParametersDic['matterDensity_initial']*avogadroNum_g_1/ ((1e-5)**3))\
                * np.exp(-dist/ init.physicalParametersDic['matterDensity_exponentialDecay'])
        # ============================
        # Polynomial Profile
        elif init.physicalParametersDic['matterDensity_profile'] == 2.0:
            nBaryon_1_km3_LAMBDA =lambda dist : (init.physicalParametersDic['matterDensity_initial']*avogadroNum_g_1/ ((1e-5)**3))\
                * ((init.rf_km- dist)/ init.rf_km)**init.physicalParametersDic['matterDensity_polynomialDecay']
        # ============================
        # Electron fraction is constant.
        Y_e_LAMBDA = lambda dist : init.physicalParametersDic['electronFraction_constant']
    # Read SN data from file.
    else:
        # Full path of data file (red with extra space)
        SN_MATTER_PROFILE_DATA_DIR=\
            osCommand.getcwd() + '/backgroundProfiles/'+init.physicalParametersDic['matterProfile_fileName'][:-1]
        # Read SN data from file
        if osCommand.path.exists(SN_MATTER_PROFILE_DATA_DIR):
            print(tLog.WARNING+'The first three columns of file '
                  + SN_MATTER_PROFILE_DATA_DIR + ' have to be [distance [cm],'
                  ' Y_e, n_baryon [g/cm^3] ].')
            # ============================
            # Read distance from file [km]
            distRead_km = (np.loadtxt(SN_MATTER_PROFILE_DATA_DIR)[:, 0]) * (1e-5)
            # Read Y_e from file
            tmp0 = np.loadtxt(SN_MATTER_PROFILE_DATA_DIR)[:, 1]
            # Read n_baryon from file [g/km^3]
            tmp1 = np.loadtxt(SN_MATTER_PROFILE_DATA_DIR)[:, 2]*avogadroNum_g_1/((1e-5)**3)
            
            # ============================
            # Interpolate n_baryon and Y_e w.r.t red distance
            nBaryon_1_km3_LAMBDA = interp1d(distRead_km, tmp0, kind='linear'); del tmp0
            Y_e_LAMBDA = interp1d(distRead_km, tmp1, kind='linear'); del tmp1
        else:
            print(tLog.ERROR+'The file '+SN_MATTER_PROFILE_DATA_DIR+'does not exist.')
            print(tLog.ERROR+'Please check the file name and path.')
            exit(tLog.EXIT)
# ============================

# ============================
# Electromagnetic Interaction Hamiltonian
# It is outside due to speed consideration.
if init.hamEM:
    # ============================
    # Default Values
    if init.physicalParametersDic['use_defaultMagneticProfile']:
        # ============================
        # Constant Profile
        if   init.physicalParametersDic['magneticField_profile'] == 0:
            magMom_extMagField_km_1_LAMBDA=\
                lambda dist : init.nuMagMom__Gauss_km* init.physicalParametersDic['magneticField_initial']
        # ============================
        # Exponential Profile
        elif init.physicalParametersDic['magneticField_profile'] == 1.0:
            magMom_extMagField_km_1_LAMBDA=\
                lambda dist : init.nuMagMom__Gauss_km* init.physicalParametersDic['magneticField_initial']\
                * np.exp(-dist/ init.physicalParametersDic['magneticField_exponentialDecay'])
        # ============================
        # Polynomial Profile
        elif init.physicalParametersDic['magneticField_profile'] == 2.0:
            magMom_extMagField_km_1_LAMBDA=\
                lambda dist : init.nuMagMom__Gauss_km* init.physicalParametersDic['magneticField_initial']\
                * ((init.physicalParametersDic['magneticField_polynomialDecayDistance']\
                /dist)**init.physicalParametersDic['magneticField_polynomialDecayPower'])
    # Read SN data from file.
    else:
        # Full path of data file (red with extra space)
        SN_MAGNETIC_FIELD_PROFILE_DATA_DIR=\
            osCommand.getcwd() + '/backgroundProfiles/'+init.physicalParametersDic['magneticField_fileName'][:-1]
        # Read SN data from file
        if osCommand.path.exists(SN_MAGNETIC_FIELD_PROFILE_DATA_DIR):
            print(tLog.WARNING+'The first two columns of file '
                  + SN_MAGNETIC_FIELD_PROFILE_DATA_DIR + ' have to be [distance [cm],'
                  ' external magnetic field [Gauss] ].')
            # ============================
            # Read distance from file [km]
            distRead_km = (np.loadtxt(SN_MAGNETIC_FIELD_PROFILE_DATA_DIR)[:, 0]) * (1e-5)
            # Read External Field from file and multiply with neutrino magnetic moment
            tmp= init.nuMagMom__Gauss_km* np.loadtxt(SN_MAGNETIC_FIELD_PROFILE_DATA_DIR)[:, 1]
                        
            # ============================
            # Interpolate magMom_extMagField_km_1_LAMBDA red distance
            magMom_extMagField_km_1_LAMBDA = interp1d(distRead_km, tmp, kind='linear'); del tmp
        else:
            print(tLog.ERROR+'The file '+SN_MAGNETIC_FIELD_PROFILE_DATA_DIR+'does not exist.')
            print(tLog.ERROR+'Please check the file name and path.')
            exit(tLog.EXIT)
# ============================

# ============================
class hamiltonian:
    # ============================
    # Variables
    def __init__(self, rhoAll, dist):
        self.dist = dist
        self.rhoAll = rhoAll
    # ============================

    # ============================
    # Matter Interaction Hamiltonian (Traceless) [km^-1]
    # Two Flavor
    def hamMat2Flav_km_1(self):
        return (np.sqrt(2)/2)* G_f_km2* Y_e_LAMBDA(self.dist)\
            * nBaryon_1_km3_LAMBDA(self.dist)* np.array([[1, 0], [0, -1]])
    # Three Flavor
    def hamMat3Flav_km_1(self):
        return (np.sqrt(2)/3)* G_f_km2* Y_e_LAMBDA(self.dist)\
            * nBaryon_1_km3_LAMBDA(self.dist)\
            * np.array([[2, 0, 0], [0, -1, 0], [0, 0, -1]])
    # Four Flavor
    def hamMat4Flav_km_1(self):
    # Hamiltonian= diag[V_NC+ V_CC, V_NC, V_NC, 0]
        tmpYe= Y_e_LAMBDA(self.dist)
        return (np.sqrt(2)/8)* G_f_km2* nBaryon_1_km3_LAMBDA(self.dist)\
            * np.array([[7* tmpYe - 1, 0, 0, 0], [0, -tmpYe- 1, 0, 0]\
            , [0, 0, -tmpYe- 1, 0], [0, 0, 0, -5* tmpYe+ 3]])
    # ============================

    # ============================
    # Self Interaction Hamiltonian (Traceless) [km^-1]
    # Two or Three Flavor
    def hamSelfSA23Flav_km_1(self):
        coeff1= (mt.sqrt(2) * G_f_km2 * mt.pi) * (1 - mt.sqrt(1 - (init.ri_km / self.dist) ** 2)) ** 2

        # Define Zero Hamiltonian
        hamSelf_Flav = np.zeros((init.totFlav, init.totFlav), dtype= np.complex128)

        # Calculate Hamiltonian
        hamSelf_Flav = coeff1\
            * (simps(init.TraceTerm[:, np.newaxis, np.newaxis]* self.rhoAll[0], init.E_MeV, axis=0)
            - simps(init.TraceTermb[:, np.newaxis, np.newaxis]* self.rhoAll[1], init.E_MeV, axis=0))
        # Remove trace and return
        return hamSelf_Flav- (np.sum(np.diag(hamSelf_Flav))* np.eye(init.totFlav))/ init.totFlav    
    # Four Flavor
    def hamSelfSA4Flav_km_1(self):
        coeff1= (mt.sqrt(2) * G_f_km2 * mt.pi) * (1 - mt.sqrt(1 - (init.ri_km / self.dist) ** 2)) ** 2

        # Define the matrix of integral
        integral = np.zeros(np.shape(self.rhoAll), dtype=np.complex128)
        # Multiply Trace term to neutrino and anti neutrino densities
        # TraceTerm [1/(MeV*km^3)]
        # \sum_i[(L_i/<E_i>)*f_i] *(1/(4*(SN_R^2)*mt.pi^2))
        # L [MeV/km], <E> [MeV], f [1/MeV], SN_R [km], rhoAll [ ]
        integral[0]= init.TraceTerm[:, np.newaxis, np.newaxis]* self.rhoAll[0]
        integral[1]= init.TraceTermb[:, np.newaxis, np.newaxis]* self.rhoAll[1]
        # Multiply Sterile (Anti)Neutrinos with G_sCoeff
        integral[:, :, init.totFlav-1, init.totFlav-1]= integral[:, :, init.totFlav-1, init.totFlav-1]* init.physicalParametersDic['couplingConstant']
        # Make Non-Diagonal Terms Zero
        # No interaction between active and sterile
        integral[:, :, 0:init.totFlav-1, init.totFlav-1]= 0
        integral[:, :, init.totFlav-1, 0:init.totFlav-1]= 0

        # Calculate Hamiltonian
        hamSelf_Flav= coeff1\
            * (simps(integral[0], init.E_MeV, axis=0)
             - simps(integral[1], init.E_MeV, axis=0))
        # Remove trace and return
        return hamSelf_Flav- (np.sum(np.diag(hamSelf_Flav))* np.eye(init.totFlav))/ init.totFlav    
    # ============================

    # ============================
    # Total Hamiltonian
    # Two Flavor
    def totalHam2Flav_km_1(self):
        totHam_km_1= np.zeros((2, init.Emod, init.totFlav, init.totFlav), dtype=np.complex128)
        # Add Oscillation Hamiltonian [km^-1]
        if init.hamOsc:
            totHam_km_1[0]= totHam_km_1[0]+ HAM_OSC_FLAV_KM_1
            totHam_km_1[1]= totHam_km_1[1]+ HAM_OSC_FLAV_KM_1
        # Add Matter Interaction Hamiltonian [km^-1]
        if init.hamMat:
            hamMat_km_1= self.hamMat2Flav_km_1()
            totHam_km_1[0]= totHam_km_1[0]+ hamMat_km_1
            totHam_km_1[1]= totHam_km_1[1]- hamMat_km_1
        # Add Self Interaction Hamiltonian [km^-1]
        if init.hamSelfSA:
            hamSelfSA_km_1= self.hamSelfSA23Flav_km_1()
            totHam_km_1[0]= totHam_km_1[0]+ hamSelfSA_km_1
            totHam_km_1[1]= totHam_km_1[1]- hamSelfSA_km_1
        return totHam_km_1
    # Three Flavor
    def totalHam3Flav_km_1(self):
        totHam_km_1= np.zeros((2, init.Emod, init.totFlav, init.totFlav), dtype=np.complex128)
        # Add Oscillation Hamiltonian [km^-1]
        if init.hamOsc:
            totHam_km_1[0]= totHam_km_1[0]+ HAM_OSC_FLAV_KM_1
            totHam_km_1[1]= totHam_km_1[1]+ HAM_OSC_FLAV_KM_1
        # Add Matter Interaction Hamiltonian [km^-1]
        if init.hamMat:
            hamMat_km_1= self.hamMat3Flav_km_1()
            totHam_km_1[0]= totHam_km_1[0]+ hamMat_km_1
            totHam_km_1[1]= totHam_km_1[1]- hamMat_km_1
        # Add Self Interaction Hamiltonian [km^-1]
        if init.hamSelfSA:
            hamSelfSA_km_1= self.hamSelfSA23Flav_km_1()
            totHam_km_1[0]= totHam_km_1[0]+ hamSelfSA_km_1
            totHam_km_1[1]= totHam_km_1[1]- hamSelfSA_km_1
        return totHam_km_1
    # Four Flavor
    def totalHam4Flav_km_1(self):
        totHam_km_1= np.zeros((2, init.Emod, init.totFlav, init.totFlav), dtype=np.complex128)
        # Add Oscillation Hamiltonian [km^-1]
        if init.hamOsc:
            totHam_km_1[0]= totHam_km_1[0]+ HAM_OSC_FLAV_KM_1
            totHam_km_1[1]= totHam_km_1[1]+ HAM_OSC_FLAV_KM_1
        # Add Matter Interaction Hamiltonian [km^-1]
        if init.hamMat:
            hamMat_km_1= self.hamMat4Flav_km_1()
            totHam_km_1[0]= totHam_km_1[0]+ hamMat_km_1
            totHam_km_1[1]= totHam_km_1[1]- hamMat_km_1
        # Add Self Interaction Hamiltonian [km^-1]
        if init.hamSelfSA:
            hamSelfSA_km_1= self.hamSelfSA4Flav_km_1()
            totHam_km_1[0]= totHam_km_1[0]+ hamSelfSA_km_1
            totHam_km_1[1]= totHam_km_1[1]- hamSelfSA_km_1
        return totHam_km_1
# ============================

# ============================
# If dimension is totFlav*2 x totFlav*2
class hamiltonian_BigRho:
    # ============================
    # Variables
    def __init__(self, rhoAll, dist):
        self.dist = dist
        self.rhoAll = rhoAll
    # ============================

    # ============================
    # Matter Interaction Hamiltonian (Traceless) [km^-1]
    # Two Flavor
    def hamMat2Flav_km_1(self):
        # V_CC+V_NC
        tempCC_NC=(np.sqrt(2)* G_f_km2/ 2)* nBaryon_1_km3_LAMBDA(self.dist)* (3 * Y_e_LAMBDA(self.dist)- 1)
        # V_NC
        tempNC=-(np.sqrt(2)* G_f_km2/ 2)* nBaryon_1_km3_LAMBDA(self.dist)* (1 - Y_e_LAMBDA(self.dist))
        return np.array([[tempCC_NC, 0 , 0, 0]\
            , [0, tempNC, 0, 0], [0, 0, -tempCC_NC, 0]\
            , [0, 0, 0, -tempNC]])
    # Three Flavor
    def hamMat3Flav_km_1(self):
        # V_CC+V_NC
        tempCC_NC=(np.sqrt(2)* G_f_km2/ 2)* nBaryon_1_km3_LAMBDA(self.dist)* (3 * Y_e_LAMBDA(self.dist)- 1)
        # V_NC
        tempNC=-(np.sqrt(2)* G_f_km2/ 2)* nBaryon_1_km3_LAMBDA(self.dist)* (1 - Y_e_LAMBDA(self.dist))
        return np.array([[tempCC_NC, 0 , 0, 0, 0, 0]\
            , [0, tempNC, 0, 0, 0, 0], [0, 0, tempNC, 0, 0, 0]\
            , [0, 0, 0, -tempCC_NC,0, 0], [0, 0, 0, 0, -tempNC, 0]\
            , [0, 0, 0, 0, 0, -tempNC]])   
    # Four Flavor
    def hamMat4Flav_km_1(self):
        # V_CC+V_NC
        tempCC_NC=(np.sqrt(2)* G_f_km2/ 2)* nBaryon_1_km3_LAMBDA(self.dist)* (3 * Y_e_LAMBDA(self.dist)- 1)
        # V_NC
        tempNC=-(np.sqrt(2)* G_f_km2/ 2)* nBaryon_1_km3_LAMBDA(self.dist)* (1 - Y_e_LAMBDA(self.dist))
        return np.array([\
            [tempCC_NC, 0 , 0, 0, 0, 0, 0, 0]\
            , [0, tempNC, 0, 0, 0, 0, 0, 0]\
            , [0, 0, tempNC, 0, 0, 0, 0, 0]\
            , [0, 0, 0, 0, 0, 0, 0, 0]\
            , [0, 0, 0, 0, -tempCC_NC,0, 0, 0]\
            , [0, 0, 0, 0, 0, -tempNC, 0, 0]\
            , [0, 0, 0, 0, 0, 0, -tempNC, 0]\
            , [0, 0, 0, 0, 0, 0, 0, 0]])
    # ============================

    # ============================
    # ElectroMagnetic Interaction Hamiltonian (Traceless) [km^-1]
    # Two Flavor
    def hamEM2Flav_km_1(self):
        return magMom_extMagField_km_1_LAMBDA(self.dist)* np.array([[0, 0, 0, 1]\
            , [0, 0, -1, 0]\
            , [0, -1, 0, 0]\
            , [1, 0, 0, 0]])
    # Three Flavor
    def hamEM3Flav_km_1(self):
        # Sasaki et. al., https://arxiv.org/pdf/2106.02181.pdf
        # Eqn. (5)
        return magMom_extMagField_km_1_LAMBDA(self.dist)* np.array([\
            [0, 0, 0,  0,  1, 1]\
            , [0, 0, 0, -1,  0, 1]\
            , [0, 0, 0, -1, -1, 0]\
            , [0, -1, -1, 0, 0, 0]\
            , [1,  0, -1, 0, 0, 0]\
            , [1,  1,  0, 0, 0, 0]])
    # Four Flavor
    def hamEM4Flav_km_1(self):
        return magMom_extMagField_km_1_LAMBDA(self.dist)* np.array([\
            [0, 0, 0, 0,  0,  1, 1, 0]\
            , [0, 0, 0, 0, -1,  0, 1, 0]\
            , [0, 0, 0, 0, -1, -1, 0, 0]\
            , [0, 0, 0, 0,  0,  0, 0, 0]\
            , [0, -1, -1, 0, 0, 0, 0, 0]\
            , [1,  0, -1, 0, 0, 0, 0, 0]\
            , [1,  1,  0, 0, 0, 0, 0, 0]\
            , [0,  0,  0, 0, 0, 0, 0, 0]])
    # ============================

    # ============================
    # Self Interaction Hamiltonian (Traceless) [km^-1]
    # Two or Three Flavor
    def hamSelfSA23Flav_km_1(self):
        coeff1= (mt.sqrt(2) * G_f_km2 * mt.pi) * (1 - mt.sqrt(1 - (init.ri_km / self.dist) ** 2)) ** 2
        # Define Zero Hamiltonian
        hamSelf_Flav = np.zeros((init.totFlav* 2, init.totFlav* 2), dtype= np.complex128)
        # Define the matrix of integral
        integral = np.zeros(np.shape(self.rhoAll), dtype=np.complex128)
        
        # G Factor, see [Abbar:2020ggq] eq (10)
        gFactorSelf=\
            np.block([[ np.eye(init.totFlav, init.totFlav)    ,  np.zeros((init.totFlav, init.totFlav))]\
                        ,[ np.zeros((init.totFlav, init.totFlav)), -np.eye(init.totFlav, init.totFlav)]])

        # Rho_C is modified rho matrix, see [Abbar:2020ggq] eq (11)
        # Rho_C=[[rho[2:4,2:4] , rho[0:2,2:4]*]
        #        [rho[0:2,2:4]T. rho[0:2,0:2] ]]
        rho_C= self.rhoAll.copy()
        rho_C[:, 0:2, 0:2], rho_C[:, 2:4, 2:4]= rho_C[:, 2:4, 2:4].copy(), rho_C[:, 0:2, 0:2].copy()
        rho_C[:, 0:2, 2:4], rho_C[:, 2:4, 0:2]= np.conjugate(rho_C[:, 0:2, 2:4].copy()), np.conjugate(rho_C[:, 2:4, 0:2].copy())

        # Integral Part
        integral= (init.TraceTerm+ init.TraceTermb)\
            * (  (np.conjugate(np.transpose(gFactorSelf))* (self.rhoAll- rho_C)* gFactorSelf)\
            + (0.5)* np.conjugate(np.transpose(gFactorSelf))* np.trace((self.rhoAll- rho_C)* gFactorSelf))

        # Self Interaction Hamiltonian [km^-1]
        # After taking integral over energy, hamSelf_Flav dim. is [km^-1]
        # integral has (energyMod,totFlav*2,totFlav*2) dimension.
        # simps(integral, init.E_MeV, axis=0) means integrate
        # integral[:,0,0] and integral[:,0,1] and so on over energy.
        hamSelf_Flav = coeff1 * simps(integral, init.E_MeV, axis=0)

        # Remove trace and return
        return hamSelf_Flav- (np.sum(np.diag(hamSelf_Flav))* np.eye(init.totFlav* 2)) / init.totFlav*2    
    # Four Flavor
    def hamSelfSA4Flav_km_1(self):
        coeff1= (mt.sqrt(2) * G_f_km2 * mt.pi) * (1 - mt.sqrt(1 - (init.ri_km / self.dist) ** 2)) ** 2
        # Define Zero Hamiltonian
        hamSelf_Flav = np.zeros((init.totFlav* 2, init.totFlav* 2), dtype= np.complex128)
        # Define the matrix of integral
        integral = np.zeros(np.shape(self.rhoAll), dtype=np.complex128)
        
        # G Factor, see [Abbar:2020ggq] eq (10)
        gFactorSelf=\
            np.block([[ np.eye(init.totFlav, init.totFlav)    ,  np.zeros((init.totFlav, init.totFlav))]\
                        ,[ np.zeros((init.totFlav, init.totFlav)), -np.eye(init.totFlav, init.totFlav)]])

        # Rho_C is modified rho matrix, see [Abbar:2020ggq] eq (11)
        # Rho_C=[[rho[2:4,2:4] , rho[0:2,2:4]*]
        #        [rho[0:2,2:4]T. rho[0:2,0:2] ]]
        rho_C= self.rhoAll.copy()
        rho_C[:, 0:2, 0:2], rho_C[:, 2:4, 2:4]= rho_C[:, 2:4, 2:4].copy(), rho_C[:, 0:2, 0:2].copy()
        rho_C[:, 0:2, 2:4], rho_C[:, 2:4, 0:2]= np.conjugate(rho_C[:, 0:2, 2:4].copy()), np.conjugate(rho_C[:, 2:4, 0:2].copy())

        # Integral Part
        integral= (init.TraceTerm+ init.TraceTermb)\
            * (  (np.conjugate(np.transpose(gFactorSelf))* (self.rhoAll- rho_C)* gFactorSelf)\
            + (0.5)* np.conjugate(np.transpose(gFactorSelf))* np.trace((self.rhoAll- rho_C)* gFactorSelf))

        # Self Interaction Hamiltonian [km^-1]
        # After taking integral over energy, hamSelf_Flav dim. is [km^-1]
        # integral has (energyMod,totFlav*2,totFlav*2) dimension.
        # simps(integral, init.E_MeV, axis=0) means integrate
        # integral[:,0,0] and integral[:,0,1] and so on over energy.
        hamSelf_Flav = coeff1 * simps(integral, init.E_MeV, axis=0)

        # Sterile-Sterile Neutrino Self Interactions might be different then active-active self interactions
        tmp1= np.zeros(np.shape(hamSelf_Flav))
        # Neutrinos
        tmp1[init.totFlav-1,init.totFlav-1]     = tmp1[init.totFlav-1,init.totFlav-1]    *init.G_sCoeff
        # Anti-Neutrinos
        tmp1[init.totFlav*2-1,init.totFlav*2-1] = tmp1[init.totFlav*2-1,init.totFlav*2-1]*init.G_sCoeff

        # There is not any interactions between active and sterile neutrinos
        # 3+1 Neutrino Model
        # hamSelfMaj_Flav = [. . . 0 . . . 0]
        #                   [. . . 0 . . . 0]
        #                   [. . . 0 . . . 0]
        #                   [0 0 0 . 0 0 0 .]
        #                   [. . . 0 . . . 0]
        #                   [. . . 0 . . . 0]
        #                   [. . . 0 . . . 0]
        #                   [0 0 0 . 0 0 0 .]
        tF= init.totFlav
        tmp1[:, 0:tF- 1, 0:tF- 1], tmp1[:, tF:tF*2- 1, tF:tF*2- 1], tmp1[:, 0 :tF- 1  , tF:tF*2- 1], tmp1[:, tF:tF*2- 1, 0 :tF- 1]=\
            hamSelf_Flav[:, 0:tF- 1, 0:tF- 1], hamSelf_Flav[:, tF:tF*2- 1, tF:tF*2- 1], hamSelf_Flav[:, 0 :tF- 1  , tF:tF*2- 1], hamSelf_Flav[:, tF:tF*2- 1, 0 :tF- 1]
        tmp1[:, tF- 1, tF- 1], tmp1[:, tF*2- 1, tF*2- 1]= hamSelf_Flav[:, tF- 1, tF- 1], hamSelf_Flav[:, tF*2- 1, tF*2- 1]       
    
        # Remove trace and return
        return tmp1- (np.sum(np.diag(tmp1))* np.eye(init.totFlav*2)) / init.totFlav*2
    # ============================

    # ============================
    # Total Hamiltonian
    # Two Flavor
    def totalHam2Flav_km_1(self):
        totHam_km_1= np.zeros((init.Emod, init.totFlav*2, init.totFlav*2), dtype=np.complex128)
        # Add Oscillation Hamiltonian [km^-1]
        if init.hamOsc:
            totHam_km_1= totHam_km_1+ HAM_OSC_FLAV_KM_1
        # Add Matter Interaction Hamiltonian [km^-1]
        if init.hamMat:
            totHam_km_1= totHam_km_1+ self.hamMat2Flav_km_1()
        # Add ElectroMagnetic Interaction Hamiltonian [km^-1]
        if init.hamEM:
            totHam_km_1= totHam_km_1+ self.hamEM2Flav_km_1()
        # Add Self Interaction Hamiltonian [km^-1]
        if init.hamSelfSA:
            totHam_km_1= totHam_km_1- self.hamSelfSA23Flav_km_1()
        return totHam_km_1
    # Three Flavor
    def totalHam3Flav_km_1(self):
        totHam_km_1= np.zeros((init.Emod, init.totFlav*2, init.totFlav*2), dtype=np.complex128)
        # Add Oscillation Hamiltonian [km^-1]
        if init.hamOsc:
            totHam_km_1= totHam_km_1+ HAM_OSC_FLAV_KM_1
        # Add Matter Interaction Hamiltonian [km^-1]
        if init.hamMat:
            totHam_km_1= totHam_km_1+ self.hamMat3Flav_km_1()
        # Add ElectroMagnetic Interaction Hamiltonian [km^-1]
        if init.hamEM:
            totHam_km_1= totHam_km_1+ self.hamEM3Flav_km_1()
        # Add Self Interaction Hamiltonian [km^-1]
        if init.hamSelfSA:
            totHam_km_1= totHam_km_1- self.hamSelfSA23Flav_km_1()
        return totHam_km_1
    # Four Flavor
    def totalHam4Flav_km_1(self):
        totHam_km_1= np.zeros((init.Emod, init.totFlav*2, init.totFlav*2), dtype=np.complex128)
        # Add Oscillation Hamiltonian [km^-1]
        if init.hamOsc:
            totHam_km_1= totHam_km_1+ HAM_OSC_FLAV_KM_1
        # Add Matter Interaction Hamiltonian [km^-1]
        if init.hamMat:
            totHam_km_1= totHam_km_1+ self.hamMat4Flav_km_1()
        # Add ElectroMagnetic Interaction Hamiltonian [km^-1]
        if init.hamEM:
            totHam_km_1= totHam_km_1+ self.hamEM4Flav_km_1()
        # Add Self Interaction Hamiltonian [km^-1]
        if init.hamSelfSA:
            totHam_km_1= totHam_km_1- self.hamSelfSA4Flav_km_1()
        return totHam_km_1
# ============================

