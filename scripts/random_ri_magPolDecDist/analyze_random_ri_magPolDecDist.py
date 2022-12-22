import matplotlib.pyplot as plt
import numpy as np
import os as osCommand # File operations
from pathlib import Path

# Current Directory
currDir = osCommand.getcwd() + "/"
# collNuPy main directory
COLLECTIVE_NU_OSC_DIR= str(Path(str(Path(currDir).parent)+ "/").parent) + '/'
# ============================

# ============================
# Constants
G_f_MeV_2      = 1.1663786999999999e-11 # [MeV^-2]
km_1__MeV      = 1.973269804593025e-16
MeV__km_1      = 5067730716156395.0
MeV_1_TO_km    = 1.973269804593025e-16
erg_s_TO_MeV_km = 2.0819433270935597
# ============================

# ============================
# Colors For Printing
class tLog:
    INFO  = '\033[94m'+'[INFO]    '+'\033[0m'
    OK      = '\033[92m'+'[OK]        '+'\033[0m'
    WARNING = '\033[93m'+'[WARNING]   '+'\033[0m'
    ERROR   = '\033[91m'+'[ERROR]     '+'\033[0m'
    EXIT    = '\033[91m'+'[EXIT]      analyze_random_ri_magPolDecDist.py' + '\033[0m'
# ============================

# ============================
# MSW Resonance Variables
def mswResonanceVariablesCreate(hierarchy, E_MeV, distAll_km\
    , deltaMsqr_MeV2, theta, Vcc_km_1, Vnc_km_1, hamTot_km_1):
    # ============================
    # Shortcuts
    #Delta_MeV= deltaMsqr_MeV2/ (4* E_MeV)
    #Vcc__MeV= Vcc_km_1* km_1__MeV
    #Vnc_MeV= Vnc_km_1* km_1__MeV
    # ============================

    # ============================
    # M and M bar [MeV] for All Distance
    #M_MeV    = np.sqrt( ((Delta_MeV* np.sin(2* theta))**2)\
    #    + ((Delta_MeV* np.cos(2* theta) - Vcc__MeV/2)**2) )
    #Mbar_MeV = np.sqrt( ((Delta_MeV* np.sin(2* theta))**2)\
    #    + ((Delta_MeV* np.cos(2* theta) + Vcc__MeV/2)**2) )
    # Eigenvalues of Matter Basis (omega's) for All Distance [MeV]
    #eigValMat_MeV = np.array([\
    #      (- M_MeV+ Vnc_MeV +(Vcc__MeV/ 2) )\
    #    , (+ M_MeV+ Vnc_MeV +(Vcc__MeV/ 2) )\
    #    , (+ Mbar_MeV- Vnc_MeV -(Vcc__MeV/ 2) )\
    #    , (- Mbar_MeV- Vnc_MeV -(Vcc__MeV/ 2) )])
    # ThetaM, ThetaBar_M, gamma and Theta_Tilde for All Distance
    #theta_M   = np.pi/2 + np.arctan2(Delta_MeV* np.sin(2* theta)\
    #    , (Delta_MeV* np.cos(2* theta)- Vcc__MeV/2))/2
    #thetabar_M=           np.arctan2(Delta_MeV* np.sin(2* theta)\
    #    , (Delta_MeV* np.cos(2* theta)+ Vcc__MeV/2))/2
    #gamma_AllDist= theta_M - thetabar_M
    # ============================

    # ============================
    # Find locations of MSW resonance
    if hierarchy == 1:
        # For Normal Hierarchy
        tmp= np.array(np.where((((np.roll(np.sign(hamTot_km_1[:,0,0,0]\
            - hamTot_km_1[:,0,1,1]), 1) - np.sign(hamTot_km_1[:,0,0,0]\
            - hamTot_km_1[:,0,1,1])) != 0)) == 1)).flatten()
    else:
        # For Inverted Hierarchy
        tmp= np.array(np.where((((np.roll(np.sign(hamTot_km_1[:,0,2,2]\
            - hamTot_km_1[:,0,3,3]), 1) - np.sign(hamTot_km_1[:,0,2,2]\
            - hamTot_km_1[:,0,3,3])) != 0)) == 1)).flatten()
    # if tmp is empty (NO MSW resonance)
    if tmp.size == 0:
        print(tLog.ERROR + "No MSW resonance found!")
        return None,None,None,None
    else:
        # Remove the element which is zero
        if tmp[0] == 0: 
            tmp= tmp[1:]
        else:
            exit(tLog.EXIT + "The first element of tmp is not zero. It is "\
                + str(tmp[0]))
        # Create Zero Arrays
        r_resMSWs_km_idx= np.zeros(len(tmp), dtype=int)
        r_resMSWs_km=np.zeros(len(tmp))
        adiabaticity_MSW= np.zeros(len(tmp))
        P_LZ_MSW= np.zeros(len(tmp)); itTmp=0
        # Assign Resonance Index and Distance
        for it1 in range(len(tmp)):
            if tmp[it1] == 0:
                continue
            r_resMSWs_km_idx[itTmp]= tmp[it1]
            r_resMSWs_km[itTmp]= distAll_km[tmp[it1]]
            itTmp= itTmp+1
        # If neutrino oscillation is not considered
        if theta == 0: 
            r_resMSWs_km_idx=0; r_resMSWs_km=0
        # ============================
        
        # ============================
        # Adiabaticity (Bulmus:2022gyz, eqn (48))
        Vcc_gradient_km_2= np.abs(np.gradient(Vcc_km_1, distAll_km[1]-distAll_km[0]))
        Vcc_gradient_MeV2= Vcc_gradient_km_2* (km_1__MeV**2)
        for it1 in range(len(adiabaticity_MSW)):
            adiabaticity_MSW[it1]=\
                (((deltaMsqr_MeV2/ (2* E_MeV))* np.sin(2* theta))**2)\
                / Vcc_gradient_MeV2[r_resMSWs_km_idx[it1]]
        # ============================

        # ============================
        # Landau Zener Transition Probability
        P_LZ_MSW= np.exp((-2)* np.pi* adiabaticity_MSW)
        # ============================

        # ============================
        # Printing
        print(tLog.INFO + "MSW resonance is at " + str(r_resMSWs_km) + " km(s).")
        print(tLog.INFO + "Adiabaticity is(are) " + str(adiabaticity_MSW) + ".")
        print(tLog.INFO + "Landau Zener Transition Probability is(are) "\
            + str(P_LZ_MSW) + ".")
        # ============================

        return r_resMSWs_km_idx, r_resMSWs_km, adiabaticity_MSW, P_LZ_MSW
# ============================

# ============================
# SFP Resonance Variables
def sfpResonanceVariables(hierarchy, E_MeV, distAll_km\
    , deltaMsqr_MeV2, Vcc_km_1, Vnc_km_1, muB_km_1, hamTot_km_1):
    # ============================
    # Shortcuts
    #Delta_MeV= deltaMsqr_MeV2/ (4* E_MeV)
    #Vcc__MeV= Vcc_km_1* km_1__MeV
    #Vnc_MeV= Vnc_km_1* km_1__MeV
    # ============================

    # ============================
    # Find locationS of SFP resonance
    if hierarchy == 1:
        # For Normal Hierarchy
        tmp= np.array(np.where((((np.roll(np.sign(hamTot_km_1[:,0,1,1]\
            - hamTot_km_1[:,0,2,2]), 1) - np.sign(hamTot_km_1[:,0,1,1]\
            - hamTot_km_1[:,0,2,2])) != 0)) == 1)).flatten()
    else:
        # For Inverted Hierarchy
        tmp= np.array(np.where((((np.roll(np.sign(hamTot_km_1[:,0,0,0]\
            - hamTot_km_1[:,0,3,3]), 1) - np.sign(hamTot_km_1[:,0,0,0]\
            - hamTot_km_1[:,0,3,3])) != 0)) == 1)).flatten()
    # if tmp is empty (NO SFP resonance)
    if tmp.size == 0:
        print(tLog.ERROR + "No SFP resonance found!")
        return None,None,None,None
    else:
        # Remove the element which is zero
        if tmp[0] == 0: 
            tmp= tmp[1:]
        else:
            exit(tLog.EXIT + "The first element of tmp is not zero. It is "\
                + str(tmp[0]))
        for it1 in range(len(tmp)):
            if tmp[it1] == 0:
                continue
            if distAll_km[tmp[it1]] > 5e+4:
                continue
            else:
                r_resSFP_km_idx= tmp[it1]
                r_resSFP_km= distAll_km[r_resSFP_km_idx]
                break
        # ============================

        # ============================
        # Adiabaticity (Bulmus:2022gyz, eqn (38))
        # We need (1-2*Ye) in the denominator
        # (-2* Vnc - Vcc) = A*(1-Ye) - A*Ye = A* (1- 2*Ye)
        denominator_km_2= np.real(np.abs(np.gradient(((-2* Vnc_km_1)- Vcc_km_1)\
            , distAll_km[1]-distAll_km[0])))
        adiabaticity_SFP = (np.real(muB_km_1[r_resSFP_km_idx])**2)\
            / denominator_km_2[r_resSFP_km_idx]
        # ============================

        # ============================
        # Landau Zener Transition Probability
        P_LZ_SFP= np.exp((-2)* np.pi* adiabaticity_SFP)
        # ============================

        # ============================
        # Printing
        print(tLog.INFO + "SFP resonance is at " + str(r_resSFP_km) + " km.")
        print(tLog.INFO + "SFP adiabaticity is " + str(adiabaticity_SFP))
        print(tLog.INFO + "SFP Landau Zener Transition Probability is " + str(P_LZ_SFP))
        # ============================
        
        return r_resSFP_km_idx, r_resSFP_km, adiabaticity_SFP, P_LZ_SFP
# ============================

# ============================
def analyzeAllResults():
    RESULTS_DIR = COLLECTIVE_NU_OSC_DIR+'results/'
    # https://stackoverflow.com/questions/29769181/count-the-number-of-folders-in-a-directory-and-subdirectories
    noOfSim=len(next(osCommand.walk(RESULTS_DIR))[1])
    print(tLog.INFO+'Number of simulations is '+str(noOfSim)+'.')
    # ============================

    # ============================
    E_MeV= np.zeros(noOfSim); TraceTerm= np.zeros(noOfSim)
    distINIT_km= np.zeros(noOfSim)
    # ============================

    # ============================
    # Similar for all simulations
    # First Data Load
    DATA_DIR= RESULTS_DIR+'simulation1/data/'
    DATA_rhoFlavAll= np.load(DATA_DIR+'rhoFlavAll.npz', allow_pickle=True)
    #DATA_hamiltonians= np.load(DATA_DIR+'hamiltonians.npz', allow_pickle=True)
    DATA_parametersDic= np.load(DATA_DIR+'parametersDic.npz', allow_pickle=True)
    #print(DATA_parametersDic.files); exit()
    # ============================
    # Variables
    totFlav                   = DATA_rhoFlavAll['totFlav']
    dim_rho_2totFlav_Bool     = DATA_rhoFlavAll['dim_rho_2totFlav_Bool']
    #hamBoolAll_OscMatEMSelfSA = DATA_hamiltonians['hamBoolAll_OscMatEMSelfSA']
    use_defaultMatterProfile  = DATA_parametersDic['use_defaultMatterProfile']
    #deltaMsqr_MeV2            = DATA_parametersDic['deltaM_square21']
    #theta                     = DATA_parametersDic['theta12']
    #hierarchy                 = DATA_parametersDic['neutrino_hierarchy']
    neutrino_MagneticMoment   = DATA_parametersDic['neutrino_MagneticMoment']
    print(tLog.INFO+'Neutrino Magnetic Moment is '+str(neutrino_MagneticMoment)+'.')
    #if not use_defaultMatterProfile:
    #    matterProfile_fileName= str(DATA_parametersDic['matterProfile_fileName'])[:-1]
    #else:
    #    #TODO Create function for default matter profile
    #    pass
    #DATA_matterProfile= np.loadtxt(DATA_DIR+ matterProfile_fileName)
    #matterProfile_distAll_km= DATA_matterProfile[:,0]* 1e-5
    #matterProfile_baryonDens_g_cm3= DATA_matterProfile[:,2]
    #matterProfile_Ye= DATA_matterProfile[:,1]
    # If all elements of matterProfile_Ye are same, then took only one of them
    #if np.all(matterProfile_Ye == matterProfile_Ye[0]):
    #    matterProfile_Ye= matterProfile_Ye[0]
    print(tLog.INFO+'============================\n')
    # ============================

    # Create Zero Arrays
    rhoFlavAll_finalAveragedDiag= np.zeros((noOfSim, 4)); avIndex= 50
    rhoFlavAll_initialDiag= np.zeros((noOfSim, 4))
    rhoFlavAll_initialMat= np.zeros((noOfSim,4,4))
    energySpec= np.zeros((2,noOfSim, 4))
    # ============================
    # Combine all the results
    for it_sim in range(1, noOfSim+1):
        # ============================
        DATA_DIR          = RESULTS_DIR+'simulation'+str(it_sim)+'/data/'
        DATA_rhoFlavAll   = np.load(DATA_DIR+'rhoFlavAll.npz', allow_pickle=True)
        #DATA_hamiltonians = np.load(DATA_DIR+'hamiltonians.npz', allow_pickle=True)
        if use_defaultMatterProfile:
            exit(tLog.ERROR+'use_defaultMatterProfile is True.')
        # ============================
        
        # ============================
        # Load Variables
        E_MeV[it_sim-1]= DATA_rhoFlavAll['E_MeV']
        print(tLog.INFO+ 'Energy='+str(E_MeV[it_sim-1]) +' MeV is starting.')
        rhoFlavAll = DATA_rhoFlavAll['rhoFlavAll']
        distAll_km = DATA_rhoFlavAll['distAll_km']
        #hamOsc_km_1= DATA_hamiltonians['hamOsc_km_1']
        #hamMat_km_1= DATA_hamiltonians['hamMat_km_1']
        #hamEM_km_1 = DATA_hamiltonians['hamEM_km_1']
        #hamTot_km_1= DATA_hamiltonians['hamTot_km_1']
        #Vcc_km_1= np.real(hamMat_km_1[:,0,0]- hamMat_km_1[:,1,1])
        #Vnc_km_1= np.real(hamMat_km_1[:,1,1])
        #muB_km_1= hamEM_km_1[:,0,-1]
        # ============================

        # ============================
        #r_resMSWs_km_idx, r_resMSWs_km, adiabaticity_MSW, P_LZ_MSW_tmpArray=\
        #    mswResonanceVariablesCreate(hierarchy, E_MeV[it_sim-1], distAll_km\
        #    , deltaMsqr_MeV2, theta, Vcc_km_1, Vnc_km_1, hamTot_km_1)
        #print(tLog.INFO+ '============================')
        #r_resSFP_km_idx, r_resSFP_km, adiabaticity_SFP, P_LZ_SFP=\
        #    sfpResonanceVariables(hierarchy, E_MeV[it_sim-1], distAll_km\
        #    , deltaMsqr_MeV2, Vcc_km_1, Vnc_km_1, muB_km_1, hamTot_km_1)
        #print(tLog.INFO+ '============================\n')
        # ============================

        # ============================
        # Last Elements of rhoFlavAll
        #rhoFlavAll_finalDiag[it_sim-1]= np.real(np.diag(rhoFlavAll[-1,0]))
        rhoFlavAll_finalAveragedDiag[it_sim-1,0]=\
            np.real(rhoFlavAll[-avIndex:,0,0,0].mean())
        rhoFlavAll_finalAveragedDiag[it_sim-1,1]=\
            np.real(rhoFlavAll[-avIndex:,0,1,1].mean())
        rhoFlavAll_finalAveragedDiag[it_sim-1,2]=\
            np.real(rhoFlavAll[-avIndex:,0,2,2].mean())
        rhoFlavAll_finalAveragedDiag[it_sim-1,3]=\
            np.real(rhoFlavAll[-avIndex:,0,3,3].mean())

        rhoFlavAll_initialDiag[it_sim-1,0]= np.real(rhoFlavAll[0,0,0,0])
        rhoFlavAll_initialDiag[it_sim-1,1]= np.real(rhoFlavAll[0,0,1,1])
        rhoFlavAll_initialDiag[it_sim-1,2]= np.real(rhoFlavAll[0,0,2,2])
        rhoFlavAll_initialDiag[it_sim-1,3]= np.real(rhoFlavAll[0,0,3,3])
        
        rhoFlavAll_initialMat[it_sim-1]= np.real(rhoFlavAll[0,0])
        # ============================

        # ============================
        # Test Averaged Results For Last
        #flavType= 0
        #plt.plot(distAll_km, np.real(rhoFlavAll[:,0,flavType,flavType]), label='rho11')
        #plt.plot(distAll_km, rhoFlavAll_finalAveragedDiag[it_sim-1,flavType]\
        #    * np.ones(len(distAll_km)), label='rho11_Final_Averaged')
        #plt.legend()
        #plt.show()
        #plt.close('all')
        # ============================
        
        # ============================
        # If Fermi Dirac Distribution is used
        # Multiply with trace terms
        '''
        # Calculate Constant Term
        if totFlav == 2:
            if dim_rho_2totFlav_Bool: # TraceTerm = TraceTermb
                if DATA_parametersDic['neutrino_distributionParameter'] == 1:
                    luminosityArray_MeV_km= np.zeros(4); tempArray_MeV= np.zeros(4)
                    luminosityArray_MeV_km[0]= DATA_parametersDic['luminosity_e']\
                        * erg_s_TO_MeV_km
                    luminosityArray_MeV_km[1]= DATA_parametersDic['luminosity_mu']\
                        * erg_s_TO_MeV_km
                    luminosityArray_MeV_km[2]= DATA_parametersDic['luminosity_eb']\
                        * erg_s_TO_MeV_km
                    luminosityArray_MeV_km[3]= DATA_parametersDic['luminosity_mub']\
                        * erg_s_TO_MeV_km
                    tempArray_MeV[0]= DATA_parametersDic['temperature_e']
                    tempArray_MeV[1]= DATA_parametersDic['temperature_mu']
                    tempArray_MeV[2]= DATA_parametersDic['temperature_eb']
                    tempArray_MeV[3]= DATA_parametersDic['temperature_mub']
                    for i3 in range(totFlav*2):
                        constant= (km_1__MeV**2)*(DATA_parametersDic['TraceTerm'][0]\
                            * (1 / (4 * (distAll_km[0]**2) * np.pi**2))\
                            * luminosityArray_MeV_km[i3] / (3.1514 * tempArray_MeV[i3]))
                        rhoFlavAll_initialDiag[it_sim-1,i3]=\
                            rhoFlavAll_initialDiag[it_sim-1,i3]/ constant
                        rhoFlavAll_finalAveragedDiag[it_sim-1,i3]=\
                            rhoFlavAll_finalAveragedDiag[it_sim-1,i3]/ constant
                        print(constant)
                        input('')
        '''
        TraceTerm[it_sim-1]= np.real(DATA_parametersDic['TraceTerm'][0])

        energySpec[0,it_sim-1,0]= np.real(rhoFlavAll[0,0,0,0])#* TraceTerm[it_sim-1]
        energySpec[0,it_sim-1,1]= np.real(rhoFlavAll[0,0,1,1])#* TraceTerm[it_sim-1]
        energySpec[0,it_sim-1,2]= np.real(rhoFlavAll[0,0,2,2])#* TraceTerm[it_sim-1]
        energySpec[0,it_sim-1,3]= np.real(rhoFlavAll[0,0,3,3])#* TraceTerm[it_sim-1]

        energySpec[1,it_sim-1,0]= np.real(rhoFlavAll[-1,0,0,0])#* TraceTerm[it_sim-1]
        energySpec[1,it_sim-1,1]= np.real(rhoFlavAll[-1,0,1,1])#* TraceTerm[it_sim-1]
        energySpec[1,it_sim-1,2]= np.real(rhoFlavAll[-1,0,2,2])#* TraceTerm[it_sim-1]
        energySpec[1,it_sim-1,3]= np.real(rhoFlavAll[-1,0,3,3])#* TraceTerm[it_sim-1]

        distINIT_km[it_sim-1]= distAll_km[0]
    if totFlav == 2:
        if dim_rho_2totFlav_Bool: # TraceTerm = TraceTermb
            if DATA_parametersDic['neutrino_distributionParameter'] == 1:
                luminosityArray_MeV_km= np.zeros(4); tempArray_MeV= np.zeros(4)
                luminosityArray_MeV_km[0]= DATA_parametersDic['luminosity_e']\
                    * erg_s_TO_MeV_km
                luminosityArray_MeV_km[1]= DATA_parametersDic['luminosity_mu']\
                    * erg_s_TO_MeV_km
                luminosityArray_MeV_km[2]= DATA_parametersDic['luminosity_eb']\
                    * erg_s_TO_MeV_km
                luminosityArray_MeV_km[3]= DATA_parametersDic['luminosity_mub']\
                    * erg_s_TO_MeV_km
                tempArray_MeV[0]= DATA_parametersDic['temperature_e']
                tempArray_MeV[1]= DATA_parametersDic['temperature_mu']
                tempArray_MeV[2]= DATA_parametersDic['temperature_eb']
                tempArray_MeV[3]= DATA_parametersDic['temperature_mub']

                distributionInit_nu_flav= np.zeros((len(E_MeV),totFlav*2, totFlav*2))
                rho= np.zeros((len(E_MeV),totFlav*2, totFlav*2))

                for i2 in range(len(E_MeV)):
                    for i3 in range(totFlav*2):
                        distributionInit_nu_flav[i2, i3, i3]=\
                            (1/ (1.803 * (tempArray_MeV[i3]**3)))* (E_MeV[i2]**2)\
                            / (1 + np.exp(E_MeV[i2] / tempArray_MeV[i3]))
                
                        rho[i2, i3, i3]  = (1/ (4 * (distINIT_km[i2]**2) * np.pi**2))\
                            * distributionInit_nu_flav[i2,i3, i3]\
                            * luminosityArray_MeV_km[i3] / (3.1514 * tempArray_MeV[i3])
                        # TraceTerm => np.trace(rho[i2])
                    energySpec[0,i2]= energySpec[0,i2]* np.trace(rho[i2])
                    energySpec[1,i2]= energySpec[1,i2]* np.trace(rho[i2])
                '''
                plt.plot(E_MeV, energySpec[0,:,0], label= 'e', linewidth= 3)
                plt.plot(E_MeV, rho[:,0,0], label= 'eRho')
                plt.plot(E_MeV, energySpec[0,:,1], label= 'mu', linewidth= 3)
                plt.plot(E_MeV, rho[:,1,1], label= 'muRho')
                plt.plot(E_MeV, energySpec[0,:,2], label= 'eb', linewidth= 3)
                plt.plot(E_MeV, rho[:,2,2], label= 'ebRho')
                plt.plot(E_MeV, energySpec[0,:,3], label= 'mub', linewidth= 3)
                plt.plot(E_MeV, rho[:,3,3], label= 'mubRho')
                plt.legend();plt.show();exit('')
                '''
                
                hold_E= np.zeros(int(noOfSim/3))
                hold_energySpec= np.zeros((int(noOfSim/3), 4, 2))
                for it in range(int(noOfSim/3)):
                    hold_E[it]= E_MeV[it*3]
                    for it2 in range(4):
                        hold_energySpec[it,it2]=\
                            ([np.min(energySpec[1,it*3:(it+1)*3,it2])\
                            , np.max(energySpec[1,it*3:(it+1)*3,it2])])
                   
                plt.plot(E_MeV, energySpec[0,:,0], color='r', linestyle='dashed'\
                    , label=r'$\rho_{i,ee}$')
                plt.plot(E_MeV, energySpec[0,:,1], color='k', linestyle='dashed'\
                    , label=r'$\rho_{i,xx}$')
                plt.plot(E_MeV, energySpec[0,:,2], color='b', linestyle='dashed'\
                    , label=r'$\rho_{i,ebeb}$')
                plt.plot(E_MeV, energySpec[0,:,3], color='g', linestyle='dashed'\
                    , label=r'$\rho_{i,xbxb}$')

                plt.scatter(E_MeV, energySpec[1,:,0], color='r', s=1\
                    , label=r'$\rho_{ee}$')
                plt.scatter(E_MeV, energySpec[1,:,1], color='k', s=1\
                    , label=r'$\rho_{xx}$')
                plt.scatter(E_MeV, energySpec[1,:,2], color='b', s=1\
                    , label=r'$\rho_{ebeb}$')
                plt.scatter(E_MeV, energySpec[1,:,3], color='g', s=1\
                    , label=r'$\rho_{xbxb}$')

                #plt.fill_between(\
                #    hold_E, hold_energySpec[:,0,0], hold_energySpec[:,0,1]\
                #    , color='r', alpha=0.6)
                #plt.fill_between(\
                #    hold_E, hold_energySpec[:,1,0], hold_energySpec[:,1,1]\
                #    , color='k', alpha=0.6)
                #plt.fill_between(\
                #    hold_E, hold_energySpec[:,2,0], hold_energySpec[:,2,1]\
                #    , color='b', alpha=0.6)
                #plt.fill_between(\
                #    hold_E, hold_energySpec[:,3,0], hold_energySpec[:,3,1]\
                #    , color='g', alpha=0.6)
                plt.legend()        
                plt.show(); plt.close('all'); exit()

    #plt.plot(E_MeV, energySpec[0,:,0], label='e');plt.show();exit('')
    # If Fermi Dirac Distribution is used
    if DATA_parametersDic['neutrino_distributionParameter'] == 1:
        plt.plot(E_MeV, rhoFlavAll_initialDiag[:,0]* TraceTerm, color='r'\
            , linestyle='dashed', label=r'$\rho_{i,ee}$')
        plt.plot(E_MeV, rhoFlavAll_initialDiag[:,1]* TraceTerm, color='k'\
            , linestyle='dashed', label=r'$\rho_{i,xx}$')
        plt.plot(E_MeV, rhoFlavAll_initialDiag[:,2]* TraceTerm, color='b'\
            , linestyle='dashed', label=r'$\rho_{i,ebeb}$')
        plt.plot(E_MeV, rhoFlavAll_initialDiag[:,3]* TraceTerm, color='g'\
            , linestyle='dashed', label=r'$\rho_{i,xbxb}$')

        plt.scatter(E_MeV, rhoFlavAll_finalAveragedDiag[:,0]* TraceTerm, color='r'\
            , s=0.3, label=r'$\rho_{ee}$')
        plt.scatter(E_MeV, rhoFlavAll_finalAveragedDiag[:,1]* TraceTerm, color='k'\
            , s=0.3, label=r'$\rho_{xx}$')
        plt.scatter(E_MeV, rhoFlavAll_finalAveragedDiag[:,2]* TraceTerm, color='b'\
            , s=0.3, label=r'$\rho_{ebeb}$')
        plt.scatter(E_MeV, rhoFlavAll_finalAveragedDiag[:,3]* TraceTerm, color='g'\
            , s=0.3, label=r'$\rho_{xbxb}$')
        plt.legend()
        plt.show()
        plt.close()
        exit()

    # Scattering Graph
    plt.scatter(E_MeV, rhoFlavAll_finalAveragedDiag[:,0], color='k', s=0.3\
        , label=r'$\rho_{ee}$')
    plt.legend()
    plt.show()
    plt.close()

    plt.scatter(E_MeV, rhoFlavAll_finalAveragedDiag[:,1], color='k', s=0.3\
        , label=r'$\rho_{xx}$')
    plt.legend()
    plt.show()
    plt.close()

    plt.scatter(E_MeV, rhoFlavAll_finalAveragedDiag[:,2], color='k', s=0.3\
        , label=r'$\rho_{ebeb}$')
    plt.legend()
    plt.show()
    plt.close()

    plt.scatter(E_MeV, rhoFlavAll_finalAveragedDiag[:,3], color='k', s=0.3\
        , label=r'$\rho_{xbxb}$')
    plt.legend()
    plt.show()
    plt.close()

if __name__ == '__main__':
    analyzeAllResults()
