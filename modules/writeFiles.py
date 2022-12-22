"""
ID      : writeFiles.py
Author  : Taygun Bulmus
E-Mail  : bulmust@gmail.com
-----
"""
import numpy as np
import hamilt as hamil
from shutil import copy2

# ============================
class tLog:
    INFO  = '\033[94m'+'[INFO]    '+'\033[0m'
    OK      = '\033[92m'+'[OK]        '+'\033[0m'
    WARNING = '\033[93m'+'[WARNING]   '+'\033[0m'
    ERROR   = '\033[91m'+'[ERROR]     '+'\033[0m'
    EXIT    = '\033[91m'+'[EXIT]      writeFiles.py'+'\033[0m'
# ============================

# ============================
# NPZ Files
# ============================
# rhoFlavAll.npz
def rhoFlavAll_npz(init, fullFilePath, rhoFlavAll):
    # Save file as npz format
    np.savez_compressed(fullFilePath, totFlav=init.totFlav, rhoFlavAll= rhoFlavAll\
        , distAll_km= init.distAll_km\
        , E_MeV=init.E_MeV, dim_rho_2totFlav_Bool=init.dim_rho_2totFlav_Bool)
# ============================
# parametersDic.npz
def parametersDic_npz(init, fullFilePath):
    # Save file as npz format
    np.savez_compressed(fullFilePath, **init.technicalParametersDic\
        , **init.physicalParametersDic, distAll_km= init.distAll_km\
        , E_MeV=init.E_MeV, TraceTerm=init.TraceTerm, TraceTermb=init.TraceTermb\
        , dim_rho_2totFlav_Bool=init.dim_rho_2totFlav_Bool)
# ============================

# ============================
# hamiltonians.npz
def hamiltonians_npz(init, fullFilePath, rhoFlavAll):
    # All possibilities
    p1=init.hamOsc; p2=init.hamMat; p3=init.hamEM; p4=init.hamSelfSA
    hamBoolAll_OscMatEMSelfSA = [p1, p2, p3, p4]

    # ============================
    # Zero Matrices
    if init.dim_rho_2totFlav_Bool:
        dim= init.totFlav* 2 
        hamTot_km_1 = np.zeros((len(init.distAll_km), init.Emod, dim, dim)\
            , dtype=np.complex128)
    else:
        dim=init.totFlav
        hamTot_km_1 = np.zeros((len(init.distAll_km), 2, init.Emod, dim, dim)\
            , dtype=np.complex128)
    if p2: hamiltMat   = np.zeros((len(init.distAll_km), dim, dim), dtype=np.complex128)
    if p3: hamiltEM    = np.zeros((len(init.distAll_km), dim, dim), dtype=np.complex128)
    if p4: hamiltSelfSA= np.zeros((len(init.distAll_km), dim, dim), dtype=np.complex128)
    # ============================

    # ============================
    # Hold all Hamiltonian's and total Hamiltonian
    for it1 in range(len(init.distAll_km)):
        # ============================
        if init.dim_rho_2totFlav_Bool:
            # ============================
            # Matter Hamiltonian Matrix
            if p2:
                if init.totFlav == 2:
                    # Two Flavor
                    hamiltMat[it1] = hamil.hamiltonian_BigRho(rhoFlavAll[it1]\
                        , init.distAll_km[it1]).hamMat2Flav_km_1()
                elif init.totFlav == 3:
                    # Three Flavor
                    hamiltMat[it1] = hamil.hamiltonian_BigRho(rhoFlavAll[it1]\
                        , init.distAll_km[it1]).hamMat3Flav_km_1()
                elif init.totFlav == 4:
                    # Four Flavor
                    hamiltMat[it1] = hamil.hamiltonian_BigRho(rhoFlavAll[it1]\
                        , init.distAll_km[it1]).hamMat4Flav_km_1()
            # ============================
            # EM Hamiltonian Matrix
            if p3:
                if init.totFlav == 2:
                    # Two Flavor
                    hamiltEM[it1] = hamil.hamiltonian_BigRho(rhoFlavAll[it1]\
                        , init.distAll_km[it1]).hamEM2Flav_km_1()
                elif init.totFlav == 3:
                    # Three Flavor
                    hamiltEM[it1] = hamil.hamiltonian_BigRho(rhoFlavAll[it1]\
                        , init.distAll_km[it1]).hamEM3Flav_km_1()
                elif init.totFlav == 4:
                    # Four Flavor
                    hamiltEM[it1] = hamil.hamiltonian_BigRho(rhoFlavAll[it1]\
                        , init.distAll_km[it1]).hamEM4Flav_km_1()
            # ============================
            # Self SA Hamiltonian Matrix
            if p4:
                if init.totFlav == 2 or init.totFlav == 3:
                    # Two or Three Flavor
                    hamiltSelfSA[it1] = hamil.hamiltonian_BigRho(rhoFlavAll[it1]\
                        , init.distAll_km[it1]).hamSelfSA23Flav_km_1()
                elif init.totFlav == 4:
                    # Four Flavor
                    hamiltSelfSA[it1] = hamil.hamiltonian_BigRho(rhoFlavAll[it1]\
                        , init.distAll_km[it1]).hamSelfSA4Flav_km_1()
            # ============================
            # Total Hamiltonian
            if init.totFlav == 2:
                # Two Flavor
                hamTot_km_1[it1]= hamil.hamiltonian_BigRho(rhoFlavAll[it1]\
                    , init.distAll_km[it1]).totalHam2Flav_km_1()
            elif init.totFlav == 3:
                # Three Flavor
                hamTot_km_1[it1]= hamil.hamiltonian_BigRho(rhoFlavAll[it1]\
                    , init.distAll_km[it1]).totalHam3Flav_km_1()
            elif init.totFlav == 4:
                # Four Flavor
                hamTot_km_1[it1]= hamil.hamiltonian_BigRho(rhoFlavAll[it1]\
                    , init.distAll_km[it1]).totalHam4Flav_km_1()
        else:
            # ============================
            # Matter Hamiltonian Matrix
            if p2:
                if init.totFlav==2:
                    # Two Flavor
                    hamiltMat[it1] = hamil.hamiltonian(rhoFlavAll[it1]\
                        , init.distAll_km[it1]).hamMat2Flav_km_1()
                elif init.totFlav==3:
                    # Three Flavor
                    hamiltMat[it1] = hamil.hamiltonian(rhoFlavAll[it1]\
                        , init.distAll_km[it1]).hamMat3Flav_km_1()
                elif init.totFlav==4:
                    # Four Flavor
                    hamiltMat[it1] = hamil.hamiltonian(rhoFlavAll[it1]\
                        , init.distAll_km[it1]).hamMat4Flav_km_1()
            # ============================
            # Self SA Hamiltonian Matrix
            if p4:
                if init.totFlav == 2 or init.totFlav == 3:
                    # Two or Three Flavor
                    hamiltSelfSA[it1] = hamil.hamiltonian(rhoFlavAll[it1]\
                        , init.distAll_km[it1]).hamSelfSA23Flav_km_1()
                elif init.totFlav == 4:
                    # Four Flavor
                    hamiltSelfSA[it1] = hamil.hamiltonian(rhoFlavAll[it1]\
                        , init.distAll_km[it1]).hamSelfSA4Flav_km_1()
            # ============================
            # Total Hamiltonian
            if init.totFlav == 2:
                # Two Flavor
                hamTot_km_1[it1]= hamil.hamiltonian(rhoFlavAll[it1]\
                    , init.distAll_km[it1]).totalHam2Flav_km_1()
            elif init.totFlav == 3:
                # Three Flavor
                hamTot_km_1[it1]= hamil.hamiltonian(rhoFlavAll[it1]\
                    , init.distAll_km[it1]).totalHam3Flav_km_1()
            elif init.totFlav == 4:
                # Four Flavor
                hamTot_km_1[it1]= hamil.hamiltonian(rhoFlavAll[it1]\
                    , init.distAll_km[it1]).totalHam4Flav_km_1()
    # ============================
    
    # ============================
    # Save file as npz format
    if p1 and p2 and p3 and p4:
        np.savez_compressed(fullFilePath, totFlav=init.totFlav\
            , hamBoolAll_OscMatEMSelfSA= hamBoolAll_OscMatEMSelfSA\
            ,hamOsc_km_1= hamil.HAM_OSC_FLAV_KM_1, hamMat_km_1= hamiltMat\
            , hamEM_km_1= hamiltEM, hamSelfSA_km_1= hamiltSelfSA\
            , hamTot_km_1= hamTot_km_1, distAll_km= init.distAll_km\
                , E_MeV= init.E_MeV, dim_rho_2totFlav_Bool=init.dim_rho_2totFlav_Bool)
    elif not p1 and p2 and p3 and p4:
        np.savez_compressed(fullFilePath, totFlav=init.totFlav\
            , hamBoolAll_OscMatEMSelfSA= hamBoolAll_OscMatEMSelfSA\
            , hamMat_km_1= hamiltMat\
            , hamEM_km_1= hamiltEM, hamSelfSA_km_1= hamiltSelfSA\
            , hamTot_km_1= hamTot_km_1, distAll_km= init.distAll_km, E_MeV= init.E_MeV\
                , dim_rho_2totFlav_Bool=init.dim_rho_2totFlav_Bool)
    elif p1 and not p2 and p3 and p4:
        np.savez_compressed(fullFilePath, totFlav=init.totFlav\
            , hamBoolAll_OscMatEMSelfSA= hamBoolAll_OscMatEMSelfSA\
            , hamOsc_km_1= hamil.HAM_OSC_FLAV_KM_1\
            , hamEM_km_1= hamiltEM, hamSelfSA_km_1= hamiltSelfSA\
            , hamTot_km_1= hamTot_km_1, distAll_km= init.distAll_km, E_MeV= init.E_MeV\
            , dim_rho_2totFlav_Bool=init.dim_rho_2totFlav_Bool)
    elif p1 and p2 and not p3 and p4:
        np.savez_compressed(fullFilePath, totFlav=init.totFlav\
            , hamBoolAll_OscMatEMSelfSA= hamBoolAll_OscMatEMSelfSA\
            , hamOsc_km_1= hamil.HAM_OSC_FLAV_KM_1, hamMat_km_1= hamiltMat\
            , hamSelfSA_km_1= hamiltSelfSA\
            , hamTot_km_1= hamTot_km_1, distAll_km= init.distAll_km, E_MeV= init.E_MeV\
            , dim_rho_2totFlav_Bool=init.dim_rho_2totFlav_Bool)
    elif p1 and p2 and p3 and not p4:
        np.savez_compressed(fullFilePath, totFlav=init.totFlav\
            , hamBoolAll_OscMatEMSelfSA= hamBoolAll_OscMatEMSelfSA\
            , hamOsc_km_1= hamil.HAM_OSC_FLAV_KM_1, hamMat_km_1= hamiltMat\
            , hamEM_km_1= hamiltEM\
            , hamTot_km_1= hamTot_km_1, distAll_km= init.distAll_km, E_MeV= init.E_MeV\
            , dim_rho_2totFlav_Bool=init.dim_rho_2totFlav_Bool)
    elif not p1 and not p2 and p3 and p4:
        np.savez_compressed(fullFilePath, totFlav=init.totFlav\
            , hamBoolAll_OscMatEMSelfSA= hamBoolAll_OscMatEMSelfSA\
            , hamEM_km_1= hamiltEM, hamSelfSA_km_1= hamiltSelfSA\
            , hamTot_km_1= hamTot_km_1, distAll_km= init.distAll_km, E_MeV= init.E_MeV\
            , dim_rho_2totFlav_Bool=init.dim_rho_2totFlav_Bool)
    elif not p1 and p2 and not p3 and p4:
        np.savez_compressed(fullFilePath, totFlav=init.totFlav\
            , hamBoolAll_OscMatEMSelfSA= hamBoolAll_OscMatEMSelfSA\
            , hamMat_km_1= hamiltMat\
            , hamSelfSA_km_1= hamiltSelfSA\
            , hamTot_km_1= hamTot_km_1, distAll_km= init.distAll_km, E_MeV= init.E_MeV\
            , dim_rho_2totFlav_Bool=init.dim_rho_2totFlav_Bool)
    elif not p1 and p2 and p3 and not p4:
        np.savez_compressed(fullFilePath, totFlav=init.totFlav\
            , hamBoolAll_OscMatEMSelfSA= hamBoolAll_OscMatEMSelfSA\
            , hamMat_km_1= hamiltMat\
            , hamEM_km_1= hamiltEM\
            , hamTot_km_1= hamTot_km_1, distAll_km= init.distAll_km, E_MeV= init.E_MeV\
            , dim_rho_2totFlav_Bool=init.dim_rho_2totFlav_Bool)
    elif p1 and not p2 and not p3 and p4:
        np.savez_compressed(fullFilePath, totFlav=init.totFlav\
            , hamBoolAll_OscMatEMSelfSA= hamBoolAll_OscMatEMSelfSA\
            , hamOsc_km_1= hamil.HAM_OSC_FLAV_KM_1\
            , hamSelfSA_km_1= hamiltSelfSA\
            , hamTot_km_1= hamTot_km_1, distAll_km= init.distAll_km, E_MeV= init.E_MeV\
            , dim_rho_2totFlav_Bool=init.dim_rho_2totFlav_Bool)
    elif p1 and not p2 and p3 and not p4:
        np.savez_compressed(fullFilePath, totFlav=init.totFlav\
            , hamBoolAll_OscMatEMSelfSA= hamBoolAll_OscMatEMSelfSA\
            , hamOsc_km_1= hamil.HAM_OSC_FLAV_KM_1\
            , hamEM_km_1= hamiltEM\
            , hamTot_km_1= hamTot_km_1, distAll_km= init.distAll_km, E_MeV= init.E_MeV\
            , dim_rho_2totFlav_Bool=init.dim_rho_2totFlav_Bool)
    elif p1 and p2 and not p3 and not p4:
        np.savez_compressed(fullFilePath, totFlav=init.totFlav\
            , hamBoolAll_OscMatEMSelfSA= hamBoolAll_OscMatEMSelfSA\
            , hamOsc_km_1= hamil.HAM_OSC_FLAV_KM_1, hamMat_km_1= hamiltMat\
            , hamTot_km_1= hamTot_km_1, distAll_km= init.distAll_km\
            , E_MeV= init.E_MeV, dim_rho_2totFlav_Bool=init.dim_rho_2totFlav_Bool)
    elif not p1 and not p2 and not p3 and p4:
        np.savez_compressed(fullFilePath, totFlav=init.totFlav\
            , hamBoolAll_OscMatEMSelfSA= hamBoolAll_OscMatEMSelfSA\
            , hamSelfSA_km_1= hamiltSelfSA\
            , hamTot_km_1= hamTot_km_1, distAll_km= init.distAll_km\
            , E_MeV= init.E_MeV, dim_rho_2totFlav_Bool=init.dim_rho_2totFlav_Bool)
    elif not p1 and not p2 and p3 and not p4:
        np.savez_compressed(fullFilePath, totFlav=init.totFlav\
            , hamBoolAll_OscMatEMSelfSA= hamBoolAll_OscMatEMSelfSA\
            , hamEM_km_1= hamiltEM\
            , hamTot_km_1= hamTot_km_1, distAll_km= init.distAll_km\
            , E_MeV= init.E_MeV, dim_rho_2totFlav_Bool=init.dim_rho_2totFlav_Bool)
    elif not p1 and p2 and not p3 and not p4:
        np.savez_compressed(fullFilePath, totFlav=init.totFlav\
            , hamBoolAll_OscMatEMSelfSA= hamBoolAll_OscMatEMSelfSA\
            , hamMat_km_1= hamiltMat\
            , hamTot_km_1= hamTot_km_1, distAll_km= init.distAll_km\
            , E_MeV= init.E_MeV, dim_rho_2totFlav_Bool=init.dim_rho_2totFlav_Bool)
    elif p1 and not p2 and not p3 and not p4:
        np.savez_compressed(fullFilePath, totFlav=init.totFlav\
            , hamBoolAll_OscMatEMSelfSA= hamBoolAll_OscMatEMSelfSA\
            , hamOsc_km_1= hamil.HAM_OSC_FLAV_KM_1\
            , hamTot_km_1= hamTot_km_1, distAll_km= init.distAll_km\
            , E_MeV= init.E_MeV, dim_rho_2totFlav_Bool=init.dim_rho_2totFlav_Bool)
    else:
        print(tLog.EXIT)
        exit('')
# ============================

# ============================
# eigVal_eigVec.npz
def eigVal_eigVec_npz(init, fullFilePath, rhoFlavAll):
    # Zero Matrices
    if init.dim_rho_2totFlav_Bool:
        hamTot_km_1 = np.zeros((len(init.distAll_km)\
            , init.Emod, init.totFlav*2, init.totFlav*2), dtype=np.complex128)
    else:
        hamTot_km_1 = np.zeros((len(init.distAll_km)\
            , 2, init.Emod, init.totFlav, init.totFlav), dtype=np.complex128)
    
    # ============================
    # Total Hamiltonian
    for it1 in range(len(init.distAll_km)):
        if init.dim_rho_2totFlav_Bool:
            if init.totFlav == 2:
                # Two Flavor
                hamTot_km_1[it1]= hamil.hamiltonian_BigRho(rhoFlavAll[it1]\
                    , init.distAll_km[it1]).totalHam2Flav_km_1()
            elif init.totFlav == 3:
                # Three Flavor
                hamTot_km_1[it1]= hamil.hamiltonian_BigRho(rhoFlavAll[it1]\
                    , init.distAll_km[it1]).totalHam3Flav_km_1()
            elif init.totFlav == 4:
                # Four Flavor
                hamTot_km_1[it1]= hamil.hamiltonian_BigRho(rhoFlavAll[it1]\
                    , init.distAll_km[it1]).totalHam4Flav_km_1()
        else:
            if init.totFlav == 2:
                # Two Flavor
                hamTot_km_1[it1]= hamil.hamiltonian(rhoFlavAll[it1]\
                    , init.distAll_km[it1]).totalHam2Flav_km_1()
            elif init.totFlav == 3:
                # Three Flavor
                hamTot_km_1[it1]= hamil.hamiltonian(rhoFlavAll[it1]\
                    , init.distAll_km[it1]).totalHam3Flav_km_1()
            elif init.totFlav == 4:
                # Four Flavor
                hamTot_km_1[it1]= hamil.hamiltonian(rhoFlavAll[it1]\
                    , init.distAll_km[it1]).totalHam4Flav_km_1()

    # Eigenvalues and Eigenvectors
    eigVal_km_1, eigVec = np.linalg.eigh(hamTot_km_1)
    # Save npz file
    np.savez_compressed(fullFilePath, totFlav=init.totFlav\
        , eigVal_Ascending_km_1=eigVal_km_1, eigVec_Ascending= eigVec\
        , distAll_km= init.distAll_km, E_MeV= init.E_MeV\
        , dim_rho_2totFlav_Bool= init.dim_rho_2totFlav_Bool)
# ============================

# ============================
# TXT Files
# ============================
# RhoFlavAll.txt
def rhoFlavAll_txt(init,fullFilePath, rhoFlavAll):
    # Save file as txt format
    np.savetxt(fullFilePath[:-4]+'distAll_km.txt', init.distAll_km, fmt='%-7.5f')
    np.savetxt(fullFilePath[:-4]\
        +'_realPart_flatten'+str(np.shape(rhoFlavAll))\
        +'.txt', np.real(rhoFlavAll.flatten()), fmt='%-7.5f')
    np.savetxt(fullFilePath[:-4]\
        +'_imagPart_flatten'+str(np.shape(rhoFlavAll))\
        +'.txt', np.imag(rhoFlavAll.flatten()), fmt='%-7.5f')
# ============================
# hamiltonians.txt
# TODO: Implement hamiltonians_txt function
def hamiltonians_txt(init, fullFilePath, rhoFlavAll):
    print(tLog.ERROR+'hamiltonians.txt file not implemented yet.')
# ============================
# eigVal_eigVec.txt
# TODO: Implement eigVal_eigVec_txt function
def eigVal_eigVec_txt(init, fullFilePath, rhoFlavAll):
    print(tLog.ERROR+'eigVal_eigVec.txt file not implemented yet.')
# ============================

# ============================
# Save Background Profile
def saveBackgroundProfiles(init, backgroundFolderPath, saveFolderPath):
    # Do you use default background profile?
    if not init.physicalParametersDic['use_defaultMatterProfile']:
        matterProfile_fileName=\
            str(init.physicalParametersDic['matterProfile_fileName'][:-1])
        # Save file
        copy2(backgroundFolderPath+matterProfile_fileName, saveFolderPath)
        print(tLog.INFO+'Matter background profile file is copied to: '+saveFolderPath)
    if not init.physicalParametersDic['use_defaultMagneticProfile']:
        magneticProfile_fileName=\
            str(init.physicalParametersDic['magneticProfile_fileName'][:-1])
        # Save file
        copy2(backgroundFolderPath+magneticProfile_fileName, saveFolderPath)
        print(tLog.INFO+'EM background profile file is copied to: '+saveFolderPath)
# ============================