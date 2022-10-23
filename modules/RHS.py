"""
-----
ID      : RHS.py
Author  : Taygun Bulmus
E-Mail  : bulmust@gmail.com
-----
"""
import numpy as np
import initial as init
import datetime
import hamilt as hamil
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
# RHS For Neutrinos f(r,rho) !!Only LSODA!!
# Two Flavor
def rhs2Flav_noInterData(rhoAll, dist):
    # ============================
    # Define the RHS
    totHamAP = hamil.hamiltonian(rhoAll, dist).totalHam2Flav_km_1()
    # rhoAll must be in the form of [2, energyMod, totFlav*2, totFlav*2]

    # Total Hamiltionian (totHam) [2, energyMod, totFlav, totFlav]
    # For neutrinos     -> totHam[0]
    # For antineutrinos -> totHam[1]

    # Note: 3d array + 2d array gives a result that add all 2d array with each element of 3d array
    # Take Commutator of density and Hamiltonian for each energy
    return np.array([(-1j)* ((totHamAP[0] @ rhoAll[0])- (rhoAll[0] @ totHamAP[0]))\
        , (-1j)*((totHamAP[1] @ rhoAll[1])- (rhoAll[1] @ totHamAP[1]))]).flatten()
def rhs2Flav(rhoAll, dist):
    # ============================
    # Show intermediate distance and hold the values
    if dist >= init.distTemp:
        # Save to npz file
        np.savez_compressed('intermediateData/intermediateRhoFlav'+ str(int(np.round(dist, decimals=0)))+ 'km.npz'\
            , dist_km= dist, rhoAll_flav= rhoAll, dim_rho_2totFlav_Bool= init.dim_rho_2totFlav_Bool\
            , technicalParametersDic = init.technicalParametersDic\
            , physicalParametersDic=init.physicalParametersDic)
        
        # Hold new distance
        init.distTemp = init.distTemp+ init.technicalParametersDic['holdData_Every']

        # To show current time
        currentTimeDate = datetime.datetime.now()
        fileObject = open('distanceAndTime.log','a+')
        # Do not color the [INFO]
        print('Current distance is %.7f' % dist, 'km and time :', currentTimeDate.strftime("%Y-%m-%d %H:%M:%S"), file=fileObject)
        fileObject.close()
        print(tLog.INFO+'Current distance is %.7f' % dist, 'km and time :', currentTimeDate.strftime("%Y-%m-%d %H:%M:%S"))
    # ============================
    
    # ============================
    # Define the RHS
    totHamAP = hamil.hamiltonian(rhoAll, dist).totalHam2Flav_km_1()
    # rhoAll must be in the form of [2, energyMod, totFlav*2, totFlav*2]

    # Total Hamiltionian (totHam) [2, energyMod, totFlav, totFlav]
    # For neutrinos     -> totHam[0]
    # For antineutrinos -> totHam[1]

    # Note: 3d array + 2d array gives a result that add all 2d array with each element of 3d array
    # Take Commutator of density and Hamiltonian for each energy
    return np.array([(-1j)* ((totHamAP[0] @ rhoAll[0])- (rhoAll[0] @ totHamAP[0]))\
        , (-1j)*((totHamAP[1] @ rhoAll[1])- (rhoAll[1] @ totHamAP[1]))]).flatten()
# Three Flavor
def rhs3Flav_noInterData(rhoAll, dist):
    # ============================
    # Define the RHS
    totHamAP = hamil.hamiltonian(rhoAll, dist).totalHam3Flav_km_1()
    # rhoAll must be in the form of [2, energyMod, totFlav*2, totFlav*2]

    # Total Hamiltionian (totHam) [2, energyMod, totFlav, totFlav]
    # For neutrinos     -> totHam[0]
    # For antineutrinos -> totHam[1]

    # Note: 3d array + 2d array gives a result that add all 2d array with each element of 3d array
    # Take Commutator of density and Hamiltonian for each energy
    return np.array([(-1j)* ((totHamAP[0] @ rhoAll[0])- (rhoAll[0] @ totHamAP[0]))\
        , (-1j)*((totHamAP[1] @ rhoAll[1])- (rhoAll[1] @ totHamAP[1]))]).flatten()
def rhs3Flav(rhoAll, dist):
    # ============================
    # Show intermediate distance and hold the values
    if dist >= init.distTemp:
        # Save to npz file
        np.savez_compressed('intermediateData/intermediateRhoFlav'+ str(int(np.round(dist, decimals=0)))+ 'km.npz'\
            , dist_km= dist, rhoAll_flav= rhoAll, dim_rho_2totFlav_Bool= init.dim_rho_2totFlav_Bool\
            , technicalParametersDic = init.technicalParametersDic\
            , physicalParametersDic=init.physicalParametersDic)
        
        # Hold new distance
        init.distTemp = init.distTemp+ init.technicalParametersDic['holdData_Every']

        # To show current time
        currentTimeDate = datetime.datetime.now()
        fileObject = open('distanceAndTime.log','a+')
        # Do not color the [INFO]
        print('Current distance is %.7f' % dist, 'km and time :', currentTimeDate.strftime("%Y-%m-%d %H:%M:%S"), file=fileObject)
        fileObject.close()
        print(tLog.INFO+'Current distance is %.7f' % dist, 'km and time :', currentTimeDate.strftime("%Y-%m-%d %H:%M:%S"))
    # ============================
    
    # ============================
    # Define the RHS
    totHamAP = hamil.hamiltonian(rhoAll, dist).totalHam3Flav_km_1()
    # rhoAll must be in the form of [2, energyMod, totFlav*2, totFlav*2]

    # Total Hamiltionian (totHam) [2, energyMod, totFlav, totFlav]
    # For neutrinos     -> totHam[0]
    # For antineutrinos -> totHam[1]

    # Note: 3d array + 2d array gives a result that add all 2d array with each element of 3d array
    # Take Commutator of density and Hamiltonian for each energy
    return np.array([(-1j)* ((totHamAP[0] @ rhoAll[0])- (rhoAll[0] @ totHamAP[0]))\
        , (-1j)*((totHamAP[1] @ rhoAll[1])- (rhoAll[1] @ totHamAP[1]))]).flatten()
# Four Flavor
def rhs4Flav_noInterData(rhoAll, dist):
    # ============================
    # Define the RHS
    totHamAP = hamil.hamiltonian(rhoAll, dist).totalHam4Flav_km_1()
    # rhoAll must be in the form of [2, energyMod, totFlav*2, totFlav*2]

    # Total Hamiltionian (totHam) [2, energyMod, totFlav, totFlav]
    # For neutrinos     -> totHam[0]
    # For antineutrinos -> totHam[1]

    # Note: 3d array + 2d array gives a result that add all 2d array with each element of 3d array
    # Take Commutator of density and Hamiltonian for each energy
    return np.array([(-1j)* ((totHamAP[0] @ rhoAll[0])- (rhoAll[0] @ totHamAP[0]))\
        , (-1j)*((totHamAP[1] @ rhoAll[1])- (rhoAll[1] @ totHamAP[1]))]).flatten()
def rhs4Flav(rhoAll, dist):
    # ============================
    # Show intermediate distance and hold the values
    if dist >= init.distTemp:
        # Save to npz file
        np.savez_compressed('intermediateData/intermediateRhoFlav'+ str(int(np.round(dist, decimals=0)))+ 'km.npz'\
            , dist_km= dist, rhoAll_flav= rhoAll, dim_rho_2totFlav_Bool= init.dim_rho_2totFlav_Bool\
            , technicalParametersDic = init.technicalParametersDic\
            , physicalParametersDic=init.physicalParametersDic)
        
        # Hold new distance
        init.distTemp = init.distTemp+ init.technicalParametersDic['holdData_Every']

        # To show current time
        currentTimeDate = datetime.datetime.now()
        fileObject = open('distanceAndTime.log','a+')
        # Do not color the [INFO]
        print('Current distance is %.7f' % dist, 'km and time :', currentTimeDate.strftime("%Y-%m-%d %H:%M:%S"), file=fileObject)
        fileObject.close()
        print(tLog.INFO+'Current distance is %.7f' % dist, 'km and time :', currentTimeDate.strftime("%Y-%m-%d %H:%M:%S"))
    # ============================
    
    # ============================
    # Define the RHS
    totHamAP = hamil.hamiltonian(rhoAll, dist).totalHam4Flav_km_1()
    # rhoAll must be in the form of [2, energyMod, totFlav*2, totFlav*2]

    # Total Hamiltionian (totHam) [2, energyMod, totFlav, totFlav]
    # For neutrinos     -> totHam[0]
    # For antineutrinos -> totHam[1]

    # Note: 3d array + 2d array gives a result that add all 2d array with each element of 3d array
    # Take Commutator of density and Hamiltonian for each energy
    return np.array([(-1j)* ((totHamAP[0] @ rhoAll[0])- (rhoAll[0] @ totHamAP[0]))\
        , (-1j)*((totHamAP[1] @ rhoAll[1])- (rhoAll[1] @ totHamAP[1]))]).flatten()
# ============================

# ============================
# RHS For Neutrinos f(r,rho) !!Only LSODA!!
# Two Flavor
def rhs2Flav_bigRho_noInterData(rhoAll, dist):
    # ============================
    # Define the RHS
    totHamAP = hamil.hamiltonian_BigRho(rhoAll, dist).totalHam2Flav_km_1()
    # rhoAll must be in the form of [energyMod, totFlav*2, totFlav*2]
    # Total Hamiltionian (totHam) [energyMod,totFlav*2,totFlav*2]
    # totHam[N,0:totFlav,0:totFlav] => Neutrinos
    # totHam[N,totFlav:totFlav*2,totFlav:totFlav*2] => Anti-Neutrinos
    # totHam[N,0:totFlav,totFlav:totFlav*2] => Neutrino to Anti-neutrino
    
    # Note: 3d array + 2d array gives a result that add all 2d array with each element of 3d array
    # Take Commutator of density and Hamiltonian for each energy

    return ((-1j)* ((totHamAP @ rhoAll)- (rhoAll @ totHamAP))).flatten()
def rhs2Flav_bigRho(rhoAll, dist):
    # ============================
    # Show intermediate distance and hold the values
    if dist >= init.distTemp:
        # Save to npz file
        np.savez_compressed('intermediateData/intermediateRhoFlav'+ str(int(np.round(dist, decimals=0)))+ 'km.npz'\
            , dist_km= dist, rhoAll_flav= rhoAll, dim_rho_2totFlav_Bool= init.dim_rho_2totFlav_Bool\
            , technicalParametersDic = init.technicalParametersDic\
            , physicalParametersDic=init.physicalParametersDic)
        
        # Hold new distance
        init.distTemp = init.distTemp+ init.technicalParametersDic['holdData_Every']

        # To show current time
        currentTimeDate = datetime.datetime.now()
        fileObject = open('distanceAndTime.log','a+')
        # Do not color the [INFO]
        print('Current distance is %.7f' % dist, 'km and time :', currentTimeDate.strftime("%Y-%m-%d %H:%M:%S"), file=fileObject)
        fileObject.close()
        print(tLog.INFO+'Current distance is %.7f' % dist, 'km and time :', currentTimeDate.strftime("%Y-%m-%d %H:%M:%S"))
    # ============================
    
    # ============================
    # Define the RHS
    totHamAP = hamil.hamiltonian_BigRho(rhoAll, dist).totalHam2Flav_km_1()
    # rhoAll must be in the form of [energyMod, totFlav*2, totFlav*2]
    # Total Hamiltionian (totHam) [energyMod,totFlav*2,totFlav*2]
    # totHam[N,0:totFlav,0:totFlav] => Neutrinos
    # totHam[N,totFlav:totFlav*2,totFlav:totFlav*2] => Anti-Neutrinos
    # totHam[N,0:totFlav,totFlav:totFlav*2] => Neutrino to Anti-neutrino
    
    # Note: 3d array + 2d array gives a result that add all 2d array with each element of 3d array
    # Take Commutator of density and Hamiltonian for each energy

    return ((-1j)* ((totHamAP @ rhoAll)- (rhoAll @ totHamAP))).flatten()
# Three Flavor
def rhs3Flav_bigRho_noInterData(rhoAll, dist):
    # ============================
    # Define the RHS
    totHamAP = hamil.hamiltonian_BigRho(rhoAll, dist).totalHam3Flav_km_1()
    # rhoAll must be in the form of [energyMod, totFlav*2, totFlav*2]
    # Total Hamiltionian (totHam) [energyMod,totFlav*2,totFlav*2]
    # totHam[N,0:totFlav,0:totFlav] => Neutrinos
    # totHam[N,totFlav:totFlav*2,totFlav:totFlav*2] => Anti-Neutrinos
    # totHam[N,0:totFlav,totFlav:totFlav*2] => Neutrino to Anti-neutrino
    
    # Note: 3d array + 2d array gives a result that add all 2d array with each element of 3d array
    # Take Commutator of density and Hamiltonian for each energy

    return ((-1j)* ((totHamAP @ rhoAll)- (rhoAll @ totHamAP))).flatten()
def rhs3Flav_bigRho(rhoAll, dist):
    # ============================
    # Show intermediate distance and hold the values
    if dist >= init.distTemp:
        # Save to npz file
        np.savez_compressed('intermediateData/intermediateRhoFlav'+ str(int(np.round(dist, decimals=0)))+ 'km.npz'\
            , dist_km= dist, rhoAll_flav= rhoAll, dim_rho_2totFlav_Bool= init.dim_rho_2totFlav_Bool\
            , technicalParametersDic = init.technicalParametersDic\
            , physicalParametersDic=init.physicalParametersDic)
        
        # Hold new distance
        init.distTemp = init.distTemp+ init.technicalParametersDic['holdData_Every']

        # To show current time
        currentTimeDate = datetime.datetime.now()
        fileObject = open('distanceAndTime.log','a+')
        # Do not color the [INFO]
        print('Current distance is %.7f' % dist, 'km and time :', currentTimeDate.strftime("%Y-%m-%d %H:%M:%S"), file=fileObject)
        fileObject.close()
        print(tLog.INFO+'Current distance is %.7f' % dist, 'km and time :', currentTimeDate.strftime("%Y-%m-%d %H:%M:%S"))
    # ============================
    
    # ============================
    # Define the RHS
    totHamAP = hamil.hamiltonian_BigRho(rhoAll, dist).totalHam3Flav_km_1()
    # rhoAll must be in the form of [energyMod, totFlav*2, totFlav*2]
    # Total Hamiltionian (totHam) [energyMod,totFlav*2,totFlav*2]
    # totHam[N,0:totFlav,0:totFlav] => Neutrinos
    # totHam[N,totFlav:totFlav*2,totFlav:totFlav*2] => Anti-Neutrinos
    # totHam[N,0:totFlav,totFlav:totFlav*2] => Neutrino to Anti-neutrino
    
    # Note: 3d array + 2d array gives a result that add all 2d array with each element of 3d array
    # Take Commutator of density and Hamiltonian for each energy

    return ((-1j)* ((totHamAP @ rhoAll)- (rhoAll @ totHamAP))).flatten()
# Four Flavor
def rhs4Flav_bigRho_noInterData(rhoAll, dist):
    # ============================
    # Define the RHS
    totHamAP = hamil.hamiltonian_BigRho(rhoAll, dist).totalHam4Flav_km_1()
    # rhoAll must be in the form of [energyMod, totFlav*2, totFlav*2]
    # Total Hamiltionian (totHam) [energyMod,totFlav*2,totFlav*2]
    # totHam[N,0:totFlav,0:totFlav] => Neutrinos
    # totHam[N,totFlav:totFlav*2,totFlav:totFlav*2] => Anti-Neutrinos
    # totHam[N,0:totFlav,totFlav:totFlav*2] => Neutrino to Anti-neutrino
    
    # Note: 3d array + 2d array gives a result that add all 2d array with each element of 3d array
    # Take Commutator of density and Hamiltonian for each energy

    return ((-1j)* ((totHamAP @ rhoAll)- (rhoAll @ totHamAP))).flatten()
def rhs4Flav_bigRho(rhoAll, dist):
    # ============================
    # Show intermediate distance and hold the values
    if dist >= init.distTemp:
        # Save to npz file
        np.savez_compressed('intermediateData/intermediateRhoFlav'+ str(int(np.round(dist, decimals=0)))+ 'km.npz'\
            , dist_km= dist, rhoAll_flav= rhoAll, dim_rho_2totFlav_Bool= init.dim_rho_2totFlav_Bool\
            , technicalParametersDic = init.technicalParametersDic\
            , physicalParametersDic=init.physicalParametersDic)
        
        # Hold new distance
        init.distTemp = init.distTemp+ init.technicalParametersDic['holdData_Every']

        # To show current time
        currentTimeDate = datetime.datetime.now()
        fileObject = open('distanceAndTime.log','a+')
        # Do not color the [INFO]
        print('Current distance is %.7f' % dist, 'km and time :', currentTimeDate.strftime("%Y-%m-%d %H:%M:%S"), file=fileObject)
        fileObject.close()
        print(tLog.INFO+'Current distance is %.7f' % dist, 'km and time :', currentTimeDate.strftime("%Y-%m-%d %H:%M:%S"))
    # ============================
    
    # ============================
    # Define the RHS
    totHamAP = hamil.hamiltonian_BigRho(rhoAll, dist).totalHam4Flav_km_1()
    # rhoAll must be in the form of [energyMod, totFlav*2, totFlav*2]
    # Total Hamiltionian (totHam) [energyMod,totFlav*2,totFlav*2]
    # totHam[N,0:totFlav,0:totFlav] => Neutrinos
    # totHam[N,totFlav:totFlav*2,totFlav:totFlav*2] => Anti-Neutrinos
    # totHam[N,0:totFlav,totFlav:totFlav*2] => Neutrino to Anti-neutrino
    
    # Note: 3d array + 2d array gives a result that add all 2d array with each element of 3d array
    # Take Commutator of density and Hamiltonian for each energy

    return ((-1j)* ((totHamAP @ rhoAll)- (rhoAll @ totHamAP))).flatten()
# ============================