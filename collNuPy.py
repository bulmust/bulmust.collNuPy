"""
--------------------------------------------------------------------------------
Copyright (c) 2019, Taygun Bulmus
All rights reserved.
See the LICENSE file for license information.
--------------------------------------------------------------------------------
-----
ID      : collNuPy.py
Author  : Taygun Bulmus
E-Mail  : bulmust@gmail.com
-----
"""
# ============================
import os as osCommand # File operations
from sys import path # Insert modules paths
# ============================
# Insert modules path
if __name__ == "__main__":
    # Collective neutrino oscillation simulation directory location
    COLLECTIVE_NU_OSC_DIR = osCommand.getcwd() + "/"
    # Insert modules path
    path.insert(0, COLLECTIVE_NU_OSC_DIR+"modules")
# ============================
# User defined modules
# Differential Eqn Solvers
from diffSolvers import odeintwSolver
# Initial Conditions
import initial as init
import writeFiles as write2file
import plotgraphs as plot
# ============================

# ============================
# Colors For Printing
class tLog:
    STATUS  = '\033[94m'+'[STATUS]    '+'\033[0m'
    OK      = '\033[92m'+'[OK]        '+'\033[0m'
    WARNING = '\033[93m'+'[WARNING]   '+'\033[0m'
    ERROR   = '\033[91m'+'[ERROR]     '+'\033[0m'
# ============================

# ============================
# Start of the program
if __name__ == "__main__":
    # TODO Add option to create user defined name. If it is not given use simulation#
    
    # ============================
    # Make Results Directory
    RESULTS_DIR = COLLECTIVE_NU_OSC_DIR+'results/'
    if not osCommand.path.exists(RESULTS_DIR):
        osCommand.makedirs(RESULTS_DIR)
        print(tLog.OK+"Results Folder Created.")
    else:
        print(tLog.STATUS+"Results Folder Already Exists.")
    # ============================

    # ============================
    # Make sim# directory.
    # It tries to make a folder that starts from simulation_1 to simulation10000
    for simNum in range(1, 10000): # 10000 is big enough
        # Make simulation# Folder
        RESULTS_SIMULATION_DIR = RESULTS_DIR+'simulation'+str(simNum)+'/'
        # Make Simulation_# Folder
        if not osCommand.path.exists(RESULTS_SIMULATION_DIR):
            osCommand.makedirs(RESULTS_SIMULATION_DIR)
            print(tLog.OK+"simulation"+str(simNum)+" Directory Created.")
            break # Exit from for loop
    # ============================

    # ============================
    # Differential equation solvers
    if init.technicalParametersDic['initialValueProblem_solverMethod'] == 'LSODA':
        # odeintw uses LSODA method with complex domain
        rhoFlavAll = odeintwSolver(init, RESULTS_SIMULATION_DIR)
    # ============================
    
    # ============================
    # Save final npz data
    DATA_FOLDER_NPZ= RESULTS_SIMULATION_DIR+ 'data/'
    # rhoFlavAll.npz
    write2file.rhoFlavAll_npz(init, DATA_FOLDER_NPZ+'rhoFlavAll.npz', rhoFlavAll)
    print(tLog.OK+"Saved: rhoFlavAll.npz.")
    write2file.parametersDic_npz(init, DATA_FOLDER_NPZ+'parametersDic.npz')
    print(tLog.OK+"Saved: parametersDic.npz.")
    # ============================
    if init.technicalParametersDic['output_distanceHamiltonianAllE']:
        # Hamiltonians hamiltonians.npz
        write2file.hamiltonians_npz(init, DATA_FOLDER_NPZ+'hamiltonians.npz', rhoFlavAll)
        print(tLog.OK+"Saved: hamiltonians.npz.")
    if init.technicalParametersDic['output_distance_eigenValuesAllE']:
        # Eigenvalues eigenvalues.npz
        write2file.eigVal_eigVec_npz(init, DATA_FOLDER_NPZ+'eigVal_eigVec.npz', rhoFlavAll)
        print(tLog.OK+"Saved: eigVal_eigVec.npz.")
    # ============================

    # ============================
    # Save final human readable data
    # TODO: Comment out print after wrote the functions
    if init.technicalParametersDic['output_humanReadable']:
        DATA_FOLDER_TXT= RESULTS_SIMULATION_DIR+ 'data/txt/'
        if not osCommand.path.exists(DATA_FOLDER_TXT):
            osCommand.makedirs(DATA_FOLDER_TXT)
            print(tLog.OK+"TXT Data Folder Created.")
        else:
            print(tLog.STATUS+"TXT Data Folder Already Exists.")
        # ============================
        # rhoFlavAll.txt
        write2file.rhoFlavAll_txt(init,DATA_FOLDER_TXT+'rhoFlavAll.txt', rhoFlavAll)
        print(tLog.OK+"Saved: rhoFlavAll.txt.")
        # ============================
        # hamiltonians.txt
        if init.technicalParametersDic['output_distanceHamiltonianAllE']:
            # Hamiltonians hamiltonians.txt
            write2file.hamiltonians_txt(init, DATA_FOLDER_TXT+'hamiltonians.txt', rhoFlavAll)
            #print(tLog.OK+"Saved: hamiltonians.txt.")
        # ============================
        # eigenvalues.txt
        if init.technicalParametersDic['output_distance_eigenValuesAllE']:
            # Eigenvalues eigenvalues.txt
            write2file.eigVal_eigVec_txt(init, DATA_FOLDER_TXT+'eigVal_eigVec.txt', rhoFlavAll)
            #print(tLog.OK+"Saved: eigVal_eigVec.txt.")
    # ============================

    # ============================
    # Save final plots
    if init.technicalParametersDic['plotGraphs']:
        if init.technicalParametersDic['plotGraphs_diagonalRhoAllEnergy_2distance']:
            plot.plotMain(RESULTS_SIMULATION_DIR).distDiag()
        if init.technicalParametersDic['plotGraphs_diagonalRhoFinal_2energy']:
            plot.plotMain(RESULTS_SIMULATION_DIR).energyDiag()
        if init.technicalParametersDic['plotGraphs_hamiltonianAllEnergy_2distance']:
            plot.plotMain(RESULTS_SIMULATION_DIR).distHamiltDiag()
    # ============================

    #! ============================
    #! WORKS ONLY FOR UNIX SYSTEMS
    # Clear Unused Files and Folders
    osCommand.chdir(COLLECTIVE_NU_OSC_DIR)
    # Find *.pyc files and remove them
    osCommand.system('find . -name "*.pyc" -exec rm -f {} \;')
    # Find all __pycache__ folders and remove them
    osCommand.system('find . | grep -E "(__pycache__|.pyc|.pyo$)" | xargs rm -rf')
    print(tLog.OK+"Unused Files and Folders Removed.")
    #! ============================
    print(tLog.STATUS+'[FINISHED] =====> simulation'+str(simNum)+' <===== [FINISHED]\n')
    exit()
# ============================