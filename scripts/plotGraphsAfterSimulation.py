"""
--------------------------------------------------------------------------------
Copyright (c) 2019, Taygun Bulmus
All rights reserved.
See the LICENSE file for license information.
--------------------------------------------------------------------------------
-----
ID      : plotGraphsAfterSimulation.py
Author  : Taygun Bulmus
E-Mail  : bulmust@gmail.com
-----
"""
# ============================
import os as osCommand # File operations
from sys import path
# Insert plotgraphs.py path
import plotgraphs as plot
from pathlib import Path
# Current Directory
currDir = osCommand.getcwd() + "/"
# collNuPy main directory
COLLECTIVE_NU_OSC_DIR= str(Path(currDir).parent)+ "/"
# Insert Modules Path
path.insert(0, COLLECTIVE_NU_OSC_DIR+'modules')

# ============================
# Colors For Printing
class tLog:
    INFO  = '\033[94m'+'[INFO]    '+'\033[0m'
    OK      = '\033[92m'+'[OK]        '+'\033[0m'
    WARNING = '\033[93m'+'[WARNING]   '+'\033[0m'
    ERROR   = '\033[91m'+'[ERROR]     '+'\033[0m'
    EXIT    = '\033[91m' + '[EXIT]      plotGraphsAfterSimulation.py' + '\033[0m'
# ============================

if __name__ == "__main__":
    RESULTS_SIMULATION_DIR= COLLECTIVE_NU_OSC_DIR+ 'results/'\
        + input('Enter the simulation name: ')+'/'
    NPZ_FILE_PATH= RESULTS_SIMULATION_DIR+ 'data/rhoFlavAll.npz'

    # Try to reach the file
    print(tLog.INFO+'Try to reach the file:', NPZ_FILE_PATH)
    if Path(NPZ_FILE_PATH).is_file():
        print(tLog.OK+'File Found.')
    else:
        print(tLog.ERROR+'File Not Found.')
        exit(tLog.EXIT)
    print('')
    print(tLog.INFO+'Choose plot(s): Available plots are following')
    print(tLog.INFO+'1 = distDiag => Diagonal elements of density matrix to distance.')
    print(tLog.INFO+'2 = energyDiag => Diagonal elements of density matrix to energy.')
    print(tLog.INFO+\
        '3 = distHamiltDiag => Diagonal elements of Hamiltonian to distance.')
    print('')
    print(tLog.INFO+'Example: To plot distDiag and distHamiltDiag, enter 13')
    plotCheck= input('Enter the plot numbers: ')
    print(tLog.INFO+'The plotting is started.')
    if '1' in plotCheck:
        print(tLog.INFO+'Plotting distDiag.')
        plot.plotMain(RESULTS_SIMULATION_DIR).distDiag()
        print(tLog.OK+'Plotting distDiag is done.')
    if '2' in plotCheck:
        print(tLog.INFO+'Plotting energyDiag.')
        plot.plotMain(RESULTS_SIMULATION_DIR).energyDiag()
        print(tLog.OK+'Plotting energyDiag is done.')
    if '3' in plotCheck:
        print(tLog.INFO+'Plotting distHamiltDiag.')
        plot.plotMain(RESULTS_SIMULATION_DIR).distHamiltDiag()
        print(tLog.OK+'Plotting distHamiltDiag is done.')
    print(tLog.OK+'The plotting is done.')

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