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
import numpy as np
import random
import os as osCommand # File operations
from pathlib import Path
# Current Directory
currDir = osCommand.getcwd() + "/"
# collNuPy main directory
COLLECTIVE_NU_OSC_DIR= str(Path(currDir).parent)+ "/"

# ============================
# Colors For Printing
class tLog:
    INFO  = '\033[94m'+'[INFO]    '+'\033[0m'
    OK      = '\033[92m'+'[OK]        '+'\033[0m'
    WARNING = '\033[93m'+'[WARNING]   '+'\033[0m'
    ERROR   = '\033[91m'+'[ERROR]     '+'\033[0m'
    EXIT    = '\033[91m' + '[EXIT]      run_random_ri_magPolDecDist.py' + '\033[0m'
# ============================

if __name__ =='__main__':
    # Energies
    Ei_MeV= 1.0; Ef_MeV= 2.0; Emod=2
    energy= np.zeros(Emod); energy[0]= Ei_MeV

    # Initial Distance (Different For Every Energy)
    ri_km= np.zeros(Emod)

    # magneticField_polynomialDecayDistance (Different For Every Energy)
    magneticField_polynomialDecayDistance= np.zeros(Emod)

    # Starting Random Points
    ri_randStart_km= 49.95
    magneticField_polynomialDecayDistance_randStart= 49.95

    # Assign First Values
    ri_km[0]= ri_randStart_km+ 0.1* random.random()
    magneticField_polynomialDecayDistance[0]= 49.95+ 0.1* random.random()

    # Create Random ri_km and magneticField_polynomialDecayDistance
    for i1 in range(1, int(Emod)):
        # Create Energy
        energy[i1]= float(Ei_MeV) + ((float(Ef_MeV)- float(Ei_MeV))/ (Emod- 1))* i1
        
        # Assign Random ri_km and magneticField_polynomialDecayDistance
        ri_km[i1]= ri_randStart_km+ 0.1* random.random()
        magneticField_polynomialDecayDistance[i1]= 49.95+ 0.1* random.random()

    # Start Doing The Simulation
    for it_do in range(Emod):
        # ============================
        # Modify New ri_km and magneticField_polynomialDecayDistance 
        # in physicalParameters.dat
        file= open(COLLECTIVE_NU_OSC_DIR+ 'parameters/physicalParameters.dat', 'r')
        fileWrite= open(COLLECTIVE_NU_OSC_DIR+'parameters/physicalParameters2.dat', 'w')
        for line in file.readlines():
            # Energy
            if (line.startswith('numberOf_energyMode')):
                line= 'numberOf_energyMode=1  # Valid for only integer\n'
                fileWrite.write(line)
            elif (line.startswith('energy_initial')):
                line= 'energy_initial='+ str(energy[it_do])+ '  # [MeV]\n'
                fileWrite.write(line)
            elif (line.startswith('energy_final')):
                line= 'energy_final='+ str(energy[it_do])+ '  # [MeV]\n'
                fileWrite.write(line)
            # ri_km
            elif (line.startswith('distance_initial')):
                line= 'distance_initial='+ str(ri_km[it_do])+ '  # [km]\n'
                fileWrite.write(line)
            # magneticField_polynomialDecayDistance
            elif (line.startswith('magneticField_polynomialDecayDistance')):
                line= 'magneticField_polynomialDecayDistance='\
                    + str(magneticField_polynomialDecayDistance[it_do])\
                    +' # [km] Valid for (use_defaultMagneticProfile=1'\
                    ', magneticField_profile=2)\n'
                fileWrite.write(line)
            else:
                fileWrite.write(line)
        file.close()
        fileWrite.close()
        # Remove and Rename
        osCommand.remove(COLLECTIVE_NU_OSC_DIR+ 'parameters/physicalParameters.dat')
        osCommand.rename(COLLECTIVE_NU_OSC_DIR+ 'parameters/physicalParameters2.dat'\
            , COLLECTIVE_NU_OSC_DIR+ 'parameters/physicalParameters.dat')
        # ============================
        osCommand.chdir(COLLECTIVE_NU_OSC_DIR)
        #! ============================
        #! WORKS ONLY FOR UNIX SYSTEMS
        osCommand.system('python3 '+COLLECTIVE_NU_OSC_DIR+'collNuPy.py')
        #! ============================






