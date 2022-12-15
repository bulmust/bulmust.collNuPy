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
import argparse

# Current Directory
currDir = osCommand.getcwd() + "/"
from pathlib import Path
# collNuPy main directory
COLLECTIVE_NU_OSC_DIR= str(Path(str(Path(currDir).parent)+ "/").parent) + '/'

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
    # ============================
     # Initialize parser
    parser = argparse.ArgumentParser()
    # Adding optional argument
    parser.add_argument("-i", "--Ei_MeV", help= "Initial Energy", type= float)
    parser.add_argument("-f", "--Ef_MeV", help= "Final Energy", type= float)
    parser.add_argument("-m", "--NumberOfMod", help= "Number of Mod", type= int)
    parser.add_argument("-b", "--NumberOfMod_EachE", help= "Number of Mod For Each Energy", type= int)
    # Read arguments from command line
    args = parser.parse_args()
    # Check for --Ei_MeV
    if args.Ei_MeV:
        Ei_MeV= float(args.Ei_MeV)
        print(tLog.INFO+'Initial energy is ', Ei_MeV, ' MeV.')
    else:
        print(tLog.INFO+'Initial energy is not given. Using default value of 1.0 MeV.')
        Ei_MeV= 1.0    
    # Check for --Ef_MeV
    if args.Ef_MeV:
        Ef_MeV= float(args.Ef_MeV)
        print(tLog.INFO+'Final energy is ', Ef_MeV, ' MeV.')
    else:
        print(tLog.INFO+'Final energy is not given. Using default value of 50.0 MeV.')
        Ef_MeV= 1.0
    # Check for --NumberOfMod
    if args.NumberOfMod:
        Emod= int(args.NumberOfMod)
        print(tLog.INFO+'Number of energy mod is ', Emod, '.')
    else:
        print(tLog.INFO+'Number of energy mod is not given. Using default value of 500.')
        Emod= 500
    if args.NumberOfMod_EachE:
        Emod_EachE= int(args.NumberOfMod_EachE)
        print(tLog.INFO+'Number of energy mod for each energy is ', Emod_EachE, '.')
    else:
        print(tLog.INFO+'Number of energy mod for each energy is not given. Using default value of 1.')
        Emod= 1
    # ============================

    # Energies
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
    magneticField_polynomialDecayDistance[0]= magneticField_polynomialDecayDistance_randStart+ 0.1* random.random()

    # Create Random ri_km and magneticField_polynomialDecayDistance
    for i1 in range(1, int(Emod)):
        # Create Energy
        energy[i1]= float(Ei_MeV) + ((float(Ef_MeV)- float(Ei_MeV))/ (Emod- 1))* i1

    # Start Doing The Simulation
    for it_do in range(Emod):
        for it_do_EachE in range(Emod_EachE):
            # ============================
            # Assign Random ri_km and magneticField_polynomialDecayDistance
            ri_km= ri_randStart_km+ 0.1* random.random()
            magneticField_polynomialDecayDistance= 49.95+ 0.1* random.random()
            # ============================

            # ============================
            # Modify New ri_km and magneticField_polynomialDecayDistance in physicalParameters.dat
            file= open(COLLECTIVE_NU_OSC_DIR+ 'parameters/physicalParameters.dat', 'r')
            fileWrite= open(COLLECTIVE_NU_OSC_DIR+ 'parameters/physicalParameters2.dat', 'w')
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
                    line= 'distance_initial='+ str(ri_km)+ '  # [km]\n'
                    fileWrite.write(line)
                # magneticField_polynomialDecayDistance
                elif (line.startswith('magneticField_polynomialDecayDistance')):
                    line= 'magneticField_polynomialDecayDistance='+ str(magneticField_polynomialDecayDistance)+ ' # [km] Valid for (use_defaultMagneticProfile=1, magneticField_profile=2)\n'
                    fileWrite.write(line)
                else:
                    fileWrite.write(line)
            file.close()
            fileWrite.close()
            # Remove and Rename
            osCommand.remove(COLLECTIVE_NU_OSC_DIR+ 'parameters/physicalParameters.dat')
            osCommand.rename(COLLECTIVE_NU_OSC_DIR+ 'parameters/physicalParameters2.dat', COLLECTIVE_NU_OSC_DIR+ 'parameters/physicalParameters.dat')
            # ============================
            osCommand.chdir(COLLECTIVE_NU_OSC_DIR)
            #! ============================
            #! WORKS ONLY FOR UNIX SYSTEMS
            osCommand.system('python3 '+COLLECTIVE_NU_OSC_DIR+'collNuPy.py')
            #! ============================
    exit()






