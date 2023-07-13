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
from random import uniform
import os as osCommand # File operations
import argparse
from pathlib import Path

# Current Directory
currDir = osCommand.getcwd() + "/"
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
    parser.add_argument("-b", "--NumberOfMod_EachE"\
        , help= "Number of Mod For Each Energy", type= int)
    parser.add_argument("-t", "--distance_initial_initial"\
        , help= "Initial value of possible distance_initial"\
            , type= float)
    parser.add_argument("-y", "--distance_initial_Final"\
        , help= "Final value of possible distance_initial\n"\
            " (must be bigger than distance_initial)"\
            , type= float)
    parser.add_argument("-u", "--magneticField_polynomialDecayDistance_initial"\
        , help= "Initial value of possible magneticField_polynomialDecayDistance"\
            , type= float)
    parser.add_argument("-k", "--magneticField_polynomialDecayDistance_final"\
        , help= "Final value of possible magneticField_polynomialDecayDistance\n"\
            " (must be bigger than magneticField_polynomialDecayDistance)"\
            , type= float)
    parser.add_argument("-r", "--randomizingRmagR"\
        , help= "Randomizing distance_initial"\
            " and magneticField_polynomialDecayDistance (0 for No, 1 for Yes)"\
            , type= int)
    # Read arguments from command line
    args = parser.parse_args()
    # Check for --Ei_MeV
    if args.Ei_MeV:
        Ei_MeV= float(args.Ei_MeV)
        print(tLog.INFO+'Initial energy is ', Ei_MeV, ' MeV.')
    else:
        print(tLog.INFO+'Initial energy is not given.')
        exit(tLog.EXIT)
    # Check for --Ef_MeV
    if args.Ef_MeV:
        Ef_MeV= float(args.Ef_MeV)
        print(tLog.INFO+'Final energy is ', Ef_MeV, ' MeV.')
    else:
        print(tLog.INFO+'Final energy is not given.')
        exit(tLog.EXIT)
    # Check for --NumberOfMod
    if args.NumberOfMod:
        Emod= int(args.NumberOfMod)
        print(tLog.INFO+'Number of energy mod is ', Emod, '.')
    else:
        print(tLog.INFO\
            +'Number of energy mod is not given.')
        exit(tLog.EXIT)
    # Check for --NumberOfMod_EachE
    if args.NumberOfMod_EachE:
        Emod_EachE= int(args.NumberOfMod_EachE)
        print(tLog.INFO+'Number of energy mod for each energy is ', Emod_EachE, '.')
    else:
        print(tLog.INFO\
            +'Number of energy mod for each energy is not given.')
        exit(tLog.EXIT)
    # Check for --distance_initial_initial
    if args.distance_initial_initial:
        distance_initial_initial= float(args.distance_initial_initial)
        print(tLog.INFO+'distance_initial_initial will be '\
            , distance_initial_initial, '.')
    else:
        print(tLog.INFO\
            +'distance_initial_initial is not given')
        exit(tLog.EXIT)
    # Check for --distance_initial_Final
    if args.distance_initial_Final:
        distance_initial_Final= float(args.distance_initial_Final)
        print(tLog.INFO+'distance_initial_Final will be '\
            , distance_initial_Final, '.')
    else:
        print(tLog.INFO\
            +'distance_initial_Final is not given.')
        exit(tLog.EXIT)
    # Check for --magneticField_polynomialDecayDistance_initial
    if args.magneticField_polynomialDecayDistance_initial:
        magneticField_polynomialDecayDistance_initial=\
            float(args.magneticField_polynomialDecayDistance_initial)
        print(tLog.INFO+'magneticField_polynomialDecayDistance_initial will be '\
            , magneticField_polynomialDecayDistance_initial, '.')
    else:
        print(tLog.INFO\
            +'magneticField_polynomialDecayDistance_initial is not given.')
        exit(tLog.EXIT)
    # Check for --magneticField_polynomialDecayDistance_final
    if args.magneticField_polynomialDecayDistance_final:
        magneticField_polynomialDecayDistance_final=\
            float(args.magneticField_polynomialDecayDistance_final)
        print(tLog.INFO+'magneticField_polynomialDecayDistance_final will be '\
            , magneticField_polynomialDecayDistance_final, '.')
    else:
        print(tLog.INFO\
            +'magneticField_polynomialDecayDistance_final is not given')
        exit(tLog.EXIT)
    # Check for --randomizingRmagR
    if args.randomizingRmagR:
        randomizingRmagR= bool(args.randomizingRmagR)
        print(tLog.INFO+'randomizingRmagR will be '\
            , randomizingRmagR, '.')
    elif args.randomizingRmagR== 0:
        randomizingRmagR= bool(args.randomizingRmagR)
        print(tLog.INFO+'randomizingRmagR will be '\
            , randomizingRmagR, '.')
    else:
        print(tLog.INFO\
            +'randomizingRmagR is not given.')
        exit(tLog.EXIT)
    # ============================

    # Energies
    energy= np.zeros(Emod)
    energy[0]= Ei_MeV
    for i1 in range(1, int(Emod)):
        # Create Energy
        energy[i1]= float(Ei_MeV) + ((float(Ef_MeV)- float(Ei_MeV))/ (Emod- 1))* i1

    # Start Doing The Simulation
    for it_do in range(Emod):
        for it_do_EachE in range(Emod_EachE):
            if randomizingRmagR:
                # ============================
                # Assign Random ri_km and magneticField_polynomialDecayDistance
                ri_km= uniform(distance_initial_initial, distance_initial_Final)
                magneticField_polynomialDecayDistance=\
                    uniform(magneticField_polynomialDecayDistance_initial\
                        , magneticField_polynomialDecayDistance_final)
                # ============================
            else:
                if Emod_EachE == 1:
                    if it_do == 0:
                        ri_km= distance_initial_initial
                        magneticField_polynomialDecayDistance=\
                            magneticField_polynomialDecayDistance_initial
                    else:
                        # ============================
                        # Assign Monoton Increasing ri_km 
                        # and magneticField_polynomialDecayDistance
                        ri_km+= (distance_initial_Final- distance_initial_initial)\
                            / (Emod- 1)
                        magneticField_polynomialDecayDistance+= \
                            (magneticField_polynomialDecayDistance_final\
                        - magneticField_polynomialDecayDistance_initial)/ (Emod- 1)
                        # ============================
            # ============================
            # Modify New ri_km and magneticField_polynomialDecayDistance
            # in physicalParameters.dat
            file= open(COLLECTIVE_NU_OSC_DIR+ 'parameters/physicalParameters.dat', 'r')
            fileWrite= open(COLLECTIVE_NU_OSC_DIR\
                + 'parameters/physicalParameters2.dat', 'w')
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
                    line= 'magneticField_polynomialDecayDistance='\
                        + str(magneticField_polynomialDecayDistance)\
                        + ' # [km] Valid for (use_defaultMagneticProfile=1, '\
                        'magneticField_profile=2)\n'
                    fileWrite.write(line)
                else:
                    fileWrite.write(line)
            file.close()
            fileWrite.close()
            # Remove and Rename
            osCommand.remove(COLLECTIVE_NU_OSC_DIR+ 'parameters/physicalParameters.dat')
            osCommand.rename(COLLECTIVE_NU_OSC_DIR+'parameters/physicalParameters2.dat'\
                , COLLECTIVE_NU_OSC_DIR+ 'parameters/physicalParameters.dat')
            # ============================
            osCommand.chdir(COLLECTIVE_NU_OSC_DIR)
            #! ============================
            #! WORKS ONLY FOR UNIX SYSTEMS
            osCommand.system('python3 '+COLLECTIVE_NU_OSC_DIR+'collNuPy.py')
            #! ============================
            
    exit()