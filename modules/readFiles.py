"""
ID      : readFiles.py
Author  : Taygun Bulmus
E-Mail  : bulmust@gmail.com
-----
"""
#import numpy as np
# Colors For Printing
class tLog:
    INFO  = '\033[94m'+'[INFO]    '+'\033[0m'
    OK      = '\033[92m'+'[OK]        '+'\033[0m'
    WARNING = '\033[93m'+'[WARNING]   '+'\033[0m'
    ERROR   = '\033[91m'+'[ERROR]     '+'\033[0m'
    EXIT    = '\033[91m'+'[EXIT]      readFiles.py'+'\033[0m'
#! ============================================================
#! WORKS ONLY FOR UNIX SYSTEMS
def technicalParameters_read():
    print(tLog.OK+'Start reading technicalParameters.dat file.')
    # Name Of File
    namOfFile='parameters/technicalParameters.dat'

    # Read File
    try:
        fileID=open(namOfFile,'r')
    except IOError:
        # ERROR te technicalParameters.dat file can not be opened
        print(tLog.ERROR+'technicalParameters.dat file can not be opened.')
        exit(tLog.EXIT)

    # Read the Lines
    fileLines = fileID.readlines()

    # ----------------------------------------------------
    # hamParam Variable Names
    # !!! Add new variables at the and of VarNam
    # With Number Comment
    # Same for two and three flavor
    technicalParametersNames= ('holdData_Every',
        'holdIntermediateData',
        'plotGraphs',
        'plotGraphs_diagonalRhoAllEnergy_2distance',
        'plotGraphs_diagonalRhoFinal_2energy',
        'plotGraphs_hamiltonianAllEnergy_2distance',
        'plot_savingFormat',
        'output_humanReadable',
        'output_distanceHamiltonianAllE',
        'output_distance_eigenValuesAllE',
        'initialValueProblem_solverMethod',
        'tolerance_relativeError',
        'tolerance_absoluteError')

    # If parameters' value must be 0 or 1, add their indexes below
    err01=('holdIntermediateData','plotGraphs','plotGraphs_diagonalRhoAllEnergy_2distance'
        , 'plotGraphs_diagonalRhoFinal_2energy'\
        , 'plotGraphs_hamiltonianAllEnergy_2distance'\
        , 'output_humanReadable','output_distanceHamiltonianAllE'
        , 'output_distance_eigenValuesAllE')

    # Name of parameters that have string values
    strParams= ('plot_savingFormat', 'initialValueProblem_solverMethod')

    # Variable Length
    technicalParametersNamesLen=len(technicalParametersNames)
    technicalParameters=[None]*technicalParametersNamesLen

    # Start Reading
    # (Does Not Read Strings)
    for line in range(len(fileLines)):
        # Do not Read Comments with started '#'
        if not fileLines[line].lstrip().startswith('#'):
            # Find the index of "=" Sign
            eqSigIND=fileLines[line].find('=')

            # Check the name in the file and technicalParametersNames
            if fileLines[line].lstrip()[0:eqSigIND] in technicalParametersNames:
                # Check there is an '#' sign in the line
                if fileLines[line].find('#')!=-1:
                    # Find the index of '#' sign
                    commentIND=fileLines[line].find('#')
                    # Read the value
                    # If RHS is not a number
                    if fileLines[line].lstrip()[0:eqSigIND] in strParams:
                        technicalParameters[technicalParametersNames.index(fileLines[line].lstrip()[0:eqSigIND])]=fileLines[line].lstrip()[eqSigIND+1:commentIND]
                        # Strip the spaces
                        technicalParameters[technicalParametersNames.index(fileLines[line].lstrip()[0:eqSigIND])]=\
                            technicalParameters[technicalParametersNames.index(fileLines[line].lstrip()[0:eqSigIND])].strip()
                    else:
                        technicalParameters[technicalParametersNames.index(fileLines[line].lstrip()[0:eqSigIND])]=float(fileLines[line].lstrip()[eqSigIND+1:commentIND])
                else:
                    # Read the value
                    # If RHS is not a number
                    if fileLines[line].lstrip()[0:eqSigIND] in strParams:
                        technicalParameters[technicalParametersNames.index(fileLines[line].lstrip()[0:eqSigIND])]=fileLines[line].lstrip()[eqSigIND+1:]
                        # Strip the spaces
                        technicalParameters[technicalParametersNames.index(fileLines[line].lstrip()[0:eqSigIND])]=\
                            technicalParameters[technicalParametersNames.index(fileLines[line].lstrip()[0:eqSigIND])].strip()
                    else:
                        technicalParameters[technicalParametersNames.index(fileLines[line].lstrip()[0:eqSigIND])]=float(fileLines[line].lstrip()[eqSigIND+1:])

    # ERROR Contradictions for technicalParameters.dat
    if None in technicalParameters:
        print(tLog.ERROR+'There is a contradiction between technicalParameters.dat file'
            'and technicalParametersNames variable.')
        exit(tLog.EXIT)

    # Create Dictionary
    technicalParametersDic={}; it=0
    for key in technicalParametersNames:
        technicalParametersDic[key] = technicalParameters[it]
        it+=1
    
    # ERROR Bool Values for technicalParameters
    for i2 in range(len(err01)):
        if technicalParametersDic[err01[i2]] != 1.0\
            and technicalParametersDic[err01[i2]] != 0.0:
            print(tLog.ERROR+'Invalid number for '+technicalParametersNames[i2]+'.')
            exit(tLog.EXIT)
        else:
            if technicalParametersDic[err01[i2]] == 0.0:
                technicalParametersDic[err01[i2]] = False
            else:
                technicalParametersDic[err01[i2]] = True

    # If plotGraphs is False, then set all other plotGraphs related parameters to False
    if not technicalParametersDic['plotGraphs']:
        del technicalParametersDic['plotGraphs_diagonalRhoAllEnergy_2distance']
        del technicalParametersDic['plotGraphs_diagonalRhoFinal_2energy']
        del technicalParametersDic['plotGraphs_hamiltonianAllEnergy_2distance']
    # Print The names and their values
    print('\n\033[94m -------- Technical Parameters -------- \033[0m')
    print("\n".join("{}={}".format(k, v) for k, v in technicalParametersDic.items()))
    print(' -------------------------------------- \n')
    
    print(tLog.OK+'Finish reading technicalParameters.dat file.')
    return technicalParametersDic
# ============================================================
def physicalParameters_read():
    
    # Name Of File
    namOfFile='parameters/physicalParameters.dat'

    print(tLog.OK+'Start reading physicalParameters.dat file.')
    # Name Of File
    namOfFile='parameters/physicalParameters.dat'

    # Read File
    try:
        fileID=open(namOfFile,'r')
    except IOError:
        # ERROR physicalParameters.dat file can not be opened
        print(tLog.ERROR+'physicalParameters.dat file can not be opened.')
        exit(tLog.EXIT)

    # Read the Lines
    fileLines = fileID.readlines()

    # ----------------------------------------------------
    # hamParam Variable Names
    # !!! Add new variables at the and of VarNam
    # With Number Comment
    # Same for two and three flavor
    physicalParametersNames= ('hamiltonian_oscillation',
        'hamiltonian_matter',
        'hamiltonian_selfInteraction_singleAngle',
        'hamiltonian_electromagnetic',
        'flavor_number',
        'neutrino_hierarchy',
        'theta12',
        'theta23',
        'theta13',
        'theta_14',
        'theta_24',
        'theta_34',
        'CP_Phase',
        'deltaM_square21',
        'deltaM_square32',
        'deltaM_square41',
        'neutrino_MagneticMoment',
        'numberOf_energyMode',
        'energy_initial',
        'energy_final',
        'distance_initial',
        'distance_final',
        'neutrino_distributionParameter',
        'luminosity_e',
        'luminosity_mu',
        'luminosity_tau',
        'luminosity_sterile',
        'luminosity_eb',
        'luminosity_mub',
        'luminosity_taub',
        'luminosity_sterileb',
        'temperature_e',
        'temperature_mu',
        'temperature_tau',
        'temperature_sterile',
        'temperature_eb',
        'temperature_mub',
        'temperature_taub',
        'temperature_sterileb',
        'use_defaultMatterProfile',
        'matterDensity_profile',
        'matterDensity_initial',
        'matterDensity_exponentialDecay',
        'matterDensity_polynomialDecay',
        'electronFraction_constant',
        'matterProfile_fileName',
        'use_defaultMagneticProfile',
        'magneticField_profile',
        'magneticField_initial',
        'magneticField_exponentialDecay',
        'magneticField_polynomialDecayPower',
        'magneticField_polynomialDecayDistance',
        'magneticField_fileName',
        'couplingConstant_sterile')

    # If parameters' value must be 0 or 1, add their indexes below
    err01=('hamiltonian_oscillation','hamiltonian_matter','hamiltonian_selfInteraction_singleAngle'\
        , 'hamiltonian_electromagnetic', 'use_defaultMatterProfile'\
        , 'use_defaultMagneticProfile')

    # Name of parameters that have string values
    strParams= ('matterProfile_fileName', 'magneticField_fileName')

    # Parameters that must be integer values
    intParams= ('flavor_number', 'neutrino_hierarchy', 'numberOf_energyMode'\
        , 'neutrino_distributionParameter'\
        , 'matterDensity_profile', 'magneticField_profile')

    # Variable Length
    physicalParametersNamesLen=len(physicalParametersNames)
    physicalParameters=[None]*physicalParametersNamesLen

    # Start Reading
    # (Does Not Read Strings)
    for line in range(len(fileLines)):
        # Do not Read Comments with started '#'
        if not fileLines[line].lstrip().startswith('#'):
            # Find the index of "=" Sign
            eqSigIND=fileLines[line].find('=')

            # Check the name in the file and physicalParametersNames
            if fileLines[line].lstrip()[0:eqSigIND] in physicalParametersNames:
                # Check there is an '#' sign in the line
                if fileLines[line].find('#')!=-1:
                    # Find the index of '#' sign
                    commentIND=fileLines[line].find('#')
                    # Read the value
                    # If RHS is not a number
                    if fileLines[line].lstrip()[0:eqSigIND] in strParams:
                        physicalParameters[physicalParametersNames.index(fileLines[line].lstrip()[0:eqSigIND])]=fileLines[line].lstrip()[eqSigIND+1:commentIND]
                    else:
                        # If RHS is a integer number
                        if fileLines[line].lstrip()[0:eqSigIND] in intParams:
                            try:
                                int(fileLines[line].lstrip()[eqSigIND+1:commentIND])
                            except ValueError:
                                # ERROR If RHS is not a integer number
                                print(tLog.ERROR+'The value of '\
                                    + fileLines[line].lstrip()[0:eqSigIND]\
                                    + ' must be an integer.')
                                exit(tLog.EXIT)
                            physicalParameters[physicalParametersNames.index(fileLines[line].lstrip()[0:eqSigIND])]=int(fileLines[line].lstrip()[eqSigIND+1:commentIND])
                            
                        else:
                            physicalParameters[physicalParametersNames.index(fileLines[line].lstrip()[0:eqSigIND])]=float(fileLines[line].lstrip()[eqSigIND+1:commentIND])
                else:
                    # Read the value
                    # If RHS is not a number
                    if fileLines[line].lstrip()[0:eqSigIND] in strParams:
                        physicalParameters[physicalParametersNames.index(fileLines[line].lstrip()[0:eqSigIND])]=fileLines[line].lstrip()[eqSigIND+1:]
                    else:
                        # If RHS is a integer number
                        if fileLines[line].lstrip()[0:eqSigIND] in intParams:
                            physicalParameters[physicalParametersNames.index(fileLines[line].lstrip()[0:eqSigIND])]=int(fileLines[line].lstrip()[eqSigIND+1:])
                        else:
                            physicalParameters[physicalParametersNames.index(fileLines[line].lstrip()[0:eqSigIND])]=float(fileLines[line].lstrip()[eqSigIND+1:])

    # ERROR Contradictions
    if None in physicalParameters:
        print(physicalParameters)
        print(tLog.ERROR+'There is a contradiction between'\
            'physicalParameters.dat file and physicalParametersNames variable.')
        exit(tLog.EXIT)

    # Create Dictionary
    physicalParametersDic={}; it=0
    for key in physicalParametersNames:
        physicalParametersDic[key] = physicalParameters[it]
        it+=1
    
    #! ============================
    # ERROR Bool Values for physicalParameters
    # 0 or 1 Values
    for i2 in range(len(err01)):
        if physicalParametersDic[err01[i2]] != 1.0\
            and physicalParametersDic[err01[i2]] != 0.0:
            print(tLog.ERROR+'Invalid number for '+physicalParametersNames[i2]+'.')
            exit(tLog.EXIT)
        else:
            if physicalParametersDic[err01[i2]] == 0.0:
                physicalParametersDic[err01[i2]] = False
            else:
                physicalParametersDic[err01[i2]] = True
    # ERROR If All Hamiltionian are False
    if physicalParametersDic['hamiltonian_oscillation']\
        and physicalParametersDic['hamiltonian_matter']\
        and physicalParametersDic['hamiltonian_selfInteraction_singleAngle']\
        and physicalParametersDic['hamiltonian_electromagnetic']:
        print(tLog.ERROR+'All Hamiltonians are False.')
        exit(tLog.EXIT)
            
    # Conditions
    if not physicalParametersDic['hamiltonian_oscillation']:
        del physicalParametersDic['theta12']
        del physicalParametersDic['theta23']
        del physicalParametersDic['theta13']
        del physicalParametersDic['theta_14']    
        del physicalParametersDic['theta_24']
        del physicalParametersDic['theta_34']
        del physicalParametersDic['CP_Phase']
    
    if not physicalParametersDic['hamiltonian_matter']:
        del physicalParametersDic['use_defaultMatterProfile']
        del physicalParametersDic['matterDensity_profile']
        del physicalParametersDic['matterDensity_initial']
        del physicalParametersDic['matterDensity_exponentialDecay']
        del physicalParametersDic['matterDensity_polynomialDecay']
        del physicalParametersDic['electronFraction_constant']
        del physicalParametersDic['matterProfile_fileName']
    else:
        if not physicalParametersDic['use_defaultMatterProfile']:
            del physicalParametersDic['matterDensity_profile']
            del physicalParametersDic['matterDensity_initial']
            del physicalParametersDic['matterDensity_exponentialDecay']
            del physicalParametersDic['matterDensity_polynomialDecay']
            del physicalParametersDic['electronFraction_constant']
        else:
            del physicalParametersDic['matterProfile_fileName']
            # ERROR Invalid value for electron fraction
            if (physicalParametersDic['electronFraction_constant'] > 1.0\
                or physicalParametersDic['electronFraction_constant'] < 0):
                print(tLog.ERROR+'Invalid value for electronFraction_constant.')
                exit(tLog.EXIT)

    if not physicalParametersDic['hamiltonian_electromagnetic']:
        del physicalParametersDic['magneticField_profile']
        del physicalParametersDic['magneticField_initial']
        del physicalParametersDic['magneticField_exponentialDecay']
        del physicalParametersDic['magneticField_polynomialDecayPower']
        del physicalParametersDic['magneticField_polynomialDecayDistance']
        del physicalParametersDic['magneticField_fileName']
    else:
        if not physicalParametersDic['use_defaultMagneticProfile']:
            del physicalParametersDic['magneticField_profile']
            del physicalParametersDic['magneticField_initial']
            del physicalParametersDic['magneticField_exponentialDecay']
            del physicalParametersDic['magneticField_polynomialDecayPower']
            del physicalParametersDic['magneticField_polynomialDecayDistance']
        else:
            del physicalParametersDic['magneticField_fileName']
    
    if (not physicalParametersDic['hamiltonian_selfInteraction_singleAngle'])\
        or  physicalParametersDic['flavor_number'] != 4:
        del physicalParametersDic['couplingConstant_sterile']   

    if physicalParametersDic['flavor_number'] != 2 and\
        physicalParametersDic['flavor_number'] != 3 and\
        physicalParametersDic['flavor_number'] != 4:
        print(physicalParametersDic['flavor_number'])
        print(tLog.ERROR+'Invalid number for flavor_number.')
        exit(tLog.EXIT)
    
    # ERROR Invalid value for neutrino hierarchy
    if physicalParametersDic['neutrino_hierarchy'] != 1 and\
        physicalParametersDic['neutrino_hierarchy'] != -1:
        print(tLog.ERROR+'Invalid value for neutrino_hierarchy.')
        exit(tLog.EXIT)
    
    if physicalParametersDic['flavor_number'] == 2:
        del physicalParametersDic['theta23']
        del physicalParametersDic['theta13']
        del physicalParametersDic['theta_14']    
        del physicalParametersDic['theta_24']
        del physicalParametersDic['theta_34']
        del physicalParametersDic['deltaM_square32']
        del physicalParametersDic['deltaM_square41']
        del physicalParametersDic['luminosity_tau']
        del physicalParametersDic['luminosity_sterile']
        del physicalParametersDic['luminosity_taub']
        del physicalParametersDic['luminosity_sterileb']
        del physicalParametersDic['temperature_tau']
        del physicalParametersDic['temperature_sterile']
        del physicalParametersDic['temperature_taub']
        del physicalParametersDic['temperature_sterileb']

    if physicalParametersDic['flavor_number'] == 3:
        del physicalParametersDic['theta_14']    
        del physicalParametersDic['theta_24']
        del physicalParametersDic['theta_34']
        del physicalParametersDic['deltaM_square41']
        del physicalParametersDic['luminosity_sterile']
        del physicalParametersDic['luminosity_sterileb']
        del physicalParametersDic['temperature_sterile']
        del physicalParametersDic['temperature_sterileb']

    # ERROR Invalid number for neutrino distribution parameter
    if physicalParametersDic['neutrino_distributionParameter'] != 1 and\
        physicalParametersDic['neutrino_distributionParameter'] != 2 and\
        physicalParametersDic['neutrino_distributionParameter'] != 3 and\
        physicalParametersDic['neutrino_distributionParameter'] != 4 and\
        physicalParametersDic['neutrino_distributionParameter'] != 5 and\
        physicalParametersDic['neutrino_distributionParameter'] != 6 and\
        physicalParametersDic['neutrino_distributionParameter'] != 7:
        print(tLog.ERROR+'Invalid number for neutrino_distributionParameter.')
        exit(tLog.EXIT)

    if physicalParametersDic['neutrino_distributionParameter'] != 1 and\
        physicalParametersDic['neutrino_distributionParameter'] != 2 and\
        physicalParametersDic['neutrino_distributionParameter'] != 5:
        del physicalParametersDic['luminosity_e']
        del physicalParametersDic['luminosity_mu']
        del physicalParametersDic['luminosity_eb']
        del physicalParametersDic['luminosity_mub']        
        if physicalParametersDic['flavor_number'] != 2\
            and physicalParametersDic['flavor_number'] != 3:
            del physicalParametersDic['luminosity_tau']
            del physicalParametersDic['luminosity_taub']
            del physicalParametersDic['luminosity_sterile']
            del physicalParametersDic['luminosity_sterileb']

    # Print The names and their values
    print('\n\033[94m -------- Technical Parameters -------- \033[0m')
    print("\n".join("{}={}".format(k, v) for k, v in physicalParametersDic.items()))
    print(' -------------------------------------- \n')
    
    print(tLog.OK+'Finish reading technicalParameters.dat file.')
    return physicalParametersDic
#! ============================================================
