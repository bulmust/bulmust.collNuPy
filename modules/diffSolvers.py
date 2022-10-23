"""
-----
ID      : diffSolvers.py
Author  : Taygun Bulmus
E-Mail  : bulmust@gmail.com
-----
"""
import time
import datetime
import os as osCommand
from shutil import rmtree
import numpy as np
import RHS as func
import odeintw as integralSolver

#from .readFiles import technicalParameters_read
# Various Scripts
#import taygunsScripts as tScripts

# ============================
# Colors For Printing
class tLog:
    INFO  = '\033[94m'+'[INFO]    '+'\033[0m'
    OK      = '\033[92m'+'[OK]        '+'\033[0m'
    WARNING = '\033[93m'+'[WARNING]   '+'\033[0m'
    ERROR   = '\033[91m'+'[ERROR]     '+'\033[0m'
    EXIT    = '\033[91m'+'[EXIT]      diffSolvers.py'+'\033[0m'
# ============================

# ============================
# Scipy Differential Equation Solver: LSODA (odeintw)
def odeintwSolver(init, RESULTS_SIMULATION_DIR):
    """
    Definitions
    -----------
        Ordinary differential equation with initial value problem
         solver using modified default numpy ivp solver odeint called
         [odeintw](https://github.com/WarrenWeckesser/odeintw).
         odeintw is the modification of odeint core
         but it considers complex value inputs. The solving method is
         LSODA, see ref: [1](https://computing.llnl.gov/casc/odepack/) and
         [2](https://doi.org/10.1137/0904010)
    """
    # ============================
    # Start calculating time
    timeInitialClock = time.time()
    timeInitialDate  = datetime.datetime.now()
    # ============================
    # Create data folder and change the path
    DATA_FOLDER= RESULTS_SIMULATION_DIR+ 'data/'
    if not osCommand.path.exists(DATA_FOLDER):
        osCommand.makedirs(DATA_FOLDER)
        print(tLog.OK+"Data Folder Created.")
    else:
        print(tLog.ERROR+"Data Folder Already Exists. Check The Folder: "+ DATA_FOLDER)
        print(tLog.EXIT)
    if init.technicalParametersDic['holdIntermediateData']:
        # Create intermediateData folder
        INTERMEDIATE_DATA_FOLDER= RESULTS_SIMULATION_DIR+ 'data/intermediateData/'
        if not osCommand.path.exists(INTERMEDIATE_DATA_FOLDER):
            osCommand.makedirs(INTERMEDIATE_DATA_FOLDER)
            print(tLog.OK+"Intermediate Data Folder Created.")
        else:
            print(tLog.ERROR+"Intermediate Data Folder Already Exists. Check The Folder: "+ INTERMEDIATE_DATA_FOLDER)
            print(tLog.EXIT)
    # Change directory to data 
    # Code saves intermediate data
    osCommand.chdir(DATA_FOLDER)
    # ============================
    # Copy initial parameters
    #copy_tree('../../../parameters/', '.')
    #np.savez('rhoInit_flav.npz', rhoInit_flav=init.rhoInit_flav)
    # ============================  
    
    # ============================
    # Solve E.o.M.
    # Starting time
    print(tLog.INFO+'Start solving E.o.M with odeintw. Current date and time :', timeInitialDate.strftime("%Y-%m-%d %H:%M:%S"))
    # Hold intermediate data
    if init.technicalParametersDic['holdIntermediateData']:
        if init.dim_rho_2totFlav_Bool:
            # Default error tolerances
            if (init.technicalParametersDic['tolerance_relativeError'] == 0 and init.technicalParametersDic['tolerance_absoluteError'] == 0):
                if init.totFlav == 2:
                    # Two Flavor
                    rhoFinalAll,infow = integralSolver.odeintw(func.rhs2Flav_bigRho, init.rhoInit_flav, init.distAll_km, full_output=True, mxstep=500000000)
                elif init.totFlav == 3:
                    # Three Flavor
                    rhoFinalAll,infow = integralSolver.odeintw(func.rhs3Flav_bigRho, init.rhoInit_flav, init.distAll_km, full_output=True, mxstep=500000000)
                elif init.totFlav == 4:
                    # Four Flavor
                    rhoFinalAll,infow = integralSolver.odeintw(func.rhs4Flav_bigRho, init.rhoInit_flav, init.distAll_km, full_output=True, mxstep=500000000)
            # User defined error tolerances
            else:
                if init.totFlav == 2:
                    # Two Flavor
                    rhoFinalAll,infow = integralSolver.odeintw(func.rhs2Flav_bigRho, init.rhoInit_flav\
                        , init.distAll_km, full_output=True, rtol=init.technicalParametersDic['tolerance_relativeError']\
                        , atol=init.technicalParametersDic['tolerance_absoluteError']\
                        , mxstep=500000000)
                elif init.totFlav == 3:
                    # Three Flavor
                    rhoFinalAll,infow = integralSolver.odeintw(func.rhs3Flav_bigRho, init.rhoInit_flav\
                        , init.distAll_km, full_output=True, rtol=init.technicalParametersDic['tolerance_relativeError']\
                        , atol=init.technicalParametersDic['tolerance_absoluteError']\
                        , mxstep=500000000)
                elif init.totFlav == 4:
                    # Four Flavor
                    rhoFinalAll,infow = integralSolver.odeintw(func.rhs4Flav_bigRho, init.rhoInit_flav\
                        , init.distAll_km, full_output=True, rtol=init.technicalParametersDic['tolerance_relativeError']\
                        , atol=init.technicalParametersDic['tolerance_absoluteError']\
                        , mxstep=500000000)
        else:
            # Default error tolerances
            if (init.technicalParametersDic['tolerance_relativeError'] == 0 and init.technicalParametersDic['tolerance_absoluteError'] == 0):
                if init.totFlav == 2:
                    # Two Flavor
                    rhoFinalAll,infow = integralSolver.odeintw(func.rhs2Flav, init.rhoInit_flav, init.distAll_km, full_output=True, mxstep=500000000)
                elif init.totFlav == 3:
                    # Three Flavor
                    rhoFinalAll,infow = integralSolver.odeintw(func.rhs3Flav, init.rhoInit_flav, init.distAll_km, full_output=True, mxstep=500000000)
                elif init.totFlav == 4:
                    # Four Flavor
                    rhoFinalAll,infow = integralSolver.odeintw(func.rhs4Flav, init.rhoInit_flav, init.distAll_km, full_output=True, mxstep=500000000)    
            # User defined error tolerances
            else:
                if init.totFlav == 2:
                    # Two Flavor
                    rhoFinalAll,infow = integralSolver.odeintw(func.rhs2Flav, init.rhoInit_flav\
                        , init.distAll_km, full_output=True, rtol=init.technicalParametersDic['tolerance_relativeError']\
                        , atol=init.technicalParametersDic['tolerance_absoluteError']\
                        , mxstep=500000000)
                elif init.totFlav == 3:
                    # Three Flavor
                    rhoFinalAll,infow = integralSolver.odeintw(func.rhs3Flav, init.rhoInit_flav\
                        , init.distAll_km, full_output=True, rtol=init.technicalParametersDic['tolerance_relativeError']\
                        , atol=init.technicalParametersDic['tolerance_absoluteError']\
                        , mxstep=500000000)
                elif init.totFlav == 4:
                    # Four Flavor
                    rhoFinalAll,infow = integralSolver.odeintw(func.rhs4Flav, init.rhoInit_flav\
                        , init.distAll_km, full_output=True, rtol=init.technicalParametersDic['tolerance_relativeError']\
                        , atol=init.technicalParametersDic['tolerance_absoluteError']\
                        , mxstep=500000000)
    else:
        if init.dim_rho_2totFlav_Bool:
            # Default error tolerances
            if (init.technicalParametersDic['tolerance_relativeError'] == 0 and init.technicalParametersDic['tolerance_absoluteError'] == 0):
                if init.totFlav == 2:
                    # Two Flavor
                    rhoFinalAll,infow = integralSolver.odeintw(func.rhs2Flav_bigRho_noInterData, init.rhoInit_flav, init.distAll_km, full_output=True, mxstep=500000000)
                elif init.totFlav == 3:
                    # Three Flavor
                    rhoFinalAll,infow = integralSolver.odeintw(func.rhs3Flav_bigRho_noInterData, init.rhoInit_flav, init.distAll_km, full_output=True, mxstep=500000000)
                elif init.totFlav == 4:
                    # Four Flavor
                    rhoFinalAll,infow = integralSolver.odeintw(func.rhs4Flav_bigRho_noInterData, init.rhoInit_flav, init.distAll_km, full_output=True, mxstep=500000000)
            # User defined error tolerances
            else:
                if init.totFlav == 2:
                    # Two Flavor
                    rhoFinalAll,infow = integralSolver.odeintw(func.rhs2Flav_bigRho_noInterData, init.rhoInit_flav\
                        , init.distAll_km, full_output=True, rtol=init.technicalParametersDic['tolerance_relativeError']\
                        , atol=init.technicalParametersDic['tolerance_absoluteError']\
                        , mxstep=500000000)
                elif init.totFlav == 3:
                    # Three Flavor
                    rhoFinalAll,infow = integralSolver.odeintw(func.rhs3Flav_bigRho_noInterData, init.rhoInit_flav\
                        , init.distAll_km, full_output=True, rtol=init.technicalParametersDic['tolerance_relativeError']\
                        , atol=init.technicalParametersDic['tolerance_absoluteError']\
                        , mxstep=500000000)
                elif init.totFlav == 4:
                    # Four Flavor
                    rhoFinalAll,infow = integralSolver.odeintw(func.rhs4Flav_bigRho_noInterData, init.rhoInit_flav\
                        , init.distAll_km, full_output=True, rtol=init.technicalParametersDic['tolerance_relativeError']\
                        , atol=init.technicalParametersDic['tolerance_absoluteError']\
                        , mxstep=500000000)
        else:
            # Default error tolerances
            if (init.technicalParametersDic['tolerance_relativeError'] == 0 and init.technicalParametersDic['tolerance_absoluteError'] == 0):
                if init.totFlav == 2:
                    # Two Flavor
                    rhoFinalAll,infow = integralSolver.odeintw(func.rhs2Flav_noInterData, init.rhoInit_flav, init.distAll_km, full_output=True, mxstep=500000000)
                elif init.totFlav == 3:
                    # Three Flavor
                    rhoFinalAll,infow = integralSolver.odeintw(func.rhs3Flav_noInterData, init.rhoInit_flav, init.distAll_km, full_output=True, mxstep=500000000)
                elif init.totFlav == 4:
                    # Four Flavor
                    rhoFinalAll,infow = integralSolver.odeintw(func.rhs4Flav_noInterData, init.rhoInit_flav, init.distAll_km, full_output=True, mxstep=500000000)    
            # User defined error tolerances
            else:
                if init.totFlav == 2:
                    # Two Flavor
                    rhoFinalAll,infow = integralSolver.odeintw(func.rhs2Flav_noInterData, init.rhoInit_flav\
                        , init.distAll_km, full_output=True, rtol=init.technicalParametersDic['tolerance_relativeError']\
                        , atol=init.technicalParametersDic['tolerance_absoluteError']\
                        , mxstep=500000000)
                elif init.totFlav == 3:
                    # Three Flavor
                    rhoFinalAll,infow = integralSolver.odeintw(func.rhs3Flav_noInterData, init.rhoInit_flav\
                        , init.distAll_km, full_output=True, rtol=init.technicalParametersDic['tolerance_relativeError']\
                        , atol=init.technicalParametersDic['tolerance_absoluteError']\
                        , mxstep=500000000)
                elif init.totFlav == 4:
                    # Four Flavor
                    rhoFinalAll,infow = integralSolver.odeintw(func.rhs4Flav_noInterData, init.rhoInit_flav\
                        , init.distAll_km, full_output=True, rtol=init.technicalParametersDic['tolerance_relativeError']\
                        , atol=init.technicalParametersDic['tolerance_absoluteError']\
                        , mxstep=500000000)
    # ============================  
    
    # ============================  
    print(tLog.OK+'Finish solving E.o.M with odeintw. Current date and time :', timeInitialDate.strftime("%Y-%m-%d %H:%M:%S"))
    # Stop Calculating time
    timeFinalClock = time.time()
    # ============================
    # Calculate Elapsed Time
    tmpTimeCalculation = timeFinalClock-timeInitialClock
    # Define temporary time elapsed in second
    tmpTimeCalculationSec = 0
    # Hold it
    tmpTimeCalculationSec = tmpTimeCalculation
    # Print the time values
    days = tmpTimeCalculation // (24 * 3600)
    tmpTimeCalculation = tmpTimeCalculation % (24 * 3600)
    hours = tmpTimeCalculation // 3600
    tmpTimeCalculation %= 3600
    minutes = tmpTimeCalculation // 60
    tmpTimeCalculation %= 60
    seconds = tmpTimeCalculation
    print(tLog.INFO+"Total Process Time, %6.2fs,  d:h:m:s : %d:%d:%d:%d"% (tmpTimeCalculationSec, days, hours, minutes, seconds))
    # ============================
    # Save Total Process Time to file
    # Do not color the [INFO]
    fileObject = open('distanceAndTime.log','a+')
    print('Total Process Time, %6.2fs,  d:h:m:s : %d:%d:%d:%d'% (tmpTimeCalculationSec, days, hours, minutes, seconds), file=fileObject)
    fileObject.close()
    # ============================
    # Number Of Evolution
    print(tLog.INFO+'Number of evaluation: %d' % infow['nfe'][-1])
    # ============================
    # Check differential equation solver method changed or not
    methodIndicators = infow['mused']
    if methodIndicators[0] == 1:
        print(tLog.INFO+'Differential equation solver started with adams (nonstiff) method.')
    else:
        print(tLog.INFO+'Differential equation solver started with bdf (stiff) method.')
    indicatorChangePrint = [tLog.INFO+'Differential equation solver method is adams (nonstiff)'\
        , tLog.INFO+'Differential equation solver method is bdf   (stiff)   ']
    # ============================
    # If tmpIndicatorCheckControl == 1, the method changed.
    tmpIndicatorCheckControl = 0
    tmpIndicator = 0
    tmpIndicator = methodIndicators[0]
    # ============================
    # Check all steps, method is changed or not
    for it_methodIndicator in range(len(methodIndicators)):
        if methodIndicators[it_methodIndicator] != tmpIndicator:
            tmpIndicatorCheckControl = 1
            # Method changed
            if methodIndicators[it_methodIndicator] == 1:
                # To nonstiff
                print(str(indicatorChangePrint[0])+' at distance '\
                      +str(init.distAll_km[it_methodIndicator])+'.')
                tmpIndicator = methodIndicators[it_methodIndicator]
            if methodIndicators[it_methodIndicator] == 2:
                # To stiff
                print(str(indicatorChangePrint[1])+' at distance '\
                      +str(init.distAll_km[it_methodIndicator])+'.')
                tmpIndicator=methodIndicators[it_methodIndicator]
    # ============================
    # Print if the method is changed or not
    if tmpIndicatorCheckControl == 0:
        print(tLog.INFO+'Differential equation solver method have not changed.')
    elif tmpIndicatorCheckControl == 1:
        print(tLog.INFO+'Differential equation solver method have changed.')
    # ============================
    '''
    ALL INFOS
     'hu' vector of step sizes successfully used for each time step.
     'tcur' vector with the value of t reached for each time step.
        (will always be at least as large as the input times).
     'tolsf' vector of tolerance scale factors, greater than 1.0,
        computed when a request for too much accuracy was detected.
     'tsw' value of t at the time of the last method switch
        (given for each time step)
     'nst' cumulative number of time steps
     'nfe' cumulative number of function evaluations for each time step
     'nje' cumulative number of jacobian evaluations for each time step
     'nqu' a vector of method orders for each successful step.
     'imxer' index of the component of largest magnitude in the weighted
        local error vector (e / ewt) on an error return, -1 otherwise.
     'lenrw' the length of the double work array required.
     'leniw' the length of integer work array required.
     'mused' a vector of method indicators for each successful time step:
        1: adams (nonstiff), 2: bdf (stiff)
    '''
    # ============================
    if init.technicalParametersDic['holdIntermediateData']:
        # Remove intermediateData folder
        rmtree(INTERMEDIATE_DATA_FOLDER)
    # Return Back To RESULTS_SIMULATION_DIR Folder
    osCommand.chdir(RESULTS_SIMULATION_DIR)

    return rhoFinalAll