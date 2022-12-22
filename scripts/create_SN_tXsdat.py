import numpy as np
#import matplotlib.pyplot as plt
import argparse
from pathlib import Path
import os as osCommand # File operations

# ============================
# Colors For Printing
class tLog:
    INFO  = '\033[94m'+'[INFO]    '+'\033[0m'
    OK      = '\033[92m'+'[OK]        '+'\033[0m'
    WARNING = '\033[93m'+'[WARNING]   '+'\033[0m'
    ERROR   = '\033[91m'+'[ERROR]     '+'\033[0m'
    EXIT    = '\033[91m'+'[EXIT]       create_SN_tXsdat.py' + '\033[0m'
# ============================

# ============================
# Create Ye as a Step Function
def stepYe_for_Athar1995cx(distAll_km, innerStep_km=695.700, outerStep_km=695700\
    , innerYe=0.35, midYe=0.4998, outerYe=0.75):
    print(tLog.INFO+'Creating Ye as a Step Function.')
    Ye= np.zeros(len(distAll_km))
    for i1 in range(len(distAll_km)):
        if (distAll_km[i1] <= innerStep_km):
            Ye[i1]= innerYe
        elif (distAll_km[i1] <= outerStep_km and distAll_km[i1] > innerStep_km):
            Ye[i1]= midYe
        elif (distAll_km[i1] > outerStep_km):
            Ye[i1]= outerYe
    print(tLog.OK+'Y_e Step Function Created.')
    return Ye
# Create Ye as a Smooth Function at Inner Area
def smoothInnerYe_for_Athar1995cx(distAll_km, innerStep_km=695.700, outerStep_km=695700\
    , innerYe=0.35, midYe=0.4998, outerYe=0.75):
    print(tLog.INFO+'Creating Ye as a Smooth Function at Inner Area.')
    Ye= np.zeros(len(distAll_km))
    for i1 in range(len(distAll_km)):
        if (distAll_km[i1] <= innerStep_km):
            Ye[i1]= ((midYe- innerYe)/ innerStep_km)* distAll_km[i1]+ innerYe
        elif (distAll_km[i1] <= outerStep_km and distAll_km[i1] > innerStep_km):
            Ye[i1]= midYe
        elif (distAll_km[i1] > outerStep_km):
            Ye[i1]= outerYe
    print(tLog.OK+'Y_e Smooth Function at Inner Area Created.')
    return Ye
# ============================

# ============================
# Create Shockwave Added Baryon Density On PreSN Data
def create_shockedNb_g_cm3(dataPreSN_PATH, tPostBounce_s):
    print(tLog.INFO+'Creating Shockwave Added Baryon Density On PreSN Data.')
    dataSN_t0s= np.loadtxt(dataPreSN_PATH)
    r_cm= dataSN_t0s[:, 0]
    r_km= r_cm* 1e-5; lenData= len(r_km)

    if tPostBounce_s < 1:
        print(tLog.INFO+'Shockwave Added Baryon Density On PreSN Data is not created.')
        print(tLog.ERROR+'Consider tPostBounce_s >= 1s.')
        exit(tLog.EXIT)
    elif tPostBounce_s >= 1:
        density_shock_g_cm3 = np.zeros((lenData))
        r_s_km= (-4.6e+3)+ ((11.3e+3)* tPostBounce_s)+ (0.5* (0.2e+3)* tPostBounce_s**2)

        for it1 in range(lenData):
            
            if r_km[it1] <= r_s_km:
                f_x_function = np.exp(0.28- (0.69* np.log(r_s_km))\
                    * ((np.arcsin(1- (r_km[it1]/ r_s_km)))**(1.1)) )
                    
                density_shock_g_cm3[it1]=  10* f_x_function* dataSN_t0s[it1, 2]
            elif r_km[it1] > r_s_km:
                density_shock_g_cm3[it1]=  dataSN_t0s[it1, 2]
        print(tLog.OK+'Shockwave Added Baryon Density On PreSN Data Created.')
        return density_shock_g_cm3
# ============================

if __name__ == '__main__':
    # Current Directory
    currDir = osCommand.getcwd() + "/"
    COLLECTIVE_NU_OSC_DIR= str(Path(currDir).parent)+ "/"

    # Initialize parser
    parser = argparse.ArgumentParser()

    # Adding optional argument
    parser.add_argument("-t", "--Time", help= "Time after bounce in seconds"\
        , type = float, required=True)
    parser.add_argument("-p", "--PathSN", help= "Path to PreSN Data", required=True)
    
    # Read arguments from command line
    args = parser.parse_args()

    # Check for --Time
    if args.Time:
        tPostBounce_s= float(args.Time)
        print(tLog.INFO+'Time after bounce is', tPostBounce_s, 'seconds.')
    else:
        print(tLog.ERROR+'Time after bounce is not provided.')
        exit(tLog.EXIT)
    
    # Check for --PathSN
    if args.PathSN:
        dataPreSN_PATH= args.PathSN
        print(tLog.INFO+'Path of PreSN Data is', dataPreSN_PATH)
    else:
        print(tLog.ERROR+'Path of PreSN Data is not provided.')
        exit(tLog.EXIT)
    
    outputPATH= COLLECTIVE_NU_OSC_DIR+"/backgroundProfiles/SN_distCM_Ye_Nb_g_cm3_t"\
        + str(tPostBounce_s)+"s.dat"

    dataSN_t0s= np.loadtxt(dataPreSN_PATH)
    r_cm= dataSN_t0s[:, 0]
    r_km= dataSN_t0s[:, 0]* 1e-5

    # Create Shockwave Added Baryon Density On PreSN Data
    density_shock_g_cm3= create_shockedNb_g_cm3(dataPreSN_PATH, tPostBounce_s)

    # Create Ye as a Step Function
    Ye= stepYe_for_Athar1995cx(r_km)

    # Create Output Data
    dataOut= np.zeros([len(r_cm),3])    
    for it1 in range(len(r_km)):
        dataOut[it1, 0]= r_cm[it1]
        dataOut[it1, 1]= Ye[it1]
        dataOut[it1, 2]= density_shock_g_cm3[it1]
    np.savetxt(outputPATH, dataOut, fmt='%1.6E')

    print(tLog.OK+'Output Data Created at '+outputPATH+'.')



