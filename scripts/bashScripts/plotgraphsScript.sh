#!/usr/bin/env bash

cd ../../
# collnu Directory
CURR_DIR=$(pwd)
if [ -z "$1" ]
  then
    read -p "Enter simulation name: " SIMU_NAME
  else
    SIMU_NAME=$1
fi

RESULTS_DIR = "$CURR_DIR/results/"
NPZ_FILE_PATH="$CURR_DIR/results/$SIMU_NAME/data/rhoFlavAll.npz"

echo "Try to reach file destination:"
echo $NPZ_FILE_PATH

### Check if a directory does not exist ###
if [ ! -f "$NPZ_FILE_PATH" ]; then
    echo "Simulation $SIMU_NAME DOES NOT exists."
    exit 9999 # die with error code 9999
fi
echo "Reached file destination:"

echo "Avaible plots are following"
echo "Diagonal elements of density matrix to energy         => plot_EnergyDiag   : a"
echo "Diagonal elements of density matrix to distance plot  => plot_DistDiag     : b"
echo ""
echo "Example: To plot plot_EnergyDiag and plot_Adiabaticity, enter ag"

read -p "Which graph do you want to plot? " plotCheck

python3 <<- EOF
from sys import path
path.insert(0, "$CURR_DIR/modules")
import plotgraphs as plot
print("The plotting is started.")

RESULTS_SIMULATION_DIR= "$CURR_DIR/results/$SIMU_NAME/"

# ======================================================================
# PLOTS
# ======================================================================
if 'a' in "$plotCheck":
  print('RUN energyDiag')
  plot.plotMain(RESULTS_SIMULATION_DIR).energyDiag()
  print('Finished plot_energyDiag')
if 'b' in "$plotCheck":
  print('RUN distDiag')
  plot.plotMain(RESULTS_SIMULATION_DIR).distDiag()
  print('Finished plot_distDiag')
# ======================================================================
EOF
