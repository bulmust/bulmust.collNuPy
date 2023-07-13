#!/bin/bash

# Add Python3.6
module load centos7.3/comp/python/3.6.5-gcc

# User Name
USR_NAM=""
# User Email
USR_EMAIL=""

# ----------------------------
# Argument Check
SIM_NAM=$1
# if no arguments are given, print usage
if [ $# -eq 0 ]; then
    echo "Argument is missing, please give initial energy <Simulation Name> as first argument."
    echo "Usage: $0 <Simulation Name> <Ei_MeV> <Ef_MeV> <Emod> <Emod_EachE> <distance_initial_initial> <distance_initial_final> <magneticField_polynomialDecayDistance_initial> <magneticField_polynomialDecayDistance_final> <randomizingRmagR>"
    exit 0
fi

# Second Argument is Initial Energy, Ei_MeV
Ei_MeV=$2
# If $2 is missing, exit
if [ -z "$Ei_MeV" ]; then
    echo "Argument is missing, please give initial energy <Ei_MeV> as second argument."
    echo "Usage: $0 <Simulation Name> <Ei_MeV> <Ef_MeV> <Emod> <Emod_EachE> <distance_initial_initial> <distance_initial_final> <magneticField_polynomialDecayDistance_initial> <magneticField_polynomialDecayDistance_final> <randomizingRmagR>"
    exit 0
fi

# Third Argument is Final Energy, Ef_MeV
Ef_MeV=$3
# If $3 is missing, exit
if [ -z "$Ef_MeV" ]; then
    echo "Argument is missing, please give final energy <Ef_MeV> as third argument."
    echo "Usage: $0 <Simulation Name> <Ei_MeV> <Ef_MeV> <Emod> <Emod_EachE> <distance_initial_initial> <distance_initial_final> <magneticField_polynomialDecayDistance_initial> <magneticField_polynomialDecayDistance_final> <randomizingRmagR>"
    exit 0
fi

# Fourth Argument is Number of Energy Points or Energy Bin Emod
Emod=$4
# If $4 is missing, exit
if [ -z "$Emod" ]; then
    echo "Argument is missing, please give number of energy points as <Emod> fourth argument."
    echo "Usage: $0 <Simulation Name> <Ei_MeV> <Ef_MeV> <Emod> <Emod_EachE> <distance_initial_initial> <distance_initial_final> <magneticField_polynomialDecayDistance_initial> <magneticField_polynomialDecayDistance_final> <randomizingRmagR>"
    exit 0
fi

# Fifth Argument is Number of energy mod for each energy Emod_EachE
Emod_EachE=$5
# If $5 is missing, exit
if [ -z "$Emod_EachE" ]; then
    echo "Argument is missing, please give number of energy mod for each energy as <Emod_EachE> fifth argument."
    echo "Usage: $0 <Simulation Name> <Ei_MeV> <Ef_MeV> <Emod> <Emod_EachE> <distance_initial_initial> <distance_initial_final> <magneticField_polynomialDecayDistance_initial> <magneticField_polynomialDecayDistance_final> <randomizingRmagR>"
    exit 0
fi

# Sixth Argument is Number of energy mod for each energy distance_initial_initial
distance_initial_initial=$6
# If $6 is missing, exit
if [ -z $distance_initial_initial]; then
    echo "Argument is missing, <distance_initial_initial> sixth argument."
    echo "Usage: $0 <Simulation Name> <Ei_MeV> <Ef_MeV> <Emod> <Emod_EachE> <distance_initial_initial> <distance_initial_final> <magneticField_polynomialDecayDistance_initial> <magneticField_polynomialDecayDistance_final> <randomizingRmagR>"
    exit 0
fi

# Seventh Argument is Number of energy mod for each energy distance_initial_final
distance_initial_final=$7
# If $7 is missing, exit
if [ -z $distance_initial_final]; then
    echo "Argument is missing, <distance_initial_final> seventh argument."
    echo "Usage: $0 <Simulation Name> <Ei_MeV> <Ef_MeV> <Emod> <Emod_EachE> <distance_initial_initial> <distance_initial_final> <magneticField_polynomialDecayDistance_initial> <magneticField_polynomialDecayDistance_final> <randomizingRmagR>"
    exit 0
fi

# Eighth Argument is Number of energy mod for each energy magneticField_polynomialDecayDistance_initial
magneticField_polynomialDecayDistance_initial=$8
# If $8 is missing, exit
if [ -z $magneticField_polynomialDecayDistance_initial]; then
    echo "Argument is missing, <magneticField_polynomialDecayDistance_initial> eighth argument."
    echo "Usage: $0 <Simulation Name> <Ei_MeV> <Ef_MeV> <Emod> <Emod_EachE> <distance_initial_initial> <distance_initial_final> <magneticField_polynomialDecayDistance_initial> <magneticField_polynomialDecayDistance_final> <randomizingRmagR>"
    exit 0
fi

# Ninth Argument is Number of energy mod for each energy magneticField_polynomialDecayDistance_final
magneticField_polynomialDecayDistance_final=$9
# If $9 is missing, exit
if [ -z $magneticField_polynomialDecayDistance_final]; then
    echo "Argument is missing, <magneticField_polynomialDecayDistance_final> ninth argument."
    echo "Usage: $0 <Simulation Name> <Ei_MeV> <Ef_MeV> <Emod> <Emod_EachE> <distance_initial_initial> <distance_initial_final> <magneticField_polynomialDecayDistance_initial> <magneticField_polynomialDecayDistance_final> <randomizingRmagR>"
    exit 0
fi

# Tenth Argument is Number of energy mod for each energy randomizingRmagR
randomizingRmagR=${10}
# If ${10} is missing, exit
if [ -z $randomizingRmagR]; then
    echo "Argument is missing, <randomizingRmagR> tenth argument."
    echo "Usage: $0 <Simulation Name> <Ei_MeV> <Ef_MeV> <Emod> <Emod_EachE> <distance_initial_initial> <distance_initial_final> <magneticField_polynomialDecayDistance_initial> <magneticField_polynomialDecayDistance_final> <randomizingRmagR>"
    exit 0
fi
# ----------------------------

# ----------------------------
# Printing All Arguments
echo "Simulation Name: $SIM_NAM"
echo "Initial Energy: $Ei_MeV MeV"
echo "Final Energy: $Ef_MeV MeV"
echo "Number of Energy Points: $Emod"
echo "Number of energy mod for each energy: $Emod_EachE"
echo "distance_initial_initial: $distance_initial_initial"
echo "distance_initial_final: $distance_initial_final"
echo "magneticField_polynomialDecayDistance_initial: $magneticField_polynomialDecayDistance_initial"
echo "magneticField_polynomialDecayDistance_final: $magneticField_polynomialDecayDistance_final"
echo "Take Random R and Rmag: $randomizingRmagR"
# ----------------------------

# ----------------------------
# Orjinal Code Directory
ORJ_DIR="/truba/home/${USR_NAM}/collNuPy"
# Code Dir in Scratch
WORK_DIR="/truba_scratch/${USR_NAM}/collNuPy_${SIM_NAM}"

# Check if a directory does exist. If so, exit
if [ -f "${WORK_DIR}" ]; then
    echo "Simulation $WORK_DIR DOES NOT exists."
    exit 0 # die with error code 0
fi

# Copy collnu Orjinal Code to Job Folder execpt git folder
rsync -av ${ORJ_DIR}/ ${WORK_DIR}/ --exclude .git

# Change folder to script
cd ${WORK_DIR}

# slurm file location
SLURM_FILE="${WORK_DIR}/collnu_${SIM_NAM}.slurm"
# ----------------------------

# ----------------------------
# Create slurm file
echo "#!/bin/bash

#SBATCH --job-name=${SIM_NAM}    # Job name
#SBATCH --partition=single       # Partition requested
#SBATCH --account=${USR_NAM}     # Charge job to specified account
#SBATCH --ntasks=1               # Number Of Task
#SBATCH --cpus-per-task=1        # Number of cpus required per task
#SBATCH --time=15-00:00:00       # Time limit, MANDATORY (Format:days-hours:minutes:seconds).
#SBATCH --workdir=${WORK_DIR}
#SBATCH --output=output.log      # Output File Name
#SBATCH --error=errors.log       # Error File Name
#SBATCH --qos=normal             # Quality of service
#SBATCH --mail-type=ALL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=${USR_EMAIL} # Where to send mail
#SBATCH --no-requeue             # Not continue, Required for non-resumable codes.

echo '======================================================'
echo HOST      : \$SLURM_SUBMIT_HOST
echo JOB ID    : \$SLURM_JOB_ID
echo JOB NAME  : \$SLURM_SIM_NAME
echo NODE LIST : \$SLURM_NODELIST
echo '======================================================'

# Go Python Script
cd ${WORK_DIR}/scripts/random_ri_magPolDecDist/

# Run code
python3.6 run_random_ri_magPolDecDist.py -i ${Ei_MeV} -f ${Ef_MeV} -m ${Emod} -b ${Emod_EachE} -t ${distance_initial_initial} -y ${distance_initial_final} -u ${magneticField_polynomialDecayDistance_initial} -k ${magneticField_polynomialDecayDistance_final} -r ${randomizingRmagR}

# Time limit for levrek cluster
# Kuyruklar (Partitions)
# Partition	Node	Cekirdek	Is Suresi	  Sunucular
# single	38	  5x16	    15-00:00:00	levrek[2-12,102-128]
# short	    38	  10x16	    0-4:00:00	levrek[2-12,102-128]
# mid1	    38	  20x16	    4-00:00:00	levrek[2-12,102-128]
# mid2	    38	  20x16	    8-00:00:00	levrek[2-12,102-128]
# long	    38	  50x16	    15-00:00:00	levrek[2-12,102-128]
# Minimum amount of real memory
# #SBATCH --mem=12000
# Maximum amount of real memory per allocated  cpu required by the job
# #SBATCH --mem-per-cpu=100

" > $SLURM_FILE

# Run the code
sbatch $SLURM_FILE

# Create Backup at Home
cp -r ${WORK_DIR} ${HOME}
# ----------------------------