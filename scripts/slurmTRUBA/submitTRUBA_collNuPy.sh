#!/bin/bash

# Add Python3.6
module load centos7.3/comp/python/3.6.5-gcc

# User Name
USR_NAM=""
# User Email
USR_EMAIL=""

# Simulation Name (Also Job Name)
# if no arguments are given, print usage
if [ $# -eq 0 ]; then
    echo "Argument is missing, please give initial energy <Simulation Name> as first argument."
    echo "Usage: $0 <Simulation Name>"
    exit 0
fi
SIM_NAM=$1

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

# Python Script
cd ${WORK_DIR}

# Run the main code
python3.6 collNuPy.py

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