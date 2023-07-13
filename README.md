# collNuPy
Collective Neutrino Oscillations Code written in Python. It works with only Linux systems.

- [collNuPy](#collnupy)
- [Python Version](#python-version)
  - [Modules](#modules)
    - [Default Modules](#default-modules)
    - [Installation Needed](#installation-needed)
- [Usage](#usage)
- [Output](#output)

# Python Version
It works with Python3. It is tested with version `3.8.10`.

## Modules
### Default Modules
- `os`, `sys`
- `time`
- `datetime`
- `shutil`
- `math`
### Installation Needed
- `numpy` , Tested Version: 1.24.2
- `odeintw` , Tested Version: 1.0.0
- `scipy` , Tested Version: 1.10.0
- `matplotlib` , Tested Version: 3.7.0

# Usage
1. Modify parameter files in [parameters/physicalParameters.dat](parameters/physicalParameters.dat) and [parameters/technicalParameters.dat](parameters/technicalParameters.dat). 
   
   See [parameters/paramExplanations.md](parameters/paramExplanations.md) file for explanations.

2. Install modules.
   
   `pip install numpy odeintw scipy matplotlib`

3. Run `python3 collNuPy.py`.
   
# Output
1. All output will be in results/Simulation# folder.
2. 
