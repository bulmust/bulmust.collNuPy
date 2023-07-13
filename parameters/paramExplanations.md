- [Parameters](#parameters)
  - [Physical Parameters](#physical-parameters)
    - [`hamiltonian_oscillation`](#hamiltonian_oscillation)
    - [`hamiltonian_selfInteraction_singleAngle`](#hamiltonian_selfinteraction_singleangle)
    - [`hamiltonian_electromagnetic`](#hamiltonian_electromagnetic)
    - [`flavor_number`](#flavor_number)
    - [`neutrino_hierarchy`](#neutrino_hierarchy)
    - [`theta12`](#theta12)
    - [`theta23`](#theta23)
    - [`theta13`](#theta13)
    - [`theta_14`](#theta_14)
    - [`theta_24`](#theta_24)
    - [`theta_34`](#theta_34)
    - [`CP_Phase`](#cp_phase)
    - [`deltaM_square21`](#deltam_square21)
    - [`deltaM_square32`](#deltam_square32)
    - [`deltaM_square41`](#deltam_square41)
    - [`neutrino_MagneticMoment`](#neutrino_magneticmoment)
    - [`numberOf_energyMode`](#numberof_energymode)
    - [`energy_initial`](#energy_initial)
    - [`energy_final`](#energy_final)
    - [`distance_initial`](#distance_initial)
    - [`distance_final`](#distance_final)
    - [`neutrino_distributionParameter`](#neutrino_distributionparameter)
    - [`luminosity_e`](#luminosity_e)
    - [`luminosity_mu`](#luminosity_mu)
    - [`luminosity_tau`](#luminosity_tau)
    - [`luminosity_sterile`](#luminosity_sterile)
    - [`luminosity_eb`](#luminosity_eb)
    - [`luminosity_mub`](#luminosity_mub)
    - [`luminosity_taub`](#luminosity_taub)
    - [`luminosity_sterileb`](#luminosity_sterileb)
    - [`temperature_e`](#temperature_e)
    - [`temperature_mu`](#temperature_mu)
    - [`temperature_tau`](#temperature_tau)
    - [`temperature_sterile`](#temperature_sterile)
    - [`temperature_eb`](#temperature_eb)
    - [`temperature_mub`](#temperature_mub)
    - [`temperature_taub`](#temperature_taub)
    - [`temperature_sterileb`](#temperature_sterileb)
    - [`use_defaultMatterProfile`](#use_defaultmatterprofile)
    - [`matterDensity_profile`](#matterdensity_profile)
    - [`matterDensity_initial`](#matterdensity_initial)
    - [`matterDensity_exponentialDecay`](#matterdensity_exponentialdecay)
    - [`matterDensity_polynomialDecay`](#matterdensity_polynomialdecay)
    - [`electronFraction_constant`](#electronfraction_constant)
    - [`matterProfile_fileName`](#matterprofile_filename)
    - [`use_defaultMagneticProfile`](#use_defaultmagneticprofile)
    - [`magneticField_profile`](#magneticfield_profile)
    - [`magneticField_initial`](#magneticfield_initial)
    - [`magneticField_exponentialDecay`](#magneticfield_exponentialdecay)
    - [`magneticField_polynomialDecayPower`](#magneticfield_polynomialdecaypower)
    - [`magneticField_polynomialDecayDistance`](#magneticfield_polynomialdecaydistance)
    - [`magneticField_fileName`](#magneticfield_filename)
    - [`couplingConstant_sterile`](#couplingconstant_sterile)
  - [Technical Parameters](#technical-parameters)
    - [`holdData_Every`](#holddata_every)
    - [`holdIntermediateData`](#holdintermediatedata)
    - [`plotGraphs`](#plotgraphs)
    - [`plotGraphs_diagonalRhoAllEnergy_2distance`](#plotgraphs_diagonalrhoallenergy_2distance)
    - [`plotGraphs_diagonalRhoFinal_2energy`](#plotgraphs_diagonalrhofinal_2energy)
    - [`plotGraphs_hamiltonianAllEnergy_2distance`](#plotgraphs_hamiltonianallenergy_2distance)
    - [`plot_savingFormat`](#plot_savingformat)
    - [`output_distanceHamiltonianAllE`](#output_distancehamiltonianalle)
    - [`output_distance_eigenValuesAllE`](#output_distance_eigenvaluesalle)
    - [`output_humanReadable`](#output_humanreadable)
    - [`initialValueProblem_solverMethod`](#initialvalueproblem_solvermethod)
    - [`tolerance_relativeError`](#tolerance_relativeerror)
    - [`tolerance_absoluteError`](#tolerance_absoluteerror)
- [References](#references)

# Parameters
In parameter files,

1. Lines starting with "#" are ignored.
2. The format is "`label`=value".
3. Some parameters are valid for only specific values. Look descriptions.
4. Use 1.2e+4 format instead of 1.2*1e+4.
5. If you want to add a new Parameter, add the parameter name into readFiles.py.
6. Code does not consider some variables in specific conditions. Look descriptions.
7. Use dot (.) for decimal, not comma (,).

## Physical Parameters
### `hamiltonian_oscillation`
- Neutrino oscillation Hamiltonian.
- Acceptable values are only (0 or 1).
###`hamiltonian_matter`
- Neutrino-Matter interaction Hamiltonian.
- Acceptable values are only (0 or 1).
### `hamiltonian_selfInteraction_singleAngle`
- Neutrino-neutrino interaction Hamiltonian. 
- Acceptable values are only (0 or 1).
### `hamiltonian_electromagnetic`
- Neutrino-electromagnetic interaction Hamiltonian.
- Acceptable values are only (0 or 1).
### `flavor_number`
- Number of flavor.
- Acceptable values are only (2, 3, 4).
- 2 is for electron and x.
- 3 is for electron, muon and tau.
  - Mixing matrix: $U_{mix}=R_{23} R_{13} R_{12}$.
  - CP Violation phase is in the $R_{13}$.
- 4 is for electron, muon, tau and sterile. Sterile neutrino is a general beyond standard model particle.
  - Mixing Matrix: $U_{mix}= R_{34} R_{24} R_{14} R_{23} R_{13} R_{12}$ , [Barry:2011wb](https://inspirehep.net/literature/900636) eqn (2).
  - CP Violation phase is in the $R_{13}$.
### `neutrino_hierarchy`
- Neutrino hierarchy. 
- Acceptable values are only (1,-1). 1 is for normal hierarchy, -1 is for inverted hierarchy. 
- For normal hierarchy, $\delta m^{2}>0$ for two flavor, $\delta m_{21}^{2},\delta m_{32}^{2}, \delta m_{31}^{2} >0$ for three and four flavor.
- For inverted hierarchy $\delta m^{2}<0$ for two flavor, $\delta m_{21}^{2} > 0,\delta m_{32}^{2}, \delta m_{31}^{2} <0$ for three and four flavor.
### `theta12` 
- 1-2 Neutrino mixing angle. (Solar Mixing Angle)
- Unit is [Radian].
- Acceptable values are float numbers. Example: `0.59027`.
- Valid only if (`flavor_number`=2,3,4) and (`hamiltonian_oscillation`=1).
- For two flavor, `theta23`, `theta13`, `theta_14`, `theta_24` and `theta_34` are not considered.
### `theta23`
- 2-3 Neutrino mixing angle. (Atmospheric Mixing Angle for Normal Hierarchy)
- Unit is [Radian].
- Acceptable values are float numbers. Example: `0.84823`.
- Valid only if (`flavor_number`=3,4) , (`hamiltonian_oscillation`=1)
- For three flavor, `theta_14`, `theta_24` and `theta_34` are not considered.
### `theta13`
- 1-3 Neutrino mixing angle, (Atmospheric Mixing Angle for Inverted Hierarchy).
- Unit is [Radian].
- Acceptable values are float numbers. Example: `0.15010`.
- Valid only if (flavor_number=3,4) , (hamiltonian_oscillation=1)
- For three flavor, `theta_14`, `theta_24` and `theta_34` are not considered.
### `theta_14`
- 1 (or electorn)-4 Neutrino mixing angle.
- Unit is [Radian].
- Acceptable values are float numbers. Example: `0.1770`.
- Valid only if (`flavor_number`=4) , (`hamiltonian_oscillation`=1)
### `theta_24`
- 2  (or muon)-4 Neutrino mixing angle.
- Unit is [Radian].
- Acceptable values are float numbers. Example: `0.1645`.
- Valid only if (`flavor_number`=4) , (`hamiltonian_oscillation`=1)
### `theta_34`
- 3 (or tau)-4 Neutrino mixing angle.
- Unit is [Radian]
- Acceptable values are float numbers. Example: `0.039`.
- Valid only if (`flavor_number`=4) , (`hamiltonian_oscillation`=1)
### `CP_Phase`
- Dirac CP Violation Phase.
- Acceptable values are float numbers. Example: `3.85718`.
- Valid only if (`hamiltonian_oscillation`=1).
### `deltaM_square21`
- Mass Squared Differences.
- Unit is [MeV$^2$].
- Acceptable values are positive float numbers. Example: `7.53e-17`. Do not put negative sign for inverted hierarchy.
- Valid only if (`flavor_number`=2,3,4).
- For two flavor, `deltaM_square32` and `deltaM_square41` are not considered.
### `deltaM_square32`
- Mass Squared Differences.
- Unit is [MeV$^2$].
- Acceptable values are positive float numbers. Example: `2.56e-15`. Do not put negative sign for inverted hierarchy.
- Valid only if (`flavor_number`=3,4).
- For three flavor, `deltaM_square41` is not considered. 
### `deltaM_square41`
- Mass Squared Differences.
- Unit is [MeV$^2$].
- Acceptable values are positive float numbers. Example: `1.7e-12`. Do not put negative sign for inverted hierarchy.
- Valid only if (`flavor_number`=4).
### `neutrino_MagneticMoment`
- Neutrino Magnetic Moment.
- Unit is [Bohr Magneton, $\mu_{B}$].
- Acceptable values are positive float numbers. Example `5e-14`.
- Valid only if (`hamiltonian_electromagnetic`=1)
### `numberOf_energyMode`
- Number of Energy Mode
- Acceptable values are positive integers. Example `50`.
### `energy_initial`
- Lowest Energy Value of Neutrinos.
- Acceptable values are positive float numbers. Example `0.5`.
- Unit is [MeV].
### `energy_final`
- Highest Energy Value of Neutrinos.
- Acceptable values are positive float numbers. Example `50`.
- Unit is [MeV].
### `distance_initial`
- Initial Distance.
- Acceptable values are positive float numbers. Example `50.0`.
- Unit is [km].
- Simulation starts from this distance.
### `distance_final`
- Final Distance.
- Acceptable values are positive float numbers. Example `250.0`.
- Unit is [km].
- Simulation finishes at this distance.
### `neutrino_distributionParameter`
- Neutrino Energy Distribution.
- Acceptable values are (1,2,3,4,5,6,7). Example `1`.
- [1] : Fermi-Dirac
  - Unit is [1/MeV].
  - $\frac{1}{T^{3}}\frac{1}{F_{2}}\frac{E^{2}}{(1+\exp(E/T))}$ , [Duan:2006an](https://inspirehep.net/literature/720004)
- [2] : Pinched 
  - Unit is [$1/MeV$].
  - $(\frac{\alpha+1}{<E>})^{\alpha}\frac{E^{\alpha}}{\Gamma(\alpha+1)}\exp(-(\alpha+1)E/<E>)$ , [Keil:2002in](https://inspirehep.net/literature/591743)
- [3] : Electron Box
  - Unit is []. 
  - Electron neutrino and electron antineutrino have box spectrum.
- [4] : muBox
  - Unit is []. 
  - Muon neutrino and muon antineutrino have box spectrum.
- [5] : Fermi, Zero
  - Unit is [1/MeV].
  - Active neutrinos and antineutrino have Fermi-Dirac Dist.
  - Sterile neutrinos: Zero.
  - Valid only if (`flavor_number`=4).
- [6] : Only eBox
  - Unit is []. 
  - Only electron neutrino has box spectrum.
- [7] Only ebarBox
  - Unit is []. 
  - Only electron antineutrino has box spectrum.
### `luminosity_e`
- Luminosity of electron neutrinos.
- Unit is [erg/s].
- Valid only if (`flavor_number`=2,3,4, `neutrino_densityParameter`=1,2,5).
- Acceptable values are positive float numbers. Example `1e+52`.
### `luminosity_mu`
- Luminosity of muon neutrinos.
- Unit is [erg/s].
- Valid only if (`flavor_number`=2,3,4, `neutrino_densityParameter`=1,2,5)
- Acceptable values are positive float numbers. Example `1e+52`.
### `luminosity_tau`
- Luminosity of tau neutrinos.
- Unit is [erg/s].
- Valid only if (`flavor_number`=3,4, `neutrino_densityParameter`=1,2,5)
- Acceptable values are positive float numbers. Example `1e+52`.
### `luminosity_sterile`
- Luminosity of sterile neutrinos.
- Unit is [erg/s].
- Valid only if (`flavor_number`=4, `neutrino_densityParameter`=1,2,5)
- Acceptable values are positive float numbers. Example `1e+52`.
### `luminosity_eb`
- Luminosity of electron antineutrinos.
- Unit is [erg/s].
- Valid only if (`flavor_number`=2,3,4, `neutrino_densityParameter`=1,2,5)
- Acceptable values are positive float numbers. Example `1e+52`.
### `luminosity_mub`
- Luminosity of muon antineutrinos.
- Unit is [erg/s].
- Valid only if (`flavor_number`=2,3,4, `neutrino_densityParameter`=1,2,5)
- Acceptable values are positive float numbers. Example `1e+52`.
### `luminosity_taub`
- Luminosity of tau antineutrinos.
- Unit is [erg/s].
- Valid only if (`flavor_number`=3,4, `neutrino_densityParameter`=1,2,5)
- Acceptable values are positive float numbers. Example `1e+52`.
### `luminosity_sterileb`
- Luminosity of sterile antineutrinos.
- Unit is [erg/s].
- Valid only if (`flavor_number`=4, `neutrino_densityParameter`=1,2,5)
- Acceptable values are positive float numbers. Example `1e+52`.
### `temperature_e`
- Temperature of electron neutrinos.
- Unit is [MeV].
- Valid only if (`flavor_number`=2,3,4, `neutrino_densityParameter`=1,2,5)
- Acceptable values are positive float numbers. Example `3`.
### `temperature_mu`
- Temperature of muon neutrinos.
- Unit is [MeV].
- Valid only if (`flavor_number`=2,3,4, `neutrino_densityParameter`=1,2,5)
- Acceptable values are positive float numbers. Example `7`.
### `temperature_tau`
- Temperature of tau neutrinos.
- Unit is [MeV].
- Valid only if (`flavor_number`=3,4, `neutrino_densityParameter`=1,2,5)
- Acceptable values are positive float numbers. Example `7`.
### `temperature_sterile`
- Temperature of sterile neutrinos.
- Unit is [MeV].
- Valid only if (`flavor_number`=4, `neutrino_densityParameter`=1,2,5)
- Acceptable values are positive float numbers. Example `0`.
### `temperature_eb`
- Temperature of electron antineutrinos.
- Unit is [MeV].
- Valid only if (`flavor_number`=2,3,4, `neutrino_densityParameter`=1,2,5)
- Acceptable values are positive float numbers. Example `5`.
### `temperature_mub`
- Temperature of muon antineutrinos.
- Unit is [MeV].
- Valid only if (`flavor_number`=2,3,4, `neutrino_densityParameter`=1,2,5)
- Acceptable values are positive float numbers. Example `7`.
### `temperature_taub`
- Temperature of tau antineutrinos.
- Unit is [MeV].
- Valid only if (`flavor_number`=3,4, `neutrino_densityParameter`=1,2,5)
- Acceptable values are positive float numbers. Example `7`.
### `temperature_sterileb`
- Temperature of Sterile antineutrinos.
- Unit is [MeV].
- Valid only if (`flavor_number`=4, `neutrino_densityParameter`=1,2,5)
- Acceptable values are positive float numbers. Example `0`.
### `use_defaultMatterProfile`
- Use Default Matter Profile.
- Acceptable values are (0 or 1). Example `0`.
- Valid only if (`hamiltonian_matter`=1).
### `matterDensity_profile`
- Default Matter Profile.
- Valid only if (`use_defaultMatterProfile`=1).
- Acceptable values are (0, 1, 2). Example `1`.
- [1] : Constant Matter Profile
  - $\rho_{matter}(r)=$`matterDensity_initial`
- [2] : Exponential Matter Profile
  - $\rho_{matter}(r)=$`matterDensity_initial` $\exp(-r/$`matterDensity_exponentialDecay`$)$
- [3] : Polynomial Matter Profile
  - $\rho_{matter}(r)=$`matterDensity_initial` $\times[($ `distance_initial`$-r)/$`distance_final`$]$^ `matterDensity_polynomialDecay`
### `matterDensity_initial`
- Initial Matter Density.
- Unit is [g/cm^3].
- Valid only if (`use_defaultMatterProfile`=1).
- Acceptable values are positive float numbers. Example `1.0e+6`.
### `matterDensity_exponentialDecay`
- Exponential Decay Parameter for Matter Density.
- Unit is [km].
- Valid only if (`use_defaultMatterProfile`=1, `matterDensity_profile`=1)
- Acceptable values are positive float numbers. Example `200`.
### `matterDensity_polynomialDecay`
- Polynomial Decay Power Parameter for Matter Density.
- Unit is [].
- Valid only if (`use_defaultMatterProfile`=1, `matterDensity_profile`=2)
- Acceptable values are positive float numbers. Example `3.0`.
### `electronFraction_constant`
- Electron Fraction of Background Matter.
- Unit is [].
- Valid only if (`use_defaultMatterProfile`=1)
- Acceptable values are positive float numbers. Example `0.5`.
### `matterProfile_fileName`
- File Name of Matter Profile.
- Valid only if (`use_defaultMatterProfile`=0)
- The columns have to be **Distance [cm]**, **Electron Fraction** and **$\rho_{matter}$ [g/cm$^{3}$]**.
- It must be under "[backgroundProfiles](backgroundProfiles/)" folder. Example `distYeRhoMat.dat`.
### `use_defaultMagneticProfile`
- Default External Magnetic Field Profile.
- Acceptable values are (0 or 1). Example `0`.
- Valid only if (`hamiltonian_electromagnetic`=1).
### `magneticField_profile`
- Default External Magnetic Field Profile.
- Valid only if (`use_defaultMagneticProfile`=1).
- Acceptable values are (0, 1, 2). Example `1`.
- [1] : Constant Matter Profile
  - $\rho_{EM}(r)=$`magneticField_initial`.
- [2] : Exponential Matter Profile
  - $\rho_{EM}(r)=$`magneticField_initial` $\exp(-r/$`magneticField_exponentialDecay`$)$.
- [3] : Polynomial Matter Profile
  - $\rho_{EM}(r)=$`magneticField_initial` $\times[($ `magneticField_polynomialDecayDistance`$/r)/$`distance_final`$]$^ `magneticField_polynomialDecayPower`.
### `magneticField_initial`
- Initial Magnetic Field.
- Unit is [Gauss].
- Valid only if (`hamiltonian_electromagnetic`=1, `use_defaultMagneticProfile`=1)
- Acceptable values are positive float numbers. Example `1e+15`.
### `magneticField_exponentialDecay`
- Exponential Decay Parameter for Magnetic Field.
- Unit is [km].
- Valid only if (`use_defaultMagneticProfile`=1, `magneticField_profile`=1)
- Acceptable values are positive float numbers. Example `200`.
### `magneticField_polynomialDecayPower`
- Polynomial Decay Power Parameter for Magnetic Field.
- Unit is [km].
- Valid only if (`use_defaultMagneticProfile`=1, `magneticField_profile`=2)
- Acceptable values are positive float numbers. Example `2.0`.
### `magneticField_polynomialDecayDistance`
- Polynomial Decay Distance Parameter for Magnetic Field.
- Unit is [km].
- Valid only if (`use_defaultMagneticProfile`=1, `magneticField_profile`=2)
- Acceptable values are positive float numbers. Example `50.0`.
### `magneticField_fileName`
- File Name of Magnetic Field Profile.
- Valid only if (`use_defaultMagneticProfile`=0)
- The columns have to be **Distance [cm]** and **$\rho_{EM}$ [Gauss$]**.
- It must be under "[backgroundProfiles](backgroundProfiles/)" folder. Example `distRhoEM.dat`.
- Valid only if (`use_defaultMagneticProfile`=0).
### `couplingConstant_sterile`
- Coupling Constant For Beyond Standard Model.
- Unit is Fermi Constant, [$G_{F}$].
- It is defined in [Sigl:1993ctk](https://inspirehep.net/literature/33842) eq (5.2).
- Valid only if (`flavor_number`=4, `hamiltonian_selfInteraction_singleAngle`=1)
- Acceptable values are positive float numbers. Example `0.1`.

## Technical Parameters
### `holdData_Every`
- Holding Data Value For Every Distance At The End Of Simulation.
- Unit is [km].
- Acceptable values are positive float numbers. Example `1`.
### `holdIntermediateData`
- Saving Intermediate Data.
- Valid only if (0 or 1). Example `1`.
### `plotGraphs`
- Plot Graphs
- Valid only if (0 or 1). Example `0`.
### `plotGraphs_diagonalRhoAllEnergy_2distance`
- Plot and Save Diagonal Elements of Neutrino Density To Distance For All Energy.
- Valid only if (0 or 1) and (`plotGraphs`= 1). Example `1`.
### `plotGraphs_diagonalRhoFinal_2energy`
- Plot and Save Diagonal Elements of Neutrino Density To Energy For Final And Initial Distance.
- Valid only if (0 or 1) and (`plotGraphs`= 1). Example `1`.
### `plotGraphs_hamiltonianAllEnergy_2distance`
- Plot and Save Diagonal Elements of Hamiltonian (For EM Interaction, One of the Non-Diagonal Element) To Distance For All Energy.
- Valid only if (0 or 1) and (`plotGraphs`= 1). Example `1`.
### `plot_savingFormat`
- Saving Format Of Graphs
- Valid only if (png,rgba,eps,svgz,pgf,ps,raw,svg,pdf), (`plotGraphs`= 1). Example `pdf`.
### `output_distanceHamiltonianAllE`
- Saving Hamiltonians in Separate File
- Valid only if (0 or 1). Example `0`.
### `output_distance_eigenValuesAllE`
- Saving Hamiltonians Eigenvalues in Separate File
- Valid only if (0 or 1). Example `0`.
### `output_humanReadable`
- Saving Neutrino Density in Human Readable Form
- Valid only if (0 or 1). Example `0`.
### `initialValueProblem_solverMethod`
- Initial Value Problem Solver Method (UNTESTED PARAMETER. USE `LSODA`)
- Valid only if (RK45, RK23, BDF, LSODA). Example `LSODA`.
### `tolerance_relativeError`
- Relative Error For Initial Value Problem Solver
- Acceptable values are positive float numbers. Example `0.0`.
- `0.0` is default value of solver.
### `tolerance_absoluteError`
- Absolute Error For Initial Value Problem Solver
- Acceptable values are positive float numbers. Example `0.0`.
- `0.0` is default value of solver.

# References
- [NuFit](http://www.nu-fit.org/)
- [PDG](https://pdg.lbl.gov/)