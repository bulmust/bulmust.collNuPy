"""
ID      : unitConv.py
Author  : Taygun Bulmus
E-Mail  : bulmust@gmail.com

Notes:
---------------
    Data are taken from https://physics.nist.gov/cuu/Constants/index.html
"""
from scipy import constants as scConst

# MeV to km^-1
eV_TO_m_1=(scConst.eV/(scConst.hbar*scConst.c)) # 5067730.716156395
eV_TO_km_1=eV_TO_m_1*1e+3  # 5067730716.156395
MeV_TO_km_1=eV_TO_m_1*1e+9 # 5067730716156395.0
km_TO_MeV_1=MeV_TO_km_1    # 5067730716156395.0 
MeV_1_TO_km=1/MeV_TO_km_1  # 1.973269804593025e-16 
km_1_TO_MeV=1/MeV_TO_km_1  # 1.973269804593025e-16

# s to km
# [m s^-1] * scConst.c = 1
# [s] = scConst.c [m]
s_TO_km=scConst.c*1e-3
# -------------------------------------------------------

# Fermi Coupling Constant [GeV^-2]
G_f_GeV_2=scConst.value('Fermi coupling constant') # 1.1663787e-05
G_f_MeV_2=G_f_GeV_2*1e-6 # 1.1663786999999999e-11
G_f_km2=G_f_MeV_2/(MeV_TO_km_1**2) #4.541637983249213e-43
# -------------------------------------------------------

# 12 g ^{12}C has N_a atom, where N_a is Avogadro Number
#  1 g ^{12}C has N_a baryon
#  1 baryon = (1/N_a) g
# (It does not equal to 2/(m_p+m_n)) due to the bounding energy)
# scConst.Avogadro = 6.022140857e+23
# (2/((scConst.proton_mass+scConst.neutron_mass)*1e+3)) = 5.974519804012485e+23
avogadroNum_g_1=scConst.Avogadro; # 6.02214076e+23
#nucleon_Mass_km_1=(938.95*(10**6))*(1/197)*(10**12)

# -------------------------------------------------------

# erg/s to MeV/km
# 1 [erg] = scConst.erg [J]
# 1 [MeV] = scConst.eV*1e+6 [J]
# 1 [erg] = scConst.erg/scConst.eV*1e+6 [MeV]
# 1 [erg/s] = scConst.erg/scConst.eV*1e+6 [MeV s^-1]
# 1 [erg/s] = (scConst.erg/scConst.eV*1e+6)/s__km [MeV km^-1]
erg_s_TO_MeV_km = (scConst.erg/ (scConst.eV* 1e+6))/ s_TO_km; # 2.0819433270935597

# Bohr Magneton [eV/Tesla]
muB_TO_eVTesla_1 = scConst.value('Bohr magneton in eV/T') # 5.788381806e-05
# 1 Tesla = 1e+4 Gauss
muB_TO_eVGauss_1 = muB_TO_eVTesla_1 * 1e-4 # 5.788381806e-09
# 1 eV = 5067730716.156395 km^-1
muB_TO_Gauss__km = muB_TO_eVGauss_1 * eV_TO_km_1 # 29.333960275107025

# -------------------------------------------------------
# Speed Of Light [m s-1]
# scConst.c
# -------------------------------------------------------
# -------------------------------------------------------
# Planck Constant [eV s] & Planck Constant Over Pi [eV s]
# scConst.h  & scConst.hbar
# -------------------------------------------------------

# -------------------------------------------------------
# Proton Mass [kg]
# scConst.proton_mass
# print(scConst.unit(u'proton mass'))
# -------------------------------------------------------

# -------------------------------------------------------
# Neutron Mass [kg]
# scConst.neutron_mass
# print(scConst.unit(u'neutron mass'))
# -------------------------------------------------------
