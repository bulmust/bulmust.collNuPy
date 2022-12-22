import numpy as np
import matplotlib.pyplot as plt
import matplotlib.pylab as pylab
from matplotlib import rc
import os as osCommand # File operations

# ============================
# Colors For Printing
class tLog:
    INFO  = '\033[94m'+'[INFO]    '+'\033[0m'
    OK      = '\033[92m'+'[OK]        '+'\033[0m'
    WARNING = '\033[93m'+'[WARNING]   '+'\033[0m'
    ERROR   = '\033[91m'+'[ERROR]     '+'\033[0m'
    EXIT    = '\033[91m' + '[EXIT]      plotgraphs.py' + '\033[0m'
# ============================

class plotMain:
    # ============================
    def __init__(self, SIMULATION_PATH):
        # Modify Plot Parameters
        rc('figure', figsize=(25.6, 14.4))
        params = {'legend.fontsize':'x-large',
                    'axes.labelsize' :'xx-large',
                    'axes.titlesize' :'xx-large',
                    'xtick.labelsize':'x-large',
                    'ytick.labelsize':'x-large'}
        pylab.rcParams.update(params)

        # ============================
        # Paths
        self.SIMULATION_PATH = SIMULATION_PATH
        self.DATA_PATH = SIMULATION_PATH+'data/'
        self.FIGURE_PATH = SIMULATION_PATH+'figures/'
        # ============================
        # rhoFlavAll_npz load
        rhoFlavAll_npz = np.load(self.DATA_PATH+'rhoFlavAll.npz','r',allow_pickle=True)
        # Assign Data
        self.totFlav= rhoFlavAll_npz['totFlav']
        self.distAll_km= rhoFlavAll_npz['distAll_km']
        self.E_MeV= rhoFlavAll_npz['E_MeV']; self.Emod= len(self.E_MeV)
        self.rhoFlavAll = rhoFlavAll_npz['rhoFlavAll']
        self.dim_rho_2totFlav_Bool = rhoFlavAll_npz['dim_rho_2totFlav_Bool']
        # ============================
        # parametersDic_npz load
        parametersDic_npz = np.load(self.DATA_PATH+\
            'parametersDic.npz','r',allow_pickle=True)
        # Assign Data
        #print(parametersDic_npz.files)
        
        # physicalParameters.dat
        self.TraceTerm = parametersDic_npz['TraceTerm']
        self.TraceTermb= parametersDic_npz['TraceTermb']

        # technicalParameters.dat
        self.pltFormat = str(parametersDic_npz['plot_savingFormat'])
        self.output_distanceHamiltonianAllE=\
            parametersDic_npz['output_distanceHamiltonianAllE']
        self.plotGraphs_hamiltonianAllEnergy_2distance=\
            parametersDic_npz['plotGraphs_hamiltonianAllEnergy_2distance']
        
        # ============================
        # Figures Folder
        if not osCommand.path.exists(self.FIGURE_PATH):
            osCommand.makedirs(self.FIGURE_PATH)
            print(tLog.OK+"figures Folder Created.")
        else:
            print(tLog.INFO+"figures Folder Already Exists.")
    # ============================

    # ============================
    def ldistDiagNu_Flav(self, totFlav):
        nu_e   = '$\\nu_{e}$';     nub_e   = '$\\overline{\\nu}_{e}$'
        nu_mu  = '$\\nu_{\\mu}$';  nub_mu  = '$\\overline{\\nu}_{\\mu}$'
        nu_x   = '$\\nu_{x}$';     nub_x   = '$\\overline{\\nu}_{x}$'
        nu_tau = '$\\nu_{\\tau}$'; nub_tau = '$\\overline{\\nu}_{\\tau}$'
        nu_ster= '$\\nu_{s}$';     nub_ster= '$\\overline{\\nu}_{s}$'
    
        if totFlav == 2:
            labelsArrayFlav= np.array([nu_e, nu_x, nub_e, nub_x])
        if totFlav == 3:
            labelsArrayFlav= np.array([nu_e, nu_mu, nu_tau, nub_e, nub_mu, nub_tau])
        if totFlav == 4:
            labelsArrayFlav= np.array(\
                [nu_e, nu_mu, nu_tau, nu_ster, nub_e, nub_mu, nub_tau, nub_ster])

        return labelsArrayFlav
    def cdistDiag_Flav(self, totFlav):
        nu_e   = 'r';  nub_e   = 'b'
        nu_mu  = 'k';  nub_mu  = 'g'
        nu_x   = 'k';  nub_x   = 'g'
        nu_tau = 'm';  nub_tau = 'purple'
        nu_ster= 'orange'; nub_ster= 'brown'
    
        if totFlav == 2:
            colorsArrayFlav= np.array([nu_e, nu_x, nub_e, nub_x])
        if totFlav == 3:
            colorsArrayFlav= np.array([nu_e, nu_mu, nu_tau, nub_e, nub_mu, nub_tau])
        if totFlav == 4:
            colorsArrayFlav= np.array(\
                [nu_e, nu_mu, nu_tau, nu_ster, nub_e, nub_mu, nub_tau, nub_ster])

        return colorsArrayFlav
    def lwdistDiag(self, totFlav):
        nu_e   = 4.5; nub_e   = 2.0
        nu_mu  = 3.5; nub_mu  = 1.5
        nu_x   = 3.5; nub_x   = 1.5
        nu_tau = 2.5; nub_tau = 0.5
        nu_ster= 1.5; nub_ster= 0.5
    
        if totFlav == 2:
            linewidthsArrayFlav= np.array([nu_e, nu_x, nub_e, nub_x])
        if totFlav == 3:
            linewidthsArrayFlav= np.array([nu_e, nu_mu, nu_tau, nub_e, nub_mu, nub_tau])
        if totFlav == 4:
            linewidthsArrayFlav= np.array(\
                [nu_e, nu_mu, nu_tau, nu_ster, nub_e, nub_mu, nub_tau, nub_ster])

        return linewidthsArrayFlav
    # ============================

    # ============================
    def distDiag(self):
        # Paths
        FIG_DISTDIAG_PATH = self.FIGURE_PATH+'distDiag/'
        if not osCommand.path.exists(FIG_DISTDIAG_PATH):
            osCommand.makedirs(FIG_DISTDIAG_PATH)
            print(tLog.OK+"figures/distDiag/ Folder Created.")
        else:
            print(tLog.INFO+"figures/distDiag/ Folder Already Exists.")

        # Create Y Axises
        yAxises = np.zeros((len(self.distAll_km), self.Emod, self.totFlav*2))
        
        # Assign Diagonal Elements of rhoFlavAll to yAxises
        for it_E in range(self.Emod):
            for it_flav in range(self.totFlav* 2):
                if self.dim_rho_2totFlav_Bool:
                    yAxises[:, it_E, it_flav]= np.real(self.rhoFlavAll\
                        [:, it_E, it_flav, it_flav])
                else:
                    if it_flav < self.totFlav:
                        yAxises[:, it_E, it_flav]= np.real(self.rhoFlavAll\
                            [:, 0, it_E, it_flav, it_flav])
                    else:
                        yAxises[:, it_E, it_flav]= np.real(self.rhoFlavAll\
                            [:, 1, it_E, it_flav% self.totFlav, it_flav% self.totFlav])
        
        # Plot
        labelsArrayFlav= self.ldistDiagNu_Flav(self.totFlav)
        colorsArrayFlav= self.cdistDiag_Flav(self.totFlav)
        linewidthsArrayFlav= self.lwdistDiag(self.totFlav)
        
        for it_E in range(self.Emod):
            for it_plt in range(self.totFlav*2):
                plt.plot(self.distAll_km, yAxises[:, it_E, it_plt]\
                    , label=labelsArrayFlav[it_plt]\
                    , color=colorsArrayFlav[it_plt]\
                    , linewidth=linewidthsArrayFlav[it_plt])
            plt.legend(loc = 'upper right')

            # Labels
            plt.ylabel('$\\rho_{\\alpha\\alpha}$')
            plt.xlabel('Distance [km]')
            # Save The Graph
            plt.savefig(FIG_DISTDIAG_PATH+'distDiag_E'\
                +str(int(np.round(self.E_MeV[it_E], decimals=0)))+'MeV.'+self.pltFormat)

            # Clear Figure
            plt.clf()
            plt.close('all')
    # ============================
    print(tLog.OK+"Diagonal elements of density matrix '\
        'to distance graph(s) for each energy are created.")
    # ============================
    def energyDiag(self):
        # Paths
        FIG_DISTDIAG_PATH = self.FIGURE_PATH+'energyDiag/'
        if not osCommand.path.exists(FIG_DISTDIAG_PATH):
            osCommand.makedirs(FIG_DISTDIAG_PATH)
            print(tLog.OK+"figures/energyDiag/ Folder Created.")
        else:
            print(tLog.INFO+"figures/energyDiag/ Folder Already Exists.")

        # Create Y Axises
        yAxisesBeg = np.zeros((self.Emod, self.totFlav*2))
        yAxisesFin = np.zeros((self.Emod, self.totFlav*2))
        
        # Assign Diagonal Elements of rhoFlavAll to yAxises
        for it_E in range(self.Emod):
            for it_flav in range(self.totFlav* 2):
                if self.dim_rho_2totFlav_Bool:
                    # TraceTerm = TraceTermb
                    yAxisesBeg[it_E, it_flav]= np.real(self.rhoFlavAll\
                        [ 0, it_E, it_flav, it_flav])* np.real(self.TraceTerm[it_E])
                    yAxisesFin[it_E, it_flav]= np.real(self.rhoFlavAll\
                        [-1, it_E, it_flav, it_flav])* np.real(self.TraceTerm[it_E])
                else:
                    if it_flav < self.totFlav:
                        # Neutrino
                        yAxisesBeg[it_E, it_flav]= np.real(self.rhoFlavAll\
                            [ 0, 0, it_E, it_flav, it_flav])\
                            * np.real(self.TraceTerm[it_E])
                        yAxisesFin[it_E, it_flav]= np.real(self.rhoFlavAll\
                            [-1, 0, it_E, it_flav, it_flav])\
                            * np.real(self.TraceTerm[it_E])
                    else:
                        # Antineutrino
                        yAxisesBeg[it_E, it_flav]= np.real(self.rhoFlavAll\
                            [ 0, 1, it_E, it_flav% self.totFlav\
                                , it_flav% self.totFlav])\
                            * np.real(self.TraceTermb[it_E])
                        yAxisesFin[it_E, it_flav]= np.real(self.rhoFlavAll\
                            [-1, 1, it_E, it_flav% self.totFlav\
                                , it_flav% self.totFlav])\
                            * np.real(self.TraceTermb[it_E])
    
        # Plot
        labelsArrayFlav= self.ldistDiagNu_Flav(self.totFlav)
        colorsArrayFlav= self.cdistDiag_Flav(self.totFlav)
        linewidthsArrayFlav= self.lwdistDiag(self.totFlav)
        
        for it_plt in range(self.totFlav* 2):
            plt.plot(self.E_MeV, yAxisesFin[:, it_plt]\
                , label=labelsArrayFlav[it_plt]\
                , color=colorsArrayFlav[it_plt]\
                , linewidth=linewidthsArrayFlav[it_plt])
            plt.plot(self.E_MeV, yAxisesBeg[:, it_plt]\
                , label=labelsArrayFlav[it_plt]\
                , color=colorsArrayFlav[it_plt]\
                , linewidth=linewidthsArrayFlav[it_plt]\
                , linestyle='dashed')
        plt.legend(loc = 'upper right')

        # Labels
        plt.ylabel('$\\rho_{\\alpha\\alpha}$')
        plt.xlabel('Distance [km]')
        # Save The Graph
        plt.savefig(FIG_DISTDIAG_PATH+'energyDiag.'+self.pltFormat)

        # Clear Figure
        plt.clf()
        plt.close('all')
        print(tLog.OK+"Diagonal elements of final density'\
            'matrix to energy graph(s) are created.")
    def distHamiltDiag(self):
        # ============================
        # hamiltonians.npz load
        try:
            hamiltonians_npz = np.load(self.DATA_PATH+'hamiltonians.npz'\
                ,'r',allow_pickle=True)
        except FileNotFoundError:
            print(tLog.ERROR+'File hamiltonians.npz not found.')
            print(tLog.ERROR+'Run hamiltonians.npz creation script.')
            print(tLog.ERROR+'Skip plotting plotGraphs_hamiltonianAllEnergy_2distance.')
            # Return to main
            return
        # ============================
        # Paths
        FIG_DISTDIAG_PATH = self.FIGURE_PATH+'distHamiltDiag/'
        if not osCommand.path.exists(FIG_DISTDIAG_PATH):
            osCommand.makedirs(FIG_DISTDIAG_PATH)
            print(tLog.OK+"figures/distHamiltDiag/ Folder Created.")
        else:
            print(tLog.INFO+"figures/distHamiltDiag/ Folder Already Exists.")
        # ============================

        # ============================
        # Assign Hamiltonians
        # Hamiltonian Bool
        hamBoolAll_OscMatEMSelfSA= hamiltonians_npz['hamBoolAll_OscMatEMSelfSA']
        # ============================
        # Oscillation Hamiltonian
        if hamBoolAll_OscMatEMSelfSA[0]:
            hamOsc_km_1= hamiltonians_npz['hamOsc_km_1']
            
            # ============================
            # Plot
            for it_E in range(self.Emod):
                plt.plot(self.distAll_km\
                    , np.ones(len(self.distAll_km))\
                        *(np.abs(np.real(hamOsc_km_1[it_E, 0, 0])))\
                    , color='r'\
                    , alpha=0.4)
            # Labels
            plt.ylabel('$H_{osc}(1,1)$ [km$^{-1}$]')
            plt.xlabel('Distance [km]')
            # Save The Graph
            plt.savefig(FIG_DISTDIAG_PATH+'distHamiltDiag_Osc.'+self.pltFormat)
            # Clear Figure
            plt.clf()
            plt.close('all')
            print(tLog.OK+"(0,0) elements of oscillation Hamiltonian '\
                '(km^-1) to distance graph(s) for each energy are created.")
        # ============================
        # Matter Hamiltonian
        if hamBoolAll_OscMatEMSelfSA[1]:
            hamMat_km_1= hamiltonians_npz['hamMat_km_1']
            # Plot
            plt.plot(self.distAll_km, (np.abs(np.real(hamMat_km_1[:, 0, 0])))\
                , color='k', label='$H_{mat}(1,1)$')
            plt.plot(self.distAll_km, (np.abs(np.real(hamMat_km_1[:, 1, 1])))\
                , color='r', label='$H_{mat}(2,2)$')
            if self.totFlav == 3:
                plt.plot(self.distAll_km, (np.abs(np.real(hamMat_km_1[:, 2, 2])))\
                    , color='b', label='$H_{mat}(3,3)$')
            plt.ylabel('$H_{mat}$ [km$^{-1}$]')
            # Labels and legend
            plt.xlabel('Distance [km]')
            plt.legend(loc = 'upper right')
            # Save The Graph
            plt.savefig(FIG_DISTDIAG_PATH+'distHamiltDiag_Mat.'+self.pltFormat)
            # Clear Figure
            plt.clf()
            plt.close('all')
            print(tLog.OK+"(0,0), (1,1) elements of matter Hamiltonian '\
                '(km^-1) to distance graph(s) are created.")
        # ============================
        # EM Hamiltonian
        if hamBoolAll_OscMatEMSelfSA[2]:
            hamEM_km_1= hamiltonians_npz['hamEM_km_1']
            # Plot
            plt.plot(self.distAll_km, (np.abs(np.real(hamEM_km_1[:, 0, -1])))\
                , color='k', label='$\\mu_{\\nu}B$')
            plt.ylabel('$\\mu_{\\nu}B$ [km$^{-1}$]')
            # Labels and legend
            plt.xlabel('Distance [km]')
            plt.legend(loc = 'upper right')
            # Save The Graph
            plt.savefig(FIG_DISTDIAG_PATH+'distHamiltDiag_EM.'+self.pltFormat)
            # Clear Figure
            plt.clf()
            plt.close('all')
            print(tLog.OK+"(0,-1) elements of electromagnetic Hamiltonian '\
                '(km^-1) or muB to distance graph(s) are created.")
        # ============================
        # SelfSA Hamiltonian
        if hamBoolAll_OscMatEMSelfSA[3]:
            hamSelfSA_km_1= hamiltonians_npz['hamSelfSA_km_1']
            # Plot
            plt.plot(self.distAll_km, (np.abs(np.real(hamSelfSA_km_1[:, 0, 0])))\
                , color='k', label='$H_{self}(1,1)$')
            plt.ylabel('$H_{self}(1,1)$ [km$^{-1}$]')
            # Labels and legend
            plt.xlabel('Distance [km]')
            plt.legend(loc = 'upper right')
            # Save The Graph
            plt.savefig(FIG_DISTDIAG_PATH+'distHamiltDiag_SelfSA.'+self.pltFormat)
            # Clear Figure
            plt.clf()
            plt.close('all')
            print(tLog.OK+"(0,0) elements of self interaction Hamiltonian '\
                '(km^-1) to distance graph(s) are created.")
    # ============================
