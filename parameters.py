#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 18 08:20:58 2018

The mathematical model of the photosynthetic electron transport chain defines methods to calculate reaction rates
and set of ten differential equations based on the model published by Ebenhoeh et al. in 2014

Copyright (C) 2014-2018  Anna Matuszyńska, Oliver Ebenhöh

This program is free software: you can redistribute and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with this program (license.txt).  If not, see <http://www.gnu.org/licenses/>.
"""
__author__ = "Anna Matuszyńska"
__copyright__ = "Copyright 2018, Heinrich-Heine University Dusseldorf"
__credits__ = ["Anna Matuszynska", "Oliver Ebenhoeh"]
__maintainer__ = "Anna Matuszynska"
__email__ = "Anna.Matuszynska@uni-duesseldorf.de"
__status__ = "Development"

from dotmap import DotMap

class Parameters(DotMap):
    def __init__(self):
        super().__init__()  

        self.convf = 3.2*10e-3

        # pool sizes
        self.PSIItot = 2.5 # [mmol/molChl] total concentration of PSII
        self.PSItot = 2.5
        self.PQtot = 17.5 # [mmol/molChl]
        self.PCtot = 4. # Bohme1987 but other sources give different values - seems to depend greatly on organism and conditions
        self.Fdtot = 5. # Bohme1987
        self.Ctot = 2.5 #source unclear (Schoettler says 0.4...?, but plausible to assume that complexes (PSII,PSI,b6f) have approx. same abundance)
        self.NADPtot = 0.8 # estimate from ~ 0.8 mM, Heineke1991
        self.APtot = 2.55 # [mmol/molChl] Bionumbers ~2.55mM (=81mmol/molChl) (FIXME: Soma had 50)
        self.PsbStot = 1 # relative pool of PsbS
        self.Xtot = 1 # relative pool of carotenoids (V+A+Z)

        # parameters associated with photosystem II
        self.kH = 5e9
        self.kH0 = 5e8 # base quenchingself. after calculation with Giovanni
        self.kF = 6.25e8 # 6.25e7 fluorescence 16ns
        self.k1 = 5e9 # excitation of Pheo / charge separation 200ps
        self.k1rev = 1e10
        self.k2 = 5e9 # original 5e9 (charge separation limiting step ~ 200ps) - made this faster for higher Fs fluorescence
        self.kdeg = 100    # rate of PSII damage corresponds to p.k2 / .5e8
        self.krep = 5.55e-4 # rate of repair fo PSII

        # parameters associated with photosystem I
        self.kStt7 = 0.0035 # [s-1] fitted to the FM dynamics
        self.kPph1 = 0.0013 # [s-1] fitted to the FM dynamics
        self.KM_ST = 0.2 # Switch point (half-activity of Stt7) for 20% PQ oxidised (80% reduced)
        self.n_ST = 2. # Hill coefficient of 4 -> 1/(2.5^4)~1/40 activity at PQox=PQred
        self.staticAntI = 0.37     # corresponds to PSI - LHCI supercomplex, when chlorophyll decreases more relative fixed antennae
        self.staticAntII = 0.1     # corresponds to PSII core
        self.prob_attach = 1            # probability of antena attaching to PSI


         # ATP and NADPH parameters
        self.kActATPase = 0.05  # on 14.09 increased from 0.01 to saturate between 1-2 min, not 10
                                # paramter relating the rate constant of activation of the ATPase in the light
        self.kDeactATPase = 0.002   # paramter relating the deactivation of the ATPase at night
        self.kATPsynth = 20.    # taken from MATLAB
        self.kATPcons = 10.     # taken from MATLAB
        self.kATPimport = 0.    # TODO possibility for ATP import at night - NOT YET IMPLEMENTED!
        self.ATPcyt = 0.5       # only relative levels are relevant (normalised to 1) to set equilibrium
        self.Pi_mol = 0.01
        self.DeltaG0_ATP = 30.6 # 30.6kJ/mol / RT
        self.HPR = 14./3.  #Vollmar et al. 2009 (after Zhu et al. 2013)
        self.kNADPHimport = 0. # TODO possibility for NADPH import - NOT YET IMPLEMENTED!
        self.kNADPHcons = 15. # taken from MATLAB
        self.NADPHcyt = 0.5 # only relatice levels

        # global conversion factor of PFD to excitation rate
        #self.cPFD = 4. # [m^2/mmol PSII]

        # pH and protons
        self.pHstroma = 7.8
        self.kLeak = 10#0.010 # [1/s] leakage rate -- inconsistency with Kathrine
        self.bH = 100. # proton buffer: ratio total / free protons

        # rate constants
        self.kPQred = 250. # [1/(s*(mmol/molChl))]
        self.kCytb6f = 2.5 # a rough estimate: transfer PQ->cytf should be ~10ms
        self.kPTOX = .01 # ~ 5 electrons / seconds. This gives a bit more (~20)
        self.kPCox = 2500. # a rough estimate: half life of PC->P700 should be ~0.2ms
        self.kFdred = 2.5e5 # a rough estimate: half life of PC->P700 should be ~2micro-s
        self.kcatFNR = 500. # Carrillo2003 (kcat~500 1/s)
        self.kcyc = 1.

        self.O2ext = 8. # corresponds to 250 microM corr. to 20%
        self.kNDH = .002 # re-introduce e- into PQ pool. Only positive for anaerobic (reducing) condition
        self.kNh = 0.05
        self.kNr = 0.004
        self.NPQsw = 5.8
        self.nH = 5.

        self.EFNR = 3. # Bohme1987
        self.KM_FNR_F = 1.56 # corresponds to 0.05 mM (Aliverti1990)
        self.KM_FNR_N = 0.22 # corresponds to 0.007 mM (Shin1971 Aliverti2004)

        # quencher fitted parameters
        self.gamma0 = 0.1          # slow quenching of Vx present despite lack of protonation
        self.gamma1 = 0.25         # fast quenching present due to the protonation
        self.gamma2 = 0.6          # slow quenching of Zx present despite lack of protonation
        self.gamma3 = 0.15         # fastest possible quenching

        # non-photochemical quenching PROTONATION
        self.kDeprotonation = 0.0096
        self.kProtonationL = 0.0096
        self.kphSatLHC = 5.8
        self.nH = 5.
        self.NPQsw = 5.8

        # non-photochemical quenching XANTOPHYLLS
        self.kDeepoxV = 0.0024
        self.kEpoxZ = 0.00024      # 6.e-4        # converted to [1/s]
        self.kphSat = 5.8          # [-] half-saturation pH value for activity de-epoxidase highest activity at ~pH 5.8
        self.kHillX = 5.     # [-] hill-coefficient for activity of de-epoxidase
        self.kHillL = 3.     # [-] hill-coefficient for activity of de-epoxidase
        self.kZSat = 0.12          # [-] half-saturation constant (relative conc. of Z) for quenching of Z

        # standard redox potentials (at pH=0) in V
        self.E0_QA = -0.140
        self.E0_PQ = 0.354
        self.E0_cytf = 0.350
        self.E0_PC = 0.380
        self.E0_P700 = 0.480
        self.E0_FA = -0.550
        self.E0_Fd = -0.430
        self.E0_NADP = -0.113

        # physical constants
        self.F = 96.485 # Faraday constant
        self.R = 8.3e-3 # universal gas constant
        self.T = 298. # Temperature in K - for now assumed to be constant at 25 C

        # CBB cycle associated parameter set according to Pettersson and Pettersson 1988
        self.CN = 0.5
        self.CO2 = 0.2
        self.Cp = 15+2.05#15.0
        self.Ca = 0.5
        self.pHmedium = 7.6
        self.pHstroma = 7.9
        self.Pext = 0.5

        #Vmaxes of Calvin cycle enzymes
        self.V1 = 0.34*8
        self.V6 = 0.2*8
        self.V9 = 0.04*8
        self.V13 = 0.9999*8
        self.Vst = 0.04*8
        self.Vx = 0.25*8
        
        #equilibrium constants of calvin cycle enzymes
        self.q2 = 3.1 * (10.0 ** (-4.0))
        self.q3 = 1.6 * (10.0**7.0)
        self.q4 = 22.0
        self.q5 = (7.1)
        self.q7 = 0.084
        self.q8 = (13.0)
        self.q10 = 0.85
        self.q11 = 0.4
        self.q12 = 0.67
        self.q14 = 2.3
        self.q15 = 0.058
        
        #michaelis constants of calvin cycle enzymes
        self.Km1 = 0.02
        self.KmCO2 = 0.0107 #millimol laut witzel
        self.Km6 = 0.03
        self.Km9 = 0.013
        self.Km131 = 0.05
        self.Km132 = 0.05
        self.Km161 = 0.014
        self.Km162 = 0.3
        self.Kmst1 = 0.08
        self.Kmst2 = 0.08
        self.Kmnadph = 0.19#ausgerechneter wert (ideal wert)
        self.Kpga = 0.25
        self.Kgap = 0.075
        self.Kdhap = 0.077
        self.Kpi = 0.63
        self.Kpxt = 0.74
        self.Ki11 = 0.04
        self.Ki12 = 0.04
        self.Ki13 = 0.075
        self.Ki14 = 0.9
        self.Ki15 = 0.07
        self.Ki61 = 0.7
        self.Ki62 = 12.0
        self.Ki9 = 12.0
        self.Ki131 = 2.0
        self.Ki132 = 0.7
        self.Ki133 = 4.0
        self.Ki134 = 2.5
        self.Ki135 = 0.4
        self.Kist = 10.0
        self.Kast1 = 0.1
        self.Kast2 = 0.02
        self.Kast3 = 0.02
        
        self.k = 10.0**8.0*8
        self.oxPPP=0.

        # Equilibrium constants


        # light
        self.PFD = 100.
        self.Ton = 360
        self.Toff = 1800
        self.dT=120

        self.ox = True

