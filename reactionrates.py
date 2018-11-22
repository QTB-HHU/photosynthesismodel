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

import numpy as np
import math
import warnings
# ================================================= #
# Photosynthetic Electron Transport Reaction Rates  #
# ================================================= #

class Reactions:

    def fluorescence(self, p, P, Pred, LHC, Psbs, Viola, **kwargs):
        cs = self.crossSection(p, LHC)
        B = self.ps2states(p, P, Pred, LHC, Psbs, Viola, **kwargs)
        fluo = cs * ((p.kF / (p.kF + p.kH0 + p.k2)) * B[0] + (p.kF / (p.kF + p.kH0)) * B[2])
        return fluo

  
    # ====================================================================== #
    # Composed parameters #
    # ====================================================================== #
    def RT(self,p):
        return p.R*p.T

    def dG_pH(self,p):
        return np.log(10)*self.RT(p)

    def Hstroma(self,p):
        return 3.2e4*10**(-p.pHstroma) 

    def kProtonation(self,p):
        return 4e-3 / self.Hstroma(p)

    def pH(self,x):
        return (-np.log(x*(2.5e-4))/np.log(10))
    
    def pHstroma(self,x):
        return (-np.log(x*(3.2e-5))/np.log(10))    
    
    def pHinv(self,x):
        return (4e3*10**(-x))

    def Keq_PQred(self,p):
        DG1 = -p.E0_QA * p.F
        DG2 = -2 * p.E0_PQ * p.F
        DG = -2 * DG1 + DG2 + 2 * p.pHstroma * self.dG_pH(p)
        K = np.exp(-DG/self.RT(p))
        return K

    def Keq_cyc(self,p):
        DG1 = -p.E0_Fd * p.F
        DG2 = -2 * p.E0_PQ * p.F
        DG = -2 * DG1 + DG2 + 2 * self.dG_pH(p) * p.pHstroma
        K = np.exp(-DG/self.RT(p))
        return K

    def Keq_cytfPC(self,p):
        DG1 = -p.E0_cytf * p.F
        DG2 = -p.E0_PC * p.F
        DG = -DG1 + DG2
        K = np.exp(-DG/self.RT(p))
        return K

    def Keq_FAFd(self,p):
        DG1 = -p.E0_FA * p.F
        DG2 = -p.E0_Fd * p.F
        DG = -DG1 + DG2
        K = np.exp(-DG/self.RT(p))
        return K

    def Keq_PCP700(self,p):
        DG1 = -p.E0_PC * p.F
        DG2 = -p.E0_P700 * p.F
        DG = -DG1 + DG2
        K = np.exp(-DG/self.RT(p))
        return K

    def Keq_NDH(self,p):
        DG1 = -2 * p.E0_NADP * p.F
        DG2 = -2 * p.E0_PQ * p.F
        DG = -DG1 + DG2 + self.dG_pH(p) * p.pHstroma
        K = np.exp(-DG/self.RT(p))
        return K

    def Keq_FNR(self,p):
        DG1 = -p.E0_Fd * p.F
        DG2 = -2 * p.E0_NADP * p.F
        DG = -2 * DG1 + DG2 + self.dG_pH(p) * p.pHstroma
        K = np.exp(-DG/self.RT(p))
        return K

    def Keq_ATP(self,p, pH, Pi):
        DG = p.DeltaG0_ATP - self.dG_pH(p) * p.HPR * (p.pHstroma - pH)
        Keq = p.Pi_mol * np.exp(-DG/self.RT(p))
        #Keq = np.exp(-DG/self.RT(p))
        return Keq

    def Keq_cytb6f(self,p, pH):
        DG1 = -2 * p.F * p.E0_PQ
        DG2 = -p.F * p.E0_PC
        DG = - (DG1 + 2*self.dG_pH(p) * pH) + 2 * DG2 + 2*self.dG_pH(p) * (p.pHstroma - pH)
        Keq = np.exp(-DG/self.RT(p))
        return Keq

    # ====================================================================== #
    # Conserved quantities -> for algebraic module #
    # ====================================================================== #
    def Pimoiety(self, p, y):
        PGA,BPGA,GAP,DHAP,FBP,F6P,G6P,G1P,SBP,S7P,E4P,X5P,R5P,RUBP,RU5P,ATP = y
        return np.array([p.Cp - (PGA + 2*BPGA + GAP + DHAP + 2*FBP + F6P + G6P + G1P + 2*SBP + S7P + E4P + X5P + R5P + 2*RUBP + RU5P + ATP)])

    def pqmoiety(self, p, PQ):
        return p.PQtot - PQ

    def pcmoiety(self, p, PC):
        return p.PCtot - PC

    def fdmoiety(self, p, Fd):
        return p.Fdtot - Fd

    def adpmoiety(self, p, ATP):
        return p.APtot - ATP

    def nadpmoiety(self, p, NADPH):
        return p.NADPtot - NADPH

    def N(self, p, y):
        """Used several times to calculate the rate of vPGA, vGAP and vDHAP"""
        Pi, PGA, GAP, DHAP = y
        return np.array([(1+(1+(p.Kpxt/p.Pext))*((Pi/p.Kpi)
                   +(PGA/p.Kpga)
                   +(GAP/p.Kgap)
                   +(DHAP/p.Kdhap)))])

    # ====================================================================== #
    # Light function
    # ====================================================================== #

    def light(self, p, **kwargs):
        '''
        :return: light intensity at certain point of time. 
        Typical PAM light function
        '''
        return p.PFD
        #return ((np.sin(kwargs['t']*0.06)+2)*p.PFD)
    
#        if kwargs['t']%p.dT <=0.8:
#            print('BLA')
#            return 5000
#        elif ((kwargs['t'] > p.Ton) and (kwargs["t"] < p.Toff) and \
#              (kwargs['t']%p.dT>0.8)):
#            #print(p.PFD)
#            return p.PFD
#        else:
#            return 0.0000001


    # ====================================================================== #
    # Reaction rates
    # ====================================================================== #
    #def ps2states(self, p, PQox, PQred, LHC, PSIItot, Psbs, Viola):
    def ps2states(self, p, PQox, PQred, LHC, Psbs, Viola, **kwargs):
        """ 
        QSSA, calculates the states of the photosystem II
        accepts values:
            Pox: oxidised fraction of the PQ pool (PQH2)
            Q: quencher
            L: light, int or array of the n x 1 dimension, that gives the light intensity

            returns:
            B: array of arrays with the states of PSII; rows: time, columns states: 1 and 3 excited
        """
        L = self.LII(p, LHC,**kwargs)

        Q = self.vQuencher4states(p, Psbs, Viola)

        k2 = p.k2
        kF = p.kF
        kH = p.kH0 + p.kH * Q
        k3p = p.kPQred * PQox
        k3m = p.kPQred * PQred / self.Keq_PQred(p)
        M = np.array([[-L-k3m, kH+kF,       k3p, 0],
                      [L,      -(kH+kF+k2), 0,   0],
                      [0,      0,           L,   -(kH+kF)],
                      [1,      1,           1,   1]])

        A = np.array([0, 0, 0, p.PSIItot])
        B = np.linalg.solve(M, A)
        return B

    def ps1states(self, p, PC, PCred, Fd, Fdred, LHC,**kwargs):
        """ 
        QSSA calculates open state of PSI
        depends on reduction states of plastocyanin and ferredoxin
        C = [PC], F = [Fd] (ox. forms)
        accepts: light, y as an array of arrays
        returns: array of PSI open
        """
        L = self.LI(p, LHC,**kwargs)

        A1 = p.PSItot / (1 + L/(p.kFdred * Fd) + (1 + Fdred/(self.Keq_FAFd(p) * Fd))
                          * (PC/(self.Keq_PCP700(p) * PCred)
                             + L/(p.kPCox * PCred))
        )
        return A1

    ###############################################################################
    # method to calculate cross sections
    ###############################################################################
    def crossSection(self, p, LHC):
        """ calculates the cross section of PSII """
        cs = p.staticAntII + (1 - p.staticAntII - p.staticAntI) * LHC
        return cs

    def LII(self,p,LHC,**kwargs):
        return self.crossSection(p, LHC) * self.light(p, **kwargs)

    def LI(self,p,LHC,**kwargs):
        return (1-self.crossSection(p, LHC)) * self.light(p, **kwargs)

    ###############################################################################
    # Reaction rates
    ###############################################################################
    #def vDegradation(self, p, P, Pred, LHC, B, Psbs, Viola):
    def vDegradation(self, p, P, Pred, LHC, B, Psbs, Viola, **kwargs):
        """ rate of damage to photosystem II, 
        proportional to the occupation of 'closed' states"""
        #B = self.ps2states(p, P, Pred, LHC, B, Psbs, Viola)
        B = self.ps2states(p, P, Pred, LHC, B, Psbs, Viola, **kwargs)
        return p.kdeg*(B[1]+B[3])
    
    def vRepair(self, p, B):
        """ rate of the repair of photosystem 
        krep as in Nikolaou et al. 2015"""
        return p.krep * (1-B/p.PSIItot)
                
    #def vPS2(self, p, P, Pred, LHC, B, Psbs, Viola):        
    def vPS2(self, p, P, Pred, LHC,  Psbs, Viola, **kwargs):
        """ reaction rate constant for photochemistry """
        #Q = self.quencher(p,Q,H)
        B = self.ps2states(p, P, Pred, LHC,  Psbs, Viola, **kwargs)
        v = p.k2 * B[1] / 2
        return v

    def vPS1(self, p, PC, PCred, Fd, Fdred, LHC, **kwargs):
        """ reaction rate constant for open PSI """
        L = self.LI(p, LHC, **kwargs)
        A = self.ps1states(p, PC, PCred, Fd, Fdred, LHC, **kwargs)
        v = L * A
        return v
    
    def oxygen(self, p, **kwargs):
        """ return oxygen and NDH concentration as a function of time
        used to simulate anoxia conditions as in the paper"""
        if p.ox == True:
            ''' by default we assume constant oxygen supply'''
            return p.O2ext, p.kNDH
        else:
            if kwargs['t']<= p.Ton or kwargs['t']>=p.Toff:
                return p.O2ext, 0
            else:
                return 0, p.kNDH

    def vPTOX(self, p, Pox, Pred, **kwargs):
        """ calculates reaction rate of PTOX """
        v = Pred * p.kPTOX * self.oxygen(p, **kwargs)[0] 
        return v

    def vNDH(self, p, Pox, **kwargs):
        """ 
        calculates reaction rate of PQ reduction under absence of oxygen
        can be mediated by NADH reductase NDH
        """
        v = self.oxygen(p, **kwargs)[1] * Pox
        return v

    def vB6f(self, p, Pox, Pred, PC, PCred, H):
        """ calculates reaction rate of cytb6f """
        ph = self.pH(H)
        Keq = self.Keq_cytb6f(p,ph)
        v = max(p.kCytb6f * (Pred * PC**2 - (Pox * PCred**2)/Keq), -p.kCytb6f)
        return v

    def vCyc(self, p, Pox, Fdred):
        """
        calculates reaction rate of cyclic electron flow
        considered as practically irreversible
        """
        v = p.kcyc * ((Fdred**2) * Pox)
        return v

    def vFNR(self, p, Fd, Fdred, NADPH, NADP):
        """
        Reaction rate mediated by the Ferredoxin—NADP(+) reductase (FNR)
        Kinetic: convenience kinetics Liebermeister and Klipp, 2006
        Compartment: lumenal side of the thylakoid membrane
        Units:
        Reaction rate: mmol/mol Chl/s
        [F], [Fdred] in mmol/mol Chl/s
        [NADPH] in mM
        """
        fdred = Fdred/p.KM_FNR_F
        fdox = Fd/p.KM_FNR_F
        nadph = (NADPH/p.convf)/p.KM_FNR_N  # NADPH requires conversion to mmol/mol of chlorophyll
        nadp = (NADP/p.convf)/p.KM_FNR_N
        v = (p.EFNR * p.kcatFNR *
            ((fdred**2) * nadp - ((fdox**2) * nadph) / self.Keq_FNR(p)) /
            ((1+fdred+fdred**2) * (1+nadp) + (1+fdox+fdox**2) * (1+nadph) - 1))
        return v

    def vLeak(self, p, Hlf):
        """ 
        rate of leak of protons through the membrane
        """
        v = p.kLeak * (Hlf - self.pHinv(p.pHstroma))
        return v

    def vSt12(self,p, Pox,Ant):
        """ 
        reaction rate of state transitions from PSII to PSI
        Ant depending on module used corresponds to non-phosphorylated antennae
        or antennae associated with PSII
        """
        kKin = p.kStt7 * ( 1 / (1 + ((Pox /p.PQtot)/p.KM_ST)**p.n_ST))
        v = kKin * Ant
        return v

    def vSt21(self,p, Ant):
        """
        reaction rate of state transitions from PSI to PSII
        """
        v = p.kPph1 * (1 - Ant)
        return v

    def vATPsynthase(self, p, Pi, ATPm, ADPm, H):
        """
        Reaction rate of ATP production
        Kinetic: simple mass action with PH dependant equilibrium
        Compartment: lumenal side of the thylakoid membrane
        Units:
        Reaction rate: mmol/mol Chl/s
        [ATP], [ADP] in mM
        """
        ATP = ATPm/p.convf  #change units from mmol/molChl to mM
        ADP = ADPm/p.convf
        
        #Vmax_fw = p.kATPsynth * p.APtot
        #Vmax_bw = p.kATPsynth * p.APtot
        
        #Km_fw = 60 #where is this number from?
        #pH = self.pH(H)
        #Km_bw = self.Keq_ATP(p, pH, Pi)*(Vmax_bw/Vmax_fw)* Km_fw
        #v = ((Vmax_fw * Pi/p.convf/1000 * ADP/Km_fw) - (Vmax_bw * ATP/Km_bw)/(1 + ADP/Km_fw + ATP/Km_bw))
        
        ph = self.pH(H)
        v = p.kATPsynth * (ADP - ATP / self.Keq_ATP(p, ph, Pi))     
        return v


    def vQuencherBase(self, p, Q, H):
        vNPQdyn = H**p.nH / (H ** p.nH + self.pHinv(p.NPQsw)**p.nH)
        v = (1-Q) * p.kNh * vNPQdyn - p.kNr * Q
        return v    
    
    
    def vDeepox(self, p, V, Hlf):
        """
        activity of xantophyll cycle: de-epoxidation of violaxanthin, modelled by Hill kinetics
        """
        nH = p.kHillX
        vf = p.kDeepoxV * ((Hlf ** nH)/ (Hlf ** nH + self.pHinv(p.kphSat) ** nH)) * V
        return vf
    
    def vEpox(self, p, V, Hlf):
        """
        activity of xantophyll cycle: epoxidation
        """
        nH = p.kHillX
        vr = p.kEpoxZ * (p.Xtot - V)
        return vr

    def vLhcprotonation(self, p, L, Hlf):
        """
        activity of PsbS protein protonation: protonation modelled by Hill kinetics
        """
        nH = p.kHillL # shall be changed, for now =5; also half saturation of pH could be change
        vf = p.kProtonationL * ((Hlf ** nH)/ (Hlf ** nH + self.pHinv(p.kphSatLHC) ** nH)) * L
        return vf

    def vLhcdeprotonation(self, p, L, Hlf):
        """
        activity of PsbS protein protonation: deprotonation
        """
        vr = p.kDeprotonation * (p.PsbStot - L)
        return vr

    def vQuencher4states(self, p, L, V):
        """ 
        co-operatiove 4-state quenching mechanism
        """
        Lp = p.PsbStot - L
        Z = p.Xtot - V
        ZAnt = Z / (Z + p.kZSat)
        Q = p.gamma0 * V * L + p.gamma1 * V * Lp + p.gamma2 * ZAnt * Lp + p.gamma3 * ZAnt * L

        return Q

    
    """ Calvin-Benson-Bassham Cycle Reaction Rates """

    def H_stroma(self, p):
        return (10.0**((-1.0)*p.pHstroma))*1000.0

    def v1(self, p, RUBP, PGA, FBP, SBP, P, NADPH):
        """ rate of RuBisCO
        3 Ribulose-1,5-bisphosphate + 3 CO2
        -- RuBisCO -->
        6 3-Phosphoglycerate

        RuBp + CO2 -> PGA
        """
        return (p.V1*RUBP*p.CO2)/((RUBP+p.Km1*(1+(PGA/p.Ki11)+(FBP/p.Ki12)+(SBP/p.Ki13)+(P/p.Ki14)+(NADPH/p.Ki15)))*(p.CO2+p.KmCO2))
    
    def v2(self, p, ATP,PGA,ADP, BPGA):
        """
        6 3-Phosphoglycerate + 6 ATP
        -- Phosphoglycerate kinase (PGK) -->
        6 1,3-Bisphosphoglycerate + 6 ADP

        PGA + ATP -> BPGA + ADP
        Assumed to be at equilibrium
        """
        return p.k*((ATP*PGA)-(1/p.q2)*(ADP*BPGA))

    def v3(self, p, NADPH, BPGA, GAP, NADP, P):
        """
        6 1,3-Bisphosphoglycerate + 6 NADPH + 6 H+
        -- Glyceraldehyde 3-phosphate dehydrogenase (GADPH)-->
        1 G3P + 5 Glyceraldehyde 3-phosphate

        BPGA + NADPH -> GAP + NADP
        Assumed to be at equilibrium
        Stroma pH is assumed to be constant
        """
        
        return p.k*((NADPH*BPGA*self.H_stroma(p))-(1/p.q3)*(GAP*NADP*P))

    def v4(self, p,GAP, DHAP):
        """ 
        (5) Glyceraldehyde 3-phosphate
        -- Triose phosphate isomerae (TPI)-->
        (?) Dihydroxyacetone phosphate

        GAP -> DHAP
        Assumed to be at equilibrium
        """
        return p.k*((GAP)-(1/p.q4)*(DHAP))

    def v5(self, p,GAP, DHAP, FBP):
        """
        (5) Glyceraldehyde 3-phosphate + (?) Dihydroxyacetone phosphate
        -- Aldolase (ALD)-->
        Frucose 1,6-bisphosphate

        GAP + DHAP -> FBP
        Assumed to be at equilibrium
        """
        return p.k*((GAP*DHAP)-(1/p.q5)*(FBP))

    def v6(self, p,FBP, F6P, P):
        """
        (?) Fructose 1,6-bisphosphate + (?) H20
        --Fructose 1,6-bisphosphatase (FBPase) -->
        (6) Fructose 6-phosphate + (Pi)
        FBP -> F6P
        """
        return (p.V6*FBP)/(FBP+p.Km6*(1+(F6P/p.Ki61)+(P/p.Ki62)))

    def v7(self, p,GAP, F6P, X5P, E4P):
        """
        (?) Fructose 6-phosphate + (?) Glyceraldehyde 3-phosphate
        -- Transketolase (TK) -->
        (?) Xylulose 5-phosphate + (?) Erythrose 4-phosphate

        GAP + F6P -> X5P + E4P
        Assumed to be at equilibrium
        """
        return p.k*((GAP*F6P)-(1/p.q7)*(X5P*E4P))

    def v8(self, p, DHAP, E4P, SBP):
        """
        (?) Dihydroxyacetone phosphate + (?) Erythrose 4-phosphate
        -- Aldolase (ALD)-->
        (?) Sedoheptulose 1,7-bisphosphate

        DHAP + E4P -> SBP
        Assumed to be at equilibrium
        """
        return p.k*((DHAP*E4P)-(1/p.q8)*(SBP))

    def v9(self, p,SBP, P):
        """
        (?) Sedoheptulose 1,7-bisphosphate + H20
        --Sedoheptulose 1,7-bisphosphatase (SBPase)-->
        (?) Sedoheptulose 7-phosphate + (?) Pi

        SBP -> S7P
        """

        return (p.V9*SBP)/(SBP+p.Km9*(1+(P/p.Ki9)))

    def v10(self, p,GAP, S7P, X5P, R5P):
        return p.k*((GAP*S7P)-(1/p.q10)*(X5P*R5P))

    def v11(self, p,R5P, RU5P):
        return p.k*((R5P)-(1/p.q11)*(RU5P))

    def v12(self, p, X5P, RU5P):
        return p.k*((X5P)-(1/p.q12)*(RU5P))

    def v13(self, p, RU5P, ATP, PGA, RUBP, P, ADP):
        return (p.V13*RU5P*ATP)/((RU5P+p.Km131*(1+(PGA/p.Ki131)+(RUBP/p.Ki132)+(P/p.Ki133)))*(ATP*(1+(ADP/p.Ki134))+p.Km132*(1+(ADP/p.Ki135))))

    def v14(self, p, F6P, G6P):
        return p.k*((F6P)-(1/p.q14)*(G6P))

    def v15(self, p, G6P, G1P):
        return p.k*((G6P)-(1/p.q15)*(G1P))

    def vpga(self, p, PGA, N):
        return (p.Vx*PGA)/(N*p.Kpga)

    def vgap(self, p, GAP, N):
        return (p.Vx*GAP)/(N*p.Kgap)

    def vdhap(self, p, DHAP, N):
        return (p.Vx*DHAP)/(N*p.Kdhap)

    def vStarch(self, p, G1P, ATP, ADP, P, PGA, F6P, FBP):
        """G1P -> Gn-1 ; Starch production"""
        return (p.Vst*G1P*ATP)/((G1P+p.Kmst1)*((1+(ADP/p.Kist))*(ATP+p.Kmst2)+((p.Kmst2*P)/(p.Kast1*PGA+p.Kast2*F6P+p.Kast3*FBP))))

    def oxPPP(self, p,Pi):
        #return (p.oxPPP)*(15/Pi)#/((Pi*p.Kpi))
        return p.oxPPP#*(15-Pi)#/((Pi*p.Kpi))