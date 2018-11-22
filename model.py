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
__credits__ = ["Anna Matuszynska", "Nima Saadat", "Oliver Ebenhoeh"]
__maintainer__ = "Anna Matuszynska"
__email__ = "Anna.Matuszynska@uni-duesseldorf.de"
__status__ = "Development"

import modelbase


class Merged2018(modelbase.Model):
    def __init__(self, p, r):
        super().__init__(p)
        
        compounds = [
                #"B",  #photosystem II protein concentration
                "PQ",  # oxidised plastoquinone
                "PC",  # oxidised plastocyan
                "Fd",  # oxidised ferrodoxin
                "ATP",  # stromal concentration of ATP
                "NADPH",  # stromal concentration of NADPH
                "H",  # lumenal protons
                #"Hstr", #stromal protons

                "LHC",  # non-phosphorylated antenna
                "Psbs", # PsBs
                "Vx", # violaxanthin (as ratio of all xantophylls)
                "PGA", 'BPGA', 'GAP', 'DHAP', 'FBP', 'F6P', 'G6P', 'G1P', 'SBP', 'S7P', 'E4P', 'X5P', 'R5P', 'RUBP', 'RU5P',
                "Fluo",
                "Light"]
        
        self.set_cpds(compounds)

        # Algebraic modules
        self.add_algebraicModule(r.pqmoiety, "pq_alm", ["PQ"],  ["PQred"])
        self.add_algebraicModule(r.pcmoiety, "pc_alm", ["PC"],  ["PCred"])
        self.add_algebraicModule(r.fdmoiety, "fd_alm", ["Fd"],  ["Fdred"])
        self.add_algebraicModule(r.adpmoiety, "adp_alm", ["ATP"],  ["ADP"])
        self.add_algebraicModule(r.nadpmoiety, "nadp_alm", ["NADPH"], ["NADP"])
        self.add_algebraicModule(r.Pimoiety, "phosphate_alm", ["PGA", "BPGA", "GAP", "DHAP", "FBP", "F6P", "G6P", "G1P", "SBP", "S7P", "E4P",
                                  "X5P", "R5P", "RUBP", "RU5P", "ATP"], ["Pi"])
        self.add_algebraicModule(r.N, "N_alm",["Pi", "PGA", "GAP", "DHAP"], ["N"]) 

        # add light driven electron transport chain reaction rates        
        #self.set_ratev("vDegradation", r.vDegradation, "PQ", "PQred", "LHC", "B", "Psbs", "Vx") #q
        #self.set_ratev("vDegradation", r.vDegradation, "PQ", "PQred", "LHC", "B") # no quencher
        #self.set_rate("vRepair", r.vRepair, "B")
        self.set_ratev('vPS2', r.vPS2, "PQ", "PQred", "LHC", "Psbs", "Vx")
        #self.set_ratev('vPS2', r.vPS2, "PQ", "PQred", "LHC", "B")
        self.set_ratev("vPS1", r.vPS1, "PC", "PCred", "Fd", "Fdred", "LHC")
        self.set_rate("vPTOX", r.vPTOX, "PQ", "PQred")
        self.set_rate("vB6f", r.vB6f, "PQ", "PQred", "PC", "PCred", "H")
        self.set_rate("vNDH", r.vNDH, "PQ")
        self.set_rate("vCyc", r.vCyc, "PQ", "Fdred")
        self.set_rate("vFNR", r.vFNR, "Fd", "Fdred", "NADPH", "NADP")
        self.set_rate("vLeak", r.vLeak, "H")
        self.set_rate("vSt12", r.vSt12, "PQ", "LHC")
        self.set_rate("vSt21", r.vSt21, "LHC")
        self.set_rate("vATPsynthase", r.vATPsynthase, "Pi", "ATP", "ADP", "H")
        self.set_rate("vLHCprotonation", r.vLhcprotonation, "Psbs", "H")
        self.set_rate("vLHCdeprotonation", r.vLhcdeprotonation, "Psbs", "H")
        self.set_rate("vDeepox", r.vDeepox, "Vx", "H")
        self.set_rate("vEpox", r.vEpox, "Vx", "H")
        self.set_rate("vQuencher4states", r.vQuencher4states, "Psbs", "Vx")
        self.set_ratev("fluorescence", r.fluorescence, "PQ", "PQred", "LHC",  "Psbs", "Vx")
        self.set_ratev("light", r.light)

        #self.set_stoichiometry_byCpd("B", {"vRepair": 1, "vDegradation": -1})
        self.set_stoichiometry_byCpd("PQ", {"vPS2": -1, "vB6f":1, "vCyc": -1, "vPTOX": 1, "vNDH": -1})
        self.set_stoichiometry_byCpd("PC", {"vB6f": -2, "vPS1": 1})
        self.set_stoichiometry_byCpd("Fd", {"vPS1": -1, "vFNR": 2, "vCyc": 2})
        self.set_stoichiometry_byCpd("ATP", {"vATPsynthase": 1*self.par.convf, "vPGA_kinase": -1, "v13": -1, "vStarch": -1})
        self.set_stoichiometry_byCpd("NADPH", {"vFNR": 1*self.par.convf, "vBPGA_dehydrogenase": -1})
        self.set_stoichiometry_byCpd("H", {"vPS2": 2/self.par.bH, "vB6f": 4/self.par.bH, "vATPsynthase": -self.par.HPR/p.bH, "vLeak": -1/self.par.bH})
        #self.set_stoichiometry_byCpd("Hstr", {"vPS2": -1*(2/self.par.bH), "vB6f": -1*(4/self.par.bH), "vATPsynthase": -1*(-self.par.HPR/p.bH), "vLeak": -1*(-1/self.par.bH)})

        self.set_stoichiometry_byCpd("LHC", {"vSt12": -1, "vSt21": 1})
        self.set_stoichiometry_byCpd("Psbs", {"vLHCprotonation": -1, "vLHCdeprotonation": 1})
        self.set_stoichiometry_byCpd("Vx", {"vDeepox": -1, "vEpox": 1})


        # add CBB reaction rates
        self.set_rate("vRuBisCO", r.v1, "RUBP", "PGA", "FBP", "SBP", "Pi", "NADPH")
        self.set_rate("vPGA_kinase", r.v2, "ATP", "PGA", "ADP", "BPGA")
        self.set_rate("vBPGA_dehydrogenase", r.v3, "NADPH", "BPGA", "GAP", "NADP", "Pi")
        self.set_rate("vTPI", r.v4, "GAP", "DHAP")
        self.set_rate("vAldolase", r.v5, "GAP", "DHAP", "FBP")
        self.set_rate("vFBPase", r.v6, "FBP", "F6P", "Pi")
        self.set_rate("vF6P_Transketolase", r.v7, "GAP", "F6P", "X5P", "E4P")
        self.set_rate("v8", r.v8, "DHAP", "E4P", "SBP")
        self.set_rate("v9", r.v9, "SBP", "Pi")
        self.set_rate("v10", r.v10, "GAP", "S7P", "X5P", "R5P")
        self.set_rate("v11", r.v11, "R5P", "RU5P")
        self.set_rate("v12", r.v12, "X5P", "RU5P")
        self.set_rate("v13", r.v13, "RU5P", "ATP", "PGA", "RUBP", "Pi", "ADP")
        self.set_rate("vG6P_isomerase", r.v14, "F6P", "G6P")
        self.set_rate("vPhosphoglucomutase", r.v15, "G6P", "G1P")
        self.set_rate("vpga", r.vpga, "PGA", "N")
        self.set_rate("vgap", r.vgap, "GAP", "N")
        self.set_rate("vDHAP", r.vdhap, "DHAP", "N")
        self.set_rate("vStarch", r.vStarch, "G1P", "ATP", "ADP", "Pi", "PGA", "F6P", "FBP")
        
        ##################
        self.set_rate("voxPPP", r.oxPPP,"Pi")
        #####################
        
        self.set_stoichiometry_byCpd("PGA", {"vRuBisCO": 2, "vPGA_kinase": -1, "vpga": -1})
        self.set_stoichiometry_byCpd("BPGA", {"vPGA_kinase": 1, "vBPGA_dehydrogenase": -1})
        self.set_stoichiometry_byCpd("GAP", {"vBPGA_dehydrogenase": 1, "vTPI": -1, "vAldolase": -1, "vF6P_Transketolase": -1, "v10": -1, "vgap": -1})
        self.set_stoichiometry_byCpd("DHAP", {"vTPI": 1, "vAldolase": -1, "v8": -1, "vDHAP": -1})
        self.set_stoichiometry_byCpd("FBP", {"vAldolase": 1, "vFBPase": -1})
        self.set_stoichiometry_byCpd("F6P", {"vFBPase": 1, "vF6P_Transketolase": -1, "vG6P_isomerase": -1})
        self.set_stoichiometry_byCpd("G6P", {"vG6P_isomerase": 1, "vPhosphoglucomutase": -1})
        self.set_stoichiometry_byCpd("G1P", {"vPhosphoglucomutase": 1, "vStarch": -1})
        self.set_stoichiometry_byCpd("SBP", {"v8": 1, "v9": -1})
        self.set_stoichiometry_byCpd("S7P", {"v9": 1, "v10": -1})
        self.set_stoichiometry_byCpd("E4P", {"vF6P_Transketolase": 1, "v8": -1})
        self.set_stoichiometry_byCpd("X5P", {"vF6P_Transketolase": 1, "v10": 1, "v12": -1})
        self.set_stoichiometry_byCpd("R5P", {"v10": 1, "v11": -1})
        self.set_stoichiometry_byCpd("RUBP", {"v13": 1, "vRuBisCO": -1})
        
        ########################
        self.set_stoichiometry_byCpd("RU5P", {"v11": 1, "v12": 1,"voxPPP":1, "v13": -1})
        ###############################
        
        # finally add fluorescence as a dynamic variable to monitor teh dynamics of photosynthesis
        self.set_stoichiometry_byCpd("Fluo", {"fluorescence": 1})
        self.set_stoichiometry_byCpd("Light", {"light":1})
