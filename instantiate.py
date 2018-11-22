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
__status__ = "Production/Stable"

# Import built-in libraries
import modelbase
import numpy as np

# Instantiate 3 modelling objects
import model
import parameters
import reactionrates


def instantiate():
    ''' returns four objects: 
        p: the parameter dictionary
        r: reaction rates
        m: model object
        s: simulator
    '''
    p = parameters.Parameters()
    r = reactionrates.Reactions()
    m = model.Merged2018(p, r)
    s = modelbase.Simulator(m)
    
    print('Created a virtual organism. Experiment ready to run')
    return p, r, m, s
    
    
if __name__=='__main__':
    p, r, m, s = instantiate()

