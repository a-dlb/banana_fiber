#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov 14 18:28:48 2026

@author: antoinedlbck
"""

import numpy as np

def max_stress(lam):
    
    for k, p in enumerate(lam.plies):
        
        S_c = [p.X_c, p.Y_c, p.S]
        S_t = [p.X_t, p.Y_t, p.S]
        
        # Top of the ply
        sigma_xyz = p.top.sigma_xyz
        SF = sigma_xyz/S_t if sigma_xyz >= 0 else -(sigma_xyz)/S_c
        
    return 