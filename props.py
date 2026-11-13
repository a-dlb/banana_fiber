#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 13 18:25:28 2026

@author: antoinedlbck
"""

def _E_x(E_m, V_m, E_f, V_f):
    return E_m * V_m + E_f * V_f

def _E_y(E_m, V_m, E_f, V_f):
    return E_m * ((1-V_f) + E_m/E_f * V_f) ** (-1)

def _G(G_m, V_m, G_f, V_f):
    return G_m * ((1-V_f) + G_m/G_f * V_f) ** (-1)

def _nu_xy(nu_m, V_m, nu_f, V_f):
    return nu_m * V_m + nu_f * V_f

def _rho(rho_m, V_m, rho_f, V_f):
    return rho_m * V_m + rho_f * V_f

def mixture(fiber, matrix):
    
    E_f, V_f, G_f, nu_f, rho_f = fiber
    E_m, V_m, G_m, nu_m, rho_m = matrix
    
    E_x = _E_x(E_m, V_m, E_f, V_f)
    E_y = _E_y(E_m, V_m, E_f, V_f)
    G = _G(G_m, V_m, G_f, V_f)
    nu_xy = _nu_xy(nu_m, V_m, nu_f, V_f)
    rho = _rho(rho_m, V_m, rho_f, V_f)
    
    return E_x, E_y, G, nu_xy, rho