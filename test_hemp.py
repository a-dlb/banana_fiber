#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov 14 15:35:08 2026

@author: antoinedlbck
"""

from props import mixture
from ply import Ply
from laminate import Laminate

E_f = 52e3 # MPa
V_f = 0.57 # -
G_f = 20e3 # MPa
nu_f = 0.39 # - 
rho_f = 1.48e3 # kg/m^3

fiber = [E_f, V_f, G_f, nu_f, rho_f]

E_m = 3e3 # MPa
V_m = 1 - V_f # -
G_m = 1.15e3 # MPa
nu_m = 0.29 # -
rho_m = 1.10e3 # kg/m^3

matrix = [E_m, V_m, G_m, nu_m, rho_m]

hemp = mixture(fiber, matrix)

ply1 = Ply(*hemp, t=0.29, theta=90)
ply2 = Ply(*hemp, t=0.29, theta=-45)
ply3 = Ply(*hemp, t=0.29, theta=45)
ply4 = Ply(*hemp, t=0.29, theta=0)
ply5 = Ply(*hemp, t=0.29, theta=0)
ply6 = Ply(*hemp, t=0.29, theta=45)
ply7 = Ply(*hemp, t=0.29, theta=-45)
ply8 = Ply(*hemp, t=0.29, theta=90)

lam = Laminate([ply1, ply2, ply3, ply4, ply5, ply6, ply7, ply8])

loading = [0,5, 0, 0, 0, 0 , 0]
lam.loading = loading