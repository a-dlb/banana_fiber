#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 13 17:55:45 2026

@author: antoinedlbck
"""

from ply import Ply
from props import mixture

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

E_x, E_y, G, nu_xy, rho = mixture(fiber, matrix)

ply1 = Ply(E_x, E_y, G, nu_xy, theta=0)

sigma_1 = 10 # MPa
length = 200 # mm
thick = 0.5 # mm
F_1 = sigma_1 * length * thick # N

ply1.sigma_123 = [10, 0, 0]

print(ply1.eps_123)
