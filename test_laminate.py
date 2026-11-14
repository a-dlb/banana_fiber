#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov 14 13:09:24 2026

@author: antoinedlbck
"""

from ply import Ply
from laminate import Laminate

E_x = 38000 # MPa
E_y = 9000 # MPa
G_xy = 3600 # MPa
nu_xy = 0.32
rho = 1500 # kg/m^3
mat = [E_x, E_y, G_xy, nu_xy, rho]

ply1 = Ply(*mat, t=1.5, theta=15)
ply2 = Ply(*mat, t=1,   theta=-30)
ply3 = Ply(*mat, t=1.5, theta=-15)
ply4 = Ply(*mat, t=1,   theta=30)

lam = Laminate([ply1, ply2, ply3, ply4])

print("A matrix : \n", lam.A)
print("B matrix : \n", lam.B)
print("D matrix : \n", lam.D)

loading = [1000, 500, 250, 0, 0, 0]

lam.loading = loading

print("eps0 : \n", lam.eps[:3])
print("k : \n", lam.eps[3:])