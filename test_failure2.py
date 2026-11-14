#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov 14 21:40:31 2026

@author: antoinedlbck
"""

# TP4 Ex. 3 (cantilever beam)

from ply import Ply
from laminate import Laminate
import numpy as np

Q = np.array([[200, 5, 0], [5, 20, 0], [0, 0, 6]]) * 1e3

nu_xy = Q[0, 1] / Q[1, 1]
nu_yx = Q[0, 1] / Q[0, 0]
m = 1 / (1 - nu_xy * nu_yx)
E_x = Q[0, 0] / m
E_y = Q[1, 1] / m
G_xy = Q[2, 2]
mat = [E_x, E_y, G_xy, nu_xy, 1500]

strength = dict(X_t=1270, X_c=2500, Y_t=60, Y_c=200, S=70) # MPa

ply1 = Ply(*mat, t=1, theta= 15, **strength)
ply2 = Ply(*mat, t=1, theta=-15, **strength)
ply3 = Ply(*mat, t=1, theta= 90, **strength)
ply4 = Ply(*mat, t=1, theta= 90, **strength)
ply5 = Ply(*mat, t=1, theta=-15, **strength)
ply6 = Ply(*mat, t=1, theta= 15, **strength)

lam = Laminate([ply1, ply2, ply3, ply4, ply5, ply6])
print("A matrix : \n", lam.A)
print("B matrix : \n", lam.B)
print("D matrix : \n", lam.D)

loading = [2500, 0, 0, 1500, 0, 0]

lam.loading = loading

for k, p in enumerate(lam.plies):
    failure = p.failure("max_stress")
    print(f"ply{k+1} : SF = {failure[0]:.3f}, loc = {failure[1]}")