#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 13 20:00:57 2026

@author: antoinedlbck
"""

from dataclasses import dataclass, field
import numpy as np

@dataclass
class Laminate:
    
    plies: list
    
    # Same idea as in class Ply, store only loading and recover eps from loading
    _loading: np.ndarray = field(default_factory = lambda: np.zeros(6), init=False)
    
    @property
    def t(self):
        """
        List of the ply thicknesses (list, dtype=float) (mm)
        """
        return [p.t for p in self.plies]
    
    @property
    def e(self):
        """
        Thickness of the whole laminate (float) (mm)
        """
        return sum(self.t)
    
    @property
    def h(self):
        """
        h for all plies (list) (mm)
        """
        return np.cumsum(self.t) - self.e / 2
    
    def h_k(self, k):
        """
        h(k) for ply k (float) (mm)
        """
        if k == -1:
            return -self.e / 2
        
        else:
            return self.h[k]
    
    def z_k(self, k):
        """
        z(k) for ply k (float) (mm)
        """
        return self.h_k(k) - self.t[k]/2
    
    @property
    def z(self):
        """
        z(k) for all plies (list) (mm)
        """
        return [self.z_k(k) for k in range(len(self.plies))]
    
    @property
    def A(self):
        """
        A matrix (3x3) (np.array) (N/mm)
        """
        return np.asarray(sum((self.h_k(k) - self.h_k(k-1)) * p.Q_123 
                              for k, p in enumerate(self.plies)))
    
    @property
    def B(self):
        """
        B matrix (3x3) (np.array) (N)
        """
        return np.asarray(1/2 * sum((self.h_k(k) ** 2 - self.h_k(k-1) ** 2) * p.Q_123 
                                    for k, p in enumerate(self.plies)))
    
    @property
    def D(self):
        """
        D matrix (3x3) (np.array) (N.mm)
        """
        return np.asarray(1/3 * sum((self.h_k(k) ** 3 - self.h_k(k-1) ** 3) * p.Q_123 
                                    for k, p in enumerate(self.plies)))
    
    # We use _loading as the local variable to recover eps, so as in class Ply,
    # we create getter and setter functions for _loading and modify _loading
    # when we cant to recover epsilon.
    
    @property
    def loading(self):
        return self._loading
    
    @loading.setter
    def loading(self, vector):
        """
        Loading of the laminate in structural axes (np.array)
        Elements 0, 1, 2 are the membrane forces N_1, N_2, N_3 (N/mm)
        Elements 3, 4, 5 are the membrane moments M_1, M_2, M_3 (N)
        """
        self._loading = np.asarray(vector)
        
        eps0 = self.eps0
        K = self.K
        
        for k, p in enumerate(self.plies):
            p.eps_123 = eps0 + self.h[k] * K
    
    @property
    def eps0_K(self):
        
        A = self.A
        B = self.B
        D = self.D
        
        ABD = np.block([[A, B], [B, D]])
        
        return np.linalg.solve(ABD, self._loading)
    
    @property
    def eps0(self):
        return self.eps0_K[:3]
    
    @property
    def K(self):
        return self.eps0_K[3:]
    
    