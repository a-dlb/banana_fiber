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
    
    plies: field(default_factory=list)
    
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
    
    def h_k(self, k):
        """
        h(k) for ply k (float) (mm)
        """
        
        t = self.t
        cum_t = np.cumsum(t)
        
        # s is the height of the midplane w.r.t. the bottom of the laminate
        s = self.e / 2
        
        # i is the index of the ply where the midplane is located
        # side="left" means it gives the left index (so the lower ply)
        # if the midplane lies at the interface of two plies
        i = np.searchsorted(cum_t, s, side="left")
        
        # If we are above the midplane, h is the cumulated thickness of the plies
        # above the midplane s + the thickness of the ply i above the midplane.
        # So we truncate t at (i+1), take its cumsum and select the (k-(i+1))th element
        # since we removed i+1 elements by truncating t. 
        # We then add cum_t[i] - s (= 0 if s is at the interface of i and i+1).
        if k > i:
            return np.cumsum(t[i+1:])[k-(i+1)] + (cum_t[i] - s)
        
        # If we are under the midplane, h is the  opposite of the cumulated 
        # thickness of the plies under the midplane s + the thickness of the ply i under
        # the midplane.
        # So we truncate t at (i-1), REVERSE IT (because we start the cumsum 
        # FROM the highest ply index under the midplane), THEN REVERSE IT AGAIN,
        # so we can take the k^th element of the cumsum. 
        # We then add t[i] - (cum_t[i] - s) (=t[i] if s is at the interface of i and i+1),
        # then take the opposite of all that. Whew.
        elif k <= i:
            
            # If k = -1, we return h(0) - t[0].
            if k == -1:
                return -(np.cumsum(t[i-1::-1])[::-1][0] + t[i] - (cum_t[i] - s)) - t[0]
            
            else:
                return -(np.cumsum(t[i-1::-1])[::-1][k] + t[i] - (cum_t[i] - s))
        
    @property
    def h(self):
        """
        h(k) for all plies (list) (mm)
        """
        return [self.h_k(k) for k in range(self.plies)]
    
    
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
        return [self.z_k(k) for k in range(self.plies)]
    
    @property
    def A(self):
        """
        A matrix (3x3) (np.array) (N/mm)
        """
        return sum((h[k] - h[k-1]) * p.Q_123 for k, p in enumerate(self.plies))
    
    @property
    def B(self):
        """
        B matrix (3x3) (np.array) (N)
        """
        return 1/2 * sum((h[k] ** 2 - h[k-1] ** 2) * p.Q_123 for k, p in enumerate(self.plies))
    
    @property
    def D(self):
        """
        D matrix (3x3) (np.array) (N.mm)
        """
        return 1/3 * sum((h[k] ** 3 - h[k-1] ** 3) * p.Q_123 for k, p in enumerate(self.plies))