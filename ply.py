#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 13 16:52:26 2026

@author: antoinedlbck
"""

from dataclasses import dataclass, field
from copy import copy
import numpy as np

@dataclass
class Ply:
    
    E_x: float # MPa
    E_y: float # MPa
    G_xy: float  # MPa
    nu_xy: float # -
    rho: float # kg/m^3
    
    t: float # mm
    theta: float = 0 # deg
    
    X_t: float = 0 # MPa
    X_c: float = 0 # MPa
    Y_t: float = 0 # MPa
    Y_c: float = 0 # MPa
    S: float   = 0 # MPa
    
    # Note: cannot initialize np.ndarray = np.zeros(3) in a class, so we need to use field()
    #       default_factory = function called when creating an instance of the class.
    #       repr = controls if the variable appears in print(ply).
    #       init = controls if the variable can be initialized when creating an instance of the class.
    # => the two latter args are set to False, since _sigma_123 is a 'private' (kinda) 
    # Only sigma123 is stored, see below

    _sigma_123: np.ndarray = field(default_factory = lambda : np.zeros(3),
                                   init = False, repr = False) # MPa
    
    _K_123: np.ndarray = field(default_factory = lambda : np.zeros(3),
                               init = False, repr = False) # -
    
    @property
    def nu_yx(self):
        return self.nu_xy * (self.E_y/self.E_x)
    
    @property
    def m(self):
        return 1/(1 - self.nu_xy * self.nu_yx)
    
    @property
    def Q_xyz(self):
        
        E_x = self.E_x
        E_y = self.E_y
        m = self.m
        nu_xy = self.nu_xy
        nu_yx = self.nu_yx
        G_xy = self.G_xy
        
        return np.array([[m * E_x, m * nu_yx * E_x, 0], [m * nu_xy * E_y, m * E_y, 0], [0, 0, G_xy]])
    
    @property
    def Q_123(self):
        return self.T_sigma(source="xyz") @ self.Q_xyz @ self.T_eps(source="123")
    
    def T_sigma(self, source="xyz"):
        """ 
        Transformation matrix (stresses): 
            from 123 -> xyz : source = "123" 
            from xyz -> 123 : source = "xyz"
        """
        
        theta = np.deg2rad(self.theta)
        
        if source == "123":
            theta = -theta
        
        c = np.cos(theta)
        s = np.sin(theta)
        
        c2 = c ** 2
        s2 = s ** 2
        
        cs = c * s
        
        c2s2 = c2 - s2
        
        return np.array([[c2, s2, -2 * cs], [s2, c2, 2 * cs], [cs, -cs, c2s2]])
    
    def T_eps(self, source="xyz"):
        """ 
        Transformation matrix (strains): 
            from 123 -> xyz : source = "123" 
            from xyz -> 123 : source = "xyz"
        """
        
        theta = np.deg2rad(self.theta)
        
        if source == "123":
            theta = -theta
            
        c = np.cos(theta)
        s = np.sin(theta)
        
        c2 = c ** 2
        s2 = s ** 2
        
        cs = c * s
        
        c2s2 = c2 - s2        
        
        return np.array([[c2, s2, -cs], [s2, c2, cs], [2 * cs, -2 * cs, c2s2]])
    
    
    # Idea: store only _sigma_123, so when we set any field, the other 3
    # are adjusted accordingly. 
    #
    # Principle: define the getter (@property), then define the setter (@field.setter)
    # -> getter: the decorator wraps the field into a property, so that ply.field is called
    #            instead of ply.field() ; it's a shortcut for field = property(field).
    # -> setter: the decorator is a shortcut for field = field.setter(_set),
    #            where _set would be the separately declared setter function.
    #            (Note: .setter() is a method that adds a custom setter function to an object).
    #
    # @property: getter function (recover the field from _sigma_123)
    # @field.setter: setter function (set the field by setting the corresponding _sigma_123)

    @property
    def sigma_123(self):
        return self._sigma_123
    
    @sigma_123.setter
    def sigma_123(self, sigma):
        """
        Stresses in structural axes (MPa)
        """
        self._sigma_123 = np.asarray(sigma)
        
    @property
    def eps_123(self):
        return np.linalg.solve(self.Q_123, self._sigma_123)
    
    @eps_123.setter
    def eps_123(self, eps):
        """
        Strains in structural axes (-)
        """
        self._sigma_123 = self.Q_123 @ eps
        
    @property
    def sigma_xyz(self):
        return self.T_sigma(source="123") @ self._sigma_123
    
    @sigma_xyz.setter
    def sigma_xyz(self, sigma):
        """
        Stresses in material axes (MPa)
        """
        self._sigma_123 = self.T_sigma(source="xyz") @ sigma
        
    @property
    def eps_xyz(self):
        return np.linalg.solve(self.Q_xyz, (self.T_sigma(source="123") @ self._sigma_123))
    
    @eps_xyz.setter
    def eps_xyz(self, eps):
        """
        Strains in material axes (-)
        """
        self._sigma_123 = self.Q_123 @ (self.T_eps(source="xyz") @ eps)
        
    
    @property
    def K_123(self):
        return self._K_123
    
    @K_123.setter
    def K_123(self, K):
        """
        Ply curvature in structural axes (-)
        """
        self._K_123 = np.asarray(K)
        
    @property
    def top(self):
        """
        Get the ply stress and strain fields at the top of the ply.
        """
        
        p = copy(self)
        p.eps_123 = self.eps_123 + self.t / 2 * self.K_123
        
        return p
    
    @property
    def bot(self):
        """
        Get the ply stress and strain fields at the bottom of the ply.
        """
        
        p = copy(self)
        p.eps_123 = self.eps_123 - self.t / 2 * self.K_123
        
        return p
    

    def _max_stress(self):
        
        sigma_xyz = self.sigma_xyz
        
        S_t = [self.X_t, self.Y_t, self.S]
        S_c = [self.X_c, self.Y_c, self.S]
        
        SF = [np.inf if sig == 0
              else S_t[i] / sig if sig > 0
              else -S_c[i] / sig
              for i, sig in enumerate(self.sigma_xyz)]
        
        return min(SF)
    
    
    def failure(self, crit="max_stress"):
        
        if crit == "max_stress":
            
            top_SF = self.top._max_stress()
            bot_SF = self.bot._max_stress()
            
            return (top_SF, "top") if top_SF <= bot_SF else (bot_SF, "bottom")
        
            
        elif crit == "tsai-wu":
            return
        
        else:
            raise ValueError(f"Failure criteria must be 'max_stress' or 'tsai-wu', not {crit}.")