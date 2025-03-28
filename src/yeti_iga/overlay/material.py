# !/usr/bin/env python
# -*- coding: utf-8 -*-

class Material:
    def __init__(self, density=0.):
        self.density = density

class ElasticMaterial(Material):
    def __init__(self, young_modulus, poisson_ratio, density=0.):
        super().__init__(density)
        self.young_modulus = young_modulus
        self.poisson_ratio = poisson_ratio

