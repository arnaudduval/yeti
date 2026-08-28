from yeti_iga.pymfiga.iga.boundary import BoundaryCondition
from yeti_iga.pymfiga.common.material.mechanical import IsotropicMat
from yeti_iga.pymfiga.iga.geometry import SinglePatch
from .mechanical import MechanicalModel
from typing import Optional
from time import time
import numpy as np
import logging

logger = logging.getLogger("SRC.IGA.MODEL")


class ExplicitDynamicsModel(MechanicalModel):


    _diagonal_mass: Optional[np.ndarray] = None

    def __init__(
        self,
        material: IsotropicMat,
        patch: SinglePatch,
        boundary: BoundaryCondition,
        mass_type: str = "diagonal_mass",
    ):
        super().__init__(material, patch, boundary)

        assert mass_type in ["diagonal_mass", "consistent_mass"]
        self.mass_type = mass_type

    @property
    def diagonal_mass(self):
        if self._diagonal_mass is None:
            return np.array([])
        return self._diagonal_mass

    def clear_properties(self):
        super().clear_properties()
        if self.material.hasnonlinearmass:
            self._diagonal_mass = None

    
    def compute_mf_mass(self, array_in: np.ndarray, **mf_args) -> np.ndarray:
        if self.mass_type == "consistent_mass":
            return MechanicalModel.compute_mf_mass(self, array_in, **mf_args)

        # sinon masse lumpée
        start = time()

        if self._mass_property is None:
            self.compute_mass_property(**mf_args)

        if self._diagonal_mass is None:
            ones = np.reshape(np.ones_like(array_in), (self.ndim, -1))
            mass = np.zeros_like(ones)

            for i in range(self.ndim):
                mass[i] = self.operator_engine.compute_mf_scalar_u_v(
                    self.part.quadrule_list,
                    self.mass_property,
                    ones[i],
                    allow_lumping=True,
                    nurbs_weights=self.part.nurbs_weights,
                )

            self._diagonal_mass = np.ravel(mass)

        return self.diagonal_mass * array_in


    
    
    def solve_linearized_system(self, array_in: np.ndarray, **kwargs) -> np.ndarray:
        if self.mass_type == "consistent_mass":
            return super().solve_linearized_system(array_in, **kwargs)

        # masse lumpée : résolution directe diagonale
        start = time()

        if self._diagonal_mass is None:
            self._diagonal_mass = self.compute_mf_mass(
                np.ones_like(array_in),
                **kwargs
            )

        array_out = array_in / self.diagonal_mass

        if self.update_manager.should_update_material:
            self.clear_properties()

        logger.debug(f"Solving linearized system in {time() - start:.2e} seconds")
        return array_out
