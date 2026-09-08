from yeti_iga.pymfiga.common.base.enum import Constants
from yeti_iga.pymfiga.common.numerics.solvers import LinearSolver
from yeti_iga.pymfiga.iga.single_model import MechanicalModel as SpMechanicalModel
from .core import MultiModel
from typing import Tuple, Dict, Union, Literal
from time import time
import numpy as np
import logging

logger = logging.getLogger("SRC.IGA.MULTI_MODEL")


class MechanicalModel(MultiModel):
    def __init__(
        self,
        patch_models_dict: Dict[Union[str, int], SpMechanicalModel],
        lagrange_type: Literal["augmented", "standard"] = "augmented",
    ):
        super().__init__(patch_models_dict, lagrange_type)

    def compute_residual(
        self, array_in: np.ndarray, **kwargs
    ) -> Tuple[np.ndarray, dict]:
        start = time()
        external_force: np.ndarray = kwargs["external_force"]
        plastic_vars: dict = kwargs.get("plastic_vars", {})
        displacement_cutted = super().cut(array_in)
        external_force_cutted = super().cut(external_force)

        residual_cutted_list, mf_args = {}, {}
        new_plastic_vars = {}
        for pid, model in zip(self.patch_ids, self.patch_models):
            kwargs_cutted = {
                "external_force": external_force_cutted[pid],
                "plastic_vars": plastic_vars.get(pid, {}),
            }
            residual_cutted, mf_args_cutted = model.compute_residual(
                displacement_cutted[pid],
                **kwargs_cutted,
            )
            residual_cutted_list.update({pid: residual_cutted})
            mf_args.update({f"mf_{pid}": mf_args_cutted})
            new_plastic_vars.update({pid: mf_args_cutted.get("new_plastic_vars", {})})

        residual = super().glue(residual_cutted_list)
        lagrange_res = self.lagrange_solver.compute_residual(residual, array_in)
        mf_args.update({"new_plastic_vars": new_plastic_vars})
        logger.debug(f"Computing residual in {time() - start:.2e} seconds")
        return lagrange_res, mf_args

    def compute_mf_tangent(self, array_in: np.ndarray, **kwargs):
        array_in_cutted = super().cut(array_in)
        array_out = {}
        for pid, model in zip(self.patch_ids, self.patch_models):
            assert isinstance(model, SpMechanicalModel)
            array_out_cutted = model.compute_mf_stiffness(
                array_in=array_in_cutted[pid],
                **kwargs[f"mf_{pid}"],
            )
            array_out.update({pid: array_out_cutted})
        return super().glue(array_out)

    def solve_linearized_system(self, array_in: np.ndarray, **kwargs) -> np.ndarray:
        start = time()

        linear_solver: LinearSolver = kwargs["linear_solver_backend"]
        inner_tolerance: float = (
            kwargs.get("inner_tolerance") or linear_solver.config.tolerance
        )
        maxiters = linear_solver.config.maxiters
        self.lagrange_solver.update(tolerance=inner_tolerance, maxiters=maxiters)

        # Update Lagrange data
        lagrange_penalty = Constants.SMALL
        for model in self.patch_models:
            assert isinstance(model, SpMechanicalModel)
            mean_mass = np.array(model.scalar_mean_mass)
            mean_stiffness = np.array(model.scalar_mean_stiffness)

            # TODO: include scalar_coefs to treat implicit dynamic problems
            # scalar_coefs: tuple = kwargs.get("scalar_coefs")
            scalar_coefs = [0.0, 1.0]
            if not model.has_DOFsBlocked_atLeastOnce():
                scalar_coefs[0] = 0.1  # To avoid division by zero

            # Update preconditioner
            if model.update_manager.should_update_preconditioner:
                model.preconditioner.add_scalar_space_time_correctors(
                    mass_corrector=mean_mass.tolist(),
                    stiffness_corrector=mean_stiffness.tolist(),
                )
                model.preconditioner.update_space_eigenvalues(scalar_coefs=scalar_coefs)

            # Update lagrange penalty
            vals = scalar_coefs[0] * mean_mass + scalar_coefs[1] * mean_stiffness
            lagrange_penalty = max([lagrange_penalty, float(np.max(vals))])

        if self.update_manager.should_update_penalty:
            self.preconditioner.update(lagrange_penalty=lagrange_penalty)
            self.lagrange_solver.update(lagrange_penalty=lagrange_penalty)

        # Solve increment
        U = kwargs["current_solution"]
        increment = self.lagrange_solver.compute_increment(
            apply_T=self.compute_mf_tangent,
            apply_P=self.preconditioner.apply_spatial_preconditioner,
            current_res=array_in,
            current_sol=U,
            **kwargs,
        )

        logger.debug(f"Solving linearized system in {time() - start:.2e} seconds")
        return increment
