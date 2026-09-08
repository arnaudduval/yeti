from yeti_iga.pymfiga.common.base.enum import Constants
from yeti_iga.pymfiga.common.numerics.solvers import LinearSolver
from yeti_iga.pymfiga.iga.single_model import ThermalModel as SpThermalModel
from .core import MultiModel
from typing import Tuple, Dict, Union, List, Literal
from time import time
import numpy as np
import logging

logger = logging.getLogger("SRC.IGA.MULTI_MODEL")


class ThermalModel(MultiModel):
    def __init__(
        self,
        patch_models_dict: Dict[Union[str, int], SpThermalModel],
        lagrange_type: Literal["augmented", "standard"] = "augmented",
    ):
        super().__init__(patch_models_dict, lagrange_type)

    def compute_residual(
        self, array_in: np.ndarray, **kwargs
    ) -> Tuple[np.ndarray, dict]:
        start = time()
        external_force: np.ndarray = kwargs["external_force"]
        scalar_coefs: tuple = kwargs["scalar_coefs"]
        arr_reshaped = np.reshape(array_in, (2, -1))
        temperature_cutted = super().cut(arr_reshaped[0])
        flux_cutted = super().cut(arr_reshaped[1])
        external_force_cutted = super().cut(external_force)

        residual_cutted_list, mf_args = {}, {}
        for pid, model in zip(self.patch_ids, self.patch_models):
            kwargs_cutted = {
                "external_force": external_force_cutted[pid],
                "scalar_coefs": scalar_coefs,
            }
            id_temp_and_flux = np.hstack((temperature_cutted[pid], flux_cutted[pid]))
            residual_cutted, mf_args_cutted = model.compute_residual(
                id_temp_and_flux,
                **kwargs_cutted,
            )
            res_temp = np.reshape(residual_cutted, (2, -1))
            residual_cutted_list.update({pid: res_temp[0]})
            mf_args.update({f"mf_{pid}": mf_args_cutted})

        residual = super().glue(residual_cutted_list)
        lagrange_res = self.lagrange_solver.compute_residual(residual, arr_reshaped[0])
        arr_flatten = np.hstack((lagrange_res, np.zeros_like(lagrange_res)))
        logger.debug(f"Computing residual in {time() - start:.2e} seconds")
        return arr_flatten, mf_args

    def compute_mf_tangent(self, array_in: np.ndarray, **kwargs):
        array_in_cutted = super().cut(array_in)
        array_out = {}
        for pid, model in zip(self.patch_ids, self.patch_models):
            kwargs_id = {"scalar_coefs": kwargs["scalar_coefs"]}
            kwargs_id.update(kwargs[f"mf_{pid}"])
            assert isinstance(model, SpThermalModel)
            array_out_cutted = model.compute_mf_tangent(
                array_in=array_in_cutted[pid], **kwargs_id
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
        for model in self._patch_models:
            assert isinstance(model, SpThermalModel)
            mean_capacity = np.array(model.scalar_mean_capacity)
            mean_conductivity = np.array(model.scalar_mean_conductivity)

            scalar_coefs: List[float] = list(kwargs["scalar_coefs"])
            if not model.has_DOFsBlocked_atLeastOnce():
                scalar_coefs[0] = 0.1  # To avoid division by zero

            ## Update preconditioner
            if model.update_manager.should_update_preconditioner:
                model.preconditioner.add_scalar_space_time_correctors(
                    mass_corrector=mean_capacity.tolist(),
                    stiffness_corrector=mean_conductivity.tolist(),
                )
                model.preconditioner.update_space_eigenvalues(scalar_coefs=scalar_coefs)

            ## Update lagrange penalty
            vals = scalar_coefs[0] * mean_capacity + scalar_coefs[1] * mean_conductivity
            lagrange_penalty = max([lagrange_penalty, float(np.max(vals))])

        if self.update_manager.should_update_penalty:
            self.preconditioner.update(lagrange_penalty=lagrange_penalty)
            self.lagrange_solver.update(lagrange_penalty=lagrange_penalty)

        # Solve increment
        flux_factor: tuple = kwargs["flux_factor"]
        UV = kwargs["current_solution"]
        UV_reshaped = np.reshape(UV, (2, -1))
        U = UV_reshaped[0]
        res_UV_reshaped = np.reshape(array_in, (2, -1))
        residual_U = res_UV_reshaped[0]

        increment = self.lagrange_solver.compute_increment(
            apply_T=self.compute_mf_tangent,
            apply_P=self.preconditioner.apply_spatial_preconditioner,
            current_res=residual_U,
            current_sol=U,
            **kwargs,
        )

        logger.debug(f"Solving linearized system in {time() - start:.2e} seconds")
        return np.hstack((increment, increment * flux_factor))
