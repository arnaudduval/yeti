from __future__ import annotations
from yeti_iga.pymfiga.common.base.cls import BaseSingleModel, BaseMultiModel
from yeti_iga.pymfiga.common.numerics.solvers import OuterToleranceSetter
from yeti_iga.pymfiga.common.numerics.solvers import LinearSolver
from .core import Physics, clear_bcs
from typing import  Literal, Union, List
from copy import deepcopy
from time import time
import numpy as np
import logging
import os 

logger = logging.getLogger("SRC.PHYSICS")


class StaticElastoPlasticity(Physics):
    def __init__(self, **solver_args):
        super().__init__(**solver_args)

    def solve(
        self,
        model: Union[BaseSingleModel, BaseMultiModel],
        displacement: np.ndarray,
        external_force: np.ndarray,
    ):

        assert displacement.ndim == 1, external_force.ndim == 1
        external_force = external_force.copy()
        clear_bcs(model, displacement, external_force)
        model.set_update_manager(self.update_manager)

        logger.info(f"Static linear elastic solver")
        start = time()
        nonlinearsolver = self.nonlinear_solver
        if self.input_args.auto_outer_tolerance:
            outer_tolerance_generator = OuterToleranceSetter(
                outer_tolerance_type="mesh_dependent",
                outer_tolerance_args={"default": self.input_args.tolerance_nonlinear},
            )
            nonlinearsolver.update(tolerance=outer_tolerance_generator.compute(model))
        output = nonlinearsolver.solve(
            displacement,
            model.compute_residual,
            model.solve_linearized_system,
            residual_args={"plastic_vars": {}, "external_force": external_force},
            increment_args={"linear_solver_backend": self.linear_solver},
        )
        self.output_args = output
        logger.info(f"Solve elastoplastic problem in {time() - start:.2e} seconds")


class IncrementalElastoPlasticity(Physics):
    def __init__(self, **solver_args):
        super().__init__(**solver_args)

    def solve(
        self,
        model: Union[BaseSingleModel, BaseMultiModel],
        displacement_list: np.ndarray,
        external_force_list: np.ndarray,
        save_plastic_vars: bool = False,
    ) -> List[dict]:

        assert displacement_list.ndim == 2, external_force_list.ndim == 2
        external_force_list = external_force_list.copy()
        clear_bcs(model, displacement_list, external_force_list)
        model.set_update_manager(self.update_manager)

        logger.info(f"Quasi-static elastoplastic solver")
        start = time()

        # Time-stepping problem
        all_plastic_vars: List[dict] = []
        plastic_vars: dict = {}
        nsteps = external_force_list.shape[0] - 1
        nonlinearsolver = self.nonlinear_solver
        if self._input_args.auto_outer_tolerance:
            outer_tolerance_generator = OuterToleranceSetter(
                outer_tolerance_type="mesh_dependent",
                outer_tolerance_args={"default": self._input_args.tolerance_nonlinear},
            )
            nonlinearsolver.update(tolerance=outer_tolerance_generator.compute(model))
        for it in range(1, nsteps + 1):
            start_step = time()
            logger.info(f"Time marching step: {it}/{nsteps}")
            self.update_manager.increment_step()

            # Predict values of new step
            Fext_n1 = np.copy(external_force_list[it])
            dj_n1 = np.copy(displacement_list[it - 1])

            residual_args = {"plastic_vars": plastic_vars, "external_force": Fext_n1}
            output = nonlinearsolver.solve(
                dj_n1,
                model.compute_residual,
                model.solve_linearized_system,
                residual_args=residual_args,
                increment_args={"linear_solver_backend": self.linear_solver},
            )

            displacement_list[it] = np.copy(dj_n1)
            plastic_vars: dict = output["extra_args"]["new_plastic_vars"]
            if save_plastic_vars:
                all_plastic_vars.append(deepcopy(plastic_vars))
            logger.info(f"Time-step {it} in {time() - start_step:.2e} seconds")

        logger.info(f"Solve elastoplastic problem in {time() - start:.2e} seconds")
        return all_plastic_vars


class ExplicitLinearDynamics(Physics):
    def __init__(self):
        super().__init__()

    def solve(
        self,
        model,
        displacement_initial: np.ndarray,
        external_force,
        time_list,
        velocity_initial: np.ndarray | None = None,
        nb_save: int = 500,
        use_preconditioner: bool = True,
        preconditioner_type: Literal["fastdiag", "scaled_mass", None] = "fastdiag",
    ):

        assert displacement_initial.ndim == 1

        model.set_update_manager(self.update_manager)

        def predict_displacement(dis, vel, acc, dt):
            return dis + dt * vel + 0.5 * dt**2 * acc

        def update_velocity(vel, acc_old, acc_new, dt):
            return vel + 0.5 * dt * (acc_old + acc_new)

        def compute_acceleration(res, linear_solver_type: str = "cg"):
            linear_solver = LinearSolver(
                maxiters=self.input_args.maxiters_linear,
                tolerance=self.input_args.tolerance_linear,
                linear_type=linear_solver_type,
            )

            return model.solve_linearized_system(
                res,
                linear_solver_backend=linear_solver,
                use_preconditioner=use_preconditioner,
                preconditioner_type=preconditioner_type,
            )

        def get_external_force(it, t):
            if callable(external_force):
                Fext = external_force(t, it)
            else:
                Fext = np.asarray(external_force)

            Fext = np.copy(Fext)

            if Fext.ndim != 1:
                raise ValueError("external_force doit être un vecteur 1D.")

            return Fext

        logger.info("Explicit dynamics solver")
        total_start = time()

        time_list = np.asarray(time_list)
        nsteps = len(time_list) - 1

        if nsteps < 1:
            raise ValueError("time_list doit contenir au moins deux instants.")

        if len(time_list) <= nb_save:
            save_indices = np.arange(len(time_list))
        else:
            save_indices = np.unique(
                np.round(np.linspace(0, nsteps, nb_save)).astype(int)
            )

        save_set = set(save_indices)

        saved_displacements = []
        saved_times = []

        d_n0 = np.copy(displacement_initial)

        if velocity_initial is None:
            v_n0 = np.zeros_like(d_n0)
        else:
            v_n0 = np.copy(velocity_initial)

        Fext = get_external_force(0, time_list[0])

        if Fext.shape != d_n0.shape:
            raise ValueError(
                f"Dimensions incompatibles : displacement_initial.shape = {d_n0.shape}, "
                f"external_force.shape = {Fext.shape}"
            )

        if v_n0.shape != d_n0.shape:
            raise ValueError(
                f"Dimensions incompatibles : velocity_initial.shape = {v_n0.shape}, "
                f"displacement_initial.shape = {d_n0.shape}"
            )

        residual0 = model.compute_residual(
            d_n0,
            external_force=Fext,
        )[0]

        a_n0 = compute_acceleration(residual0)

        if 0 in save_set:
            saved_displacements.append(np.copy(d_n0))
            saved_times.append(time_list[0])

        for it in range(1, nsteps + 1):

            start = time()
            logger.info(f"Time marching step: {it}/{nsteps}")
            self.update_manager.increment_step()

            dt = time_list[it] - time_list[it - 1]

            Fext = get_external_force(it, time_list[it])

            d_n1 = predict_displacement(d_n0, v_n0, a_n0, dt)

            residual = model.compute_residual(
                d_n1,
                external_force=Fext,
            )[0]

            a_n1 = compute_acceleration(residual)
            v_n1 = update_velocity(v_n0, a_n0, a_n1, dt)

            if it in save_set:
                saved_displacements.append(np.copy(d_n1))
                saved_times.append(time_list[it])

            d_n0 = np.copy(d_n1)
            v_n0 = np.copy(v_n1)
            a_n0 = np.copy(a_n1)

            logger.info(f"Time-step {it} in {time() - start:.2e} seconds")

        logger.info(
            f"Solve explicit dynamics problem in {time() - total_start:.2e} seconds"
        )

        return np.array(saved_displacements), np.array(saved_times) 
        