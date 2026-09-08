from yeti_iga.pymfiga.common.base.cls import BaseSingleModel, BaseMultiModel
from yeti_iga.pymfiga.common.numerics.solvers import OuterToleranceSetter
from .core import Physics, clear_bcs
from typing import List, Union, Literal
from time import time
import numpy as np
import logging

logger = logging.getLogger("SRC.PHYSICS")


class SteadyHeatTransfer(Physics):
    def __init__(self, **solver_args):
        super().__init__(**solver_args)

    def solve(
        self,
        model: Union[BaseSingleModel, BaseMultiModel],
        temperature: np.ndarray,
        external_force: np.ndarray,
    ):
        assert temperature.ndim == 1 and external_force.ndim == 1
        external_force = external_force.copy()
        clear_bcs(model, temperature, external_force)
        model.set_update_manager(self.update_manager)

        logger.info(f"Steady heat transfer solver")
        start = time()
        temp_and_flux = np.hstack((temperature, np.zeros_like(temperature)))
        nonlinearsolver = self.nonlinear_solver
        if self.input_args.auto_outer_tolerance:
            outer_tolerance_generator = OuterToleranceSetter(
                outer_tolerance_type="mesh_dependent",
                outer_tolerance_args={"default": self.input_args.tolerance_nonlinear},
            )
            nonlinearsolver.update(tolerance=outer_tolerance_generator.compute(model))
        output = nonlinearsolver.solve(
            temp_and_flux,
            model.compute_residual,
            model.solve_linearized_system,
            residual_args={
                "external_force": external_force,
                "scalar_coefs": (0, 1),
                "flux_factor": 0.0,
            },
            increment_args={
                "scalar_coefs": (0, 1),
                "flux_factor": 0.0,
                "linear_solver_backend": self.linear_solver,
            },
        )
        temp_and_flux = np.reshape(temp_and_flux, (2, -1))
        temperature[:] = temp_and_flux[0]  # [:] us used to replace old values
        self.output_args = output
        logger.info(f"Solve heat transfer problem in {time() - start:.2e} seconds")


class TransientHeatTransfer(Physics):
    def __init__(self, **solver_args):
        super().__init__(**solver_args)

    def solve(
        self,
        model: Union[BaseSingleModel, BaseMultiModel],
        temperature_list: np.ndarray,
        external_force_list: np.ndarray,
        incremental_type: Literal["alpha", "bdf"],
        **kwargs,
    ):
        def predict_temperature(
            histtemp: List[np.ndarray],
            norder: int,
        ):
            assert len(histtemp) == norder, "Size problem."
            if norder == 1:
                return np.copy(histtemp[0])
            elif norder == 2:
                [y1, y2] = histtemp
                return (4 * y2 - y1) / 3
            elif norder == 3:
                [y1, y2, y3] = histtemp
                return (18 * y3 - 9 * y2 + 2 * y1) / 11
            elif norder == 4:
                [y1, y2, y3, y4] = histtemp
                return (48 * y4 - 36 * y3 + 16 * y2 - 3 * y1) / 25
            else:
                raise ValueError("Order not supported")

        def select_snapshots(i, y, norder):
            "Creates a list of previous solution values."
            return [y[i - k] for k in range(norder, 0, -1)]

        def select_parameter(norder):
            assert norder in [1, 2, 3, 4], "Order not supported"
            parameters = [1.0, 2.0 / 3.0, 6.0 / 11.0, 12.0 / 25.0]
            return parameters[norder - 1]

        assert temperature_list.ndim == 2 and external_force_list.ndim == 2
        external_force_list = external_force_list.copy()
        clear_bcs(model, temperature_list, external_force_list)
        model.set_update_manager(self.update_manager)

        logger.info(f"Transient heat transfer solver")
        start = time()

        # Default values
        alpha, norder = 1.0, 1

        if incremental_type == "alpha":
            # Get variables for alpha (or theta) method
            time_list: Union[np.ndarray, list] = kwargs["time_list"]
            alpha: float = kwargs.get("alpha", alpha)
            assert all(
                x is not None for x in [time_list, alpha]
            ), "time_list or alpha are not defined"
            nsteps = int(len(time_list) - 1)
        elif incremental_type == "bdf":
            # Get variables for BDF (or theta) method
            tspan: List[float] = kwargs["tspan"]
            nsteps: int = kwargs["nsteps"]
            norder: int = kwargs.get("norder", norder)
            assert all(
                x is not None for x in [tspan, nsteps, norder]
            ), "tspan or nsteps or norder are not defined"
            time_list = np.linspace(tspan[0], tspan[1], nsteps + 1)
        else:
            raise NotImplementedError()

        nonlinearsolver = self.nonlinear_solver
        if self.input_args.auto_outer_tolerance:
            outer_tolerance_generator = OuterToleranceSetter(
                outer_tolerance_type="mesh_dependent",
                outer_tolerance_args={"default": self.input_args.tolerance_nonlinear},
            )
            nonlinearsolver.update(tolerance=outer_tolerance_generator.compute(model))

        # Initialize time and solution arrays
        dt = time_list[1] - time_list[0]
        dj_n1 = np.copy(temperature_list[0])
        vj_n1 = np.zeros_like(temperature_list[0])

        logger.info(f"Static step")
        residual_args = {
            "external_force": external_force_list[0],
            "scalar_coefs": (1.0, 1.0),
        }
        increment_args = {
            "scalar_coefs": (1.0, 0.0),
            "flux_factor": 0.0,
            "linear_solver_backend": self.linear_solver,
        }
        dj_vj_n1 = np.hstack((dj_n1, np.zeros_like(vj_n1)))
        nonlinearsolver.solve(
            dj_vj_n1,
            model.compute_residual,
            model.solve_linearized_system,
            residual_args=residual_args,
            increment_args=increment_args,
        )

        dj_vj_n1 = np.reshape(dj_vj_n1, (2, -1))
        vj_n1 += dj_vj_n1[1]

        # Main loop to solve the ODE using the BDF method
        for it in range(1, nsteps + 1):

            start_step = time()
            logger.info(f"Time marching step: {it}/{nsteps}")
            self.update_manager.increment_step()

            if incremental_type == "alpha":
                # Get delta time
                dt = time_list[it] - time_list[it - 1]
                # Predict current temperature
                dj_n1 = np.copy(temperature_list[it - 1])
            elif incremental_type == "bdf":
                # Get current order and ensure it does not exceed norder
                currorder = min(it, norder)
                alpha = select_parameter(currorder)
                # Get values of last steps
                d_list = select_snapshots(it, temperature_list, currorder)
                # Predict current temperature
                dj_n1 = predict_temperature(d_list, currorder)

            # Update
            dj_n1 += (1 - alpha) * dt * vj_n1

            # Predict values of new step
            Fext = external_force_list[it]
            vj_n1 = np.zeros_like(Fext)

            residual_args = {
                "external_force": Fext,
                "scalar_coefs": (1.0, 1.0),
            }
            increment_args = {
                "scalar_coefs": (1.0 / (alpha * dt), 1),
                "flux_factor": 1.0 / (alpha * dt),
                "linear_solver_backend": self.linear_solver,
            }
            dj_vj_n1 = np.hstack((dj_n1, vj_n1))
            nonlinearsolver.solve(
                dj_vj_n1,
                model.compute_residual,
                model.solve_linearized_system,
                residual_args=residual_args,
                increment_args=increment_args,
            )
            dj_vj_n1 = np.reshape(dj_vj_n1, (2, -1))
            temperature_list[it] = np.copy(dj_vj_n1[0])
            vj_n1 = np.copy(dj_vj_n1[1])

            logger.info(f"Time-step {it} in {time() - start_step:.2e} seconds")

        logger.info(f"Solve heat transfer problem in {time() - start:.2e} seconds")


class SpaceTimeHeatTransfer(Physics):
    def __init__(self, **solver_args):
        super().__init__(**solver_args)

    def solve(
        self,
        model: Union[BaseSingleModel, BaseMultiModel],
        temperature: np.ndarray,
        external_force: np.ndarray,
    ):
        assert temperature.ndim == 1 and external_force.ndim == 1
        external_force = external_force.copy()
        clear_bcs(model, temperature, external_force)
        model.set_update_manager(self.update_manager)

        logger.info(f"Space-time heat transfer solver")
        start = time()
        nonlinearsolver = self.nonlinear_solver
        if self.input_args.auto_outer_tolerance:
            outer_tolerance_generator = OuterToleranceSetter(
                outer_tolerance_type="mesh_dependent",
                outer_tolerance_args={"default": self.input_args.tolerance_nonlinear},
            )
            nonlinearsolver.update(tolerance=outer_tolerance_generator.compute(model))
        spacetime_type: str = self.input_args.inner_tolerance_type.split()[-1]
        output = nonlinearsolver.solve(
            temperature,
            model.compute_residual,
            model.solve_linearized_system,
            residual_args={"external_force": external_force},
            increment_args={
                "spacetime_type": spacetime_type,
                "linear_solver_backend": self.linear_solver,
            },
        )

        self.output_args = output
        logger.info(f"Solve heat transfer problem in {time() - start:.2e} seconds")
