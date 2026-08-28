from yeti_iga.pymfiga.common.base.cls import BaseSingleModel, BaseMultiModel
from yeti_iga.pymfiga.common.numerics.solvers import (
    LinearSolver,
    NonLinearSolver,
    InnerToleranceSetter,
    UpdateManager,
)
from typing import List, Dict, Union, Literal
from dataclasses import dataclass
import numpy as np


def clear_bcs(model: Union[BaseSingleModel, BaseMultiModel], *arrays: np.ndarray):
    constraint_nodes = model.get_free_and_constraint_nodes()[-1]
    for arr in arrays:
        arr[..., constraint_nodes] = 0.0


@dataclass
class PhysicsArgs:
    # Default values
    linear_solver_type: Literal["cg", "bicgstab", "gmres"] = "gmres"
    maxiters_linear: int = 100
    tolerance_linear: float = 1e-8
    maxiters_nonlinear: int = 10
    tolerance_nonlinear: float = 1e-6
    allow_acceleration: bool = False
    allow_line_search: bool = False
    auto_outer_tolerance: bool = False
    inner_tolerance_type: Literal[
        "exact picard",
        "exact newton",
        "inexact picard",
        "inexact newton",
    ] = "exact picard"


class Physics:

    # Internal variable
    _linear_solver = None
    _nonlinear_solver = None
    _inner_tolerance_setter = None
    _update_manager = None

    def __init__(self, **solver_args):
        # Public variables
        self._input_args = PhysicsArgs(**solver_args)
        self._output_args = {}

    def clear(self):
        self._linear_solver = None
        self._nonlinear_solver = None
        self._inner_tolerance_setter = None

    @property
    def input_args(self):
        return self._input_args

    @property
    def output_args(self):
        return self._output_args

    @output_args.setter
    def output_args(self, args):
        assert isinstance(args, dict)
        linear_tolerance_list: List[float] = args.get("linear_tolerance_list", [])
        nonlinear_rate_list: List[float] = args.get("nonlinear_rate_list", [])
        nonlinear_residual_list: List[float] = args.get("nonlinear_residual_list", [])
        nonlinear_time_list: List[float] = args.get("nonlinear_time_list", [])
        solution_history_list: Dict[str, np.ndarray] = args.get(
            "solution_history_list", {}
        )
        self._output_args = {
            "linear_tolerance_list": linear_tolerance_list,
            "nonlinear_rate_list": nonlinear_rate_list,
            "nonlinear_residual_list": nonlinear_residual_list,
            "nonlinear_time_list": nonlinear_time_list,
            "solution_history_list": solution_history_list,
        }

    @property
    def linear_solver(self) -> LinearSolver:
        if not isinstance(self._linear_solver, LinearSolver):
            linear_type = self.input_args.linear_solver_type
            maxiters = self.input_args.maxiters_linear
            tolerance = self.input_args.tolerance_linear
            linear_solver = LinearSolver(
                maxiters=maxiters, tolerance=tolerance, linear_type=linear_type
            )
            self._linear_solver = linear_solver
        return self._linear_solver

    @property
    def nonlinear_solver(self):
        if not isinstance(self._nonlinear_solver, NonLinearSolver):
            # Initialize stopping criteria parameters
            maxiters = self.input_args.maxiters_nonlinear
            tolerance = self.input_args.tolerance_nonlinear
            allow_acceleration = self.input_args.allow_acceleration
            allow_linesearch = self.input_args.allow_line_search

            # Solve using a nonlinear solver
            nonlinearsolver = NonLinearSolver(
                maxiters=maxiters,
                tolerance=tolerance,
                allow_acceleration=allow_acceleration,
                allow_line_search=allow_linesearch,
            )
            nonlinearsolver.update(
                inner_tolerance_setter=self.inner_tolerance_setter,
                update_manager=self.update_manager,
            )
            self._nonlinear_solver = nonlinearsolver
        return self._nonlinear_solver

    @property
    def inner_tolerance_setter(self):
        if not isinstance(self._inner_tolerance_setter, InnerToleranceSetter):
            inner_tolerance_type = self.input_args.inner_tolerance_type
            default_tolerance = self.input_args.tolerance_linear
            inner_tolerance = InnerToleranceSetter(
                inner_tolerance_type=inner_tolerance_type,
                inner_tolerance_args={
                    "default": default_tolerance,
                },
            )
            self._inner_tolerance_setter = inner_tolerance
        return self._inner_tolerance_setter

    @property
    def update_manager(self):
        if not isinstance(self._update_manager, UpdateManager):
            self._update_manager = UpdateManager()
        return self._update_manager
