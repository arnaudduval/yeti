from yeti_iga.pymfiga.common.base.enum import Constants
from yeti_iga.pymfiga.common.base.cls import BaseSingleModel, BaseMultiModel
from typing import Union, Optional, Literal, Dict
from dataclasses import dataclass
import numpy as np
import logging

logger = logging.getLogger("SRC.SOLVER_PARAMETERS")


@dataclass
class OuterToleranceArgs:
    """Configuration for outer tolerance computation."""

    default: float = Constants.SMALL
    scaling_factor: float = 0.1

    def __post_init__(self):
        """Validate configuration after initialization."""
        if self.default <= 0:
            raise ValueError(f"default must be > 0, got {self.default}")
        if self.scaling_factor <= 0:
            raise ValueError(f"scaling_factor must be > 0, got {self.scaling_factor}")


class OuterToleranceSetter:
    """
    Manages outer tolerance computation for iterative solvers.

    Supports two tolerance strategies:
    - mesh free: Fixed tolerance at default value
    - mesh dependent: Adaptive tolerance based on mesh parameters and polynomial degree

    The mesh-dependent tolerance is computed as:
        tolerance = max(default, factor * h^p)
    where h is the mesh parameter and p is the polynomial degree.
    """

    VALID_TYPES = {"mesh_free", "mesh_dependent"}

    def __init__(
        self,
        outer_tolerance_type: Literal["mesh_free", "mesh_dependent"],
        outer_tolerance_args: Dict[str, float],
    ):
        """
        Initialize OuterTolerance.

        Args:
            outer_tolerance_type: Strategy for tolerance computation
            outer_tolerance_args: Configuration dictionary (must include 'default')

        Raises:
            ValueError: If tolerance type is invalid or args missing 'default'
            TypeError: If arguments are of wrong type
        """
        self._tolerance_type = self._validate_tolerance_type(outer_tolerance_type)
        self._config = self._initialize_config(outer_tolerance_args)

    @property
    def tolerance_type(self) -> str:
        """Get current tolerance type."""
        return self._tolerance_type

    @property
    def config(self) -> OuterToleranceArgs:
        """Get current tolerance configuration."""
        return self._config

    def compute(
        self,
        model: Optional[Union[BaseSingleModel, BaseMultiModel]] = None,
    ) -> float:
        """
        Compute outer tolerance based on current strategy.

        Args:
            model: Model instance (required for mesh-dependent tolerance)

        Returns:
            Computed outer tolerance value

        Raises:
            ValueError: If model is required but not provided
            NotImplementedError: If tolerance type is unknown
        """
        if self.tolerance_type == "mesh_free":
            return self._compute_mesh_free()
        elif self.tolerance_type == "mesh_dependent":
            return self._compute_mesh_dependent(model)
        else:
            raise NotImplementedError(f"Unknown tolerance type: {self.tolerance_type}")

    def _validate_tolerance_type(self, value: str) -> str:
        """Validate and normalize tolerance type."""
        if not isinstance(value, str):
            raise TypeError(f"Tolerance type must be string, got {type(value)}")

        normalized = value.lower()
        if normalized not in self.VALID_TYPES:
            raise ValueError(
                f"Invalid tolerance type '{value}'. "
                f"Must be one of {self.VALID_TYPES}"
            )
        return normalized

    def _initialize_config(self, args: Dict[str, float]) -> OuterToleranceArgs:
        """Initialize tolerance configuration from dictionary."""
        if not isinstance(args, dict):
            raise TypeError(f"Arguments must be dict, got {type(args)}")
        if "default" not in args:
            raise ValueError("Arguments must include 'default' value")

        # Extract known parameters with defaults
        config_params = {
            "default": args["default"],
            "scaling_factor": args.get("scaling_factor", 0.25),
        }

        return OuterToleranceArgs(**config_params)

    def _compute_mesh_free(self) -> float:
        """
        Compute mesh-free tolerance (fixed at default value).

        Returns:
            Default tolerance value
        """
        tolerance = self.config.default
        logger.debug(f"Mesh-free tolerance: {tolerance:.2e}")
        return tolerance

    def _compute_mesh_dependent(
        self,
        model: Optional[Union[BaseSingleModel, BaseMultiModel]],
    ) -> float:
        """
        Compute mesh-dependent tolerance based on mesh parameters.

        Args:
            model: Model instance containing mesh information

        Returns:
            Computed tolerance based on mesh size and polynomial degree

        Raises:
            ValueError: If model is None or of wrong type
        """
        if model is None:
            raise ValueError(
                "Model must be provided for mesh-dependent tolerance computation"
            )

        if not isinstance(model, (BaseSingleModel, BaseMultiModel)):
            raise TypeError(
                f"Model must be BaseSingleModel or BaseMultiModel, "
                f"got {type(model)}"
            )

        if isinstance(model, BaseSingleModel):
            tolerance = self._compute_single_model_tolerance(model)
        else:
            tolerance = self._compute_multi_model_tolerance(model)

        return tolerance

    def _compute_single_model_tolerance(
        self,
        model: BaseSingleModel,
    ) -> float:
        """
        Compute tolerance for a single model.

        Args:
            model: Single model instance

        Returns:
            Tolerance based on global mesh parameter and polynomial degree
        """
        mesh_param, degree = self._extract_mesh_parameters(model)

        tolerance = self.config.scaling_factor * np.power(mesh_param, degree)
        tolerance = max(self.config.default, tolerance)

        logger.debug(
            f"Single model tolerance: {tolerance:.2e} "
            f"(h={mesh_param:.2e}, p={degree})"
        )

        return float(tolerance)

    def _compute_multi_model_tolerance(
        self,
        model: BaseMultiModel,
    ) -> float:
        """
        Compute tolerance for a multi-model system.

        Takes the maximum tolerance across all sub-models to ensure
        convergence in the coarsest/highest-degree model.

        Args:
            model: Multi-model instance

        Returns:
            Maximum tolerance across all sub-models
        """
        sub_models = model.recover_models()

        if not sub_models:
            logger.warning("Multi-model has no sub-models, using default tolerance")
            return self.config.default

        tolerances = [
            self._compute_single_model_tolerance(sub_model) for sub_model in sub_models
        ]

        max_tolerance = float(np.max(tolerances))

        logger.debug(
            f"Multi-model tolerance: {max_tolerance:.2e} "
            f"(max of {len(tolerances)} sub-models)"
        )

        return max_tolerance

    def _extract_mesh_parameters(
        self,
        model: BaseSingleModel,
    ) -> tuple[float, int]:
        """
        Extract global mesh parameter and polynomial degree from model.

        Computes mesh parameters for both spatial (part) and temporal (time)
        domains, taking the maximum mesh parameter and minimum degree.

        Args:
            model: Single model instance

        Returns:
            Tuple of (global_mesh_parameter, global_degree)
        """
        part_patch = model.part
        time_patch = model.time

        # Extract spatial mesh parameters
        mp_part, deg_part = part_patch.compute_global_mesh_parameter()

        # Extract temporal mesh parameters (if available)
        if time_patch is not None:
            mp_time, deg_time = time_patch.compute_global_mesh_parameter()
        else:
            mp_time, deg_time = mp_part, deg_part

        # Use maximum mesh parameter (coarsest mesh) and minimum degree
        global_mesh_param = max(mp_part, mp_time)
        global_degree = min(deg_part, deg_time)

        logger.debug(
            f"Mesh parameters: spatial=(h={mp_part:.2e}, p={deg_part}), "
            f"temporal=(h={mp_time:.2e}, p={deg_time}), "
            f"global=(h={global_mesh_param:.2e}, p={global_degree})"
        )

        return float(global_mesh_param), int(global_degree)

    def __repr__(self) -> str:
        """String representation of OuterTolerance."""
        return (
            f"OuterTolerance("
            f"type='{self.tolerance_type}', "
            f"config={self.config})"
        )
