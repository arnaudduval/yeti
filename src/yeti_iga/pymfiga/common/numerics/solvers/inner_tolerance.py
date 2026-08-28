from yeti_iga.pymfiga.common.base.enum import Constants
from typing import Optional, Literal, Dict
from dataclasses import dataclass
import numpy as np
import logging

logger = logging.getLogger("SRC.SOLVER_PARAMETERS")


@dataclass
class InnerToleranceArgs:
    """Arguments for inner tolerance computation."""

    initial_tolerance: float = 0.5
    static: float = 0.1
    coefficient: float = 0.9
    exponential: float = 1.5
    default: float = Constants.TINY

    def __post_init__(self):
        """Validate arguments after initialization."""
        for key, value in self.__dict__.items():
            if not np.isscalar(value):
                raise ValueError(f"{key} must be a scalar value")
            if key == "default" and np.real(value) >= 1.0:
                raise ValueError(f"default must be < 1.0, got {value}")
            if key == "initial_tolerance" and np.real(value) >= 1.0:
                raise ValueError(f"Initial value must be < 1.0, got {value}")
            if key == "static" and np.real(value) >= 1.0:
                raise ValueError(f"Static value must be < 1.0, got {value}")


class InnerToleranceSetter:
    """
    Manages inner tolerance computation for iterative solvers.

    Supports four tolerance strategies:
    - exact picard: Fixed tolerance at default value
    - exact newton: Fixed tolerance at default value
    - inexact picard: Static threshold-based tolerance
    - inexact newton: Adaptive tolerance based on residual reduction
    """

    _VALID_TYPES = {
        "exact picard",
        "exact newton",
        "inexact picard",
        "inexact newton",
    }

    def __init__(
        self,
        inner_tolerance_type: Literal[
            "exact picard",
            "exact newton",
            "inexact picard",
            "inexact newton",
        ],
        inner_tolerance_args: Dict[str, float],
    ):
        """
        Initialize InnerTolerance.

        Args:
            inner_tolerance_type: Strategy for tolerance computation
            inner_tolerance_args: Dictionary of tolerance parameters (must include 'default')

        Raises:
            AssertionError: If tolerance type is invalid or args missing 'default'
        """
        self._tolerance_type = self._validate_tolerance_type(inner_tolerance_type)
        self._config = self._initialize_config(inner_tolerance_args)

    @property
    def tolerance_type(self) -> str:
        """Get current tolerance type."""
        return self._tolerance_type

    @property
    def config(self) -> InnerToleranceArgs:
        """Get current tolerance configuration."""
        return self._config

    def update(
        self,
        inner_tolerance_type: Optional[str] = None,
        **kwargs: float,
    ) -> None:
        """
        Update tolerance configuration.

        Args:
            inner_tolerance_type: New tolerance type (optional)
            **kwargs: Additional arguments to update (e.g., default, coefficient)
        """
        if inner_tolerance_type is not None:
            self._tolerance_type = self._validate_tolerance_type(inner_tolerance_type)

        for key, value in kwargs.items():
            if key in self.config.__dict__:
                if not np.isscalar(value):
                    raise ValueError(f"{key} must be a scalar value")
                if key == "default" and np.real(value) >= 1.0:
                    raise ValueError(f"default must be < 1.0, got {value}")
                setattr(self.config, key, float(np.real(value)))
            else:
                logger.warning(f"Unknown argument '{key}' ignored")

    def compute(
        self,
        norm_residual_new: Optional[float] = None,
        norm_residual_old: Optional[float] = None,
        inner_tolerance_old: Optional[float] = None,
    ) -> float:
        """
        Compute inner tolerance based on current strategy.

        Args:
            norm_residual_new: Current iteration residual norm
            norm_residual_old: Previous iteration residual norm
            inner_tolerance_old: Previous inner tolerance

        Returns:
            Computed inner tolerance value
        """
        if self.tolerance_type in {"exact picard", "exact newton"}:
            return self._compute_exact()
        elif self.tolerance_type == "inexact picard":
            return self._compute_inexact_picard()
        elif self.tolerance_type == "inexact newton":
            return self._compute_inexact_newton(
                norm_residual_new, norm_residual_old, inner_tolerance_old
            )
        else:
            raise NotImplementedError(f"Unknown tolerance type: {self.tolerance_type}")

    def _validate_tolerance_type(self, value: str) -> str:
        """Validate and normalize tolerance type."""
        if not isinstance(value, str):
            raise TypeError(f"Tolerance type must be string, got {type(value)}")

        normalized = value.lower()
        if normalized not in self._VALID_TYPES:
            raise ValueError(
                f"Invalid tolerance type '{value}'. "
                f"Must be one of {self._VALID_TYPES}"
            )
        return normalized

    def _initialize_config(self, args: Dict[str, float]) -> InnerToleranceArgs:
        """Initialize tolerance arguments from dictionary."""
        if not isinstance(args, dict):
            raise TypeError(f"Arguments must be dict, got {type(args)}")
        if "default" not in args:
            raise ValueError("Arguments must include 'default' value")

        return InnerToleranceArgs(**args)

    def _compute_exact(self) -> float:
        """Compute tolerance for exact methods (fixed at default)."""
        tolerance = self.config.default
        logger.debug(f"Exact tolerance: {tolerance:.2e}")
        return tolerance

    def _compute_inexact_picard(self) -> float:
        """Compute tolerance for inexact Picard method."""
        tolerance = max(
            self.config.default, min(self.config.initial_tolerance, self.config.static)
        )
        logger.debug(f"Inexact Picard tolerance: {tolerance:.2e}")
        return tolerance

    def _compute_inexact_newton(
        self,
        norm_residual_new: Optional[float],
        norm_residual_old: Optional[float],
        inner_tolerance_old: Optional[float],
    ) -> float:
        """
        Compute tolerance for inexact Newton method using Eisenstat-Walker.

        Uses adaptive strategy based on residual reduction rate.
        """
        threshold = self.config.initial_tolerance

        if (
            isinstance(norm_residual_new, (float, np.floating))
            and isinstance(norm_residual_old, (float, np.floating))
            and isinstance(inner_tolerance_old, (float, np.floating))
        ):
            gamma = self.config.coefficient
            omega = self.config.exponential

            ratio = norm_residual_new / max(norm_residual_old, Constants.TINY)
            eps_choice1 = gamma * np.power(ratio, omega)
            eps_choice2 = gamma * np.power(inner_tolerance_old, omega)

            # Use safeguarded choice
            if eps_choice2 > 0.1:
                threshold = min(eps_choice1, eps_choice2)
            else:
                threshold = eps_choice1

        tolerance = np.clip(
            threshold, self.config.default, self.config.initial_tolerance
        )
        logger.debug(f"Inexact Newton tolerance: {tolerance:.2e}")
        return tolerance
