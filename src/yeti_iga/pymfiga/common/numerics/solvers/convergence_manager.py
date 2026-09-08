from yeti_iga.pymfiga.common.base.enum import Constants
from typing import Dict, Optional, Union, Literal
from abc import ABC, abstractmethod
from dataclasses import dataclass
import numpy as np
import logging

logger = logging.getLogger("SRC.SOLVER_PARAMETERS")


@dataclass
class ConvergenceCriteria:
    """Configuration for convergence criteria."""

    relative_tolerance: float
    max_iterations: int
    absolute_tolerance: float = Constants.SAFEGUARD

    def __post_init__(self):
        """Validate criteria after initialization."""
        if self.relative_tolerance <= 0 or self.relative_tolerance >= 1:
            raise ValueError(
                f"Relative tolerance must be > 0 and < 1, got {self.relative_tolerance}"
            )
        if self.max_iterations <= 0:
            raise ValueError(f"Max iterations must be > 0, got {self.max_iterations}")
        if self.absolute_tolerance < 0:
            raise ValueError(
                f"Absolute tolerance must be >= 0, got {self.absolute_tolerance}"
            )


class ConvergenceManager:
    """
    Manages convergence criteria for iterative algorithms.

    Tracks multiple convergence parameters (absolute error, relative error, iterations)
    and determines when an iterative process has converged.
    """

    class Parameter(ABC):
        """Abstract base class for convergence parameters."""

        def __init__(self, name: str = ""):
            self._name = name
            self._reference: float = 0.0
            self._current: Optional[float] = None

        @property
        def name(self) -> str:
            """Get parameter name."""
            return self._name

        @property
        def reference(self) -> float:
            """Get reference (tolerance/threshold) value."""
            return self._reference

        @reference.setter
        def reference(self, value: Union[int, float, np.floating]) -> None:
            """Set reference value (must be non-negative)."""
            if not isinstance(value, (int, float, np.floating)):
                raise TypeError(f"Reference must be numeric, got {type(value)}")
            if value < 0:
                raise ValueError(f"Reference must be non-negative, got {value}")
            self._reference = float(value)

        @property
        def current(self) -> Optional[float]:
            """Get current value."""
            return self._current

        @current.setter
        def current(self, value: Union[int, float, np.floating]) -> None:
            """Set current value (must be non-negative)."""
            if not isinstance(value, (int, float, np.floating)):
                raise TypeError(f"Current must be numeric, got {type(value)}")
            if value < 0:
                raise ValueError(f"Current must be non-negative, got {value}")
            self._current = float(value)

        def is_defined(self) -> bool:
            """Check if current value has been set."""
            return self._current is not None

        @abstractmethod
        def has_converged(self) -> bool:
            """Check if convergence criterion is satisfied."""
            pass

        def reset(self) -> None:
            """Reset current value to None."""
            self._current = None

        def __repr__(self) -> str:
            return (
                f"{self.__class__.__name__}("
                f"name='{self._name}', "
                f"reference={self._reference:.2e}, "
                f"current={self._current:.2e if self._current is not None else None})"
            )

    class AbsoluteError(Parameter):
        """Absolute error convergence criterion."""

        def has_converged(self) -> bool:
            """Check if absolute error is below threshold."""
            if self.current is None:
                return False
            threshold = max(self.reference, Constants.SAFEGUARD)
            converged = self.current < threshold
            if converged:
                logger.debug(
                    f"Absolute error converged: {self.current:.2e} < {threshold:.2e}"
                )
            return converged

    class RelativeError(Parameter):
        """Relative error convergence criterion."""

        def has_converged(self) -> bool:
            """Check if relative error is below threshold."""
            if self.current is None:
                return False
            threshold = max(self.reference, Constants.SAFEGUARD)
            converged = self.current < threshold
            if converged:
                logger.debug(
                    f"Relative error converged: {self.current:.2e} < {threshold:.2e}"
                )
            return converged

    class Iteration(Parameter):
        """Iteration count criterion (stops when max iterations reached)."""

        def has_converged(self) -> bool:
            """
            Check if max iterations exceeded.

            Raises:
                RuntimeError: If current iteration exceeds maximum

            Returns:
                Always False (iterations don't indicate convergence)
            """
            if self.current is None:
                return False

            if self.current > self.reference:
                logger.warning(
                    f"Maximum iterations ({self.reference:.0f}) reached at "
                    f"iteration {self.current:.0f}. Convergence not achieved."
                )
                raise RuntimeError(
                    f"Maximum iterations ({int(self.reference)}) exceeded. "
                    "Algorithm did not converge."
                )
            return False

    # Parameter type mapping
    _PARAMETER_TYPES = {
        "absolute_error": AbsoluteError,
        "relative_error": RelativeError,
        "iteration": Iteration,
    }

    def __init__(
        self,
        refe_relative: float,
        refe_iteration: int,
        refe_absolute: float = 0.0,
    ):
        """
        Initialize ConvergenceManager.

        Args:
            refe_relative: Relative error tolerance
            refe_iteration: Maximum number of iterations
            refe_absolute: Absolute error tolerance (default: 0.0)
        """
        criteria = ConvergenceCriteria(
            relative_tolerance=refe_relative,
            max_iterations=refe_iteration,
            absolute_tolerance=refe_absolute,
        )

        self._parameters: Dict[str, ConvergenceManager.Parameter] = {
            "absolute_error": self.AbsoluteError("absolute_error"),
            "relative_error": self.RelativeError("relative_error"),
            "iteration": self.Iteration("iteration"),
        }

        # Set references
        self._parameters["absolute_error"].reference = max(
            criteria.absolute_tolerance, Constants.SAFEGUARD
        )
        self._parameters["relative_error"].reference = criteria.relative_tolerance
        self._parameters["iteration"].reference = criteria.max_iterations

    @property
    def parameters(self) -> Dict[str, Parameter]:
        """Get all convergence parameters."""
        return self._parameters

    def add_criterion(
        self,
        name: str,
        criterion_type: Literal["absolute_error", "relative_error"],
        reference: Optional[float] = None,
    ) -> None:
        """
        Add a custom convergence criterion.

        Args:
            name: Unique name for the new criterion
            criterion_type: Type of criterion ("absolute_error" or "relative_error")
            reference: Custom reference value (uses existing if None)

        Raises:
            ValueError: If criterion type is invalid or name already exists
        """
        if name in self.parameters:
            raise ValueError(f"Criterion '{name}' already exists")

        if criterion_type not in self._PARAMETER_TYPES:
            raise ValueError(
                f"Invalid criterion type '{criterion_type}'. "
                f"Must be one of {list(self._PARAMETER_TYPES.keys())}"
            )

        # Create new parameter
        param_class = self._PARAMETER_TYPES[criterion_type]
        new_param: ConvergenceManager.Parameter = param_class(name)

        # Set reference value
        if reference is not None:
            new_param.reference = reference
        else:
            new_param.reference = self.parameters[criterion_type].reference

        self._parameters[name] = new_param
        logger.debug(f"Added criterion '{name}' of type '{criterion_type}'")

    def update(
        self,
        curr_absolute: Optional[float] = None,
        curr_relative: Optional[float] = None,
        curr_iteration: Optional[int] = None,
        **kwargs,
    ):
        """
        Update current values of convergence parameters.

        Args:
            curr_absolute: Current absolute error
            curr_relative: Current relative error
            curr_iteration: Current iteration number
            **kwargs: Additional custom parameter values
        """
        updates = {
            "absolute_error": curr_absolute,
            "relative_error": curr_relative,
            "iteration": curr_iteration,
        }

        # Update standard parameters
        for name, value in updates.items():
            if value is not None:
                self._parameters[name].current = value

        # Update custom parameters
        for name, value in kwargs.items():
            if name in self._parameters:
                if isinstance(value, (int, float, np.floating)) and value >= 0:
                    self._parameters[name].current = value
                else:
                    logger.warning(
                        f"Invalid value for parameter '{name}': {value}. Skipping."
                    )
            else:
                logger.warning(f"Unknown parameter '{name}'. Skipping.")

    def clear(self) -> None:
        """Reset all current values to None."""
        for param in self.parameters.values():
            param.reset()
        logger.debug("All convergence parameters cleared")

    def has_converged(self) -> bool:
        """
        Check if any convergence criterion is satisfied.

        Returns:
            True if at least one non-iteration criterion has converged

        Raises:
            RuntimeError: If maximum iterations exceeded
            AssertionError: If no parameters are defined
        """
        defined_params = [p for p in self.parameters.values() if p.is_defined()]

        if not defined_params:
            raise AssertionError("No convergence parameters have been set")

        current_iteration = self.parameters["iteration"].current

        # Check all criteria
        for name, param in self.parameters.items():
            if not param.is_defined():
                continue

            # Skip iteration criterion (it doesn't indicate convergence)
            if name == "iteration":
                # This will raise RuntimeError if max iterations exceeded
                param.has_converged()
                continue

            # Check if this criterion has converged
            if param.has_converged():
                iter_str = (
                    f" after {int(current_iteration)} iterations"
                    if current_iteration
                    else ""
                )
                logger.debug(
                    f"Convergence achieved: {name} = {param.current:.2e}{iter_str}"
                )
                return True

        return False

    def get_status(self) -> Dict[str, Dict[str, Optional[float]]]:
        """
        Get current status of all parameters.

        Returns:
            Dictionary with parameter names and their reference/current values
        """
        return {
            name: {
                "reference": param.reference,
                "current": param.current,
                "converged": param.has_converged() if param.is_defined() else None,
            }
            for name, param in self.parameters.items()
        }

    def __repr__(self) -> str:
        """String representation of ConvergenceManager."""
        status = self.get_status()
        return f"ConvergenceManager({status})"
