from yeti_iga.pymfiga.common.base.enum import Constants
from dataclasses import dataclass


@dataclass
class SolverArgs:
    tolerance: float
    maxiters: int

    def __post_init__(self):
        """Validate arguments after initialization."""
        if self.tolerance < Constants.SAFEGUARD or self.tolerance >= 1:
            # It should not be lees than the machine error (in practice around 1e-14)
            # and it could not be greater than 1 (an error of 100%)
            raise ValueError(
                f"Relative tolerance must be > 0 and < 1, got {self.tolerance}"
            )
        if self.maxiters <= 0:
            # It should a be a positve number of iterations
            raise ValueError(f"Max iterations must be > 0, got {self.maxiters}")
