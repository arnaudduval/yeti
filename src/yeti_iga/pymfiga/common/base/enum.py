from enum import Enum
from typing import List


class TargetPriority(Enum):
    MASTER = 0
    SLAVE = 1


class ParametricDirection(Enum):
    XI = 0
    ETA = 1
    NU = 2
    TAU = 3
    ALL = -1

    @staticmethod
    def get_integration_dirs(
        idx_direction: "ParametricDirection", ndim: int
    ) -> List["ParametricDirection"]:
        all_dirs = [
            ParametricDirection.XI,
            ParametricDirection.ETA,
            ParametricDirection.NU,
        ]
        fixed_dir = all_dirs[idx_direction.value]
        available_dirs = all_dirs[:ndim]
        integration_dirs = [d for d in available_dirs if d != fixed_dir]
        return integration_dirs


class BoundarySide(Enum):
    MIN = 0  # face 0 in [0,1]
    MAX = 1  # face 1 in [0,1]
    BOTH = -1


class DOF(Enum):
    T = 0  # Temperature
    UX = 1  # Displacement X
    UY = 2  # Displacement Y
    UZ = 3  # Displacement Z
    ROTX = 4  # Rotation X
    ROTY = 5  # Rotation Y
    ROTZ = 6  # Rotation Z
    LAG1 = 7  # Lagrange Multiplier 1
    LAG2 = 8  # Lagrange Multiplier 2
    ALL = -1


class Constants:
    INFTY = 1e18
    HIGH = 1e8
    PENALTY = 1e2
    SMALL = 1e-4
    TINY = 1e-8
    SAFEGUARD = 1e-14
