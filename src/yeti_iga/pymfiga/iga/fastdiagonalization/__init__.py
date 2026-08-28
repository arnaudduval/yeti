from .single_diagonalization import SingleFastDiagonalization as SingleFD
from .multi_diagonalization import MultiFastDiagonalization as MultiFD

# Choice of Lagrange preconditioner

# from .lagdense_diagonalization import LagrangeFastDiagonalization as LagrangeFD

from .lagmf_diagonalization import LagrangeFastDiagonalization as LagrangeFD

__all__ = ["SingleFD", "MultiFD", "LagrangeFD"]
