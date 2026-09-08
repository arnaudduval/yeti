from yeti_iga.pymfiga.common.base.cls import BaseSingleModel
from yeti_iga.pymfiga.common.material import Material
from yeti_iga.pymfiga.common.numerics.solvers import UpdateManager
from yeti_iga.pymfiga.iga.boundary import BoundaryCondition
from yeti_iga.pymfiga.iga.geometry import SinglePatch
from yeti_iga.pymfiga.iga.fastdiagonalization import SingleFD
from yeti_iga.pymfiga.iga.mass_preconditioner import ScaledMassPreconditioner
import numpy as np
import logging

logger = logging.getLogger("SRC.IGA.MODEL")


class SingleModel(BaseSingleModel):

    _update_manager = _preconditioner = _scaled_mass_preconditioner = None

    def __init__(
        self, material: Material, part: SinglePatch, boundary: BoundaryCondition
    ):
        self.TypeProblemToBeApplied = "IGA"
        self._set_material(material)
        self._set_part(part)
        self._set_boundary(boundary)
        logger.info("Single model initialized")

    @property
    def ndim(self):
        "Parametric dimension"
        return self.part.ndim  # parametric dimension

    @property
    def part(self):
        return self._part

    @property
    def boundary(self):
        return self._boundary

    @property
    def material(self):
        return self._material

    @property
    def preconditioner(self):
        if not isinstance(self._preconditioner, SingleFD):
            fd = SingleFD()
            fd.compute_space_eigendecomposition(
                self.part.quadrule_list, self.boundary.table_dirichlet
            )
            if hasattr(self, "time"):
                if isinstance(self.time, SinglePatch):
                    fd.compute_time_schurdecomposition(self.time.quadrule_list[0])
            logger.info(repr(fd))
            self._preconditioner = fd
        return self._preconditioner
    
    @property
    def scaled_mass_preconditioner(self):
        """
        Return the scaled-mass preconditioner.

        The preconditioner is built only once and then reused for subsequent
        calls. It is initialized from the current model.
        """
        if not isinstance(self._scaled_mass_preconditioner, ScaledMassPreconditioner):
            self._scaled_mass_preconditioner = ScaledMassPreconditioner()
            self._scaled_mass_preconditioner.compute(self)

        return self._scaled_mass_preconditioner


    @property
    def constraint_nodes(self):
        return self._constraint_nodes

    @property
    def free_nodes(self):
        return self._free_nodes

    @property
    def update_manager(self):
        if not isinstance(self._update_manager, UpdateManager):
            self._update_manager = UpdateManager()
        return self._update_manager

    def set_update_manager(self, manager):
        if isinstance(manager, UpdateManager):
            self._update_manager = manager

    def _set_material(self, material: Material):
        assert isinstance(material, Material)
        self._material: Material = material

    def _set_part(self, part: SinglePatch):
        assert isinstance(part, SinglePatch)
        self._part: SinglePatch = part
        self.operator_engine = part.operator_engine

    def _set_boundary(self, boundary: BoundaryCondition):
        assert isinstance(boundary, BoundaryCondition)
        self._boundary: BoundaryCondition = boundary
        output = boundary.select_nodes_for_solving()
        self._free_nodes, self._constraint_nodes = output

    def get_size_of_arrays(self):
        """Return the size of the arrays needed for the solution."""
        return self.nbvars * self.part.nbctrlpts_total

    def clear_bcs(self, array_in: np.ndarray):
        """Apply the boundary conditions by zeroing out the constrained degrees of freedom."""
        assert isinstance(self.constraint_nodes, list)
        array_in[self.constraint_nodes] = 0.0

    def get_free_and_constraint_nodes(self):
        """Return the indices of free and constrained nodes."""
        return self.free_nodes, self.constraint_nodes

    def has_DOFsBlocked_atLeastOnce(self):
        "Check if there are any degrees of freedom that are blocked (constrained) at least once."
        output = True
        for vars in range(self.nbvars):
            output *= np.any(self.boundary.table_dirichlet[vars])
            if output == False:
                return output
        return output

    def compute_residual(self, array_in, **kwargs):
        raise NotImplementedError()

    def solve_linearized_system(self, array_in, **kwargs):
        raise NotImplementedError()
