from yeti_iga.pymfiga.iga.fastdiagonalization import SingleFD
from yeti_iga.pymfiga.iga.single_model.cls.cspace import SingleSpatialModel as iga_model
from .core import Physics
from time import time
import numpy as np
import logging

logger = logging.getLogger("SRC.PHYSICS")


class L2projection(Physics):
    def __init__(self, **solver_args):
        super().__init__(**solver_args)

    def solve(self, model: iga_model, u_at_quadpts: np.ndarray) -> np.ndarray:

        assert isinstance(model, iga_model) and model.TypeProblemToBeApplied == "IGA"
        assert isinstance(u_at_quadpts, np.ndarray) and 0 < u_at_quadpts.ndim < 3

        start = time()

        def mass(x_in):
            x_out = model.part.operator_engine.compute_mf_scalar_u_v(
                model.part.quadrule_list,
                model.part.det_jac,
                x_in,
                allow_lumping=False,
                nurbs_weights=model.part.nurbs_weights,
            )
            return x_out

        has_to_ravel = False
        if u_at_quadpts.ndim == 1:
            has_to_ravel = True
            u_at_quadpts = np.atleast_2d(u_at_quadpts)
        nr = np.shape(u_at_quadpts)[0]

        prop: np.ndarray = u_at_quadpts * model.part.det_jac
        array_in = model.operator_engine.assemble_scalar_u_force(
            model.part.quadrule_list, prop, nurbs_weights=model.part.nurbs_weights
        )
        array_in = np.reshape(array_in, (nr, -1))

        fastdiag = SingleFD()
        fastdiag.compute_space_eigendecomposition(
            model.part.quadrule_list, np.zeros((1, model.ndim, 2))
        )
        fastdiag.update_space_eigenvalues(scalar_coefs=[1, 0])

        self.clear()
        self.linear_solver.verbose = False
        array_out = np.zeros_like(array_in)
        for i in range(nr):
            output = self.linear_solver.solve(
                mass,
                array_in[i],
                Pfun=fastdiag.apply_spatial_preconditioner,
            )
            array_out[i] = output["sol"]
        logger.info(f"Solve L2 projection problem in {time() - start:.2e} seconds")
        self.clear()
        return np.ravel(array_out) if has_to_ravel else array_out
