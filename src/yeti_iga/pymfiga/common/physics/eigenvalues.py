from yeti_iga.pymfiga.common.base.enum import Constants
from yeti_iga.pymfiga.common.numerics.solvers import LinearSolver
from yeti_iga.pymfiga.iga.single_model.cls.cspace import SingleSpatialModel as iga_model
from yeti_iga.pymfiga.fem.model.core import SingleModel as fem_model


from .core import Physics
from scipy.sparse import linalg as scsplin
from scipy import sparse as sp
from typing import Tuple, Literal, Union
from copy import deepcopy
from time import time
import numpy as np
import logging

logger = logging.getLogger("SRC.PHYSICS")


class EigenProblem(Physics):
    # TODO: refactorize this class to be more general
    def __init__(self, **solver_args):
        super().__init__(**solver_args)

    @staticmethod
    def _verify_model(model: Union[fem_model, iga_model]):
        if not isinstance(model, (fem_model, iga_model)):
            raise ValueError("Model is not FEM or IGA")
        if not hasattr(model, "compute_mf_mass"):
            raise ValueError("Model does not have compute_mf_mass")
        if not hasattr(model, "compute_mf_stiffness"):
            raise ValueError("Model does not have compute_mf_stiffness")

    @staticmethod
    def _select_case(model: Union[fem_model, iga_model]):
        fem_or_iga = model.TypeProblemToBeApplied
        if fem_or_iga == "FEA":
            return 0
        
        elif fem_or_iga == "IGA":
            if hasattr(model, "mass_type"):
                if model.mass_type == "diagonal_mass":
                    return 2
                elif model.mass_type == "consistent_mass":
                    return 1

            return 2 if hasattr(model, "diagonal_mass") else 1

            

    @staticmethod
    def _preconditioner_case_0(model: fem_model):

        if not hasattr(model, "assemble_mass"):
            raise Warning("Problem does not have assemble_mass")

        FREE = model.get_free_and_constraint_nodes()[0]
        mat: sp.csr_array = model.__getattribute__("assemble_mass")()
        ilu = scsplin.spilu(mat[np.ix_(FREE, FREE)])

        def Pfun(x_in):
            return ilu.solve(x_in)

        return lambda x: Pfun(x)

    @staticmethod
    def _preconditioner_case_1(model: iga_model):

        # Extract data from model
        N = model.get_size_of_arrays()
        FREE, CLEANDOD = model.get_free_and_constraint_nodes()

        # We define a 'good' preconditioners
        fastdiag = deepcopy(model.preconditioner)
        mass_corrector = [1.0] * fastdiag.space_preconditioner.nbDoFsPerNode
        if hasattr(model, "scalar_mean_mass"):

            mass_corrector = model.__getattribute__("scalar_mean_mass")
            
        fastdiag.add_scalar_space_time_correctors(mass_corrector=mass_corrector)
        fastdiag.update_space_eigenvalues(scalar_coefs=[1.0, 0.0])

        linear_solver = LinearSolver(
            maxiters=100,
            tolerance=Constants.TINY,
            cleandod=CLEANDOD,
            linear_type="cg",
            verbose=False,
        )

        def Pfun(x_in):
            # NOTE: Why just dont apply the preconditioner ?
            # eigs solver needs a good approximation of M^-1
            # and just the preconditioner is not enough
            # This is closer to ILU solver
            x_tmp = np.zeros(N)
            x_tmp[FREE] = x_in
            output = linear_solver.solve(
                model.__getattribute__("compute_mf_mass"),
                x_tmp,
                Pfun=fastdiag.apply_spatial_preconditioner,
            )
            x_out = output["sol"]
            return x_out[FREE]

        return lambda x: Pfun(x)

    @staticmethod
    def _preconditioner_case_scaled_mass(model: iga_model):

        N = model.get_size_of_arrays()
        FREE, CLEANDOD = model.get_free_and_constraint_nodes()

        scaled_mass = model.scaled_mass_preconditioner

        linear_solver = LinearSolver(
            maxiters=100,
            tolerance=Constants.TINY,
            cleandod=CLEANDOD,
            linear_type="cg",
            verbose=False,
        )

        def Pfun(x_in):
            x_tmp = np.zeros(N)
            x_tmp[FREE] = x_in

            output = linear_solver.solve(
                model.__getattribute__("compute_mf_mass"),
                x_tmp,
                Pfun=scaled_mass.apply_spatial_preconditioner,
            )

            x_out = output["sol"]
            return x_out[FREE]

        return lambda x: Pfun(x)
    

    @staticmethod
    def _preconditioner_case_2(model: iga_model):

        FREE = model.get_free_and_constraint_nodes()[0]
        DIAG = model.__getattribute__("compute_mf_mass")(
            np.ones(model.get_size_of_arrays())
        )[FREE]

        def Pfun(x_in):
            x_out = x_in / DIAG
            return x_out

        return lambda x: Pfun(x)
    

    def solve(
        self,
        model: Union[fem_model, iga_model],
        use_preconditioner: bool = True,
        preconditioner_type: Literal["fastdiag", "scaled_mass", None] = "fastdiag",
        solver: Literal["arpack", "power"] = "arpack",
        which: Literal["SM", "LM"] = "SM",
        k: int = 2,
        **kwargs,
    ) -> Tuple[np.ndarray, np.ndarray]:    

        logger.info("Eigen value problem")
        start = time()

        EigenProblem._verify_model(model)
        case_problem = EigenProblem._select_case(model)

        # Extract data from model
        N = model.get_size_of_arrays()
        FREE = model.get_free_and_constraint_nodes()[0]

        def mass(x_in):
            x_tmp = np.zeros(N)
            x_tmp[FREE] = x_in
            x_out = model.__getattribute__("compute_mf_mass")(x_tmp)
            return x_out[FREE]

        def stiff(x_in):
            x_tmp = np.zeros(N)
            x_tmp[FREE] = x_in
            x_out = model.__getattribute__("compute_mf_stiffness")(x_tmp)
            return x_out[FREE]

        if not use_preconditioner:
            preconditioner = None
        else:
            if case_problem == 0:
                if not isinstance(model, fem_model):
                    raise ValueError()
                preconditioner = EigenProblem._preconditioner_case_0(model)

            elif case_problem == 1:
                if not isinstance(model, iga_model):
                    raise ValueError()

                if preconditioner_type == "fastdiag":
                    preconditioner = EigenProblem._preconditioner_case_1(model)

                elif preconditioner_type == "scaled_mass":
                    preconditioner = EigenProblem._preconditioner_case_scaled_mass(model)

                elif preconditioner_type is None:
                    preconditioner = None

                else:
                    raise ValueError(f"Unknown preconditioner_type: {preconditioner_type}")

            elif case_problem == 2:
                if not isinstance(model, iga_model):
                    raise ValueError()
                # NOTE the diagonal is computed outside
                preconditioner = EigenProblem._preconditioner_case_2(model)

            else:
                raise NotImplementedError()
                
                
        if solver == "arpack":

            eigenvalues, eigenvectors = LinearSolver.eigs(
                N=len(FREE),
                Afun=stiff,
                Bfun=mass,
                Pfun=preconditioner,
                k=k,
                which=which,
                **kwargs,
            )

        elif solver == "power":

            eigenvalues, eigenvectors = LinearSolver.power_iteration(
                N=len(FREE),
                Afun=stiff,
                Bfun=mass,
                Pfun=preconditioner,
                **kwargs,
            )

        else:
            raise ValueError(f"Unknown solver: {solver}")    
                
        
        # Reshape eigenvectors:
        # NOTE: for the indices where the DOF if blocked, the eigenvector is zero
        eigvec = np.zeros((N, eigenvectors.shape[1]))
        eigvec[FREE] = eigenvectors

        logger.info(f"Eigen value solution in {time() - start:.2e} seconds")
        self.clear()
        return eigenvalues, eigvec
