from yeti_iga.pymfiga.common.numerics.solvers import LinearSolver
from yeti_iga.pymfiga.common.material import LinearElasticity
from yeti_iga.pymfiga.common.numerics.quadrature_rules import IGAQuadratureRule
from yeti_iga.pymfiga.common.numerics.quadrature_rules.operations import Operations
from yeti_iga.pymfiga.iga.boundary import BoundaryCondition
from yeti_iga.pymfiga.iga.geometry import SinglePatch
from yeti_iga.pymfiga.iga.single_model.cls.cspace import SingleSpatialModel
from typing import Callable, List, Tuple, Union, Sequence, Optional
from time import time
import numpy as np
import logging

logger = logging.getLogger("SRC.IGA.MODEL")


class TimoshenkoModel(SingleSpatialModel):

    # Only for in-plane beams
    FIXEDNBDOFS = 3  # UX, UY, RZ

    def __init__(
        self,
        material: LinearElasticity,
        axis_patch: SinglePatch,
        boundary: BoundaryCondition,
    ):
        if not axis_patch.ndim == 1:
            raise ValueError("Axis patch should be univariate")
        if not len(boundary.dofs_index) == self.FIXEDNBDOFS:
            raise ValueError(
                f"Boundary condition should have {self.FIXEDNBDOFS} dofs (UX, UY, RZ)"
            )

        super().__init__(material, axis_patch, boundary)
        self._compute_inplane_curvature()
        self.add_area_section(1.0, is_uniform=True)
        self.add_inertia_section(1.0, is_uniform=True)
        self._scalar_mean_mass: Optional[List[float]] = None
        self._scalar_mean_stiffness: Optional[List[np.ndarray]] = None
        self.compute_mass_property()
        self.compute_stiffness_property()

    @property
    def material(self):
        assert isinstance(self._material, LinearElasticity)
        return self._material

    @property
    def scalar_mean_mass(self):
        return self._scalar_mean_mass

    @property
    def scalar_mean_stiffness(self):
        return self._scalar_mean_stiffness

    def clear_properties(self):
        if self.material.hasnonlinearmass:
            self._scalar_mean_mass = None
        if self.material.hasnonlinearstiffness:
            self._scalar_mean_stiffness = None

    def verify_fun_args(self, args: dict):
        super().verify_fun_args(args)
        args.update({"quadpts": self.part.quadrule_list[0].quadpts})

    def add_area_section(self, inpt: Union[Callable, float], is_uniform: bool):
        self.cross_section: Callable = self.material.set_scalar_property(
            inpt, is_uniform=is_uniform
        )

    def add_inertia_section(self, inpt: Union[Callable, float], is_uniform: bool):
        self.section_inertia: Callable = self.material.set_scalar_property(
            inpt, is_uniform=is_uniform
        )

    def _compute_inplane_curvature(self):
        part: SinglePatch = self.part
        quadrule: IGAQuadratureRule = part.quadrule_list[0]
        degree: int = quadrule.degree
        knotvector: np.ndarray = quadrule.knotvector
        knots: np.ndarray = quadrule.quadpts
        ctrlpts: np.ndarray = part.ctrlpts
        nurbs_weights: np.ndarray = part.nurbs_weights
        self.inv_radius = np.zeros(part.nbqp_total)
        if degree == 1 or ctrlpts.shape[0] == 1:
            return

        # Compute the derivatives of BSplines
        basis_list = Operations.eval_ders_basis_sparse(
            degree, knotvector, knots, nders=2, is_periodic=False
        )

        def d_dxi_bspline(array_in: np.ndarray) -> Sequence[np.ndarray]:
            return [basis_list[i] @ array_in for i in range(3)]

        if nurbs_weights is None:

            _, dx_dxi, d2x_dxi2 = d_dxi_bspline(ctrlpts[0])
            _, dy_dxi, d2y_dxi2 = d_dxi_bspline(ctrlpts[1])

        else:

            def d_dxi_nurbs(array_in):
                new_array_in = nurbs_weights * array_in
                Bb, Bbdot, Bbdotdot = d_dxi_bspline(new_array_in)
                Ww, Wwdot, Wwdotdot = d_dxi_bspline(nurbs_weights)
                d_dxi = Bbdot / Ww - Bb * Wwdot / Ww**2
                d2_dxi2 = (
                    Bbdotdot / Ww
                    - Bbdot * Wwdot / Ww**2
                    - Bbdot * Wwdot / Ww**2
                    - Bb * Wwdotdot / Ww**2
                    + 2 * Bb * Wwdot**2 / Ww**3
                )
                return d_dxi, d2_dxi2

            dx_dxi, d2x_dxi2 = d_dxi_nurbs(ctrlpts[0])
            dy_dxi, d2y_dxi2 = d_dxi_nurbs(ctrlpts[1])

        self.inv_radius: np.ndarray = (
            dx_dxi * d2y_dxi2 - dy_dxi * d2x_dxi2
        ) / self.part.det_jac**3

    def compute_mass_property(self, **mf_args):
        self.verify_fun_args(mf_args)
        scalar_mean_mass = np.zeros(self.nbvars)
        p = self.material.density(mf_args)
        S = self.cross_section(mf_args)
        I = self.section_inertia(mf_args)
        detJ = self.part.det_jac

        prop = p * S * detJ
        mean_prop = np.mean(prop)
        scalar_mean_mass[0] = mean_prop
        scalar_mean_mass[1] = mean_prop

        prop = p * I * detJ
        scalar_mean_mass[2] = np.mean(prop)

        self._scalar_mean_mass = scalar_mean_mass.tolist()

    def compute_mf_mass(self, array_in: np.ndarray, **mf_args):
        
        self.verify_fun_args(mf_args)
        if self._scalar_mean_mass is None:
            self.compute_mass_property(**mf_args)

        nurbs_weights: np.ndarray = self.part.nurbs_weights
        mass: Callable = self.operator_engine.compute_mf_scalar_u_v
        p = self.material.density(mf_args)
        S = self.cross_section(mf_args)
        I = self.section_inertia(mf_args)
        detJ = self.part.det_jac
        quadlist = self.part.quadrule_list

        def m11_and_22(arr_in):
            prop = p * S * detJ
            arr_out = mass(quadlist, prop, arr_in, False, nurbs_weights=nurbs_weights)
            return arr_out

        def m33(arr_in):
            prop = p * I * detJ
            arr_out = mass(quadlist, prop, arr_in, False, nurbs_weights=nurbs_weights)
            return arr_out

        array_in = np.reshape(array_in, (self.FIXEDNBDOFS, -1))
        matrix_product = [m11_and_22, m11_and_22, m33]
        array_out = np.zeros_like(array_in)
        for i in range(self.nbvars):
            array_out[i] = matrix_product[i](array_in[i])
        return np.ravel(array_out)

    def compute_stiffness_property(self, **mf_args):
        self.verify_fun_args(mf_args)
        scalar_mean_stiffness = []
        S = self.cross_section(mf_args)
        I = self.section_inertia(mf_args)
        E = self.material.elastic_modulus
        G = self.material.lame_mu
        detJ = self.part.det_jac

        prop = E * S / detJ
        scalar_mean_stiffness.append(np.array([np.mean(prop)]))

        prop = G * S / detJ
        scalar_mean_stiffness.append(np.array([np.mean(prop)]))

        prop = E * I / detJ
        scalar_mean_stiffness.append(np.array([np.mean(prop)]))

        self._scalar_mean_stiffness = scalar_mean_stiffness

    def compute_mf_stiffness(self, array_in: np.ndarray, **mf_args):
		
        self.verify_fun_args(mf_args)
        if self._scalar_mean_stiffness is None:
            self.compute_stiffness_property(**mf_args)

        nurbs_weights: np.ndarray = self.part.nurbs_weights
        # Simplify notation
        stif: Callable = self.operator_engine.compute_mf_scalar_gradu_gradv
        mass: Callable = self.operator_engine.compute_mf_scalar_u_v
        adv_gn: Callable = self.operator_engine.compute_mf_scalar_gradu_v
        adv_ng: Callable = self.operator_engine.compute_mf_scalar_u_gradv
        lumping = False
        S = self.cross_section(mf_args)
        I = self.section_inertia(mf_args)
        E = self.material.elastic_modulus
        G = self.material.lame_mu
        detJ = self.part.det_jac
        invR = self.inv_radius
        quadlist = self.part.quadrule_list

        def k11(arr_in):
            prop = G * S * invR**2 * detJ
            arr_out = mass(quadlist, prop, arr_in, lumping, nurbs_weights=nurbs_weights)

            prop = E * S / detJ
            arr_out += stif(
                quadlist,
                np.reshape(prop, (1, 1, -1)),
                arr_in,
                nurbs_weights=nurbs_weights,
            )
            return arr_out

        def k22(arr_in):
            prop = E * S * invR**2 * detJ
            arr_out = mass(quadlist, prop, arr_in, lumping, nurbs_weights=nurbs_weights)

            prop = G * S / detJ
            arr_out += stif(
                quadlist,
                np.reshape(prop, (1, 1, -1)),
                arr_in,
                nurbs_weights=nurbs_weights,
            )
            return arr_out

        def k33(arr_in):
            prop = G * S * detJ
            arr_out = mass(quadlist, prop, arr_in, lumping, nurbs_weights=nurbs_weights)

            prop = E * I / detJ
            arr_out += stif(
                quadlist,
                np.reshape(prop, (1, 1, -1)),
                arr_in,
                nurbs_weights=nurbs_weights,
            )
            return arr_out

        def k12(arr_in):
            prop = -E * S * invR
            arr_out = adv_gn(
                quadlist, np.reshape(prop, (1, -1)), arr_in, nurbs_weights=nurbs_weights
            )

            prop = G * S * invR
            arr_out += adv_ng(
                quadlist, np.reshape(prop, (1, -1)), arr_in, nurbs_weights=nurbs_weights
            )
            return arr_out

        def k21(arr_in):
            prop = -E * S * invR
            arr_out = adv_ng(
                quadlist, np.reshape(prop, (1, -1)), arr_in, nurbs_weights=nurbs_weights
            )

            prop = G * S * invR
            arr_out += adv_gn(
                quadlist, np.reshape(prop, (1, -1)), arr_in, nurbs_weights=nurbs_weights
            )
            return arr_out

        def k13(arr_in):
            prop = -G * S * invR * detJ
            return mass(quadlist, prop, arr_in, lumping, nurbs_weights=nurbs_weights)

        def k31(arr_in):
            return k13(arr_in)

        def k23(arr_in):
            prop = -G * S
            return adv_gn(
                quadlist, np.reshape(prop, (1, -1)), arr_in, nurbs_weights=nurbs_weights
            )

        def k32(arr_in):
            prop = -G * S
            return adv_ng(
                quadlist, np.reshape(prop, (1, -1)), arr_in, nurbs_weights=nurbs_weights
            )

        array_in = np.reshape(array_in, (self.FIXEDNBDOFS, -1))
        matrix_product = [[k11, k12, k13], [k21, k22, k23], [k31, k32, k33]]
        array_out = np.zeros_like(array_in)
        for i in range(self.nbvars):
            for j in range(self.nbvars):
                array_out[i] += matrix_product[i][j](array_in[j])

        return np.ravel(array_out)

    def assemble_surface_force(self, fun_info, **external_args):
        raise NotImplementedError()

    def compute_residual(
        self, array_in: np.ndarray, **kwargs
    ) -> Tuple[np.ndarray, dict]:
        external_force: np.ndarray = kwargs.get("external_force", 0.0)
        residual = external_force - self.compute_mf_stiffness(array_in, **kwargs)
        self.clear_bcs(residual)
        # TODO: change {} for plastic variables, i.e., how to include non linear materials ?
        return residual, {}

    def solve_linearized_system(self, array_in: np.ndarray, **kwargs) -> np.ndarray:
        start = time()

        linear_solver: LinearSolver = kwargs["linear_solver_backend"]
        inner_tolerance: float = (
            kwargs.get("inner_tolerance") or linear_solver.config.tolerance
        )
        linear_solver.update(tolerance=inner_tolerance)
        linear_solver.set_cleandod(self._constraint_nodes)

        if self.update_manager.should_update_preconditioner:
            self.preconditioner.add_scalar_space_time_correctors(
                mass_corrector=self.scalar_mean_mass,
                stiffness_corrector=self.scalar_mean_stiffness,
            )
        self.preconditioner.update_space_eigenvalues(scalar_coefs=(0, 1))
        sol = linear_solver.solve(
            self.compute_mf_stiffness,
            array_in,
            Pfun=self.preconditioner.apply_spatial_preconditioner,
            **kwargs,
        )["sol"]

        if self.update_manager.should_update_material:
            self.clear_properties()

        logger.debug(f"Solving linearized system in {time() - start:.2e} seconds")
        return sol
