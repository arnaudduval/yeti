from yeti_iga.pymfiga.common.numerics.quadrature_rules import StandardGauss
from yeti_iga.pymfiga.iga.geometry import SinglePatch
from yeti_iga.pymfiga.iga.single_model.cls.core import SingleModel
from yeti_iga.pymfiga.iga.single_model.cls.cspace import SingleSpatialModel
from yeti_iga.pymfiga.iga.single_model.cls.ctime import SingleSpaceTimeModel
from typing import Union, List, Tuple, Callable, Optional, Sequence, Literal
import numpy as np
import logging

logger = logging.getLogger("SRC.IGA.NORM")


class SinglePatchNorm:
    def __init__(
        self,
        model: Union[SingleSpatialModel, SingleSpaceTimeModel],
        norm_type: Literal["l2", "h1", "semih1"],
        norm_args: dict,
    ):
        assert isinstance(model, SingleModel)
        if not norm_type in ["l2", "h1", "semih1"]:
            raise NotImplementedError()
        self._model = model
        self._norm_type = norm_type
        self._norm_args = norm_args

    @property
    def model(self):
        return self._model

    @property
    def norm_type(self):
        return self._norm_type

    @property
    def ndim(self):
        if isinstance(self.model, SingleSpatialModel):
            return self.model.ndim
        elif isinstance(self.model, SingleSpaceTimeModel):
            return self.model.ndim + 1
        else:
            raise NotImplementedError()

    def _calculate_norm(
        self,
        u: np.ndarray,
        uders: Optional[np.ndarray],
        det_jac: np.ndarray,
        weights: List[np.ndarray],
    ) -> float:

        u2_l2, u2_sh1 = 0.0, 0.0
        if self.norm_type in ["l2", "h1"]:
            u2_l2 += np.sum(u**2, axis=0)
        if self.norm_type in ["h1", "semih1"]:
            assert isinstance(uders, np.ndarray)
            u2_sh1 += np.sum(uders**2, axis=(0, 1))

        shape = [len(wgt) for wgt in weights]
        u2_det_jac = np.reshape((u2_l2 + u2_sh1) * det_jac, tuple(shape), order="F")
        text = f"{','.join('abcdef'[:self.ndim])},{'abcdef'[:self.ndim]}->"
        return np.sqrt(np.einsum(text, *weights, u2_det_jac))

    def _compute_exact_values(
        self,
        quadpts_phy: List[np.ndarray],
        quadpts_param: List[np.ndarray],
        prepare_norm_computation: Callable,
        allow_spacetime: bool,
    ) -> Sequence[Optional[np.ndarray]]:

        norm_type, norm_args = self.norm_type, self._norm_args
        u_exact, uders_exact = np.array([]), np.array([])

        # First method
        exact_fun: Optional[Callable[[dict], np.ndarray]] = norm_args.get(
            "exact_function", None
        )
        exact_dersfun: Optional[Callable[[dict], np.ndarray]] = norm_args.get(
            "exact_function_ders", None
        )
        if exact_fun is not None:
            exact_args: dict = norm_args.get("exact_args", {})
            assert isinstance(exact_args, dict), "Error type of extra args"
            assert isinstance(quadpts_phy, list)

            exact_args.update({"position": quadpts_phy[0]})
            if allow_spacetime:
                exact_args.update({"time": np.ravel(quadpts_phy[1])})

            u_exact: np.ndarray = exact_fun(exact_args)
            assert u_exact.ndim < 3
            if u_exact.ndim < 2:
                u_exact = np.atleast_2d(u_exact)

            if exact_dersfun is not None:
                uders_exact: np.ndarray = exact_dersfun(exact_args)
                assert uders_exact.ndim < 4
                if uders_exact.ndim < 3:
                    uders_exact = np.atleast_3d(uders_exact)

        # Second method
        part_ref: Optional[SinglePatch] = norm_args.get("part_ref", None)
        time_ref: Optional[SinglePatch] = norm_args.get("time_ref", None)
        u_ref: Optional[np.ndarray] = norm_args.get("u_ref", None)
        if (
            allow_spacetime
            and isinstance(part_ref, SinglePatch)
            and isinstance(time_ref, SinglePatch)
        ):
            assert isinstance(u_ref, np.ndarray), "Solution should be numpy array"
            u_exact, uders_exact = prepare_norm_computation(
                part_ref,
                time_ref,
                u_ref,
                knots_ref=quadpts_param,
                norm_type=norm_type,
            )[:2]
        elif isinstance(part_ref, SinglePatch):
            assert isinstance(u_ref, np.ndarray), "Solution should be numpy array"
            u_exact, uders_exact = prepare_norm_computation(
                part_ref, u_ref, knots_ref=quadpts_param, norm_type=norm_type
            )[:2]

        if u_exact.size == 0 or (norm_type != "l2" and uders_exact.size == 0):
            raise ValueError(
                "Exact solution or its derivatives are not properly defined."
            )
        return u_exact, uders_exact

    def _eval(
        self,
        u_ctrlpts: np.ndarray,
        prepare_norm_computation: Callable,
        allow_spacetime: bool,
    ) -> Tuple[float, float]:

        norm_type = self.norm_type

        # Compute u interp
        if allow_spacetime:
            output = prepare_norm_computation(
                self.model.part,
                self.model.time,
                u_ctrlpts,
                knots_ref=[None] * self.ndim,
                norm_type=norm_type,
            )
        else:
            output = prepare_norm_computation(
                self.model.part,
                u_ctrlpts,
                knots_ref=[None] * self.ndim,
                norm_type=norm_type,
            )
        (
            u_interp,
            uders_interp,
            parametric_position,
            quadpts_phy,
            det_jac,
            parametric_weights,
        ) = output

        # Compute u exact
        u_exact, uders_exact = self._compute_exact_values(
            quadpts_phy,
            parametric_position,
            prepare_norm_computation,
            allow_spacetime=allow_spacetime,
        )

        assert isinstance(u_exact, np.ndarray)
        assert isinstance(u_interp, np.ndarray)
        if norm_type == "l2":
            abserror = self._calculate_norm(
                u_exact - u_interp,
                None,
                det_jac,
                parametric_weights,
            )
            tmp2 = self._calculate_norm(u_exact, None, det_jac, parametric_weights)
        else:
            assert isinstance(uders_exact, np.ndarray)
            assert isinstance(uders_interp, np.ndarray)
            abserror = self._calculate_norm(
                u_exact - u_interp,
                uders_exact - uders_interp,
                det_jac,
                parametric_weights,
            )
            tmp2 = self._calculate_norm(
                u_exact, uders_exact, det_jac, parametric_weights
            )

        relerror = abserror / tmp2 if tmp2 != 0 else abserror
        if tmp2 == 0:
            logger.warning("Warning: Dividing by zero")

        return abserror, relerror, tmp2


class SpaceNormSinglePatch(SinglePatchNorm):
    def __init__(
        self,
        model: SingleSpatialModel,
        norm_type: Literal["l2", "h1", "semih1"],
        norm_args: dict,
    ):
        super().__init__(model, norm_type, norm_args)

    def _prepare_norm_computation(
        self,
        part_ref: SinglePatch,
        u_ref: np.ndarray,
        knots_ref: List[Optional[np.ndarray]],
        norm_type: str = "l2",
    ) -> Tuple[
        np.ndarray,
        Optional[np.ndarray],
        List[np.ndarray],
        List[np.ndarray],
        np.ndarray,
        List[np.ndarray],
    ]:
        assert isinstance(self.model, SingleSpatialModel)
        parametric_position, parametric_weights = [], []
        for i in range(part_ref.ndim):
            quadrule = StandardGauss(
                part_ref.degree[i],
                part_ref.knotvector[i],
                quadtype="legendre",
                is_periodic=part_ref.quadrule_list[i].is_periodic,
            )
            quadrule.export_quadrature_rules()
            knots_to_copy = knots_ref[i]
            position = (
                knots_to_copy.copy() if knots_to_copy is not None else quadrule.quadpts
            )
            parametric_position.append(position)
            parametric_weights.append(quadrule.parametric_weights)

        for quadrule, knots in zip(part_ref.quadrule_list, parametric_position):
            quadrule.knots_to_sample = knots

        det_jac, _, inv_jac, quadpts_phy_space = SinglePatch.eval_transformation(
            part_ref.quadrule_list,
            part_ref.ctrlpts,
            nurbs_weights=self.model.part.nurbs_weights,
        )

        quadpts_phy = [quadpts_phy_space]

        u_interp = self.model.operator_engine.interpolate_meshgrid(
            part_ref.quadrule_list,
            np.atleast_2d(u_ref),
            nurbs_weights=self.model.part.nurbs_weights,
        )
        uders_interp = None
        if norm_type.lower() != "l2":
            derstemp = self.model.operator_engine.eval_jacobien(
                part_ref.quadrule_list,
                np.atleast_2d(u_ref),
                nurbs_weights=self.model.part.nurbs_weights,
            )
            uders_interp = np.atleast_3d(
                np.einsum("ijl,jkl->ikl", derstemp, inv_jac, optimize=True)
            )

        for quadrule in part_ref.quadrule_list:
            quadrule.clear_sample()

        return (
            u_interp,
            uders_interp,
            parametric_position,
            quadpts_phy,
            det_jac,
            parametric_weights,
        )

    def eval(self, u_ctrlpts):
        return super()._eval(u_ctrlpts, self._prepare_norm_computation, False)


class SpaceTimeNormSinglePatch(SinglePatchNorm):
    def __init__(
        self,
        model: SingleSpaceTimeModel,
        norm_type: Literal["l2", "h1", "semih1"],
        norm_args: dict,
    ):
        super().__init__(model, norm_type, norm_args)

    def _prepare_norm_computation(
        self,
        part_ref: SinglePatch,
        time_ref: SinglePatch,
        u_ref: np.ndarray,
        knots_ref: List[Union[None, np.ndarray]],
        norm_type: str = "l2",
    ) -> Tuple[
        np.ndarray,
        Union[np.ndarray, None],
        List[np.ndarray],
        List[np.ndarray],
        np.ndarray,
        List[np.ndarray],
    ]:
        assert isinstance(self.model, SingleSpaceTimeModel)
        parametric_position, parametric_weights = [], []
        # For space variables
        for i in range(part_ref.ndim):
            quadrule = StandardGauss(
                part_ref.degree[i],
                part_ref.knotvector[i],
                quadtype="legendre",
                is_periodic=part_ref.quadrule_list[i].is_periodic,
            )
            quadrule.export_quadrature_rules()
            knots_to_copy = knots_ref[i]
            position = (
                knots_to_copy.copy() if knots_to_copy is not None else quadrule.quadpts
            )
            parametric_position.append(position)
            parametric_weights.append(quadrule.parametric_weights)

        # For time variables
        quadrule = StandardGauss(
            time_ref.degree[0],
            time_ref.knotvector[0],
            quadtype="legendre",
            is_periodic=time_ref.quadrule_list[0].is_periodic,
        )
        quadrule.export_quadrature_rules()
        parametric_position.append(
            quadrule.quadpts
            if not isinstance(knots_ref[-1], np.ndarray)
            else np.copy(knots_ref[-1])
        )
        parametric_weights.append(quadrule.parametric_weights)

        # Update space and time
        all_quadrule_list = part_ref.quadrule_list + time_ref.quadrule_list
        for quadrule, knots in zip(all_quadrule_list, parametric_position):
            quadrule.knots_to_sample = knots

        (
            det_jac_space,
            _,
            inv_jac_space,
            quadpts_phy_space,
        ) = SinglePatch.eval_transformation(
            part_ref.quadrule_list,
            part_ref.ctrlpts,
            nurbs_weights=self.model.sptm_nurbs_weights,
        )

        (
            det_jac_time,
            _,
            inv_jac_time,
            quadpts_phy_time,
        ) = SinglePatch.eval_transformation(
            time_ref.quadrule_list,
            time_ref.ctrlpts,
            nurbs_weights=self.model.sptm_nurbs_weights,
        )

        u_interp = self.model.operator_engine.interpolate_meshgrid(
            all_quadrule_list,
            np.atleast_2d(u_ref),
            nurbs_weights=self.model.sptm_nurbs_weights,
        )
        uders_interp = None
        if norm_type != "l2":
            nm = np.shape(u_interp)[0]
            derstmp = self.model.operator_engine.eval_jacobien(
                all_quadrule_list,
                np.atleast_2d(u_ref),
                nurbs_weights=self.model.sptm_nurbs_weights,
            )
            derstmp_reshaped = np.reshape(
                derstmp,
                (nm, part_ref.ndim + 1, len(parametric_position[-1]), -1),
            )
            uders_interp = np.zeros_like(derstmp_reshaped)
            uders_interp[:, :-1, :, :] = np.einsum(
                "mipk,ijk->mjpk",
                derstmp_reshaped[:, :-1, :, :],
                inv_jac_space,
                optimize=True,
            )
            uders_interp[:, -1, :, :] = np.einsum(
                "mpk,p->mpk",
                derstmp_reshaped[:, -1, :, :],
                np.ravel(inv_jac_time),
                optimize=True,
            )
            uders_interp = np.reshape(uders_interp, (nm, self.model.ndim, -1))

        quadpts_phy = [quadpts_phy_space, quadpts_phy_time]
        det_jac = np.kron(det_jac_time, det_jac_space)

        for quadrule in all_quadrule_list:
            quadrule.clear_sample()

        return (
            u_interp,
            uders_interp,
            parametric_position,
            quadpts_phy,
            det_jac,
            parametric_weights,
        )

    def eval(self, u_ctrlpts: np.ndarray) -> Tuple[float, float]:
        return super()._eval(u_ctrlpts, self._prepare_norm_computation, True)
