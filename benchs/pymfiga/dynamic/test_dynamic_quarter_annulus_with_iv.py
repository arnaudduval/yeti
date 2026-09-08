# -*- coding: utf-8 -*-
"""
Created on Wed Aug 26 19:32:54 2026
Test de validation de la dynamique explicite
Erreur L2 espace-temps
Quart d'anneau avec vitesse initiale
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import sympy as sp

from yeti_iga.pymfiga.iga.boundary import (
    BoundaryCondition,
    DOF,
    ParametricDirection,
    BoundarySide,
)
from yeti_iga.pymfiga.common.physics import ExplicitLinearDynamics
from yeti_iga.pymfiga.iga.geometry import GeomdlGenerator, SinglePatch
from yeti_iga.pymfiga.iga.single_model import ExplicitDynamicsModel
from yeti_iga.pymfiga.common.material import LinearElasticity
from yeti_iga.pymfiga.iga.norms import SpaceNormSinglePatch
from yeti_iga.pymfiga.common.numerics.solvers import LinearSolver


# ============================================================
# PARAMÈTRES DU CAS TEST
# ============================================================

DEGREE = 2
NBEL = 4

TIMESPAN = 3.0e-4

# Pas de temps imposés directement
DT_MC = 7.350345723891005e-06
DT_ML = 1.850507202627122e-05

QUADCLASS = "wq"
QUADTYPE = "2"

RIN = 0.25
REX = 1.0


# ============================================================
# VALEURS DE RÉFÉRENCE
# ============================================================

REFERENCE_L2_SPACETIME = {

    "MC_NONE":
        1.926432208417065e-02,

    "MC_FASTDIAG":
        1.926432208060571e-02,

    "MC_SCALED_MASS":
        1.926432205979960e-02,

    "ML":
        5.774078310585155e-01,
}


# ============================================================
# MATÉRIAU
# ============================================================

E = 2.1e11
NU = 0.3
RHO = 7850.0

LAM = (
    NU * E
    / (
        (1.0 + NU)
        * (1.0 - 2.0 * NU)
    )
)

MU = (
    E
    / (
        2.0
        * (1.0 + NU)
    )
)


# ============================================================
# SOLUTION EXACTE MANUFACTURÉE
# ============================================================

U0 = 1.0e-6

OMEGA = (
    np.pi
    / (
        2.0
        * TIMESPAN
    )
)


def spatial_amplitude(position):

    position = np.asarray(
        position,
        dtype=float,
    )

    x = position[0]
    y = position[1]

    r2 = (
        x**2
        + y**2
    )

    return (
        U0
        * (r2 - RIN**2)
        * (r2 - REX**2)
        * x
        * y
    )


def exact_displacement(position, t):

    amplitude = spatial_amplitude(
        position
    )

    value = (
        amplitude
        * np.sin(
            OMEGA * t
        )
    )

    return np.vstack(
        (
            value,
            value,
        )
    )


# ============================================================
# CONSTRUCTION SYMBOLIQUE DE LA FORCE VOLUMIQUE
# ============================================================

x_sym, y_sym, t_sym = sp.symbols(
    "x y t",
    real=True,
)

r2_sym = (
    x_sym**2
    + y_sym**2
)

A_sym = (
    U0
    * (r2_sym - RIN**2)
    * (r2_sym - REX**2)
    * x_sym
    * y_sym
)

ux_sym = (
    A_sym
    * sp.sin(
        OMEGA * t_sym
    )
)

uy_sym = (
    A_sym
    * sp.sin(
        OMEGA * t_sym
    )
)


# ============================================================
# ACCÉLÉRATIONS
# ============================================================

ux_tt_sym = sp.diff(
    ux_sym,
    t_sym,
    2,
)

uy_tt_sym = sp.diff(
    uy_sym,
    t_sym,
    2,
)


# ============================================================
# DÉFORMATIONS
# ============================================================

eps_xx_sym = sp.diff(
    ux_sym,
    x_sym,
)

eps_yy_sym = sp.diff(
    uy_sym,
    y_sym,
)

eps_xy_sym = (
    sp.Rational(1, 2)
    * (
        sp.diff(
            ux_sym,
            y_sym,
        )
        +
        sp.diff(
            uy_sym,
            x_sym,
        )
    )
)

trace_eps_sym = (
    eps_xx_sym
    + eps_yy_sym
)


# ============================================================
# CONTRAINTES
# ============================================================

sigma_xx_sym = (
    2.0
    * MU
    * eps_xx_sym
    +
    LAM
    * trace_eps_sym
)

sigma_yy_sym = (
    2.0
    * MU
    * eps_yy_sym
    +
    LAM
    * trace_eps_sym
)

sigma_xy_sym = (
    2.0
    * MU
    * eps_xy_sym
)


# ============================================================
# DIVERGENCE DES CONTRAINTES
# ============================================================

div_sigma_x_sym = (
    sp.diff(
        sigma_xx_sym,
        x_sym,
    )
    +
    sp.diff(
        sigma_xy_sym,
        y_sym,
    )
)

div_sigma_y_sym = (
    sp.diff(
        sigma_xy_sym,
        x_sym,
    )
    +
    sp.diff(
        sigma_yy_sym,
        y_sym,
    )
)


# ============================================================
# FORCE VOLUMIQUE MANUFACTURÉE
# ============================================================

fx_sym = sp.simplify(
    RHO
    * ux_tt_sym
    -
    div_sigma_x_sym
)

fy_sym = sp.simplify(
    RHO
    * uy_tt_sym
    -
    div_sigma_y_sym
)

fx_numpy = sp.lambdify(
    (
        x_sym,
        y_sym,
        t_sym,
    ),
    fx_sym,
    modules="numpy",
)

fy_numpy = sp.lambdify(
    (
        x_sym,
        y_sym,
        t_sym,
    ),
    fy_sym,
    modules="numpy",
)


def body_force_x(args):

    position = np.asarray(
        args["position"],
        dtype=float,
    )

    t = float(
        args["time"]
    )

    return np.asarray(
        fx_numpy(
            position[0],
            position[1],
            t,
        ),
        dtype=float,
    )


def body_force_y(args):

    position = np.asarray(
        args["position"],
        dtype=float,
    )

    t = float(
        args["time"]
    )

    return np.asarray(
        fy_numpy(
            position[0],
            position[1],
            t,
        ),
        dtype=float,
    )


# ============================================================
# FONCTION EXACTE POUR LA NORME
# ============================================================

def make_exact_function(t):

    t_value = float(t)

    def exact_function(args):

        position = np.asarray(
            args["position"],
            dtype=float,
        )

        return exact_displacement(
            position,
            t_value,
        )

    return exact_function


# ============================================================
# PROJECTION L2 DE LA VITESSE INITIALE
# ============================================================

def project_initial_velocity(
    model,
    use_preconditioner,
    preconditioner_type,
):

    def velocity_rhs_x(args):

        position = np.asarray(
            args["position"],
            dtype=float,
        )

        amplitude = spatial_amplitude(
            position
        )

        return (
            RHO
            * OMEGA
            * amplitude
        )

    def velocity_rhs_y(args):

        position = np.asarray(
            args["position"],
            dtype=float,
        )

        amplitude = spatial_amplitude(
            position
        )

        return (
            RHO
            * OMEGA
            * amplitude
        )

    rhs_velocity = model.assemble_volumetric_force(
        {
            DOF.UX: velocity_rhs_x,
            DOF.UY: velocity_rhs_y,
        }
    )

    rhs_velocity = np.asarray(
        rhs_velocity,
        dtype=float,
    ).reshape(-1)

    linear_solver = LinearSolver(
        maxiters=10000,
        tolerance=1.0e-12,
        linear_type="cg",
    )

    velocity_initial = model.solve_linearized_system(
        rhs_velocity,
        linear_solver_backend=linear_solver,
        use_preconditioner=use_preconditioner,
        preconditioner_type=preconditioner_type,
    )

    velocity_initial = np.asarray(
        velocity_initial,
        dtype=float,
    ).reshape(-1)

    if (
        velocity_initial.shape
        != rhs_velocity.shape
    ):
        raise ValueError(
            "Dimension incorrecte après projection : "
            f"velocity_initial.shape="
            f"{velocity_initial.shape}, "
            f"rhs_velocity.shape="
            f"{rhs_velocity.shape}"
        )

    return velocity_initial


# ============================================================
# CONSTRUCTION DE LA LISTE TEMPORELLE
# ============================================================

def build_time_list(timespan, dt):

    if timespan <= 0.0:
        raise ValueError(
            "timespan doit être strictement positif."
        )

    if dt <= 0.0:
        raise ValueError(
            "dt doit être strictement positif."
        )

    nbsteps = int(
        np.floor(
            timespan / dt
        )
    )

    time_list = (
        np.arange(
            nbsteps + 1,
            dtype=float,
        )
        * dt
    )

    # Le dernier pas est éventuellement raccourci
    # afin d'atteindre exactement TIMESPAN.
    if time_list[-1] < timespan:

        time_list = np.append(
            time_list,
            timespan,
        )

    return time_list


# ============================================================
# NORME L2 ESPACE-TEMPS
# ============================================================

def NormeSpaceTimeL2Exact(
    model,
    time_list,
    displacement,
):

    time_eval = np.asarray(
        time_list,
        dtype=float,
    )

    displacement = np.asarray(
        displacement,
        dtype=float,
    )

    if (
        len(time_eval)
        != displacement.shape[0]
    ):
        raise ValueError(
            "Nombre de temps incompatible avec "
            "le nombre de déplacements."
        )

    if len(time_eval) < 2:
        raise ValueError(
            "Il faut au moins deux instants."
        )

    nbctrlpts = (
        model.part.nbctrlpts_total
    )

    ndim = model.ndim

    error2_time = 0.0
    exact2_time = 0.0

    for k, t in enumerate(
        time_eval
    ):

        u_num = displacement[k].reshape(
            (
                ndim,
                nbctrlpts,
            )
        )

        norm_error, _, norm_exact = (
            SpaceNormSinglePatch(
                model,
                "l2",
                {
                    "exact_function":
                        make_exact_function(t)
                },
            ).eval(u_num)
        )

        if k == 0:

            weight = (
                0.5
                * (
                    time_eval[1]
                    - time_eval[0]
                )
            )

        elif k == len(time_eval) - 1:

            weight = (
                0.5
                * (
                    time_eval[-1]
                    - time_eval[-2]
                )
            )

        else:

            weight = (
                0.5
                * (
                    time_eval[k + 1]
                    - time_eval[k - 1]
                )
            )

        error2_time += (
            weight
            * norm_error**2
        )

        exact2_time += (
            weight
            * norm_exact**2
        )

    if exact2_time <= 0.0:
        raise ZeroDivisionError(
            "La norme exacte espace-temps est nulle."
        )

    return np.sqrt(
        error2_time
        / exact2_time
    )


# ============================================================
# CAS À TESTER
# ============================================================

CALC_CASES = [

    {
        "case_name": "MC_NONE",
        "mass_type": "consistent_mass",
        "dt": DT_MC,
        "use_preconditioner": False,
        "preconditioner_type": None,
    },

    {
        "case_name": "MC_FASTDIAG",
        "mass_type": "consistent_mass",
        "dt": DT_MC,
        "use_preconditioner": True,
        "preconditioner_type": "fastdiag",
    },

    {
        "case_name": "MC_SCALED_MASS",
        "mass_type": "consistent_mass",
        "dt": DT_MC,
        "use_preconditioner": True,
        "preconditioner_type": "scaled_mass",
    },

    {
        "case_name": "ML",
        "mass_type": "diagonal_mass",
        "dt": DT_ML,
        "use_preconditioner": False,
        "preconditioner_type": None,
    },
]


# ============================================================
# SIMULATION
# ============================================================

def simulate(calc_case):

    # --------------------------------------------------------
    # Géométrie
    # --------------------------------------------------------

    geometry = GeomdlGenerator(
        filename="quarter_annulus",
        geo_args={
            "degree": DEGREE,
            "nbel": NBEL,
        },
    ).export_geometry()

    patch = SinglePatch(
        geometry,
        quadclass=QUADCLASS,
        quadtype=QUADTYPE,
    )

    patch.generate()

    # --------------------------------------------------------
    # Matériau
    # --------------------------------------------------------

    material = LinearElasticity(
        {
            "elastic_modulus": E,
            "poisson_ratio": NU,
        }
    )

    material.add_density(
        RHO,
        is_uniform=True,
    )

    # --------------------------------------------------------
    # Conditions aux limites
    # --------------------------------------------------------

    boundary = BoundaryCondition(
        nbctrlpts=patch.nbctrlpts,
        dofs=(
            DOF.UX,
            DOF.UY,
        ),
    )

    boundary.add_constraint(
        constraint_info=[
            {
                "direction":
                    ParametricDirection.XI,
                "face":
                    BoundarySide.MIN,
                "dofs":
                    (DOF.UX, DOF.UY),
            },
            {
                "direction":
                    ParametricDirection.XI,
                "face":
                    BoundarySide.MAX,
                "dofs":
                    (DOF.UX, DOF.UY),
            },
            {
                "direction":
                    ParametricDirection.ETA,
                "face":
                    BoundarySide.MIN,
                "dofs":
                    (DOF.UX, DOF.UY),
            },
            {
                "direction":
                    ParametricDirection.ETA,
                "face":
                    BoundarySide.MAX,
                "dofs":
                    (DOF.UX, DOF.UY),
            },
        ],
        constraint_type="dirichlet",
    )

    # --------------------------------------------------------
    # Modèle
    # --------------------------------------------------------

    model = ExplicitDynamicsModel(
        material,
        patch,
        boundary,
        mass_type=calc_case["mass_type"],
    )

    # --------------------------------------------------------
    # Temps
    # --------------------------------------------------------

    time_list = build_time_list(
        timespan=TIMESPAN,
        dt=calc_case["dt"],
    )

    # --------------------------------------------------------
    # Force volumique
    # --------------------------------------------------------

    def external_force_time(t, it):

        return model.assemble_volumetric_force(
            {
                DOF.UX:
                    body_force_x,
                DOF.UY:
                    body_force_y,
            },
            time=t,
        )

    initial_force = external_force_time(
        time_list[0],
        0,
    )

    # --------------------------------------------------------
    # Déplacement initial
    # --------------------------------------------------------

    displacement_initial = np.zeros_like(
        initial_force,
        dtype=float,
    )

    # --------------------------------------------------------
    # Vitesse initiale non nulle
    # --------------------------------------------------------

    velocity_initial = project_initial_velocity(
        model=model,
        use_preconditioner=
            calc_case["use_preconditioner"],
        preconditioner_type=
            calc_case["preconditioner_type"],
    )

    if (
        velocity_initial.shape
        != displacement_initial.shape
    ):
        raise ValueError(
            "Dimensions incompatibles : "
            f"velocity_initial.shape="
            f"{velocity_initial.shape}, "
            f"displacement_initial.shape="
            f"{displacement_initial.shape}"
        )

    # --------------------------------------------------------
    # Résolution dynamique
    # --------------------------------------------------------

    displacement, saved_time_list = (
        ExplicitLinearDynamics().solve(
            model,
            displacement_initial,
            external_force_time,
            time_list,
            velocity_initial=velocity_initial,
            nb_save=500,
            use_preconditioner=
                calc_case["use_preconditioner"],
            preconditioner_type=
                calc_case["preconditioner_type"],
        )
    )

    displacement = np.asarray(
        displacement,
        dtype=float,
    )

    if not np.all(
        np.isfinite(displacement)
    ):
        raise FloatingPointError(
            f"Solution invalide pour "
            f"{calc_case['case_name']}."
        )

    # --------------------------------------------------------
    # Erreur L2 espace-temps
    # --------------------------------------------------------

    error_spacetime = (
        NormeSpaceTimeL2Exact(
            model=model,
            time_list=saved_time_list,
            displacement=displacement,
        )
    )

    return error_spacetime


# ============================================================
# TEST DE VALIDATION
# ============================================================

def test_quarter_annulus_dynamic_spacetime_l2():

    TOLERANCE = 1e-6

    results = []

    for calc_case in CALC_CASES:

        case_name = calc_case[
            "case_name"
        ]

        # ----------------------------------------------------
        # Calcul de l'erreur L2 espace-temps
        # ----------------------------------------------------

        error_spacetime = simulate(
            calc_case
        )

        # ----------------------------------------------------
        # Valeur de référence
        # ----------------------------------------------------

        error_ref = (
            REFERENCE_L2_SPACETIME[
                case_name
            ]
        )

        # ----------------------------------------------------
        # Écart relatif
        # ----------------------------------------------------

        relative_error = abs(
            (
                error_spacetime
                - error_ref
            )
            / error_ref
        )

        # ----------------------------------------------------
        # Validation
        # ----------------------------------------------------

        validated = (
            relative_error
            < TOLERANCE
        )

        results.append(
            (
                case_name,
                error_spacetime,
                error_ref,
                relative_error,
                validated,
            )
        )

    # ========================================================
    # TABLEAU FINAL
    # ========================================================

    print("\n")
    print("=" * 90)

    print(
        "RÉSULTATS DE VALIDATION DYNAMIQUE "
        "- ERREUR L2 ESPACE-TEMPS "
        "- QUART D'ANNEAU AVEC VITESSE INITIALE"
    )

    print("=" * 90)

    print(
        f"degree = {DEGREE}, "
        f"nbel = {NBEL}, "
        f"tolérance relative = {TOLERANCE:.1e}"
    )

    print(
        f"dt MC = {DT_MC:.15e}, "
        f"dt ML = {DT_ML:.15e}"
    )

    print("-" * 90)

    print(
        f"{'Cas':20s} "
        f"{'Valeur obtenue':23s} "
        f"{'Valeur de référence':23s} "
        f"{'Écart relatif':17s} "
        f"{'Validation':12s}"
    )

    print("-" * 90)

    for (
        case_name,
        error_spacetime,
        error_ref,
        relative_error,
        validated,
    ) in results:

        status = (
            "VALIDÉ"
            if validated
            else "ÉCHEC"
        )

        print(
            f"{case_name:20s} "
            f"{error_spacetime:<23.15e} "
            f"{error_ref:<23.15e} "
            f"{relative_error:<17.3e} "
            f"{status:12s}"
        )

    print("=" * 90)

    # ========================================================
    # ASSERTIONS
    # ========================================================

    for (
        case_name,
        error_spacetime,
        error_ref,
        relative_error,
        validated,
    ) in results:

        assert validated, (
            f"Échec du test {case_name} : "
            f"erreur relative = "
            f"{relative_error:.3e} "
            f"> tolérance = "
            f"{TOLERANCE:.3e}"
        )

    # ========================================================
    # VALIDATION GLOBALE
    # ========================================================

    print()
    print("=" * 115)
    print("TOUS LES CAS SONT VALIDÉS")
    print("=" * 115)


# ============================================================
# EXÉCUTION DIRECTE
# ============================================================

if __name__ == "__main__":
    test_quarter_annulus_dynamic_spacetime_l2()