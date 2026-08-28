# -*- coding: utf-8 -*-
"""
Created on Wed Aug 26 20:10:31 2026
Test de validation de la dynamique explicite
Erreur L2 espace-temps
Carré étiré
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


# ============================================================
# PARAMÈTRES DU CAS TEST
# ============================================================

DEGREE = 2
NBEL = 4

TIMESPAN = 2.0e-3

DT_MC = 1.346476119719685e-05
DT_ML = 3.146029665953533e-05

QUADCLASS = "wq"
QUADTYPE = "2"


# ============================================================
# GÉOMÉTRIE DU CARRÉ ÉTIRÉ
# ============================================================

STRETCHED_SQUARE_XY = np.array(
    [
        [0.0, 0.0],
        [1.0, 0.0],
        [4.0, 4.0],
        [0.0, 1.0],
    ],
    dtype=float,
)


# ============================================================
# VALEURS DE RÉFÉRENCE
# ============================================================

REFERENCE_L2_SPACETIME = {

    "MC_NONE":
        9.937100367360484e-01,

    "MC_FASTDIAG":
        9.937100363998788e-01,

    "MC_SCALED_MASS":
        9.937100360492571e-01,

    "ML":
        6.495661707984407e-01,
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


def temporal_function(t):

    return (
        1.0
        - np.cos(
            OMEGA * t
        )
    )


def spatial_amplitude(position):

    position = np.asarray(
        position,
        dtype=float,
    )

    x = position[0]
    y = position[1]

    left = (
        x
        / 4.0
    )

    bottom = (
        y
        / 4.0
    )

    right = (
        4.0
        - 4.0 * x
        + 3.0 * y
    ) / 5.0

    top = (
        3.0 * x
        - 4.0 * y
        + 4.0
    ) / 5.0

    return (
        U0
        * left**2
        * bottom**2
        * right**2
        * top**2
    )


def exact_displacement(position, t):

    amplitude = spatial_amplitude(
        position
    )

    g = temporal_function(
        t
    )

    ux = (
        amplitude
        * g
    )

    uy = (
        0.5
        * amplitude
        * g
    )

    return np.vstack(
        (
            ux,
            uy,
        )
    )


# ============================================================
# CONSTRUCTION SYMBOLIQUE DE LA SOLUTION
# ============================================================

x_sym, y_sym, t_sym = sp.symbols(
    "x y t",
    real=True,
)

left_sym = (
    x_sym
    / 4
)

bottom_sym = (
    y_sym
    / 4
)

right_sym = (
    4
    - 4 * x_sym
    + 3 * y_sym
) / 5

top_sym = (
    3 * x_sym
    - 4 * y_sym
    + 4
) / 5

A_sym = (
    U0
    * left_sym**2
    * bottom_sym**2
    * right_sym**2
    * top_sym**2
)

g_sym = (
    1
    - sp.cos(
        OMEGA * t_sym
    )
)

ux_sym = (
    A_sym
    * g_sym
)

uy_sym = (
    sp.Rational(1, 2)
    * A_sym
    * g_sym
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
    sp.Rational(
        1,
        2,
    )
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

    t_value = float(
        t
    )

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
        filename="trapezium",
        geo_args={
            "degree": DEGREE,
            "nbel": NBEL,
            "parameters": {
                "XY": STRETCHED_SQUARE_XY,
            },
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
    # Conditions aux limites BC1
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
                    ParametricDirection.ETA,
                "face":
                    BoundarySide.MIN,
                "dofs":
                    (DOF.UY,),
            },
            {
                "direction":
                    ParametricDirection.XI,
                "face":
                    BoundarySide.MIN,
                "dofs":
                    (DOF.UX,),
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
    # Conditions initiales
    # --------------------------------------------------------

    displacement_initial = np.zeros_like(
        initial_force,
        dtype=float,
    )

    velocity_initial = np.zeros_like(
        initial_force,
        dtype=float,
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

def test_stretched_square_dynamic_spacetime_l2():

    TOLERANCE = 1e-6

    results = []

    for calc_case in CALC_CASES:

        case_name = calc_case[
            "case_name"
        ]

        # Calcul de l'erreur L2 espace-temps
        error_spacetime = simulate(
            calc_case
        )

        # Valeur de référence
        error_ref = (
            REFERENCE_L2_SPACETIME[
                case_name
            ]
        )

        # Écart relatif
        relative_error = abs(
            (
                error_spacetime
                - error_ref
            )
            / error_ref
        )

        # Validation
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
        "- CARRÉ ÉTIRÉ"
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
    print("=" * 90)
    print("TOUS LES CAS SONT VALIDÉS")
    print("=" * 90)


# ============================================================
# EXÉCUTION DIRECTE
# ============================================================

if __name__ == "__main__":
    test_stretched_square_dynamic_spacetime_l2()