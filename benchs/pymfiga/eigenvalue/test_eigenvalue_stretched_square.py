# -*- coding: utf-8 -*-
"""
Created on Wed Aug 26 18:41:14 2026
Test de validation du pas de temps critique
Carré étiré
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np

from yeti_iga.pymfiga.iga.boundary import (
    BoundaryCondition,
    DOF,
)
from yeti_iga.pymfiga.iga.geometry import GeomdlGenerator, SinglePatch
from yeti_iga.pymfiga.iga.single_model import ExplicitDynamicsModel
from yeti_iga.pymfiga.common.material import LinearElasticity
from yeti_iga.pymfiga.common.physics import EigenProblem


# ============================================================
# PARAMÈTRES DU CAS TEST
# ============================================================

DEGREE = 1
NBEL = 2

E = 2.1e11
NU = 0.3
RHO = 7850.0

STRETCHED_SQUARE_XY = np.array(
    [
        [0.0, 0.0],
        [1.0, 0.0],
        [4.0, 4.0],
        [0.0, 1.0],
    ],
    dtype=float,
)

QUADCLASS = "wq"
QUADTYPE = "2"

WHICH = "LM"
K_EIGEN = 1


# ============================================================
# VALEURS DE RÉFÉRENCE
# ============================================================

REFERENCE_DT_CRIT = {

    "MC_ARPACK_NONE":
        4.165491045170976e-05,

    "MC_ARPACK_FASTDIAG":
        4.165491045170970e-05,

    "MC_ARPACK_SCALED_MASS":
        4.165491044982996e-05,

    "MC_POWER_NONE":
        4.165491170929749e-05,

    "MC_POWER_FASTDIAG":
        4.165491170929749e-05,

    "MC_POWER_SCALED_MASS":
        4.165491170929755e-05,

    "ML_ARPACK":
        8.687866974858425e-05,

    "ML_POWER":
        8.687867554007668e-05,
}


# ============================================================
# CAS À TESTER
# ============================================================

CALC_CASES = [

    # Masse consistante - ARPACK

    {
        "mass_type": "consistent_mass",
        "solver": "arpack",
        "use_preconditioner": False,
        "preconditioner_type": None,
        "case_name": "MC_ARPACK_NONE",
    },

    {
        "mass_type": "consistent_mass",
        "solver": "arpack",
        "use_preconditioner": True,
        "preconditioner_type": "fastdiag",
        "case_name": "MC_ARPACK_FASTDIAG",
    },

    {
        "mass_type": "consistent_mass",
        "solver": "arpack",
        "use_preconditioner": True,
        "preconditioner_type": "scaled_mass",
        "case_name": "MC_ARPACK_SCALED_MASS",
    },

    # Masse consistante - puissance itérée

    {
        "mass_type": "consistent_mass",
        "solver": "power",
        "use_preconditioner": False,
        "preconditioner_type": None,
        "case_name": "MC_POWER_NONE",
    },

    {
        "mass_type": "consistent_mass",
        "solver": "power",
        "use_preconditioner": True,
        "preconditioner_type": "fastdiag",
        "case_name": "MC_POWER_FASTDIAG",
    },

    {
        "mass_type": "consistent_mass",
        "solver": "power",
        "use_preconditioner": True,
        "preconditioner_type": "scaled_mass",
        "case_name": "MC_POWER_SCALED_MASS",
    },

    # Masse diagonalisée

    {
        "mass_type": "diagonal_mass",
        "solver": "arpack",
        "use_preconditioner": False,
        "preconditioner_type": None,
        "case_name": "ML_ARPACK",
    },

    {
        "mass_type": "diagonal_mass",
        "solver": "power",
        "use_preconditioner": False,
        "preconditioner_type": None,
        "case_name": "ML_POWER",
    },
]


# ============================================================
# CONSTRUCTION DU MODÈLE
# ============================================================

def build_model(degree, nbel, mass_type):

    geometry = GeomdlGenerator(
        filename="trapezium",
        geo_args={
            "degree": degree,
            "nbel": nbel,
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

    boundary = BoundaryCondition(
        nbctrlpts=patch.nbctrlpts,
        dofs=(DOF.UX, DOF.UY),
    )

    model = ExplicitDynamicsModel(
        material,
        patch,
        boundary,
        mass_type=mass_type,
    )

    return model


# ============================================================
# CALCUL DU PAS DE TEMPS CRITIQUE
# ============================================================

def compute_dt_crit(calc_case):

    model = build_model(
        degree=DEGREE,
        nbel=NBEL,
        mass_type=calc_case["mass_type"],
    )

    eigenvalues, _ = EigenProblem().solve(
        model=model,
        use_preconditioner=calc_case["use_preconditioner"],
        preconditioner_type=calc_case["preconditioner_type"],
        solver=calc_case["solver"],
        which=WHICH,
        k=K_EIGEN,
    )

    eigenvalues = np.asarray(eigenvalues)
    eigenvalues = np.real(eigenvalues)

    if eigenvalues.size == 0:
        raise RuntimeError(
            "Aucune valeur propre calculée."
        )

    if not np.all(np.isfinite(eigenvalues)):
        raise RuntimeError(
            "Valeur propre contenant NaN ou Inf."
        )

    lambda_max = float(
        np.max(eigenvalues)
    )

    if lambda_max <= 0.0:
        raise RuntimeError(
            "La valeur propre maximale doit être positive."
        )

    dt_crit = 2.0 / np.sqrt(
        lambda_max
    )

    return dt_crit


# ============================================================
# TEST DE VALIDATION
# ============================================================

def test_stretched_square_eigenvalues():

    TOLERANCE = 1e-6

    results = []

    for calc_case in CALC_CASES:

        case_name = calc_case["case_name"]

        # Calcul du pas de temps critique
        dt_crit = compute_dt_crit(calc_case)

        # Valeur de référence
        dt_crit_ref = REFERENCE_DT_CRIT[case_name]

        # Écart relatif
        relative_error = abs(
            (dt_crit - dt_crit_ref)
            / dt_crit_ref
        )

        # Validation
        validated = relative_error < TOLERANCE

        results.append(
            (
                case_name,
                dt_crit,
                dt_crit_ref,
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
        "RÉSULTATS DE VALIDATION DU PAS DE TEMPS CRITIQUE "
        "- CARRÉ ÉTIRÉ"
    )
    print("=" * 90)

    print(
        f"degree = {DEGREE}, "
        f"nbel = {NBEL}, "
        f"tolérance relative = {TOLERANCE:.1e}"
    )

    print("-" *90)

    print(
        f"{'Cas':25s} "
        f"{'Valeur obtenue':23s} "
        f"{'Valeur de référence':23s} "
        f"{'Écart relatif':17s} "
        f"{'Validation':12s}"
    )

    print("-" *90)

    for (
        case_name,
        dt_crit,
        dt_crit_ref,
        relative_error,
        validated,
    ) in results:

        status = "VALIDÉ" if validated else "ÉCHEC"

        print(
            f"{case_name:25s} "
            f"{dt_crit:<23.15e} "
            f"{dt_crit_ref:<23.15e} "
            f"{relative_error:<17.3e} "
            f"{status:12s}"
        )

    print("=" * 90)

    # ========================================================
    # ASSERTIONS
    # ========================================================

    for (
        case_name,
        dt_crit,
        dt_crit_ref,
        relative_error,
        validated,
    ) in results:

        assert validated, (
            f"Échec du test {case_name} : "
            f"erreur relative = {relative_error:.3e} "
            f"> tolérance = {TOLERANCE:.3e}"
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
    test_stretched_square_eigenvalues()