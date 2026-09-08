from yeti_iga.pymfiga.common.base.enum import ParametricDirection
from yeti_iga.pymfiga.common.numerics.quadrature_rules import IGAQuadratureRule, StandardGauss
from typing import List, Dict, Tuple
from abc import ABC, abstractmethod
import logging

logger = logging.getLogger("SRC.FASTDIAG")


class Template(ABC):
    """
    This class implements the fast-diagonalization preconditioner for a single patch.
    """

    @property
    @abstractmethod
    def nonzeros_by_dir(self) -> Dict[ParametricDirection, int]:
        pass

    @property
    @abstractmethod
    def indices_by_dir(self) -> Dict[Tuple[int, ParametricDirection], List[int]]:
        pass

    @staticmethod
    def rewrite_quadrature(
        quadrule_list: List[IGAQuadratureRule],
    ) -> List[StandardGauss]:
        def verify(q: IGAQuadratureRule) -> StandardGauss:
            if isinstance(q, StandardGauss) and q.quadrature_type == "legendre":
                return q

            logger.info("Rewrite quadrature rule to standard Gauss-Legendre")
            q = StandardGauss(
                degree=q.degree,
                knotvector=q.knotvector,
                quadtype="legendre",
                is_periodic=q.is_periodic,
            )
            q.export_quadrature_rules()
            return q

        assert isinstance(quadrule_list, list)
        return [verify(q) for q in quadrule_list]
