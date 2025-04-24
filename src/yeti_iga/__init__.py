from .overlay.iga_model import IgaModel
from .overlay.patch import Patch
from .overlay.material import Material, ElasticMaterial
from .overlay.boundary_condition import BoundaryCondition
from .overlay.distributed_load import DistributedLoad
from .overlay.iga_optimization import IgaOptimization
from .overlay.refinement import Refinement

try:
    from importlib.metadata import version
except ImportError:
    from importlib_metadata import version  # compatibilty < Python 3.8

__version__ = version("yeti_iga")