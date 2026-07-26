from HamiltonIO.bandstructure import (
    BandPathData,
    ElectronicBandStructure,
    calculate_band_structure,
    make_band_path,
    plot_band_structure,
)
from HamiltonIO.builder import OrbitalSpec, TightBindingBuilder, TightBindingModelData
from HamiltonIO.hamiltonian import Hamiltonian

try:
    from HamiltonIO.siesta import SiestaHam
except ImportError:
    SiestaHam = None

try:
    from HamiltonIO.wannier import WannierHam
except ImportError:
    WannierHam = None

__version__ = "0.3.8"

# Build __all__ dynamically based on available modules
__all__ = [
    "Hamiltonian",
    "BandPathData",
    "ElectronicBandStructure",
    "calculate_band_structure",
    "make_band_path",
    "plot_band_structure",
    "TightBindingBuilder",
    "TightBindingModelData",
    "OrbitalSpec",
    "__version__",
]
if SiestaHam is not None:
    __all__.append("SiestaHam")
if WannierHam is not None:
    __all__.append("WannierHam")
