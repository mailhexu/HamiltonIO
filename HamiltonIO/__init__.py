from HamiltonIO.bandstructure import (
    BandPathData,
    ElectronicBandStructure,
    calculate_band_structure,
    make_band_path,
    plot_band_structure,
)
from HamiltonIO.hamiltonian import Hamiltonian

try:
    from HamiltonIO.siesta import SiestaHam
except ImportError:
    SiestaHam = None

try:
    from HamiltonIO.wannier import WannierHam
except ImportError:
    WannierHam = None

__version__ = "0.2.5"

# Build __all__ dynamically based on available modules
__all__ = [
    "Hamiltonian",
    "BandPathData",
    "ElectronicBandStructure",
    "calculate_band_structure",
    "make_band_path",
    "plot_band_structure",
    "__version__",
]
if SiestaHam is not None:
    __all__.append("SiestaHam")
if WannierHam is not None:
    __all__.append("WannierHam")
