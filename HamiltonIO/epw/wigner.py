"""
Wigner-Seitz data structure for EPW calculations.

This module provides a dataclass to handle Wigner-Seitz data for electrons,
phonons, and electron-phonon interactions as used in EPW calculations.
"""

from dataclasses import dataclass
from typing import List, Tuple
import numpy as np


@dataclass
class WignerData:
    """
    Wigner-Seitz data container for EPW calculations.

    This dataclass stores the Wigner-Seitz vectors, degeneracies, and lengths
    for electrons (k), phonons (q), and electron-phonon (g) interactions.

    Attributes:
        dims: Number of Wannier functions (if use_ws, otherwise 1)
        dims2: Number of atoms (if use_ws, otherwise 1)

        # Electron data (k)
        nrr_k: Number of WS vectors for electrons
        irvec_k: WS vectors for electrons, shape (nrr_k, 3)
        ndegen_k: WS degeneracies for electrons, shape (nrr_k, dims, dims)
        wslen_k: WS vector lengths for electrons, shape (nrr_k,)

        # Phonon data (q)
        nrr_q: Number of WS vectors for phonons
        irvec_q: WS vectors for phonons, shape (nrr_q, 3)
        ndegen_q: WS degeneracies for phonons, shape (nrr_q, dims2, dims2)
        wslen_q: WS vector lengths for phonons, shape (nrr_q,)

        # Electron-phonon data (g)
        nrr_g: Number of WS vectors for electron-phonon
        irvec_g: WS vectors for electron-phonon, shape (nrr_g, 3)
        ndegen_g: WS degeneracies for electron-phonon, shape (nrr_g, dims, dims2)
        wslen_g: WS vector lengths for electron-phonon, shape (nrr_g,)
    """

    # Dimensions
    dims: int
    dims2: int

    # Electron data (k)
    nrr_k: int
    irvec_k: np.ndarray  # shape (nrr_k, 3)
    ndegen_k: np.ndarray  # shape (nrr_k, dims, dims)
    wslen_k: np.ndarray  # shape (nrr_k,)

    # Phonon data (q)
    nrr_q: int
    irvec_q: np.ndarray  # shape (nrr_q, 3)
    ndegen_q: np.ndarray  # shape (nrr_q, dims2, dims2)
    wslen_q: np.ndarray  # shape (nrr_q,)

    # Electron-phonon data (g)
    nrr_g: int
    irvec_g: np.ndarray  # shape (nrr_g, 3)
    ndegen_g: np.ndarray  # shape (nrr_g, dims, dims2)
    wslen_g: np.ndarray  # shape (nrr_g,)

    def __post_init__(self):
        """Validate array shapes after initialization."""
        # Convert to numpy arrays if needed
        self.irvec_k = np.asarray(self.irvec_k)
        self.ndegen_k = np.asarray(self.ndegen_k)
        self.wslen_k = np.asarray(self.wslen_k)

        self.irvec_q = np.asarray(self.irvec_q)
        self.ndegen_q = np.asarray(self.ndegen_q)
        self.wslen_q = np.asarray(self.wslen_q)

        self.irvec_g = np.asarray(self.irvec_g)
        self.ndegen_g = np.asarray(self.ndegen_g)
        self.wslen_g = np.asarray(self.wslen_g)

        # Validate shapes
        assert self.irvec_k.shape == (
            self.nrr_k,
            3,
        ), f"irvec_k shape mismatch: {self.irvec_k.shape} != ({self.nrr_k}, 3)"
        assert self.ndegen_k.shape == (
            self.nrr_k,
            self.dims,
            self.dims,
        ), f"ndegen_k shape mismatch"
        assert self.wslen_k.shape == (self.nrr_k,), f"wslen_k shape mismatch"

        assert self.irvec_q.shape == (self.nrr_q, 3), f"irvec_q shape mismatch"
        assert self.ndegen_q.shape == (
            self.nrr_q,
            self.dims2,
            self.dims2,
        ), f"ndegen_q shape mismatch"
        assert self.wslen_q.shape == (self.nrr_q,), f"wslen_q shape mismatch"

        assert self.irvec_g.shape == (self.nrr_g, 3), f"irvec_g shape mismatch"
        assert self.ndegen_g.shape == (
            self.nrr_g,
            self.dims,
            self.dims2,
        ), f"ndegen_g shape mismatch"
        assert self.wslen_g.shape == (self.nrr_g,), f"wslen_g shape mismatch"

    @classmethod
    def from_file(cls, filename: str = "wigner.fmt") -> "WignerData":
        """
        Read Wigner-Seitz data from file format written by EPW.

        Args:
            filename: Path to the wigner data file

        Returns:
            WignerData instance
        """
        with open(filename, "r") as f:
            # Read header
            header = f.readline().strip().split()
            nrr_k, nrr_q, nrr_g, dims, dims2 = map(int, header)

            # Initialize arrays
            irvec_k = np.zeros((nrr_k, 3), dtype=int)
            wslen_k = np.zeros(nrr_k, dtype=float)
            ndegen_k = np.zeros((nrr_k, dims, dims), dtype=int)

            irvec_q = np.zeros((nrr_q, 3), dtype=int)
            wslen_q = np.zeros(nrr_q, dtype=float)
            ndegen_q = np.zeros((nrr_q, dims2, dims2), dtype=int)

            irvec_g = np.zeros((nrr_g, 3), dtype=int)
            wslen_g = np.zeros(nrr_g, dtype=float)
            ndegen_g = np.zeros((nrr_g, dims, dims2), dtype=int)

            # Read electron data (k)
            for ir in range(nrr_k):
                line = f.readline().strip().split()
                irvec_k[ir, :] = [int(line[0]), int(line[1]), int(line[2])]
                wslen_k[ir] = float(line[3])

                for iw in range(dims):
                    line = f.readline().strip().split()
                    ndegen_k[ir, iw, :] = [int(x) for x in line]

            # Read phonon data (q)
            for ir in range(nrr_q):
                line = f.readline().strip().split()
                irvec_q[ir, :] = [int(line[0]), int(line[1]), int(line[2])]
                wslen_q[ir] = float(line[3])

                for na in range(dims2):
                    line = f.readline().strip().split()
                    ndegen_q[ir, na, :] = [int(x) for x in line]

            # Read electron-phonon data (g)
            for ir in range(nrr_g):
                line = f.readline().strip().split()
                irvec_g[ir, :] = [int(line[0]), int(line[1]), int(line[2])]
                wslen_g[ir] = float(line[3])

                for iw in range(dims):
                    line = f.readline().strip().split()
                    ndegen_g[ir, iw, :] = [int(x) for x in line]

        return cls(
            dims=dims,
            dims2=dims2,
            nrr_k=nrr_k,
            irvec_k=irvec_k,
            ndegen_k=ndegen_k,
            wslen_k=wslen_k,
            nrr_q=nrr_q,
            irvec_q=irvec_q,
            ndegen_q=ndegen_q,
            wslen_q=wslen_q,
            nrr_g=nrr_g,
            irvec_g=irvec_g,
            ndegen_g=ndegen_g,
            wslen_g=wslen_g,
        )

    def to_file(self, filename: str = "wigner.fmt") -> None:
        """
        Write Wigner-Seitz data to file in EPW format.

        Args:
            filename: Output file path
        """
        with open(filename, "w") as f:
            # Write header
            f.write(
                f"{self.nrr_k} {self.nrr_q} {self.nrr_g} {self.dims} {self.dims2}\n"
            )

            # Write electron data (k)
            for ir in range(self.nrr_k):
                f.write(
                    f"{self.irvec_k[ir, 0]:6d}{self.irvec_k[ir, 1]:6d}{self.irvec_k[ir, 2]:6d}"
                )
                f.write(f"{self.wslen_k[ir]:26.17E}\n")

                for iw in range(self.dims):
                    degeneracies = " ".join(
                        str(self.ndegen_k[ir, iw, j]) for j in range(self.dims)
                    )
                    f.write(f"{degeneracies}\n")

            # Write phonon data (q)
            for ir in range(self.nrr_q):
                f.write(
                    f"{self.irvec_q[ir, 0]:6d}{self.irvec_q[ir, 1]:6d}{self.irvec_q[ir, 2]:6d}"
                )
                f.write(f"{self.wslen_q[ir]:26.17E}\n")

                for na in range(self.dims2):
                    degeneracies = " ".join(
                        str(self.ndegen_q[ir, na, j]) for j in range(self.dims2)
                    )
                    f.write(f"{degeneracies}\n")

            # Write electron-phonon data (g)
            for ir in range(self.nrr_g):
                f.write(
                    f"{self.irvec_g[ir, 0]:6d}{self.irvec_g[ir, 1]:6d}{self.irvec_g[ir, 2]:6d}"
                )
                f.write(f"{self.wslen_g[ir]:26.17E}\n")

                for iw in range(self.dims):
                    degeneracies = " ".join(
                        str(self.ndegen_g[ir, iw, j]) for j in range(self.dims2)
                    )
                    f.write(f"{degeneracies}\n")

    def summary(self) -> str:
        """Return a summary string of the Wigner data."""
        return (
            f"WignerData Summary:\n"
            f"  Dimensions: {self.dims} Wannier functions, {self.dims2} atoms\n"
            f"  Electrons (k): {self.nrr_k} WS vectors\n"
            f"  Phonons (q): {self.nrr_q} WS vectors\n"
            f"  Electron-phonon (g): {self.nrr_g} WS vectors"
        )

    def find_origin_indices(self) -> Tuple[int, int, int]:
        """
        Find the indices of the origin vectors (0,0,0) for k, q, and g.

        Returns:
            Tuple of (k_origin_idx, q_origin_idx, g_origin_idx)
            Returns -1 if origin vector not found
        """
        k_origin = np.where(np.all(self.irvec_k == 0, axis=1))[0]
        q_origin = np.where(np.all(self.irvec_q == 0, axis=1))[0]
        g_origin = np.where(np.all(self.irvec_g == 0, axis=1))[0]

        return (
            k_origin[0] if len(k_origin) > 0 else -1,
            q_origin[0] if len(q_origin) > 0 else -1,
            g_origin[0] if len(g_origin) > 0 else -1,
        )

