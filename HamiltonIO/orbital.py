from dataclasses import dataclass
import numpy as np


@dataclass
class Orbital:
    name: str = "Unamed"
    n: int = 0
    l : int = 0
    m : int = 0
    atom : Union[int, str] = None
    xred : np.ndarray = np.array([0.0, 0.0, 0.0])




