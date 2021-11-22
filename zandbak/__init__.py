"""A sandbox python repository for testing GitHub features, CI pipelines, and python
packages."""
from typing import Iterable, Union

import numpy as np

__version__ = "0.2.11"

REAL_SCALAR = Union[int, float]
REAL_VECT = Union[REAL_SCALAR, Iterable, np.ndarray]
