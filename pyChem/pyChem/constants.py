"""File containing all constants used for pyChem"""

import numpy as np

from ndsl.dsl.typing import Float


# MAPL_UNDEF is set to 1E15 in the Fortran
# We keep it as is for now to match 11.5.2 GEOS
MAPL_UNDEF = Float(1e15)
