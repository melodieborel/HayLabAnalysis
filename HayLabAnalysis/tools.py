import os
import configparser
import pickle
import ast

import numpy as np

from ipyfilechooser import FileChooser
import ipywidgets as widgets
from IPython.display import display
from IPython import get_ipython

from .localConfigurations import localConf

"""Adapter shim for `HayLabAnalysis.tools`.

This module keeps the historical public API while delegating implementations to
`HayLabAnalysis.core.tools`. It is intentionally minimal to make the migration
invisible to consumers. New code should import from `HayLabAnalysis.core`.
"""

from importlib import import_module
import warnings

_core = import_module("HayLabAnalysis.core.tools")

# Re-export commonly used names from the core implementation
__all__ = [
    "superCleanPlot",
    "convertTheoricIndex2realTime",
    "find_nearest",
    "getPathComponent",
    "color",
]

superCleanPlot = _core.superCleanPlot
convertTheoricIndex2realTime = _core.convertTheoricIndex2realTime
find_nearest = _core.find_nearest
getPathComponent = _core.getPathComponent
color = _core.color

def _deprecation_notice():
    warnings.warn(
        "Module-level code has moved to HayLabAnalysis.core.tools; this shim will be removed in a future release.",
        DeprecationWarning,
        stacklevel=2,
    )

# Emit a single deprecation when the module is first imported
_deprecation_notice()
