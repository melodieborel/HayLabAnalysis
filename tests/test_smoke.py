import math
from pathlib import Path

def test_import_and_find_nearest():
    # Basic smoke test to ensure adapter import works
    from HayLabAnalysis import tools
    arr = [0.0, 1.0, 2.0, 3.0]
    idx = tools.find_nearest(arr, 1.9)
    assert idx in (1, 2)
