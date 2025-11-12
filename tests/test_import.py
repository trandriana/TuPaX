import imex, tupax

# ---------------------
# Test functions
# ---------------------
"""
simple test to check that the imex and tupax modules can be imported
"""


def test_import():
    assert hasattr(imex, "run_simulation")
    assert hasattr(tupax, "make_video")
