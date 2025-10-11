import gray_scott

# ---------------------
# Test functions
# ---------------------
"""
simple test to check that the gray_scott module can be imported
"""


def test_import():
    assert hasattr(gray_scott, "run_simulation")
