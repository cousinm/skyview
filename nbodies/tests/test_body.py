from nbodies.body import Body
from numpy import sum


def test_body():
    """
    Test the body initialisation
    """
    b = Body()
    assert b.mass == 1
    assert sum(b.X) == 0.
    assert sum(b.V) == 0.
