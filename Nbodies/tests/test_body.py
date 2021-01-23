from Nbodies.Body import Body
from numpy import max, min, abs, sum


def test_Body():
    """
    Test the body initialisation
    """
    b = Body()
    assert b.mass == 1
    assert ((max(abs(b.X)) <= 0.5) and (min(abs(b.X)) >= 0.))
    assert ((max(abs(b.V)) <= 0.5) and (min(abs(b.V)) >= 0.))
