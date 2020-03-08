from nbodies.chain import Chain
from numpy import sum


def test_chain():
    """
    Test the body-chain initialisation
    """
    #
    # Create a system with 200 bodies
    N = 200
    bsys = Chain(N)
    #
    #  test total mass
    assert bsys.tot_mass == float(N)
    #
    # test mass center
    assert sum(bsys.C) == 0.
    #
    #  test position and velocity bondaries
    for b in bsys.bodies:
        assert b.mass == 1
        assert sum(bsys.C) == 0.
        assert sum(bsys.C) == 0.
