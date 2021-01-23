from Nbodies.Nbodies import Nbodies


def test_Nbodies():
    """
    Test the body-chain initialisation
    """
    #
    # Create a system with 200 bodies
    N = 200
    sys = Nbodies(N)
    #
    #  test total mass
    assert sys.tot_mass == float(N)
    #
    #  test position and velocity bondaries
    for b in sys.bodies:
        assert b.mass == 1
        assert ((max(abs(b.X)) <= 0.5) and (min(abs(b.X)) >= 0.))
        assert ((max(abs(b.V)) <= 0.5) and (min(abs(b.V)) >= 0.))
