from numpy import zeros, sqrt, sum


class Body:
    """
    A class for the individual body description
    """
    def __init__(self):
        #
        # mass of the body
        self.mass = 1.
        #
        # The position vecteor in the 3D space
        self.X = zeros(3)
        #
        # The velocity vector
        self.V = zeros(3)

    def distance(self, bj):
        """
        return the euclidian distance between two bodies
        """
        return sqrt(sum((self.X-bj.X)**2))

    def acceleration(self, bj):
        """
        Compute the impact onto the acceleration of the current body (self)
        induce by the body bj
        """
        rij = self.distance(bj)
        aij = - bj.mass*(self.X - bj.X) / max(rij**3, 1.e-8)
        aji = - aij*self.mass/bj.mass
        return (aij, aji)

    def evolve_X(self, dt):
        """
        Evolve the position vector
        """
        self.X += self.V * dt

    def evolve_V(self, dt, acc):
        """
        Evolve the velocity vector
        """
        self.V += acc * dt / 2.
