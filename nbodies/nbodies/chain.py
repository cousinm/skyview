from nbodies.body import Body
from multiprocessing import Pool, cpu_count, sharedctypes
from numpy import zeros, array, ctypeslib
from math import sqrt


class Chain:
    """
    A class to run through all bodies and compute evolution
    """

    def __init__(self, N):
        """
        Initialisation of a body chain
        """
        # number of cpu detected
        self.ncpus = int(cpu_count() / int(2))
        #
        # initialize bodies
        with Pool(processes=self.ncpus) as p:
            self.bodies = [p.apply_async(Body, args=()).get()
                           for i in range(N)]
        #
        # acceleration
        self.acc = zeros((N, 3))
        #
        # total mass of the system
        self.tot_mass = sum(array([b.mass for b in self.bodies]))
        #
        self.C = zeros(3)

    def _acceleration(self, i):
        """
        Compute acceleration
        """
        #
        Acc_tmp = ctypeslib.as_array(data)
        for j in range(i+1, len(self.bodies)):
            aij, aji = self.bodies[i].acceleration(self.bodies[j])
            Acc_tmp[i] += aij
            Acc_tmp[j] += aji

    def evaluate(self):
        """
        Evaluate the sys
        """
        #
        # acceleration
        global data
        #
        # initialisation of accelerations
        tmp = ctypeslib.as_ctypes(zeros((len(self.bodies), 3)))
        data = sharedctypes.RawArray(tmp._type_, tmp)
        #
        with Pool(processes=min(self.ncpus, len(self.bodies)-1)) as p:
            p.map(self._acceleration, range(len(self.bodies)))
        self.acc = ctypeslib.as_array(data)
        data = None  # Reset data shared array
        #
        # Equilibrium
        self.Ek = self.kinetic_energy()
        #
        #  self.Ep = self.potential_energy()
        #
        self.C = self.mass_center()

    def kinetic_energy(self):
        """
        Return the total kinetic enrgy of the system
        """
        Ek = 0.
        for b in self.bodies:
            Ek = Ek + 5.e-1*b.mass * sqrt(sum(b.V**2))
        return Ek

    def mass_center(self):
        """
        Return the center of mass of the system
        """
        C = zeros(3)  # Init
        for b in self.bodies:
            C += b.mass * b.X
        C /= self.tot_mass
        return C

    def evolve(self, dt):
        """
        Evole the sys assuming a leap-frog scheme
        """
        #
        N = len(self.bodies)
        #
        # evolve velocities until half time step and position on full time-step
        with Pool(processes=self.ncpus) as p:
            self.bodies = [p.apply_async(self._evolve_i,
                                         args=(i, dt)).get() for i in range(N)]
        #
        # New computation of accelerations and systeme parameters
        self.evaluate()
        #
        # Update of the velocities
        with Pool(processes=self.ncpus) as p:
            self.bodies = [p.apply_async(self._evolve_ii,
                                         args=(i, dt)).get() for i in range(N)]

    def _evolve_i(self, i, dt):
        """
        Evole the ith body assuming the first leap-frog step
        """
        self.bodies[i].evolve_V(dt, self.acc[i])
        self.bodies[i].evolve_X(dt)
        return self.bodies[i]

    def _evolve_ii(self, i, dt):
        """
        Evole the ith body according to the second leap-frog step
        """
        #
        self.bodies[i].evolve_V(dt, self.acc[i])
        return self.bodies[i]
