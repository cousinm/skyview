from Nbodies.Nbodies import Nbodies
from copy import deepcopy
from numpy import array, sqrt, cos, sin
from math import pi
import matplotlib.pyplot as plt


def test_2bodies_problem():
    """
    The classical solved 2 bodies problem
    """
    #
    Nstep = 1000
    dt = 1.e-1
    #
    # Create a system with 2 bodies in the Oxy plan (z = 0.)
    N = 2
    sys = Nbodies(N)
    #
    # body 0
    # The mass of the body A is 1.
    mA = 1.
    sys.bodies[0].mass = mA
    # The mass of the body B is 1.
    mB = 1.
    sys.bodies[1].mass = mB
    #
    # Conics parameters
    r0 = 1.
    e = 0.5
    p = (e + 1.)*r0
    C = sqrt(p*(mA + mB))
    theta0 = pi/5.
    #
    # coordinates of body A are:
    rA = -mB / (mA + mB) * r0
    xA = rA * cos(theta0)
    yA = rA * sin(theta0)
    zA = 0.
    sys.bodies[0].X = array([xA, yA, zA])
    # velocities of body A are:
    # dr0_dt = 0 !!
    dtheta0A_dt = C / r0**2.
    VxA = - rA * dtheta0A_dt * sin(theta0)
    VyA = rA * dtheta0A_dt * cos(theta0)
    VzA = 0.
    sys.bodies[0].V = array([VxA, VyA, VzA])
    # coordinates of body B are:
    rB = mA / (mA + mB) * r0
    xB = rB * cos(theta0)
    yB = rB * sin(theta0)
    zB = 0.
    sys.bodies[1].X = array([xB, yB, zB])
    # velocities of body B are:
    # dr0_dt = 0 !!
    dtheta0B_dt = C / r0**2.
    VxB = - rB * dtheta0B_dt * sin(theta0)
    VyB = rB * dtheta0B_dt * cos(theta0)
    VzB = 0.
    sys.bodies[1].V = array([VxB, VyB, VzB])
    #
    # save initial positions
    X = [[deepcopy(b.X) for b in sys.bodies]]
    #
    # compute mass centre of the system
    sys.mass_center()
    C = [deepcopy(sys.C)]
    #
    # Plot theorical trajectories and initial positions
    fig = plt.figure()
    ax = plt.axes()
    #
    Npts = 350
    theta = array(range(Npts))*2*pi/(Npts-1)
    # M
    rM_theo = p / (1 + e*cos(theta-theta0))
    xM = rM_theo * cos(theta)
    yM = rM_theo * sin(theta)
    # A
    rA_theo = - mB / (mA + mB) * p / (1. + e*cos(theta-theta0))
    xA = rA_theo * cos(theta)
    yA = rA_theo * sin(theta)
    #
    # B
    rB_theo = mA / (mA + mB) * p / (1. + e*cos(theta-theta0))
    xB = rB_theo * cos(theta)
    yB = rB_theo * sin(theta)
    #
    ax.plot(xA, yA, linestyle='--', color='blue', label='A')
    ax.plot(xB, yB, linestyle='--', color='green', label='B')
    ax.plot(xM, yM, linestyle='--', color='grey', label='M')
    #
    # initial positions
    ax.plot(X[0][0][0], X[0][0][1],
            linestyle='None', marker='o', color='blue')
    ax.plot(X[0][1][0], X[0][1][1],
            linestyle='None', marker='o', color='green')
    ax.arrow(X[0][0][0], X[0][0][1], VxA, VyA,
             head_width=0.5, head_length=0.1,
             fc='blue', ec='blue')
    ax.arrow(X[0][1][0], X[0][1][1], VxB, VyB,
             head_width=0.5, head_length=0.1,
             fc='green', ec='green')
    #
    ax.axis('equal')
    #
    for i in range(Nstep):
        #
        sys.evolve(dt)
        #
        X.append([deepcopy(b.X) for b in sys.bodies])
        C.append(deepcopy(sys.C))

    # color
    color = ['blue', 'green', 'gray']
    for i in range(Nstep):
        #
        for b in range(len(X[i+1])):
            ax.plot(X[i+1][b][0], X[i+1][b][1],
                    linestyle='None',
                    marker='o',
                    color=color[b],
                    markersize=3.)
            ax.plot(C[i+1][0], C[i+1][1],
                    linestyle='None',
                    marker='x',
                    color='grey',
                    markersize=2.)
    #
    plt.savefig('trajectories.pdf')
    plt.close()
