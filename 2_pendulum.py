#!/usr/bin/env python
#############################################################################
# course:   Numerische Methoden D-PHYS
# exercise: assignment 2
# author:   Thomas Diggelmann <thomas.diggelmann@student.ethz.ch>
# date:     15.03.2015
#############################################################################
from numpy import array, sin, pi, shape, sqrt, linspace
import matplotlib.pyplot as plt
from ode45 import ode45
from scipy.optimize import fsolve
from sys import exit


def PendulumODE(y, l, g):
    """PendulumODE return the right-hand side of the math. pendulum"""
    dydt = array([y[1], -g * sin(y[0]) / l])
    return dydt


def IntegratePendulum(phi0, tEnd=1.8, l=0.6, g=9.81, flag=False):
    """IntegratePendulum solve the mathematical pendulum with ode45

    Input: phi0   ... initial condition
           tEnd   ... end time
           flag   ... flag == False return complete solution: (phi, phi', t)
                      flag == True  return solution at endtime only: phi(tEnd)
    """
    ##########################################
    #                                        #
    # Benutzen Sie ode45 zur Zeitintegration #
    #                                        #
    ##########################################
    #pass
    tspan = array([0, tEnd])
    y0 = array([phi0, 0])
    t, y = ode45(lambda t, y: PendulumODE(y, l, g), tspan, y0)

    if flag:
        return t[-1], y[-1][:]
    else:
        return t, y

def simpson(f, a, b, N):
    r"""
    Simpson quadrature of function f from a to b with N subintervals.

    f:     Function f(x)
    a, b:  Bounds of the integration interval
    N:     Number of subintervals
    """
    ##########################################################
    #                                                        #
    # Implementieren Sie eine zusammengesetzte Simpson Regel #
    #                                                        #
    ##########################################################
    x, h = linspace(a, b, 2*N+1, retstep=True)
    I = h/3.0 * sum(f(x[:-2:2]) + 4.0*f(x[1:-1:2]) + f(x[2::2]))
    return I

def elliptic_k_quad(k, N=10):
    """Compute the elliptic integral K(k)
    via composite Simpson qudrature.

    Input: k      ... function argument
           N      ... number of intervals for Simpson
    Output: K(k)  ... function value
    """
    #############################################################
    #                                                           #
    # Berechnen Sie das elliptische Integral mittels Quadrartur #
    #                                                           #
    #############################################################

    K = lambda phi: 1 / sqrt(1-k**2*sin(phi)**2)
    return simpson(K, 0, pi/2, N)


def agm(x, y, nit=5):
    """Computes the arithmetic-geometric mean of two numbers x and y.

    Input: x, y  ... two numbers
           nit   ... number of iteration
    Output: agm(x,y)
    """
    ########################################
    #                                      #
    # Implementieren Sie die AGM Iteration #
    #                                      #
    ########################################
    x0, y0 = x, y
    for _ in xrange(5):
        x = (x0 + y0) / 2
        y = sqrt(x0*y0)
        x0, y0 = x, y

    return x

def elliptic_k_agm(k, nit=5):
    """Compute the elliptic integral K(k)
    via arithmetic-geometric mean iteration.

    Input: k      ... function argument
            nit   ... number of iteration
    Output: K(k)  ... function value
    """
    ################################################################
    #                                                              #
    # Berechnen Sie das elliptische Integral mittels AGM Iteration #
    #                                                              #
    ################################################################
    return pi/2 * 1 / agm(1-k, 1+k)




if __name__ == '__main__':
    from time import time

    # Parameter values
    tEnd = 1.8
    l = 0.6
    g = 9.81


    # Unteraufgabe c)
    a1 = 0.8
    a2 = 0.99

    ####################################################
    #                                                  #
    # Loesen Sie das Anfangswertproblem fuer a1 und a2 #
    #                                                  #
    ####################################################

    # Solve ODE for extreme values to check the period length:
    t1, y1 = IntegratePendulum(a1*pi/2)
    t2, y2 = IntegratePendulum(a2*pi/2)

    plt.figure()
    ax1 = plt.subplot(211)
    plt.plot(t1, y1.transpose()[0], 'b', label=r'$f \, | \, f_0 = %1.2f \, \frac{\pi}{2}$' % a1)
    plt.plot(t1, y1.transpose()[1], 'b--', label=r'$\frac{df}{dt} \, | \, f_0 = %1.2f \, \frac{\pi}{2}$' % a1)
    plt.plot(t2, y2.transpose()[0], 'r', label=r'$f \, | \, f_0 = %1.2f \, \frac{\pi}{2}$' % a2)
    plt.plot(t2, y2.transpose()[1], 'r--', label=r'$\frac{df}{dt} \, | \, f_0 = %1.2f \, \frac{\pi}{2}$' % a2)
    plt.title('Solutions to the pendulum ODE for different initial amplitudes')
    plt.xlabel(r'$t$')
    plt.legend(loc='upper left')
    plt.grid(True)

    # Unteraufgabe d)

    # zero finding:
    print("Initial value via time-integration")

    #################################################
    #                                               #
    # Berechnen Sie hier die Anfangsauslenkung phiT #
    #                                               #
    #################################################
    T = 1.8 # the period to solve for
    F = lambda phi0: IntegratePendulum(phi0, T/4, l, g, True)[1][0]
    starttime = time()
    phiT = fsolve(F, a1*pi/2) # use extreme 1 as starting point
    endtime = time()
    print('needed %f seconds' % (endtime-starttime))
    print(phiT)

    # compute complete solution with period = tEnd for plotting
    t3, y3 = IntegratePendulum(phiT, tEnd, l, g, False)
    y3 = y3.transpose()

    ax1 = plt.subplot(212)
    plt.plot(t3, y3[0], '-', label=r'$f$')
    plt.plot(t3, y3[1], 'b--', label=r'$\frac{df}{dt}$')
    plt.title(r'Solution to the pendulum ODE for $f_0 = %1.3f \, \frac{\pi}{2} $' % (phiT*2/pi))
    plt.xlabel(r'$t$')
    plt.legend(loc='upper left')
    plt.grid(True)
    plt.tight_layout()
    plt.savefig('Pendulum.pdf')
    plt.show()


    # Unteraufgabe e)

    print("Initial value via exakt formula for T")
    print("Elliptic Integral by Quadrature")

    #################################################
    #                                               #
    # Berechnen Sie hier die Anfangsauslenkung phiT #
    #                                               #
    #################################################
    F = lambda phi0: sqrt(l/g)*elliptic_k_quad(sin(phi0/2)) - T/4
    starttime = time()
    phiT = fsolve(F, a1*pi/2)
    endtime = time()
    print(phiT)
    print('needed %f seconds' % (endtime-starttime))


    # Unteraufgabe f)

    print("Initial value via exakt formula for T")
    print("Elliptic Integral by Arithmetic-Geometric Mean")

    #################################################
    #                                               #
    # Berechnen Sie hier die Anfangsauslenkung phiT #
    #                                               #
    #################################################
    F = lambda phi0: sqrt(l/g)*elliptic_k_agm(sin(phi0/2)) - T/4
    starttime = time()
    phiT = fsolve(F, a1*pi/2)
    endtime = time()
    print(phiT)
    print('needed %f seconds' % (endtime-starttime))
