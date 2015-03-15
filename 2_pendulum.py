from numpy import array, sin, pi, shape, sqrt, linspace
import matplotlib.pyplot as plt
from ode45 import ode45
from scipy.optimize import fsolve


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
    I = 0.0
    ##########################################################
    #                                                        #
    # Implementieren Sie eine zusammengesetzte Simpson Regel #
    #                                                        #
    ##########################################################
    return I

    """
    h=(b-a)/n
    k=0.0
    x=a + h
    for i in range(1,n/2 + 1):
        k += 4*f(x)
        x += 2*h

    x = a + 2*h
    for i in range(1,n/2):
        k += 2*f(x)
        x += 2*h
    return (h/3)*(f(a)+f(b)+k)
    """


def elliptic_k_quad(k, N=10):
    """Compute the elliptic integral K(k)
    via composite Simpson qudrature.

    Input: k      ... function argument
           N      ... number of intervals for Simpson
    Output: K(k)  ... function value
    """
    K = 0.0
    #############################################################
    #                                                           #
    # Berechnen Sie das elliptische Integral mittels Quadrartur #
    #                                                           #
    #############################################################
    return K


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
    pass


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
    pass




if __name__ == '__main__':
    from time import time

    # Parameter values
    tEnd = 1.8
    l = 0.6
    g = 9.81


    # Unteraufgabe c)

    a1 = 0.8
    a2 = 0.99

    # Solve ODE for extreme values to check the period length:
    t1 = 0.0
    y1 = 0.0
    t2 = 0.0
    y2 = 0.0
    ####################################################
    #                                                  #
    # Loesen Sie das Anfangswertproblem fuer a1 und a2 #
    #                                                  #
    ####################################################
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

    #plt.show()
    #from sys import exit
    #exit()

    # Unteraufgabe d)

    # zero finding:
    print("Initial value via time-integration")

    #################################################
    #                                               #
    # Berechnen Sie hier die Anfangsauslenkung phiT #
    #                                               #
    #################################################
    T = 1.8
    F = lambda phi0: IntegratePendulum(phi0, T/4, l, g, True)[1][0]
    starttime = time()
    phiT = fsolve(F, a1*pi/2)
    endtime = time()
    print('needed %f seconds' % (endtime-starttime))
    print(phiT*2/pi)

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

    starttime = time()
    phiT = 0.0
    #################################################
    #                                               #
    # Berechnen Sie hier die Anfangsauslenkung phiT #
    #                                               #
    #################################################
    print(phiT)
    endtime = time()
    print('needed %f seconds' % (endtime-starttime))


    # Unteraufgabe f)

    print("Initial value via exakt formula for T")
    print("Elliptic Integral by Arithmetic-Geometric Mean")

    starttime = time()
    phiT = 0.0
    #################################################
    #                                               #
    # Berechnen Sie hier die Anfangsauslenkung phiT #
    #                                               #
    #################################################
    print(phiT)
    endtime = time()
    print('needed %f seconds' % (endtime-starttime))
