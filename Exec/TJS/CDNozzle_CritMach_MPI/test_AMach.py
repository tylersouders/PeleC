import numpy as np

def igas(P, T):

    return P / 188.9 / T

def func(M, gamma, Arat):

    return 1 / (M**2.0) * (2.0 / (gamma + 1.0) * \
        (1 + (gamma - 1.0)/2.0 * M**2.0))**((gamma+1.0)/(gamma-1.0))\
    - (Arat**2.0)

def main():

    # Constants
    pi = 3.14159

    # All cgs units
    P0 = 100 * 68947.6 / 10; #Pa
    T0 = 298;

    # Channel Geo
    rt = 0.18 / 100 #m
    r0 = 0.2286 / 100 #m

    # Assume const gamma for now
    gamma = 1.288

    # Compute Areas
    A0 = pi * (r0**2.0)
    At = pi * (rt**2.0)

    # Compute EOS things (assume air)
    rho0 = igas(P0, T0)
    mdot = At * np.sqrt((gamma * rho0 * P0 * (2.0 / (gamma + 1))**((gamma+1.0)/(gamma-1.0))))

    # Compute Mach number at inlet
    M_Lo = 0.01
    M_Hi = 1.0

    # Solve Mach using bisection because I am lazy
    it = 0
    err = 100
    errlim = 1e-8
    while (err > errlim) and (it < 100):

        # Eval guess
        M = 0.5 * (M_Lo + M_Hi)
        fa = func(M_Lo, gamma, A0/At)
        fb = func(M, gamma, A0/At)
        fc = func(M_Hi, gamma, A0/At)

        if fb == 0 or (M_Hi - M_Lo)/2.0 < errlim:

            break

        elif np.sign(fb) == np.sign(fa):
            M_Lo = M
        else:
            M_Hi = M

        it += 1

    # Use this Mach number to perform TTS transformations
    P = P0 / (1 + (gamma - 1.0)/2.0 * M**2.0)**((gamma)/(gamma-1.0))
    T = T0 / (1 + (gamma - 1.0)/2.0 * M**2.0)

    # Get new density from the static values
    rho = igas(P, T)

    # Use mdot to solve inlet velocity
    Vx = mdot / rho / A0

    # Print Outputs
    print('Computed Total Density = %.4f kg/m^3' % (rho0))
    print('Computed Area Ratio = %.4f' % (A0 / At))
    print('Computed Mass Flow Rate = %.4f kg/s\n' % (mdot))
    print('Computed Inlet Mach Number = %.4f' % (M))
    print('Computed Inlet Pressure Ratio = %.4f' % (P/P0))
    print('Computed Inlet Temperature Ratio = %.4f' % (T/T0))
    print('Computed Inlet Density Ratio = %.4f' % (rho/rho0))
    print('Computed Inlet Velocity = %.4f m/s' % (Vx))





if __name__ == "__main__": main()