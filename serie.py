#!/bin/python3
# Lab SS - Fourier Series
# 01-2022
# Authors: João G. & Daniel D.
from cmath import phase, sin, pi, e;

def fourier_serie(): # {{
    """
    This function calculates the Fourier Serie of a specific 
    triangle wave shifted downwards, up to the `nth` order.
    --^
    """
    # Tri@ngle wave properties
    X0 = int(input("X0: "));   # min value and
    X1 = int(input("X1: "));   # max value on the vertical axis
    t1 = float(input("t1: ")); # stationary point
    T = float(input("Period(T): "));
    w0 = (2*pi)/T;
    n = int(input("Iterations(n): "));
    # ~~~@~~~
    A0 = 1/2 * (X1-X0) + X0; # X0 < 0 (vertical shift down of the wave)
    final_expr = str(A0);
    print('\nA0 = ' + str(A0) + '\n');
    for k in range(1, n+1):
        print('k = ' + str(k));
        exp = ( sin(k*w0*(T-t1)/2) * ( (X0-X1)/(T-t1) ) * (e**(1.j*k*w0*( (T-t1)/2 ))) ); # coeficients of the 1st square wave
        exp += ( sin(k*w0*(t1/2)) * ( (X1-X0)/t1 ) * (e**(-1.j*k*w0*( t1/2 ))) );         # coeficients of the 2nd square wave
        exp /= (1.j * w0 * (k**2) *pi);                                                   # integration property of the FS
        Ak = 2*abs(exp);
        thetak = phase(exp);
        print('A' + str(k) + ' = ' + str(Ak));
        print('𝛉' + str(k) + ' = ' + str(thetak) + '\n');
        final_expr += ' + ' + str(Ak) + '*cos(' + str(k) + '*' + str(w0) + '*t + ' + str(thetak) + ')';

    return final_expr;
# }}

f = fourier_serie();
print('Fn(t): ' + f); # plot this expr :)
