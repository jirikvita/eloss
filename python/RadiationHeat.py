#!/usr/bin/python3

# jk 24.6.2021

from math import *

kYear = 365.25*24*3600 # s
kNA = 6.022e23
ke = 1.602e-19 # eV


name='Plutonium-238'
# https://en.wikipedia.org/wiki/Plutonium-238
A = 238.
mA = 1.*A # g!
rho =  19.8 # g/cm^3
# material spehere radius
Thalf = 87.7 * kYear
Ealfa = 5.59e6 # MeV

name='Americium-241'
A = 241.
mA = 1.*A # g!
rho = 13.67 # g/cm^3
# material spehere radius
Thalf = 432.2 * kYear
Ealfa = 5e6 # MeV


R = 2.5 # cm!

tau = Thalf / log(2.)
dconst = 1./tau
V = 4./3. * pi * pow(R,3)
m = rho*V
n = m / mA
N = n*kNA
E = Ealfa * ke # Jouls
# power
P = E * N * dconst # Watts

print('Material: {}'.format(name))
print('Power for a sphere of radius {} cm: {:4.2f} W'.format(R,P))
print('Power/gram: {:1.4} W/g'.format(P/m))
