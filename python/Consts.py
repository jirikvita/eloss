#!/usr/bin/python

homega0 = 28.816e-6 # eV

K = 0.307075 # MeV/mol/cm^2 for A in g/mol
betaMin = 0.05
gkjterm = 0.2000

me = 0.51100 # MeV
c = 2.99792458e10 # cm/s
h = 6.626e-34 # J*s
hc = 197.326972e-13 # MeV*cm
alpha = 1./137.036
echarge = 1.60217662e-19 # C

# re = e^2 / (4pi epsilon0 me*c*c)
re = alpha*hc / me
#re = 2.817940e-15*1e2 # cm 

NA = 6.02214e23

# for Brehms:
Konst = 4*alpha*NA


epsilon = 1e-6

#print re
#print alpha*hc / me
