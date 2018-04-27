#!/usr/bin/python

# jk 26.4.2018

from Bethe import *
from Brems import *

#print '*** Supported particles:'
#PrintParticles()
#print '*** Supported materials:'
#PrintMaterials()

particle = gParticles['Proton']
material = gMaterials['Si']

M = particle.GetM()
beta = 0.5
gamma = GetGamma(beta)
bg = beta*gamma
p = bg*M
E = gamma*M

print '*** Particle: {:} in {:} '.format(particle.GetName(), material.GetName())
# print bg, beta, gamma
print '    p={:3.1f}, beta={:1.4f}, gamma = {:3.3f}'.format(p, beta,gamma,)
print '    Ionization losses : {:1.3f} MeV/cm'.format(dEdX(beta, particle, material), )
print '    Radiation losses  : {:1.3f} MeV/cm'.format( dEdXBrems(E, particle, material), )
