#!/usr/bin/python

# jk 26.4.2018

from Bethe import *
from Brems import *

print '*** Supported particles:'
PrintParticles()
print '*** Supported materials:'
PrintMaterials()


print '***************************************************'

pairs = [ ['Alpha', 'Si'],
          ['Deuteron', 'Si'],
          ['Proton', 'Si'],
          ['Kaon', 'Si'],
          ['Pion', 'Si'],
          ['Muon', 'Si'],
          ['Electron', 'Si'],
          ['Positron', 'Si'],
]

KeepSameMomentum = False

for pair in pairs:
    pname = pair[0]
    mname = pair[1]
    particle = gParticles[pname]
    material = gMaterials[mname]

    beta = 0.
    gamma = 0.
    bg = 0.
    p = 0.
    M = particle.GetM()
    if KeepSameMomentum:
        p = 20. # MeV
        # beta*gamma:
        bg = p/M
        beta = GetBetaFromBg(bg)
        gamma = GetGammaFromBg(bg)
    else:
        beta = 0.5
        gamma = GetGamma(beta)
        bg = beta*gamma
        p = bg*M
        
    print '*** Particle: {:} in {:} '.format(particle.GetName(), material.GetName())
    # print bg, beta, gamma
    print '    p={:3.1f}, beta={:1.4f}, gamma = {:3.3f}'.format(p, beta,gamma,)
    print '    Ionization losses : {:1.3f} MeV/cm'.format(dEdX(beta, particle, material), )
    E = gamma*M
    print '    Radiation losses  : {:1.3f} MeV/cm'.format( dEdXBrems(E, particle, material), )
