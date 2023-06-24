#!/usr/bin/python

import math

Particles = { 'e': 0.511, '\\mu' : 105.7, '\\pi' : 139.6, 'p' : 938, 'K' : 493.7} # MeV

# momentum:
Momenta = [0.5, 1., 2., 10., 20., 180] # GeV



print('Beta')

line = ''
for momentum in Momenta:
    line= line + ' %4.1f & ' % (momentum,)
line = line[0:-3] + ' \\\\' 
print(line)


for particle in Particles:
    mass = Particles[particle]
    #print(particle,mass)
    line = '$%s$ $m=%4.1f$ & ' % (particle,mass,)
    for momentum in Momenta:
        alpha = momentum*1.e3 / mass
        beta = alpha / math.sqrt(1+alpha*alpha)
        line = line + ' %1.5f &' % (beta)
    line = line[0:-3] + ' \\\\' 
    print(line)
        
    

#########
line = ''
for momentum in Momenta:
    line= line + ' %4.1f & ' % (momentum,)
line = line[0:-3] + ' \\\\' 
print(line)


print('N')
for particle in Particles:
    mass = Particles[particle]

 

    line = '$%s$ $m=%4.1f$ & ' % (particle,mass,)
    for momentum in Momenta:
        alpha = momentum*1.e3 / mass
        beta = alpha / math.sqrt(1+alpha*alpha)
        line = line + ' %1.5f &' % (1/beta)
    line = line[0:-3] + ' \\\\' 
    print(line)
        
    
