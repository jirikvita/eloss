#!/usr/bin/python

#############################################################
#
# jiri kvita, 24th Feb 2015, 1.3.2018
# mass units in MeV
#
#############################################################

from math import *
from ROOT import *
from myAll import *

from Bethe import *
betaMin = 0.05

#############################################################

particle = gParticles['Proton']
M = particle.GetM()
T = 20   # MeV
beta = sqrt(T*(T + 2*M)) / M

gamma = GetGamma(beta)
E=M*gamma
p=beta*gamma*M

print('Initial particle [MeV] M=%3.1f, p=%3.1f, E=%3.1f, beta=%3.3f, gamma=%3.3f, beta*gamma=%3.3f' % (M,p,E,beta,gamma,beta*gamma,))

# silicon:
material = gMaterials['Si']

# test:
lossStep = dEdX(beta, particle, material)
print("E=%f, dE/dX=%f" % (E,lossStep,))

# prepare for the loop:
x=0
ip=0
# tracing in time:
dt = 0.5*1e-12 # 1/2 ps

g_dEdX = TGraph()#Errors()
g_dEdX.SetName("g_dEdX")
g_dEdX.SetTitle("dEdX;x [cm];|dE/dx| [MeV/cm]")
debug = 0

print("*** Loop")
count = 0
while beta > epsilon:
    count = count + 1
    if count % 10000 == 0:
        print('Processing step %i' % (count,))
    lossStep = dEdX(beta, particle, material)
    
    xStep = beta*c*dt
    x = x + xStep
    if debug or count == 1:
        print('E=%f, beta=%3.3f, xStep=%f, lossStep=%f MeV/mm, xStep*lossStep=%f MeV' % (E,beta,xStep,lossStep*1.e-1,xStep*lossStep,))
    E = E - xStep*lossStep
    p2 = E*E - M*M
    if p2 > 0: 
        p = sqrt(E*E - M*M)
    else:
        p=0
    beta = p/E
    # finer step towards the end?
    #if beta < 0.35:
    #    dt = dt / 1.1
    if debug > 1:
        print('  new beta=%3.3f, x=%f' % (beta,x,))
    g_dEdX.SetPoint(ip, x, lossStep)
    ip = ip+1

print("*** End of loop")

can = nextCan.nextTCanvas("Bragg", "Bragg", 0, 0, 1000, 800)
g_dEdX.SetMarkerStyle(20)
g_dEdX.SetMarkerSize(0.4)
g_dEdX.SetMarkerColor(2)
g_dEdX.Draw("APC")
#ROOT.gPad.SetLogy(1)
can.Print(can.GetName() + '.png')
#can.Print(can.GetName() + '.eps')

gApplication.Run()
