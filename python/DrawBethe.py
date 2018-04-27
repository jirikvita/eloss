#!/usr/bin/python

#############################################################
#
# jiri kvita, 24th Feb 2015, 1.3.2018, 21.2.2018
# mass units in MeV
#
#############################################################

from math import *
from ROOT import *
from myAll import *

from Bethe import *
betaMin = 0.1

#############################################################

def SetStyle(gr, mst, msz, mc):
    gr.SetMarkerStyle(mst)
    gr.SetMarkerSize(msz)
    gr.SetMarkerColor(mc)
    gr.SetLineColor(mc)




#############################################################
#############################################################
#############################################################
    
z = 1

particle = gParticles['Muon']
# particle = 'protons'
# particle = gParticles['Alpha']
# particle = gParticles['Deuteron']

# silicon:
material = gMaterials['Si']
rho = material.GetRho()

# beta
g_dEdX_beta = TGraph()#Errors()
g_dEdX_beta.SetName("g_dEdX")
g_dEdX_beta.SetTitle("dEdX; #beta;[MeV/cm]")

g_Delta_beta = TGraph()#Errors()
g_Delta_beta.SetName("g_Delta")
g_Delta_beta.SetTitle("Delta; #beta;[MeV/cm]")

g_ratio_beta = TGraph()#Errors()
g_ratio_beta.SetName("g_ratio")
g_ratio_beta.SetTitle("Ratio; #beta;dE/#Delta_{p}")


# beta*gamma
g_dEdX_betagamma = TGraph()#Errors()
g_dEdX_betagamma.SetName("g_dEdX")
g_dEdX_betagamma.SetTitle("dEdX; #beta#gamma;[MeV/cm]")

g_Delta_betagamma = TGraph()#Errors()
g_Delta_betagamma.SetName("g_Delta")
g_Delta_betagamma.SetTitle("Delta; #beta#gamma;[MeV/cm]")

g_ratio_betagamma = TGraph()#Errors()
g_ratio_betagamma.SetName("g_ratio")
g_ratio_betagamma.SetTitle("Ratio; #beta#gamma;dE/#Delta_{p} ")


# T
g_dEdX_T = TGraph()#Errors()
g_dEdX_T.SetName("g_dEdX")
g_dEdX_T.SetTitle("dEdX; T [GeV];[MeV/cm]")

g_Delta_T = TGraph()#Errors()
g_Delta_T.SetName("g_Delta")
g_Delta_T.SetTitle("Delta; T [GeV];[MeV/cm]")

g_ratio_T = TGraph()#Errors()
g_ratio_T.SetName("g_ratio")
g_ratio_T.SetTitle("Ratio; T [GeV];dE/#Delta_{p} ")


debug = 0

ip = 0
# thickness of the detector layer:
x = 1.4 # cm, length of particle camera chip
# x = 80e-6*100 # 80 mum -> cm
# x = 1.e-6*100 # 1mum

nsteps = 1000

bmin = betaMin
bmax = 1.
step = (bmax-bmin)/nsteps
betas = [ bmin + n*step for n in range(0, nsteps) ]

gmin = GetGamma(bmin)
gmax = 1000.
step = (gmax-gmin)/nsteps
gammas = [ gmin + n*step for n in range(0, nsteps) ] 


print betas[:20]
print gammas[:20]

for i in range(0,len(betas)):
    beta = betas[i]

    if beta >= 1: continue
    
    lossStep = dEdX(beta, particle, material) #* 1./rho
    Delta = GetMostProbableLoss(x, beta, particle, material) #*  1./rho
    ratio = x*lossStep / Delta
    g_dEdX_beta.SetPoint(ip, beta, lossStep)
    g_Delta_beta.SetPoint(ip, beta, Delta/x)
    g_ratio_beta.SetPoint(ip, beta, ratio)

    gamma = gammas[i]
    bg = gamma*GetBeta(gamma)
    T = (gamma-1)*particle.GetM() / 1000. # GeV
    
    lossStep = dEdX(GetBeta(gamma), particle, material) #* 1./rho
    Delta = GetMostProbableLoss(x, GetBeta(gamma), particle, material) #*  1./rho
    g_dEdX_betagamma.SetPoint(ip, bg, lossStep)
    g_Delta_betagamma.SetPoint(ip, bg, Delta/x)
    g_ratio_betagamma.SetPoint(ip, bg, ratio)

    g_dEdX_T.SetPoint(ip, T, lossStep)
    g_Delta_T.SetPoint(ip, T, Delta/x)
    g_ratio_T.SetPoint(ip, T, ratio)
    
    ip = ip+1

cw = 800
ch = 600

ctags = ['beta', 'betagamma', 'T']
grs = [ [g_dEdX_beta, g_Delta_beta, g_ratio_beta],
        [g_dEdX_betagamma, g_Delta_betagamma, g_ratio_betagamma],
        [g_dEdX_T, g_Delta_T, g_ratio_T],
        ]
        
cans = []
legs = []

for ctag, gr in zip(ctags,grs):
    can = nextCan.nextTCanvas("Bethe_" + ctag + '_' + particle.GetName() + '_in_' + material.GetName(), "Bethe_" + ctag + '_'  + particle.GetName() + '_in_' + material.GetName(), 0, 0, cw, ch)

    if ctag != 'beta':
        ROOT.gPad.SetLogx(1)
    ROOT.gPad.SetLogy(1)
    ROOT.gPad.SetGridx(1)
    ROOT.gPad.SetGridy(1)
    g_dEdX = gr[0]
    g_Delta = gr[1]
    g_ratio = gr[2]

    SetStyle(g_dEdX, 20, 0.4, kRed)
    SetStyle(g_Delta, 20, 0.4, kBlue)
    
    #g_Delta.Draw("APC")
    #g_dEdX.Draw("PC")
    #can.Update()
    #can.Print(can.GetName() + '_dEdX_Delta.png')

    g_dEdX.Draw("APC")
    g_dEdX.GetXaxis().SetMoreLogLabels(1)
    g_Delta.Draw("PC")
    can.Update()

    leg = TLegend(0.70, 0.70, 0.85, 0.88)
    leg.SetHeader(particle.GetName() + ' in ' + material.GetName())
    leg.AddEntry(g_dEdX, '|dE/dx|', 'PL')
    leg.AddEntry(g_Delta, '#Delta_{p}/x', 'PL')
    leg.SetBorderSize(0)
    leg.Draw()
    legs.append(leg)

    can.Update()
    can.Print(can.GetName() + '_dEdX.png')
    
    canratio = nextCan.nextTCanvas("BetheRatio_" + ctag + '_'  + particle.GetName() + '_in_' + material.GetName(), "BetheRatio_" + ctag + '_'  + particle.GetName() + '_in_' + material.GetName(), 100, 100, cw, ch)
    canratio.cd()
    #ROOT.gPad.SetLogy(1)
    ROOT.gPad.SetLogx(1)
    SetStyle(g_ratio, 20, 0.4, kRed)
    g_ratio.Draw("APC")
    g_ratio.GetXaxis().SetMoreLogLabels()
    can.Print(can.GetName() + '.png')
    canratio.Print(canratio.GetName() + '.png')
    
    cans.append(can)
    cans.append(canratio)
    
gApplication.Run()
