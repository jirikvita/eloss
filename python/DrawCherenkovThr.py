#!/snap/bin/pyroot

# jk 24.6.2023

import ROOT
from math import *
from collections import OrderedDict

Particles = OrderedDict({ #'e': 0.511,
                          'mu' : 105.7,
                          'pi' : 139.6,
                          #'K' : 493.7,
                          #'p' : 938,
                         }
                        ) # MeV

cols = OrderedDict({ 'e': ROOT.kRed,
                     'mu' : ROOT.kBlue,
                     'pi' : ROOT.kGreen+1,
                     'K' : ROOT.kMagenta,
                     'p' : ROOT.kBlack,
                    }
                   )

plabels = OrderedDict({ 'e': 'e',
                        'mu' : '#mu',
                        'pi' : '#pi',
                        'K' : 'K',                
                        'p' : 'p',
                       }
                      )


ns = [#1.006, 1.01,
    #1.01,
    1.06, 1.10, 1.15]

cans = []

# momentum:
Momenta = [220, 240, 260, 280, 300, 320, 340, 360, 
           #400, 500, 600, 700, 800, 900, 1000,
           ] # MeV
p0 = Momenta[0]
p1 = Momenta[-1]

line = ''
for momentum in Momenta:
    print(f'p={momentum}')
    betas = OrderedDict()
    for particle in Particles:
        mass = Particles[particle]
        #print(f'  particle: {particle}, mass={mass}')
        betagamma = momentum / mass
        beta = betagamma / sqrt(1+pow(betagamma,2))
        betas[particle] = beta
    print(betas)


# beta(p|m)
funs_betap = OrderedDict()

canname = 'CherenkovBetaPlot'
cw, ch = 800, 600
can1 = ROOT.TCanvas(canname, canname, 0, 0, cw, ch)
cans.append(can1)
for particle in Particles:
    mass = Particles[particle]
    name = f'fun_beta_{particle}'
    # here x/[0] = p/m = beta*gamma
    form = 'x/[0] / sqrt( 1 + (x/[0])^2 )' # beta
    fun = ROOT.TF1(name, form, p0, p1)
    fun.SetNpx(1000)
    fun.SetParameters(mass)
    funs_betap[particle] = fun
leg = ROOT.TLegend(0.75, 0.15, 0.90, 0.35)
leg.SetBorderSize(0)

opt = 'same'
nbx, nby = 100, 100
h2 = ROOT.TH2D('tmp', 'tmp;p [MeV/c];#beta', nbx, p0, p1, nby, 0., 1.1)
h2.SetStats(0)
h3 = h2.DrawCopy()
ROOT.gStyle.SetPalette(ROOT.kDeepSea)
ROOT.gStyle.SetOptTitle(0)
for particle,fun in funs_betap.items():
    fun.SetLineColor(cols[particle])
    fun.Draw('same')
    leg.AddEntry(fun, f'{plabels[particle]}', 'L')
leg.Draw()
ROOT.gPad.Update()


Funs_yieldsfp = OrderedDict()
canname = 'CherenkovIntensityPlot'
cw, ch = 1000, 800
can2 = ROOT.TCanvas(canname, canname, 200, 200, cw, ch)
can2.Divide(2,2)
jn = 0
labels = OrderedDict()
legs = OrderedDict()
for n in ns:
    labels[n] = ROOT.TLatex(0.76, 0.84, f'n={n}')
    labels[n].SetNDC()
    legs[n] = ROOT.TLegend(0.15, 0.65, 0.35, 0.85)
    legs[n].SetBorderSize(0)
    print(f'Processing refraction index n={n}')
    # Cherenkov yield proportionality factor
    funs_yieldsfp = OrderedDict()
    for particle in Particles:
        mass = Particles[particle]
        name = f'fun_yieldf_{particle}_n{n}'
        form = '(1/sqrt(1 + ([0]/x)^2 ) > 1/[1])*( 1 - 1/[1]^2*(1 + ([0]/x)^2 ) )'
        fun = ROOT.TF1(name, form, p0, p1)
        fun.SetNpx(1000)
        fun.SetParameters(mass, n)
        funs_yieldsfp[particle] = fun
    Funs_yieldsfp[n] = funs_yieldsfp
    jn = jn + 1
    can2.cd(jn)
    h4 = h2.DrawCopy()
    h4.GetYaxis().SetTitle('~ I ~ 1 - 1/(#betan)^{2}')
    h4.GetYaxis().SetRangeUser(0, 0.2)
    for particle,fun in funs_yieldsfp.items():
        fun.SetLineColor(cols[particle])
        fun.Draw('same')
        legs[n].AddEntry(fun, f'{plabels[particle]}', 'L')
    labels[n].Draw()
    legs[n].Draw()
    ROOT.gPad.SetGridx()
    ROOT.gPad.SetGridy()
can2.Update()
cans.append(can2)

for can in cans:
    can.Print(can.GetName() + '.png')
    can.Print(can.GetName() + '.pdf')

ROOT.gApplication.Run()
