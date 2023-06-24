#!/snap/bin/pyroot

# jk 24.6.2023

import ROOT

from math import *

from collections import OrderedDict

Particles = OrderedDict({ #'e': 0.511,
                          'mu' : 105.7, 'pi' : 139.6,
                          #'p' : 938, 'K' : 493.7
                         }
                        ) # MeV

cols = OrderedDict({ 'e': ROOT.kRed,
                     'mu' : ROOT.kBlue,
                     'pi' : ROOT.kGreen+1,
                     'p' : ROOT.kBlack,
                     'K' : ROOT.kMagenta
                    }
                   )

plabels = OrderedDict({ 'e': 'e',
                        'mu' : '#mu',
                        'pi' : '#pi',
                        'p' : 'p',
                        'K' : 'K'
                       }
                      )


# momentum:
Momenta = [220, 240, 260, 280, 300, 320, 340, 360, 
           #400, 500, 600, 700, 800, 900, 1000,
           ] # MeV

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


p0 = Momenta[0]
p1 = Momenta[-1]

# beta(p|m)
funs_betap = OrderedDict()


for particle in Particles:
    mass = Particles[particle]
    name = f'fun_beta_{particle}'
    # here x/[0] = p/m = beta*gamma
    form = 'x/[0] / sqrt( 1 + (x/[0])^2 )' # beta
    fun = ROOT.TF1(name, form, p0, p1)
    fun.SetNpx(1000)
    fun.SetParameters(mass)
    funs_betap[particle] = fun
    
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
ROOT.gPad.Update()


ns = [#1.006, 1.01,
    1.06, 1.10, 1.15]
Funs_yieldsfp = OrderedDict()
can = ROOT.TCanvas()
can.Divide(2,2)
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
    can.cd(jn)
    h4 = h2.DrawCopy()
    h4.GetYaxis().SetTitle('I ~ 1 - 1/(#betan)^{2}')
    h4.GetYaxis().SetRangeUser(0, 0.2)
    for particle,fun in funs_yieldsfp.items():
        fun.SetLineColor(cols[particle])
        fun.Draw('same')
        legs[n].AddEntry(fun, f'{plabels[particle]}', 'L')
    labels[n].Draw()
    legs[n].Draw()
can.Update()


ROOT.gApplication.Run()
