#!/usr/bin/env python
from __future__ import print_function
# Loads ROOT for opening files
import ROOT
from ROOT import gROOT, gBenchmark
from ROOT.TMath import Sqrt as sqrt
# Load matplotlib for plotting
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys

from physics import *

# Get input and output file names
if len(sys.argv) == 2:
    input_files = sys.argv[1]
    output_file = "Momentum.pdf"
elif len(sys.argv) == 3:
    input_files = sys.argv[1]
    output_file = sys.argv[2]
else:
    input_files = "../*.root"
    output_file = "Momentum.pdf"

# Start root and load in input_files
gROOT.SetBatch(True)
chain = ROOT.TChain('clas12')
num_files = chain.Add(input_files)

# Initialize empty array
p = []
W = []
Q2 = []

# For every event that was loaded in
for evnt in chain:
    if len(evnt.REC_Particle_pid) is not 0:
        # If the events ID is 11==Electron
        if evnt.REC_Particle_pid[0] == 11:
            # Add the Momentum to the p array
            p2 = evnt.REC_Particle_px[0]*evnt.REC_Particle_px[0] + evnt.REC_Particle_py[0]*evnt.REC_Particle_py[0] + evnt.REC_Particle_pz[0]*evnt.REC_Particle_pz[0]
            p.append(sqrt(p2))
            e_mu_p = fvec(evnt.REC_Particle_px[0], evnt.REC_Particle_py[0], evnt.REC_Particle_pz[0], get_mass('ELECTRON'))
            Q2.append(Q2_calc(e_mu, e_mu_p))
            W.append(W_calc(e_mu, e_mu_p))


### Momentum
output_file = "Momentum.pdf"
# Make the figure for plotting
fig = plt.figure(num=None, figsize=(16, 9), dpi=200, facecolor='w', edgecolor='k')
# Fill the histogram
plt.hist(p, 500, normed=1, histtype='stepfilled', alpha=0.75, range=[0, 10])
# Add labels
plt.title("Electron Momentum")
plt.xlabel("P (GeV)")
plt.ylabel("Count")
# Save to output_file name
plt.savefig(output_file)


### W
output_file = "W.pdf"
# Make the figure for plotting
fig = plt.figure(num=None, figsize=(16, 9), dpi=200, facecolor='w', edgecolor='k')
# Fill the histogram
plt.hist(W, 500, normed=1, histtype='stepfilled', alpha=0.75, range=[0, 5])
# Add labels
plt.title("W")
plt.xlabel("W (GeV)")
plt.ylabel("Count")
# Save to output_file name
plt.savefig(output_file)


### Q2
output_file = "Q2.pdf"
# Make the figure for plotting
fig = plt.figure(num=None, figsize=(16, 9), dpi=200, facecolor='w', edgecolor='k')
# Fill the histogram
plt.hist(Q2, 500, normed=1, histtype='stepfilled', alpha=0.75, range=[0, 10])
# Add labels
plt.title("$Q^2$")
plt.xlabel("$Q^2$ ($GeV^2$)")
plt.ylabel("Count")
# Save to output_file name
plt.savefig(output_file)


### WvsQ2
output_file = "WvsQ2.pdf"
# Make the figure for plotting
fig = plt.figure(num=None, figsize=(16, 9), dpi=200, facecolor='w', edgecolor='k')
# Fill the histogram
plt.hist2d(W,Q2,bins=500, normed=True, range=[[0, 5],[0, 10]])
# Add labels
plt.title("W vs $Q^2$")
plt.xlabel("W (GeV)")
plt.ylabel("$Q^2$ ($Gev^2$)")
# Save to output_file name
plt.savefig(output_file)
