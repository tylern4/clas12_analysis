#!/usr/bin/env python
from __future__ import print_function
import uproot
# Loads ROOT for opening files
import ROOT
from ROOT import gROOT, gBenchmark
from ROOT.TMath import Sqrt as sqrt
# Load matplotlib for plotting
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys
import numpy as np

from physics import *
from delta_t import delta_t, vertex_time

# Get input and output file names
if len(sys.argv) == 2:
    input_files = sys.argv[1]
else:
    sys.exit('Error opening files!')

# Start root and load in input_files
gROOT.SetBatch(True)
chain = ROOT.TChain('clas12')
num_files = chain.Add(input_files)

clas12 = uproot.open("/Users/tylern/data/sample_3050.root")['clas12']

# Initialize empty array
mom = np.array([])
W = []
Q2 = []
beta = np.array([])

p_pos = []
beta_pos = []

p_dt = []
deltat_proton = []

# For every event that was loaded in
for pid, b, px, py, pz, charge in zip(
        clas12.array('REC_Particle_pid'), clas12.array('REC_Particle_beta'),
        clas12.array('REC_Particle_px'), clas12.array('REC_Particle_py'),
        clas12.array('REC_Particle_pz'), clas12.array('REC_Particle_charge')):

    np.append(b, beta)

    # Add the Momentum to the p array
    p2 = px**2 + py**2 + pz**2
    p2 = np.abs(p2)
    np.append(mom, np.sqrt(p2))

    print(beta, mom)
    [(beta, mom) for b, p in zip(beta, mom) if b.any() != 0.0]

    for x in range(len(pid)):
        e_mu_p = fvec(px[x], py[x], pz[x], get_mass('ELECTRON'))
        Q2.append(Q2_calc(e_mu, e_mu_p))
        W.append(W_calc(e_mu, e_mu_p))
        # If the particle is positive
        #if charge[x] > 0:
        #p_pos.append(sqrt(p2[x]))
        #beta_pos.append(beta[x])

# Loop over length of REC_Scintillator
#for j in range(len(evnt.REC_Scintillator_pindex)):
# Get Electron vertex from first particle
#    electron_vertex = vertex_time(evnt.REC_Scintillator_time[0],
#                                  evnt.REC_Scintillator_path[0], 1.0)
# Get index for REC_Particle from REC_Scintillator
#    index = evnt.REC_Scintillator_pindex[j]

# Calculate momentum and fill array
#    p2 = evnt.REC_Particle_px[index]**2 + evnt.REC_Particle_py[index]**2 + evnt.REC_Particle_pz[index]**2
#    p2 = abs(p2)
# This array might be different from the p array so we are calculating and refilling it
#    p_dt.append(sqrt(p2))

# Calulating delta_t assuming the mass of a proton and append to the array
#    dt = delta_t(electron_vertex, MASS_P, sqrt(p2),
#                 evnt.REC_Scintillator_time[j],
#                 evnt.REC_Scintillator_path[j])
#    deltat_proton.append(dt)

# Momentum
output_file = "Momentum.pdf"
# Make the figure for plotting
fig = plt.figure(
    num=None, figsize=(16, 9), dpi=200, facecolor='w', edgecolor='k')
# Fill the histogram
plt.hist(mom, 500, normed=1, histtype='stepfilled', alpha=0.75, range=[0, 10])
# Add labels
plt.title("Electron Momentum")
plt.xlabel("P (GeV)")
plt.ylabel("Count")
# Save to output_file name
plt.savefig(output_file)

# W
output_file = "W.pdf"
# Make the figure for plotting
fig = plt.figure(
    num=None, figsize=(16, 9), dpi=200, facecolor='w', edgecolor='k')
# Fill the histogram
plt.hist(W, 500, normed=1, histtype='stepfilled', alpha=0.75, range=[0, 5])
# Add labels
plt.title("W")
plt.xlabel("W (GeV)")
plt.ylabel("Count")
# Save to output_file name
plt.savefig(output_file)

# Q2
output_file = "Q2.pdf"
# Make the figure for plotting
fig = plt.figure(
    num=None, figsize=(16, 9), dpi=200, facecolor='w', edgecolor='k')
# Fill the histogram
plt.hist(Q2, 500, normed=1, histtype='stepfilled', alpha=0.75, range=[0, 10])
# Add labels
plt.title("$Q^2$")
plt.xlabel("$Q^2$ ($GeV^2$)")
plt.ylabel("Count")
# Save to output_file name
plt.savefig(output_file)

# WvsQ2
output_file = "WvsQ2.pdf"
# Make the figure for plotting
fig = plt.figure(
    num=None, figsize=(16, 9), dpi=200, facecolor='w', edgecolor='k')
# Fill the histogram
plt.hist2d(W, Q2, bins=500, normed=True, range=[[0, 5], [0, 10]])
# Add labels
plt.title("W vs $Q^2$")
plt.xlabel("W (GeV)")
plt.ylabel("$Q^2$ ($Gev^2$)")
# Save to output_file name
plt.savefig(output_file)

# MomVsBeta
output_file = "MomVsBeta.pdf"
# Make the figure for plotting
fig = plt.figure(
    num=None, figsize=(16, 9), dpi=200, facecolor='w', edgecolor='k')
# Fill the histogram
plt.hist2d(mom, beta, bins=500, normed=True, range=[[0, 5], [0, 1.2]])
# Add labels
plt.title("W vs $Q^2$")
plt.xlabel("W (GeV)")
plt.ylabel("$Q^2$ ($Gev^2$)")
plt.colorbar()
# Save to output_file name
plt.savefig(output_file)
