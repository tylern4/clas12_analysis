#!/usr/bin/env python
from __future__ import print_function
from asyncio import events
import uproot
# Load matplotlib for plotting
import matplotlib.pyplot as plt
import sys
import numpy as np

from reaction import *
from physics import *
from delta_t import delta_t, vertex_time
from nicks_plot_utils import Hist1D, Hist2D
import time

clas12 = uproot.concatenate(sys.argv[1],
                            file_handler=uproot.MultithreadedFileSource,
                            num_workers=4)
num_events = len(clas12)
print(num_events, num_events/2000)
# Initialize empty array
mom = np.ones(num_events)*np.nan
W = np.ones(num_events)*np.nan
Q2 = np.ones(num_events)*np.nan
beta = np.ones(num_events)*np.nan

p_pos = []
beta_pos = []

p_dt = []
deltat_proton = []

start = time.time()
# For every event that was loaded in
for i, event in enumerate(clas12):
    if(len(event['pid']) <= 0):
        continue
    for pidi, pxi, pyi, pzi, bi, qi in zip(event['pid'],
                                           event['px'],
                                           event['py'],
                                           event['pz'],
                                           event['beta'],
                                           event['charge']):

        # Add the Momentum to the p array
        p2 = pxi**2 + pyi**2 + pzi**2
        p2 = np.abs(p2)
        mom[i] = np.sqrt(p2)

        if(pidi == 11):
            e_mu_p = fvec(pxi, pyi, pzi, get_mass('ELECTRON'))
            Q2[i] = Q2_calc(e_mu, e_mu_p)
            W[i] = W_calc(e_mu, e_mu_p)


end = time.time()

print(f"{end-start} Seconds, {num_events/(end-start)} Hz")

# Momentum
output_file = "Momentum.pdf"
# Make the figure for plotting
fig = plt.figure(
    num=None, figsize=(16, 9), dpi=200, facecolor='w', edgecolor='k')

# Fill the histogram
mom = mom[~np.isnan(mom)]
h1 = Hist1D(mom, 500)
_ = h1.histogram()
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
bad = np.isnan(W) | np.isnan(Q2)
W = W[~bad]
Q2 = Q2[~bad]
h2 = Hist1D(W, bins=500, xrange=[0, 10])
_ = h2.histogram()
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
h3 = Hist1D(Q2, bins=500, xrange=[0, 20])
_ = h3.histogram()

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
w_vs_q2 = Hist2D(W, Q2, xbins=500, ybins=500, xrange=[0, 10], yrange=[0, 20])
_ = w_vs_q2.plot()
# Add labels
plt.title("W vs $Q^2$")
plt.xlabel("W (GeV)")
plt.ylabel("$Q^2$ ($Gev^2$)")
# Save to output_file name
plt.savefig(output_file)

# WvsQ2
output_file = "WvsQ2_3D.pdf"
# Make the figure for plotting
fig = plt.figure(
    num=None, figsize=(16, 9), dpi=200, facecolor='w', edgecolor='k')
# Fill the histogram
w_vs_q2 = Hist2D(W, Q2, xbins=500, ybins=500, xrange=[0, 10], yrange=[0, 20])
_ = w_vs_q2.plot(zeros=False)
# Add labels
plt.title("W vs $Q^2$")
plt.xlabel("W (GeV)")
plt.ylabel("$Q^2$ ($Gev^2$)")
# Save to output_file name
plt.savefig(output_file)
