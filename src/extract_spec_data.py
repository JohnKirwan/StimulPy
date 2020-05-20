# -*- coding: utf-8 -*-
"""
Created on Mon Nov 25 22:25:27 2019

@author: John D. Kirwan
"""

import pandas as pd
I = pd.read_csv('spec_data.txt',sep='\t')
I.describe()

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

# Plot panels of luminance curves for the measurements

# Calculate photon irradiance in relevant range for urchins
# Remove wvs that would be quickly attenuated at depth.

# a line plot comparing wv to iradiance divided by matrial type
I.groupby('material').plot(kind='line',x='wv',y='irrad',xlim=(375,725))
#plt.subplot(1, 2, 1)
#plt.tight_layout()
plt.show()

Ihi = I.groupby(['material','dark'])
Ihi.describe()

# Typical lambda max for clear marine visual system is c. 500 nm. Opsin half
# width is c. 100 nm. -> half senstivity at 450 or 550 nm.

# If, for instance, used same (or even larger) photo flux of 500 nm light:
# Not clear if sensitivity high enough at all, even in shallow water
# In deeper water, long wv light rapidly attentuated, meaning sensitivity ought
# be extremmely high.
# Intensity of light reflected from walls and transmitted through water would
# change moving through the arena. For homogenous wv, will affect senstivity
# but not contrast.

# If used two wavelengths, e.g. 500 and 600 (need to be high contrast between
# max and min reflectance in each case), slight differential in attenuation,
# resulting in slight different min and max reflectance while moving. Could
# model by estimating change in intensity over arena radius. (Would keep the
# min totally black at distance, less so when approaching. Could move the
# max from below to above SNR threshold.)

# Probably not that big a problem if short/long wavelengths applied which are
# absorped in water. Just makes calculating the irradiance more involved and
# probably doesn't help much. However, inks should have similar and proportionate
# reflectance for whole range of wavelenghts used.

# Try a version with all wavelengths normalised to the highest value of white.
# Check if proportions at each increment is similar throughout.
import numpy as np
Ihi = I.groupby(['material','dark']).apply['irrad'](lambda x: x - np.mean(x))

# Compare two weight of material. Prefer consistent colours and the blacker min.

# Bring LED clusters down to two colours nearest torquise. Compare to sets of
# LEDs Ina has.

# Simple test - compare intensity of aquarium lights to LEds on two heat sinks
# in arena at present.

Ihi.plot(kind='line',x='wv',y='irrad',xlim=(375,725))
