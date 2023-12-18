import h5py
from matplotlib import pyplot as plt
import numpy as np

# The HDF5 file containing the data.
filename = 'dca_sp.hdf5'

# Open the file.
data = h5py.File(filename,'r')
wn = data['domains/frequency-domain/elements'] [:]

# Cluster self-energy
# shape: [wn, K, spin2, band2, spin1, band1, real/imag]
self_energy = data['functions/Self_Energy/data']

# Plot the cluster self-energy vs. Matsubara frequency for the first cluster momentum ([0.,0.]).
# All band and spin indices are set to zero.
K_ind = 0
s1 = s2 = 0
b1 = 1
b2 = 1

Nwn = len(wn)//2



select =  self_energy[Nwn:, K_ind, s2, b2, s1, b1, 0:2]
#select_re =  self_energy[Nwn:, K_ind, s2, b2, s1, b1, 0]

np.savetxt('sigma%i%i.txt'%(b1,b2),select,fmt='%2.4f')
