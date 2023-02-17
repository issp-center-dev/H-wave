import os, sys
import numpy as np
import itertools

file_chi0q = 'output/chi0q.npz'
file_chiq  = 'output/chiq.npz'

chi0q = np.load(file_chi0q)['chi0q']
chiq  = np.load(file_chiq)['chiq']

nx,ny = 32,32
nmat = 1024
norb = 1
ns = 2
nd = ns * norb

kx = np.arange(nx) * np.pi * 2 / nx
ky = np.arange(ny) * np.pi * 2 / ny

# chiq[l][k][a,ap,b,bp]
chi00 = (chiq[nmat//2, :, 0, 0, 0, 0]).reshape(nx,ny)
chi01 = (chiq[nmat//2, :, 0, 0, 1, 1]).reshape(nx,ny)

chis = chi00 - chi01
chic = chi00 + chi01

chi0 = (chi0q[nmat//2, :, 0, 0, 0, 0]).reshape(nx,ny)

for i,j in itertools.product(range(nx),range(ny)):
    print(kx[i], ky[j], chi0[i,j].real, chic[i,j].real, chis[i,j].real)

