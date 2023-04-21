import os, sys
import numpy as np
import itertools
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D



file_chi0q = 'output/chi0q.npz'
file_chiq  = 'output/chiq.npz'

chi0q = np.load(file_chi0q)['chi0q']
chiq  = np.load(file_chiq)['chiq']

nx,ny = 32,32
nmat = 1024
T = 1.0
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

X, Y = np.meshgrid(kx, ky)


dict_chi = {"chi0":chi0.real, "chic": chic.real, "chis":chis.real}

for name in ["chi0", "chic", "chis"]:
    fig = plt.figure(figsize = (8, 8))
    ax = fig.add_subplot(111, projection="3d")
    ax.set_xlabel("kx", size = 16)
    ax.set_ylabel("ky", size = 16)

    ax.set_zlabel(name, size = 16)

    # 曲面を描画
    ax.plot_surface(X, Y, dict_chi[name], cmap = "summer")

    plt.savefig(name+".png".format(T), format="png", dpi=300)
    plt.close()
    
    
