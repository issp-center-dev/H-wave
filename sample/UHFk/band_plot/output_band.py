import numpy as np
import tomli
import hwave.qlms
import matplotlib.pyplot as plt

#read toml file
with open("input.toml", "rb") as fp:
  toml_dict = tomli.load(fp)

#UHFk runs 
hwave.qlms.run(input_dict=toml_dict)

#get eigenvalues
data = np.load("output/eigen.dat.npz") 
eigenvalue = data["eigenvalue"]
wavevec_index   = data["wavevector_index"]
wavevec_unit    = data["wavevector_unit"]

#the total number of the eigen states
neigen = int(eigenvalue.shape[1]) 

#calc wavevectors from index & reciprocal lattice vectors
wavevec = np.dot(wavevec_index, wavevec_unit.T)/np.pi

#concatenate wavevector & eigenvalue
Ene_k = np.concatenate([wavevec, eigenvalue], axis=1)

#output wavenumber & eigenvalues to dat file
header_ = "# kx/pi  ky/pi  kz/pi  ene_spin0  ene_spin1"
np.savetxt("band.dat", Ene_k, header = header_)

# plot dispersion by matplotlib
size_list = np.linspace(100,10,neigen)
kmin = -1
kmax =  1
krange = np.linspace(kmin,kmax,100)
ene_exact = -2.0*np.cos(np.pi*krange)
plt.plot(krange, ene_exact, color="black", zorder=1, label="exact")

plt.xlim(kmin, kmax)
[plt.scatter(wavevec[:, 1], eigenvalue[:,i], marker="o", s=size_list[i], label="H-wave: spin{}".format(i)) for i in range(neigen) ]

plt.legend()
plt.xlabel("k/pi")
plt.ylabel("E(k)")
plt.savefig("band.png")
#plt.show()
