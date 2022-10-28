import numpy as np
import tomli

#get eigenvalues
data = np.load("output/eigen.dat.npz") 
eigenvalue = data["eigenvalue"]

nk = int(eigenvalue.shape[0]) 
neigen = int(eigenvalue.shape[1]) 

#output idx of wavenumber & eigenvalues to dat file
with open("band.dat", "w") as fp:
  for eidx in range(neigen):
    for kidx in range(nk):
      fp.write("{}  {}\n".format(kidx,eigenvalue[kidx,eidx]))
    fp.write("\n\n")
