import numpy as np
import tomli
import hwave.qlms

#read toml file
with open("input.toml", "rb") as fp:
  toml_dict = tomli.load(fp)

#UHFk runs 
hwave.qlms.run(input_dict=toml_dict)

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
