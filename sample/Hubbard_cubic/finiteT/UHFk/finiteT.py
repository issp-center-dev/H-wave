import tomli, tomli_w
import subprocess
import numpy as np
import argparse
import os


parser = argparse.ArgumentParser(
        prog='finiteT_uhfk.py',
        description='T-dependence of magnetic properties for Hubbard model on cubic lattice by using UHFk',
        epilog='end',
        add_help=True,
)

parser.add_argument('-u', action='store', dest='u',
                        type=float, choices=None,
                        required=True,
                        help=('Hubbard U'),
                        metavar=None)

parser.add_argument('-g', '--grid', action='store', dest='ngrid',
                        type=int, choices=None,
                        required=True,
                        help=('# of grid'),
                        metavar=None)

parser.add_argument('--min', action='store', dest='tmin',
                        default=0.0, 
                        type=float, choices=None,
                        help=('min of temperature/default:0.0'),
                        metavar=None)

parser.add_argument('--max', action='store', dest='tmax',
                        default=5.0, 
                        type=float, choices=None,
                        help=('max of temperature/defalut:5.0'),
                        metavar=None)


args = parser.parse_args()
U = args.u
Ngrid = args.ngrid
Tmin = args.tmin
Tmax = args.tmax

#read toml file
with open("input.toml.org", "rb") as fp:
  toml_dict = tomli.load(fp)

info_mode = toml_dict.get("mode", {})
info_file = toml_dict.get("file", {})
g_file = info_file["output"]["green"]+".npz" 

#make coulombintra.dat 
with open("coulombintra.dat.org", "r") as fp:
  lines = fp.readlines() 
  lines.append("0    0    0    1    1   {:2f}   0.000".format(U))
with open("coulombintra.dat", "w") as fp:
  fp.writelines(lines) 

#set temperature range
T_list = np.linspace(Tmin, Tmax, Ngrid)
Tround_list = np.round(T_list, 3)

#initialize outputs
mag_list = np.empty(2*Ngrid)
Tneel=-1.0

#finite-T UHFk calc
for Tidx, T in enumerate(Tround_list):
  output_T = "output_T{}".format(T)
  info_mode["param"]["T"] = T
  info_file["output"]["path_to_output"] = output_T
  if Tidx > 0:
    info_file["input"]["initial"] = "green_init.dat.npz"

  ifile = "input_T{}.toml".format(T)
  with open(ifile, "wb") as fp:
    tomli_w.dump(toml_dict, fp)

  #UHFk runs 
  uhfk_cmd = ["python3", "qlms.py", ifile]
  subprocess.run(uhfk_cmd)

  #calc magnetic moment 
  path2g_file = os.path.join(output_T, g_file)
  data = np.load(path2g_file)["green"] 
  mag  = data[0][0][0][0][0]
  mag -= data[0][1][0][1][0]
  mag *= 0.5
  
  mag_list[2*Tidx]   = T
  mag_list[2*Tidx+1] = abs(mag)
  
  #mv current input.toml to output_T directory 
  subprocess.run(["mv", ifile, output_T])
  #green.dat.npz is used as "initial" for next step 
  subprocess.run(["cp", path2g_file, "./green_init.dat.npz"])

  #Neel temperature 
  if abs(mag) < 1e-4 and Tneel==-1.0:
    Tneel = T

#output results
with open("mag_U{}.dat".format(U), "w") as fp:
  np.savetxt(fp,mag_list.reshape((Ngrid,2)))

with open("Tneel_U{}.dat".format(U), "w") as fp:
  fp.write("{}".format(Tneel))
