# Sample for calculating temperature dependence of magnetic properties

## What's this sample?

This is the sample for calculating the temperature dependence of the magnetic moment $m_z$ for the Hubbard model on cubic lattice at half filling.

## Preparation

Make sure that both `hwave` package (this project) and `uhf_dry.out` of `StdFace` are installed.

## How to run

```bash
uhf_dry.out stan.in
python3 finiteT.py -u 4 -g 31 --max 1
python3 finiteT.py -u 8 -g 31 --max 3
python3 finiteT.py -u 12 -g 31 --max 5
python3 finiteT.py -u 24 -g 31 --max 7
```

After running the commands above, the result files `mag_U***.dat` and `Tneel_U***.dat` are outputted.
In `mag_U***.dat`, the first and second columns denote the temperature $T/t$ and the magnetic moment $m_z$, respectively.
The N\acute{e}el temperature $T_N$ is written in `Tneel_U***.dat`, which is defined as the lowest temperature where $m_z$ is smaller than $10^{-4}$ in the script.

`plot.plt` is a gnuplot script to plot $T$-dependence of $m_z$:

```bash
gnuplot plot.plt
```

![temperature dependence of magnetic momonet mz](./mz_HubbardCubic_L12.png)


If you want to show the help message for `finiteT.py`, run the following command. 
```bash
$python3 finiteT.py -h
usage: finiteT_uhfk.py [-h] -u U -g NGRID [--min TMIN] [--max TMAX]

T-dependence of magnetic properties for Hubbard model on cubic lattice by using UHFk

optional arguments:
  -h, --help            show this help message and exit
  -u U                  Hubbard U
  -g NGRID, --grid NGRID
                        # of grid
  --min TMIN            min of temperature/default:0.0
  --max TMAX            max of temperature/defalut:5.0

end
```
