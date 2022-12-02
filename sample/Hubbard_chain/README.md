# Sample for plotting band

## What's this sample?

This is the sample for calculating the band of a tight-binding model on chain.

## Preparation

Make sure that both `hwave` package (this project) and `uhf_dry.out` of `StdFace` are installed.

## How to run

```bash
uhf_dry.out stan.in
python3 output_band.py  
```

After running the commands above, the result file `band.dat` is outputted.
In `band.dat`, the first and second columns denote the index of wavenumber and the band energies, respectively.

`plot.plt` is a gnuplot script to plot the band:

```bash
gnuplot plot.plt
```

![band of the tight-binding chain](./band_L16.png)
