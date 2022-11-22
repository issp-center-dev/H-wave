#!/usr/bin/env python3

import os,sys

tol = 1.0e-12

if len(sys.argv) < 3:
    print("usage: {} file1 file2 [tolerance]".format(os.path.basename(sys.argv[0])))
    exit(1)

def readfile(filename):
    tbl = {}
    with open(filename, "r") as fh:
        lines = fh.read().splitlines()
    for line in lines:
        k,v = line.split('=')
        tbl[k.strip()] = float(v)
    return tbl
        
w1 = readfile(sys.argv[1])
w2 = readfile(sys.argv[2])
if len(sys.argv) >= 4:
    tol = float(sys.argv[3])

err = 0
for k,v1 in w1.items():
    v2 = w2.get(k, None)
    if v2:
        if abs(v1 - v2) > tol:
            print("key={}, value={}, {}, diff={:.6e}".format(k,v1,v2,abs(v1-v2)))
            err += 1

if err > 0:
    exit(1)
else:
    print("ok.")
    exit(0)
