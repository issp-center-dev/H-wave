from __future__ import print_function

import itertools

import numpy as np


def read_geom(name_in):
    with open(name_in, 'r') as f:
        # skip header
        l_strip = [s.strip() for s in f.readlines()]
    norb = int(l_strip[3])

    rvec = np.zeros((3, 3), dtype=float)
    center = np.zeros((norb, 3), dtype=float)

    for count, line in enumerate(l_strip[:3]):
        values = line.split()
        rvec[count] = (float(values[0]), float(values[1]), float(values[2]))

    for count, line in enumerate(l_strip[4:]):
        values = line.split()
        center[count] = (float(values[0]), float(values[1]), float(values[2]))

    data = {"norb": norb, "rvec": rvec, "center": center}
    return data


def read_w90(name_in):

    with open(name_in, 'r') as f:
        # skip header
        l_strip = [s.strip() for s in f.readlines()[1:]]

    nr = int(l_strip[1])

    nints_per_line = 15
    skip_line = nr // nints_per_line
    if nr % nints_per_line != 0:
        skip_line += 1

    mat_size = len(l_strip[2 + skip_line:])
    data = {}
    for idx, line in enumerate(l_strip[2 + skip_line:]):
        values = line.split()
        # if data is empty, break
        if len(values) == 0:
            break
        irvec = (int(values[0]), int(values[1]), int(values[2]))
        orbvec = (int(values[3]) - 1, int(values[4]) - 1)
        data[(irvec, orbvec)] = float(values[5]) + 1J * float(values[6])

    return data

if __name__ == "__main__":
    path_to_sample = "../../sample/dir-model"
    import os
    data_geom = read_geom(os.path.join(path_to_sample, "zvo_geom.dat"))
    data_hr = read_w90(os.path.join(path_to_sample,"zvo_hr.dat"))
    print(data_geom)
    print(data_hr)