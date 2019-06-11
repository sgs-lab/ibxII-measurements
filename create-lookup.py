##%%
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import Divider, LocatableAxes, Size
import numpy as np
import sys
import pandas as pd
import os
import copy

from mpl_toolkits.axes_grid1 import Divider, LocatableAxes, Size

import intcal

binlist = [0x11, 0x27, 0x3d, 0x53, 0x69, 0x7f, 0x95, 0xab, 0xc1, 0xd7, 0xed, 0xff]
binlist0 = [0] + binlist
binwidth = [binlist0[i] - binlist0[i - 1] for i in range(1, len(binlist0))]

def getunmodifieddata(folder, length, eliminatelower = 0):
    filename = folder
    data = []
    for i in range(1, length + 1):
        f = open(prefix + filename + "/M{:03d}_original.dat".format(i), 'r')
        tdata = f.readlines()
        fdata = [int(n) for n in tdata]
        f.close()
        if eliminatelower > 0:
            for j in range(eliminatelower):
                fdata[j] = 0
        data.append(fdata)
    return data
                                    
calibrationpeak = [2614.511, 226, 6, 218, 240, 5]
lookuppeaks = [[2614.511, 226, 6, 218, 240, 5],
               [238.632, 15, 1, 12, 18, 3],
               [1274.5, 108, 3, 100, 116, 5], # Na-22
               [511, 40, 2, 33, 47, 5], # Na-22
               [661.64, 53, 2, 47, 59, 5], # Cs-137
               ]

prefix = "./"

# Read Data that is used for lookup table
nacslookup1 = getunmodifieddata('lookup-35mm-na-cs-01-data', 75)
nacslookup2 = getunmodifieddata('lookup-35mm-na-cs-02-data', 75)

# combine data to one list
lookupdata = nacslookup1 + nacslookup2

# calibrate data
calibratedlookupdata = []
peake = calibrationpeak[0]
peakpos = calibrationpeak[1]
width = calibrationpeak[2]
search = [calibrationpeak[3], calibrationpeak[4]]
for d in lookupdata:
    calibratedlookupdata.append(intcal.calibrate(d, peakpos, width, search))

# create lookup directory    
fromchannel = intcal.createlookupdirect(calibratedlookupdata, peakpos, peake, lookuppeaks)

for i in range(len(fromchannel)):
    print(i, fromchannel[i])

lowb  = "LOOKUPL:       .BYTE "
highb = "LOOKUPH:       .BYTE "
for i in range(256):
    if fromchannel[i] > 0:
        lowb += "${:02x},".format(int(fromchannel[i]) & 0xFF)
        highb += "${:02x},".format(int(fromchannel[i]) >> 8)
    else:
        lowb += "$ff,"
        highb += "$ff,"
lowb = lowb.rstrip(",")
highb = highb.rstrip(",")
print(lowb)
print(highb)
