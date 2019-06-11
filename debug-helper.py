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
                                    

prefix = "./"

# Read Data that is used for lookup table
valid1 = getunmodifieddata('valid-35mm-coco-01-data', 75)

lowb  = "TESTL:        .BYTE "
highb = "TESTH:        .BYTE "
for i in range(256):
    lowb += "${:02x},".format(int(valid1[0][i]) &  0xFF)
    highb += "${:02x},".format(int(valid1[0][i]) >> 8)
lowb = lowb.rstrip(",")
highb = highb.rstrip(",")
print(lowb)
print(highb)
