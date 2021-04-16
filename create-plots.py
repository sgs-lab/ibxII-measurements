##%%
import matplotlib.pyplot as plt
import numpy as np
import sys
import pandas as pd
import os
import copy

import intcal

################################################################################
# Some helper functions
# (Could be moved to module, but were left here for now)
################################################################################

def getdata_multiple(files, key):
    data = []
    if key in files:
        d = files[key]
        filename = d[0]
        for i in range(1, d[2] + 1):
            f = open(prefix + filename + "/M{:03d}.dat".format(i), 'r')
            tdata = f.readlines()
            fdata = [int(n) for n in tdata]
            f.close()
            data.append(fdata)
    return data

def getunmodifieddata_multiple(files, key, eliminatelower = 0):
    data = []
    if key in files:
        d = files[key]
        filename = d[0]
        for i in range(1, d[2] + 1):
            f = open(prefix + filename + "/M{:03d}_original.dat".format(i), 'r')
            tdata = f.readlines()
            fdata = [int(n) for n in tdata]
            f.close()
            if eliminatelower > 0:
                for j in range(eliminatelower):
                    fdata[j] = 0
            data.append(fdata)
    return data

def getospreydata_multiple(files, key):
    data = []
    if key in files:
        d = files[key]
        filename = d[0]
        for i in range(d[2]):
            f = open(prefix + filename + "/m{:03d}.csv".format(i), 'r')
            tdata = f.readlines()
            ttdata = [int(n) for n in tdata]
            fdata = [0] * 256
            # reduce number of channels to 256
            div = 8
            for i in range(0, 2048, div):
                fdata[i // div] = sum(ttdata[i:i+div])
            f.close()
            data.append(fdata)
    return data

def loadrawdata(metadata):
    data = {}
    for x in metadata:
        data[x] = {}
        if metadata[x][META_TYPE] == DTYPE_APPLE_CAL:
            if(metadata[x][META_ELIMINATELOWER] > 0):
                data[x]['raw'] = getunmodifieddata_multiple(metadata, x, metadata[x][META_ELIMINATELOWER])
            else:
                data[x]['raw'] = getunmodifieddata_multiple(metadata, x)
            data[x]['apple'] = getdata_multiple(metadata, x)
        if metadata[x][META_TYPE] == DTYPE_NOAPPLE_CAL:
            data[x]['raw'] = getunmodifieddata_multiple(metadata, x)
        if metadata[x][META_TYPE] == DTYPE_OSPREY:
            data[x]['raw'] = getospreydata_multiple(metadata, x)
        if metadata[x][META_DROPLAST]:
            for i in range(metadata[x][META_RUNS]):
                data[x]['raw'][i][-1] = 0
    return data

def addcalibration(metadata, data, title, mul = False):
    for runname in metadata:
        data[runname][title] = []
        for i in range(metadata[runname][META_RUNS]):
            peakpos = metadata[runname][META_CALP][1]
            width = metadata[runname][META_CALP][2]
            search = [metadata[runname][META_CALP][3], metadata[runname][META_CALP][4]]
            tdata = intcal.calibrate(data[runname]['raw'][i], peakpos, width, search)
            if mul:
                s = sum(tdata)
                data[runname][title].append([int(x * metadata[runname][META_COUNTS] / s) for x in tdata])
            else:
                data[runname][title].append(tdata)
    return data

def addlookup(metadata, data, title, calibrationtitle, mul = False):
    for runname in metadata:
        data[runname][title] = []
        lookupdata = []
        for y in metadata[runname][META_LOKD]:
            lookupdata += data[y][calibrationtitle]
        peaks = metadata[runname][META_LOKP]
        peakpos = metadata[runname][META_CALP][1]
        peake = metadata[runname][META_CALP][0]
        print("Adding lookup for {:s}".format(runname))
        lookup = intcal.createlookupdirect(lookupdata, peakpos, peake, peaks)
        data[runname][title + "_lookup"] = lookup
        data[runname][title] = []
        for i in range(metadata[runname][META_RUNS]):
            tdata = intcal.adjustlookup(data[runname][calibrationtitle][i], lookup)
            if mul:
                s = sum(tdata)
                data[runname][title].append([int(x * metadata[runname][META_COUNTS] / s) for x in tdata])
            else:
                data[runname][title].append(tdata)
    return data

def selfchisquare(metadata, data):
    dfdict = {'runname': [], 'title': [], 'chi2sum': [], 'chi2avg': [], 'chi2min': [], 'chi2max': [], 'chi2std': []}
    for runname in metadata:
        for t in data[runname]:
            if not "_lookup" in t:
                dfdict['runname'].append(runname)
                dfdict['title'].append(t)
                # print(data[runname][t], runname, t)
                template = intcal.bindata(intcal.average(data[runname][t]), binlist)
                templateChiSList = []
                for i in range(metadata[runname][META_RUNS]):
                    test = intcal.bindata(data[runname][t][i], binlist)
                    templateChiSList.append(intcal.chisquare(template, test))
                dfdict['chi2sum'].append(sum(templateChiSList))
                dfdict['chi2avg'].append(np.mean(templateChiSList))
                dfdict['chi2min'].append(min(templateChiSList))
                dfdict['chi2max'].append(max(templateChiSList))
                dfdict['chi2std'].append(np.std(templateChiSList))
    return pd.DataFrame(dfdict)

################################################################################
# Some definitions to make live easier
################################################################################

prefix = "./"
DTYPE_APPLE_CAL = 1
DTYPE_NOAPPLE_CAL = 2
DTYPE_OSPREY = 3

META_FILENAME = 0
META_TITLE = 1
META_RUNS = 2
META_TYPE = 3
META_COUNTS = 4
META_ELIMINATELOWER = 5
META_DROPLAST = 6
META_CALP = 7
META_LOKP = 8
META_LOKD = 9

################################################################################
# Settings for bins, calibration and lookup table
################################################################################

binlist = [0x11, 0x27, 0x3d, 0x53, 0x69, 0x7f, 0x95, 0xab, 0xc1, 0xd7, 0xed, 0xff]
binlist0 = [0] + binlist
binwidth = [binlist0[i] - binlist0[i - 1] for i in range(1, len(binlist0))]

calibrationpeak = [2614.511, 226, 6, 214, 236, 5]
lookuppeaks = [[2614.511, 226, 6, 214, 236, 5],
               [238.632, 15, 1, 12, 18, 3],
               [1274.5, 108, 3, 100, 116, 5],
               [511, 40, 2, 33, 47, 5],
               [661.64, 53, 2, 47, 59, 5]]

ospcalibrationpeak = [2614.511, 226, 6, 155, 175, 5]
osplookuppeaks = [[2614.511, 226, 6, 220, 235, 5],
                  [238.632, 22, 1, 18, 26, 3],
                  [1274.5, 111, 3, 101, 121, 5],
                  [511, 47, 2, 41, 53, 5],
                  [661.64, 60, 2, 53, 67, 5]]


################################################################################
# Data to load
################################################################################

metadata = {'lookup1': ['lookup-35mm-na-cs-01-data',
                        'Na-22 and Cs-137, 35mm', 75, DTYPE_NOAPPLE_CAL, 4 * 2 ** 16, 0, False,
                        calibrationpeak, lookuppeaks, ['lookup1', 'lookup2']],
            'lookup2': ['lookup-35mm-na-cs-02-data',
                        'Na-22 and Cs-137, 35mm', 75, DTYPE_NOAPPLE_CAL, 4 * 2 ** 16, 0, False,
                        calibrationpeak, lookuppeaks, ['lookup1', 'lookup2']],
            'valid1': ['valid-35mm-coco-01-data',
                      '2 Co-60 sources, 35mm', 75, DTYPE_APPLE_CAL, 4 * 2 ** 16, 0, False,
                       calibrationpeak, lookuppeaks, ['lookup1', 'lookup2']],
            'valid2': ['valid-35mm-coco-02-data',
                      '2 Co-60 sources, 35mm', 75, DTYPE_APPLE_CAL, 4 * 2 ** 16, 0, False,
                       calibrationpeak, lookuppeaks, ['lookup1', 'lookup2']],
            'invalid1': ['invalid-35mm-coco-140mm-cs-01-data',
                      '2 Co-60 sources, 35mm, Cs-137, 140mm', 75, DTYPE_APPLE_CAL, 4 * 2 ** 16, 0, False,
                       calibrationpeak, lookuppeaks, ['lookup1', 'lookup2']],
            'invalid2': ['invalid-35mm-coco-140mm-cs-01-data',
                      '2 Co-60 sources, 35mm, Cs-137, 140mm', 75, DTYPE_APPLE_CAL, 4 * 2 ** 16, 0, False,
                       calibrationpeak, lookuppeaks, ['lookup1', 'lookup2']],
            'osp-nacs': ['osprey-35mm-na-cs',
                        'Na-22 and Cs-137, 35mm (osprey)', 150, DTYPE_OSPREY, 4 * 2 ** 16, 0, False,
                       ospcalibrationpeak, osplookuppeaks, ['osp-nacs']],
}

data = loadrawdata(metadata)#
data = addcalibration(metadata, data, 'calibrated', False)
data = addlookup(metadata, data, 'lookup', 'calibrated',  False)

# Set all channels below 150keV to zero
# for lookup and apple-lookup
# and renormalize to initial count rate
lowerchannel = int(np.ceil(150 / 2614.511 * 226))
for runname in metadata:
    if metadata[runname][META_TYPE] == DTYPE_APPLE_CAL:
        for i in range(metadata[runname][META_RUNS]):
            totbefore = sum(data[runname]['apple'][i])
            for j in range(lowerchannel):
                data[runname]['apple'][i][j] = 0
            totafter = sum(data[runname]['apple'][i])
            # for j in range(256):
            #     data[runname]['apple'][i][j] = data[runname]['apple'][i][j] * totbefore / totafter
    for i in range(metadata[runname][META_RUNS]):
        totbefore = sum(data[runname]['lookup'][i])
        for j in range(lowerchannel):
            data[runname]['lookup'][i][j] = 0
        totafter = sum(data[runname]['lookup'][i])
        # for j in range(256):
        #     data[runname]['lookup'][i][j] = data[runname]['lookup'][i][j] * totbefore / totafter

################################################################################
# For all plots
plt.style.use('default')
no = 0
alphaselection = [0.2, 0.6]
bincolor = '0.5'
binlinecolor = 'white'
histcolor = '0.8'
passthreshold = 0x30
xw = 8 * 0.6
yw = 6 * 0.6


energyticks = [0, 500, 1000, 1500, 2000, 2500, 3000]
energychannels = [x / 2614.511 * 226 for x in energyticks]

energyx = [x / 226 * 2614.511 for x in range(256)]

def plotdefault():
    plt.grid(linestyle=":")
    plt.yscale("log")
    plt.xticks(energychannels, energyticks)
    plt.xlim(-3, 259)
    # plt.xlim(energychannels[0] - 0.03 * energychannels[-1], energychannels[-1] + 0.03 * energychannels[-1])
    plt.xlim(energychannels[0] - 0.02 * energychannels[-1], energychannels[-1] + 0.02 * energychannels[-1])
    # plt.xlim(energychannels[0], energychannels[-1])
    # plt.ylim(5, 1e4)
    plt.ylim(0.77, 1.3e4)
    plt.yticks([1, 10, 100, 1000, 10000], [1, 10, 100, 1000, 10000])
    
nacslines = [238.632, 511, 661.64, 1274.5, 2614.511]
nacschannels = [x / 2614.511 * 226  for x in nacslines]
def nacspeaklines():
    for p in nacschannels:
        ax.axvline(p)

colines = [238.632, 1173.228, 1332.492, 2614.511]
cochannels = [x / 2614.511 * 226  for x in colines]
def copeaklines():
    for p in cochannels:
        ax.axvline(p)

    
################################################################################
# See difference of apple calibration / lookup vs python calibration / lookup
# will use python calibration / lookup in the following

i = 0
fig, ax = plt.subplots()
ax.plot(data['valid1']['apple'][i], label="Apple Lookup")
ax.plot(data['valid1']['lookup'][i], label="Python Lookup")
plotdefault()
copeaklines()
plt.legend()
plt.savefig("figures/unused-figure-lookup-test.pdf")
plt.show()

################################################################################
# Create data set of valid / invalid items

# apple calibration/lookup
valid = copy.deepcopy(data['valid1']['apple'])
valid += data['valid2']['apple']

invalid = copy.deepcopy(data['invalid1']['apple'])
invalid += data['invalid2']['apple']

template = intcal.average(valid)

# python calibration/lookup
validpython = copy.deepcopy(data['valid1']['lookup'])
validpython += data['valid2']['lookup'] 

invalidpython = copy.deepcopy(data['invalid1']['lookup'])
invalidpython += data['invalid2']['lookup'] 

templatepython = intcal.average(validpython)

################################################################################
# Check Apple Lookup Template against Python Lookup Template
good = True
for i in range(256):
    if template[i] != templatepython[i]:
        print("Inconsistency, channel {:d}".format(i))
        good = False
if good:
    print("Template Lookup check: Good!")

###############################################################################
# Osprey lookup vs apple lookup
# Here, we need to use the python calibration and lookup, because lookup data
# for the apple was only created with the lookup measurement
# (but Python and Apple lookup routines yield exactly same results)

osp = intcal.average(data['osp-nacs']['lookup'])
ibx = intcal.average(data['lookup1']['lookup'] + data['lookup2']['lookup'])

ibxzeros = 0
while(ibx[ibxzeros] == 0.0):
    ibxzeros += 1
ibxendzeros = 255
while(ibx[ibxendzeros] == 0.0):
    ibxendzeros -= 1
ibxendzeros += 1    

ospzeros = 0
while(osp[ospzeros] == 0.0):
    ospzeros += 1
ospendzeros = 255
while(osp[ospendzeros] == 0.0):
    ospendzeros -= 1
ospendzeros += 1

peakthoriumratio = sum(osp[224:229]) / sum(ibx[224:229]) # thorium peak in 226
peakco1ratio = sum(osp[101:103]) / sum(ibx[101:103]) # lower co-60 peak in 101.39
peakco2ratio = sum(osp[115:117]) / sum(ibx[115:117]) # higher co-60 peak in 115.22
print(peakco1ratio, peakco2ratio, peakthoriumratio)


fig, ax = plt.subplots()

print("Total Counts Osprey: {:f}".format(sum(osp)))
print("Total Counts IBX: {:f}".format(sum(ibx)))

ax.plot(range(ospzeros, ospendzeros), osp[ospzeros:ospendzeros], label = "Osprey")
ax.plot(range(ibxzeros, ibxendzeros), ibx[ibxzeros:ibxendzeros], label = "IBX II")

plotdefault()
nacspeaklines()
plt.legend()
plt.savefig("figures/unused-figure-osprey-vs-ibxII.pdf")
plt.show()


ratio = peakthoriumratio
ospadjusted = [x / ratio for x in osp]

fig, ax = plt.subplots()

print("Total Counts Osprey Adjusted by Thorium Ratio: {:f}".format(sum(ospadjusted)))
print("Total Counts IBX: {:f}".format(sum(ibx)))

#osprey data adjusted so that 2.614 MeV peak heights match
ax.plot(range(ospzeros, ospendzeros), ospadjusted[ospzeros:ospendzeros], label = "Osprey")
ax.plot(range(ibxzeros, ibxendzeros), ibx[ibxzeros:ibxendzeros], label = "IBX II")

plotdefault()
nacspeaklines()
plt.legend()
plt.savefig("figures/unused-osprey-vs-ibxII-thpeak-adjusted.pdf")
plt.show()

ratio = sum(osp) / sum(ibx)
ospadjusted = [x / ratio for x in osp]

fig, ax = plt.subplots()

print("Total Counts Osprey Adjusted by Total countrate: {:f}".format(sum(ospadjusted)))
print("Total Counts IBX: {:f}".format(sum(ibx)))

#osprey data adjusted so that 2.614 MeV peak heights match
ax.plot(range(ospzeros, ospendzeros), ospadjusted[ospzeros:ospendzeros], label = "Osprey")
ax.plot(range(ibxzeros, ibxendzeros), ibx[ibxzeros:ibxendzeros], label = "IBX II")

plotdefault()
nacspeaklines()
plt.legend()
plt.savefig("figures/fig3-osprey-vs-ibxII-counts-adjusted.pdf")
plt.show()


ratio = 1.5
ospadjusted = [x / ratio for x in osp]
fig, ax = plt.subplots()

ax.plot(range(ospzeros, ospendzeros), ospadjusted[ospzeros:ospendzeros], label = "Osprey")
ax.plot(range(ibxzeros, ibxendzeros), ibx[ibxzeros:ibxendzeros], label = "IBX II")

plotdefault()
nacspeaklines()
plt.legend()
plt.savefig("figures/unused-figure-osprey-vs-ibxII-space.pdf")
plt.show()

## measure resolution
# average
osp = intcal.average(data['osp-nacs']['lookup'])
ibx = intcal.average(data['lookup1']['lookup'] + data['lookup2']['lookup'])

peaks = [[170, 330],
         [440, 600],
         [550, 780],
         [1200, 1350],
         [2450, 2800]]

fig, ax = plt.subplots(1, 5)
fig.set_size_inches(25, 8)

for i in range(5):
    s = int(peaks[i][0] * 226 / 2614.511)
    e = int(peaks[i][1] * 226 / 2614.511)
    ax[i].plot(energyx[s:e], ibx[s:e])
    ax[i].grid()
    ax[i].set_ylim(ymin=0)
plt.savefig("figures/unused-figure-resolution-ibxii.pdf")
plt.show()

fig, ax = plt.subplots(1, 5)
fig.set_size_inches(25, 8)

for i in range(5):
    s = int(peaks[i][0] * 226 / 2614.511)
    e = int(peaks[i][1] * 226 / 2614.511)
    ax[i].plot(energyx[s:e], osp[s:e])
    ax[i].grid()
    ax[i].set_ylim(ymin=0)
plt.savefig("figures/unused-figure-resolution-osp.pdf")
plt.show()

################################################################################
# Compare one invalid against template (average)
i = 0
pdata = invalid[i]
templatezeros = 0
while(template[templatezeros] == 0.0):
    templatezeros += 1
templateendzeros = 255
while(template[templateendzeros] == 0.0):
    templateendzeros -= 1
templateendzeros += 1

pdatazeros = 0
while(pdata[pdatazeros] == 0.0):
    pdatazeros += 1
pdataendzeros = 255
while(pdata[pdataendzeros] == 0.0):
    pdataendzeros -= 1
pdataendzeros += 1

fig, ax = plt.subplots()
ax.plot(range(templatezeros, templateendzeros), template[templatezeros:templateendzeros])
ax.plot(range(pdatazeros, pdataendzeros), pdata[pdatazeros:pdataendzeros])
plotdefault()
copeaklines()
plt.savefig("figures/fig4-invalid-vs-template.pdf", transparent=True)
plt.show()

################################################################################
# Compare one invalid against one valid
i = 0
pdata = invalid[i]
pdata2 = valid[i]

pdatazeros = 0
while(pdata[pdatazeros] == 0.0):
    pdatazeros += 1
pdataendzeros = 255
while(pdata[pdataendzeros] == 0.0):
    pdataendzeros -= 1
pdataendzeros += 1

pdata2zeros = 0
while(pdata2[pdata2zeros] == 0.0):
    pdata2zeros += 1
pdata2endzeros = 255
while(pdata2[pdata2endzeros] == 0.0):
    pdata2endzeros -= 1
pdata2endzeros += 1

fig, ax = plt.subplots()
ax.plot(range(pdata2zeros, pdata2endzeros), pdata2[pdata2zeros:pdata2endzeros])
ax.plot(range(pdatazeros, pdataendzeros), pdata[pdatazeros:pdataendzeros])
plotdefault()
copeaklines()
plt.savefig("figures/unused-invalid-vs-valid.pdf", transparent=True)
plt.show()

################################################################################
# Make bins for one invalid
# starting from 1 (for better use in Illustrator)

idx = 0
pdata = invalid[idx]
pdatazeros = 0
while(pdata[pdatazeros] == 0.0):
    pdatazeros += 1
pdataendzeros = 255
while(pdata[pdataendzeros] == 0.0):
    pdataendzeros -= 1
pdataendzeros += 1

fig, ax = plt.subplots()
last = 1
count = 0
for i in binlist:
    ax.axvspan(last, i, facecolor=bincolor, alpha=alphaselection[count % 2],zorder=0)
    last = i
    count += 1

x = [0] + binlist
for i in range(len(x) - 1):
    ax.axvline(x[i], color=binlinecolor, linewidth=0.7)
ax.axvline(x[-1], color=binlinecolor, linewidth=0.7)
    
ax.plot(range(pdatazeros, pdataendzeros), pdata[pdatazeros:pdataendzeros])
ax.set_yscale("log")
ax.set_ylim(0.7, 40000)

bins = intcal.bindata(invalid[idx], binlist)
binlist1 = [1] + binlist

avbins = [bins[i] / (binlist0[i+1] - binlist0[i]) for i in range(len(binlist))]

for i in range(len(binlist)):
    ax.set_prop_cycle(None)

    ax.bar(binlist1[i], avbins[i], binwidth[i], align='edge', alpha=0.5)
    ax.text(binlist0[i] + binwidth[i] / 2, avbins[i] / 2, "{:.1f} ({:.1f})".format(avbins[i], bins[i]), rotation="vertical", horizontalalignment="center")
    

plotdefault()
plt.savefig("figures/unused-figure-invalid-bins-bars.pdf")
plt.show()

################################################################################
# Make bins for one invalid
# lines instead of bars - even better for illustrator

barchartenergy = [0] + binlist

idx = 0
pdata = invalid[idx]
pdatazeros = 0
while(pdata[pdatazeros] == 0.0):
    pdatazeros += 1
pdataendzeros = 255
while(pdata[pdataendzeros] == 0.0):
    pdataendzeros -= 1
pdataendzeros += 1

fig, ax = plt.subplots()
last = 0
count = 0
for i in binlist:
    ax.axvspan(last, i, facecolor=bincolor, alpha=alphaselection[count % 2],zorder=0)
    last = i
    count += 1

x = [0] + binlist
for i in range(len(x) - 1):
    ax.axvline(x[i], color=binlinecolor, linewidth=0.7)
ax.axvline(x[-1], color=binlinecolor, linewidth=0.7)
    
ax.plot(range(pdatazeros, pdataendzeros), pdata[pdatazeros:pdataendzeros])

bins = intcal.bindata(pdata, binlist)

avbins = [bins[i] / (binlist0[i+1] - binlist0[i]) for i in range(len(binlist))]

ax.step(barchartenergy, [avbins[0]] + avbins)
for i in range(len(binlist)):
    ax.text(binlist0[i] + binwidth[i] / 2, avbins[i] / 2, "0x{:06X}".format(int(bins[i])), rotation="vertical", horizontalalignment="center")

plotdefault()
plt.savefig("figures/unused-invalid-bins.pdf")
plt.show()

out = ""
for i in intcal.bindata(pdata, binlist):
    out += "0x{:06X} ".format(int(i))
print(out)

################################################################################
# Make bins for template
# lines instead of bars - even better for illustrator

barchartenergy = [0] + binlist

pdata = template
pdatazeros = 0
while(pdata[pdatazeros] == 0.0):
    pdatazeros += 1
pdataendzeros = 255
while(pdata[pdataendzeros] == 0.0):
    pdataendzeros -= 1
pdataendzeros += 1

fig, ax = plt.subplots()
last = 0
count = 0
for i in binlist:
    ax.axvspan(last, i, facecolor=bincolor, alpha=alphaselection[count % 2],zorder=0)
    last = i
    count += 1

x = [0] + binlist
for i in range(len(x) - 1):
    ax.axvline(x[i], color=binlinecolor, linewidth=0.7)
ax.axvline(x[-1], color=binlinecolor, linewidth=0.7)
    
ax.plot(range(pdatazeros, pdataendzeros), pdata[pdatazeros:pdataendzeros])

bins = intcal.bindata(pdata, binlist)

avbins = [bins[i] / (binlist0[i+1] - binlist0[i]) for i in range(len(binlist))]

ax.step(barchartenergy, [avbins[0]] + avbins)
for i in range(len(binlist)):
    ax.text(binlist0[i] + binwidth[i] / 2, avbins[i] / 2, "0x{:06X}".format(int(bins[i])), rotation="vertical", horizontalalignment="center")

plotdefault()
plt.savefig("figures/unused-template-bins.pdf")
plt.show()

out = ""
for i in intcal.bindata(template, binlist):
    out += "0x{:06X} ".format(int(i))
print(out)

################################################################################
# Make bins for one invalid and the template
# lines instead of bars - even better for illustrator

barchartenergy = [0] + binlist

idx = 0
pdata = invalid[idx]
pdatazeros = 0
while(pdata[pdatazeros] == 0.0):
    pdatazeros += 1
pdataendzeros = 255
while(pdata[pdataendzeros] == 0.0):
    pdataendzeros -= 1
pdataendzeros += 1

pdata2 = template
pdata2zeros = 0
while(pdata2[pdata2zeros] == 0.0):
    pdata2zeros += 1
pdata2endzeros = 255
while(pdata2[pdata2endzeros] == 0.0):
    pdata2endzeros -= 1
pdata2endzeros += 1

fig, ax = plt.subplots()
last = 0
count = 0
for i in binlist:
    ax.axvspan(last, i, facecolor=bincolor, alpha=alphaselection[count % 2],zorder=0)
    last = i
    count += 1

x = [0] + binlist
for i in range(len(x) - 1):
    ax.axvline(x[i], color=binlinecolor, linewidth=0.7)
ax.axvline(x[-1], color=binlinecolor, linewidth=0.7)
    
ax.plot(range(pdatazeros, pdataendzeros), pdata[pdatazeros:pdataendzeros])

bins = intcal.bindata(pdata, binlist)
avbins = [bins[i] / (binlist0[i+1] - binlist0[i]) for i in range(len(binlist))]

ax.step(barchartenergy, [avbins[0]] + avbins, label="Invalid Item")
for i in range(len(binlist)):
    ax.text(binlist0[i] + binwidth[i] / 2, avbins[i] / 2, "0x{:06X}".format(int(bins[i])), rotation="vertical", horizontalalignment="center")


templatebins = intcal.bindata(pdata2, binlist)
templateavbins = [templatebins[i] / (binlist0[i+1] - binlist0[i]) for i in range(len(binlist))]
ax.step(barchartenergy, [templateavbins[0]] + templateavbins, label="Template")

plotdefault()
plt.legend()
plt.savefig("figures/fig5-invalid-template-bins.pdf")
plt.show()

################################################################################
# Compare items to template

validChiSList = []
invalidChiSList = []
templatebin = intcal.bindata(template, binlist)
for i in range(100):
    testbin = intcal.bindata(valid[i], binlist)
    validChiSList.append(intcal.chisquare(templatebin, testbin))
    testbin = intcal.bindata(invalid[i], binlist)
    invalidChiSList.append(intcal.chisquare(templatebin, testbin))

print("Average Valid Chi2 {:f} +- {:f}".format(np.mean(validChiSList), np.std(validChiSList)))
print("Average Invalid Chi2 {:f} +- {:f}".format(np.mean(invalidChiSList), np.std(invalidChiSList)))

fig, ax = plt.subplots()

ax.axvspan(0, passthreshold, facecolor='#4aac0f', zorder=0)
ax.axvspan(passthreshold, 800, facecolor='orangered', zorder=0)

ax.hist(validChiSList, range = (0, 200), bins = 30, label="Valid Item", color=histcolor)
ax.hist(invalidChiSList, range = (0, 800), bins = 120, label="Invalid Item", color=histcolor)

plt.grid(linestyle=":")

# ylim = 55
ylim = 60
xlim = 800
plt.ylim(0 - 0.03 * ylim, ylim + 0.03 * ylim)
plt.xlim(0 - 0.02 * xlim, xlim + 0.02 * xlim)
plt.savefig("figures/fig6-valid-invalid-margins.pdf")
plt.show()

fig, ax = plt.subplots()

ax.axvspan(0, passthreshold, facecolor='#4aac0f', zorder=0)
ax.axvspan(passthreshold, 800, facecolor='orangered', zorder=0)

ax.hist(validChiSList, range = (0, 200), bins = 30, label="Valid Item", color=histcolor)
ax.hist(invalidChiSList, range = (0, 800), bins = 120, label="Invalid Item", color=histcolor)

plt.grid(linestyle=":")

ylim = 60
xlim = 800
plt.ylim(0 - 0.03 * ylim, ylim + 0.03 * ylim)
plt.xlim(0 - 0.02 * xlim, xlim + 0.02 * xlim)

plt.xticks(range(0, xlim, 96), ["0x{:X}".format(i) for i in range(0, xlim, 96)])
plt.ylim(0 - 0.03 * ylim, ylim + 0.03 * ylim)
plt.xlim(0 - 0.02 * xlim, xlim + 0.02 * xlim)
plt.savefig("figures/unused-figure-valid-invalid-margins-hex.pdf")
plt.show()


