import sys
import os

if len(sys.argv) != 3:
    print("Please give a file name and number of measurements as argument")
    exit(-1)

filename = sys.argv[1]
measurements = int(sys.argv[2])

if not os.path.exists(filename):
    print("File {:s} does not exist".format(filename))
    exit(-1)

if measurements < 1 or measurements > 75:
    print("{:d} is an unreasonable number (0 < m <= 75)".format(measurements))
    exit(-1)

plainfile = filename[:-4]  # assumes .nib
rawfoldername = plainfile + "-raw"
datafoldername = plainfile + "-data"

if not os.path.exists(rawfoldername):
    os.mkdir(rawfoldername)
if not os.path.exists(datafoldername):
    os.mkdir(datafoldername)

dskfilename = filename.replace("nib", "dsk")    
os.system('nib2dsk {:s} {:s}'.format(filename, dskfilename))


for i in range(1, measurements + 1):
    applefilename = "M{:03d}".format(i)
    localfilename = os.path.join(rawfoldername, applefilename + ".bin")
    os.system('dos33 {:s} LOAD {:s} {:s}'.format(dskfilename,
                                                 applefilename,
                                                 localfilename))
    f = open(localfilename, 'rb')
    fb = f.read()
    if type(fb[0]) == int: # python3 treats bytes different from python2
        bdata = fb[4:] 
    else:
        bdata = [ord(a) for a in fb][4:]
    print(len(bdata))
    f.close()
    ldata = bdata[:256]
    hdata = bdata[256:512]
    lodata = bdata[512:768]
    hodata = bdata[768:]

    datfile = os.path.join(datafoldername, "M{:03d}.dat".format(i))
    data = [0] * 256
    for j in range(256):
        data[j] = hdata[j] * 256 + ldata[j]
    f = open(datfile, 'w')
    for j in range(256):
        f.write(str(data[j]) + '\n')
    f.close()

    originaldatfile = os.path.join(datafoldername, "M{:03d}_original.dat".format(i))
    data = [0] * 256
    for j in range(256):
        data[j] = hodata[j] * 256 + lodata[j]
    f = open(originaldatfile, 'w')
    for j in range(256):
        f.write(str(data[j]) + '\n')
    f.close()

