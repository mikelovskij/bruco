import os
import fnmatch
import numpy as np
import subprocess
# ---------------------------------Ligo import functions -----------------------------------------------------
try:
    import nds2

    from glue import lal
    from pylal import frutils, Fr

    # wrapper around the LIGO function to find where data is, returns a list of files
    def find_LIGO_data(observatory, gpsb, gpse):
        o = subprocess.Popen(["/usr/bin/gw_data_find", "-o", observatory[0],
                                "-t", observatory[0] + "1_R", "-s", str(gpsb), "-e", str(gpse), "-u", "file"],
                                stdout=subprocess.PIPE).communicate()[0]
        return o.splitlines()

except ImportError:
    print "Cannot import tne nds2, glue or pylal module. It should not be a problem unless the selected interferometer is Ligo H or L"


def ligo_list_channels(opt, gpsb, dt, scratchdir, minfs):
    # find the location of the GWF files
    files = find_LIGO_data(opt.ifo, gpsb, gpsb+dt)
    # create the cache file
    c = lal.Cache.from_urls(files)
    d = frutils.FrameCache(c, scratchdir=scratchdir, verbose=True)

    ###### Extract the list of channels and remove undesired ones ######################

    print ">>>>> Creating chache and extracting list of channels...."
    # read the list of channels and sampling rates from the first file
    firstfile = c[0].path
    os.system('/usr/bin/FrChannels ' + firstfile + ' > bruco.channels')
    f = open('bruco.channels')
    lines = f.readlines()
    channels = []
    sample_rate = []
    for l in lines:
        ll = l.split()
        if ll[0][1] != '0':
            # remove all L0/H0 channels
            channels.append(ll[0])
            sample_rate.append(int(ll[1]))

    # keep only channels with high enough sampling rate
    idx = find(sample_rate >= minfs)
    channels = channels[idx]
    sample_rate = sample_rate[idx]
    return np.array(channels), np.array(sample_rate), d

# ------------------------------ Virgo import functions ---------------------------------------------------------

try:
    import nap

    # Function to get data from raw
    def getRawData(ch, gps, dt):
        """Read data from RAW file:
        ch  = channel name
        gps = gps time
        dt  = duration
        """
        try:
            str1 = nap.FrameIChannel('/virgoData/ffl/raw.ffl', ch, dt, gps)
            d = nap.SeqViewDouble()
            str1.GetData(d)
            fs = 1.0 / d.GetSampling()
            return d.row(0), fs
        except:
            print "Error reading channel %s at gps %d + %d" % (ch, gps, dt)
            return -1

except ImportError:
    print "Cannot import tne NAP module. It should not be a problem unless the selected interferometer is Virgo"


def virgo_list_channels(opt, minfs, gpsb): #the exclusion part could be in common?
    if opt.channelfile is None:  # get the auxiliary channel list from the raw.ffl file
        command = "/virgoApp/Fr/v8r22/Linux-x86_64-SL6/FrDump.exe " \
                  + "-i /virgoData/ffl/raw.ffl -d 4 -f %d -l %d | grep ndata | grep -v V1:x=" % (gpsb, gpsb)
        p, g = os.popen4(command)
        L = g.readlines()
        channels = []
        for c in L:
            ch = c.split()[0].split('Vector:')[1]
            f = int(c.split('ndata=')[1].split()[0])/10
            # TODO: ndata conta quanti sample ci sono in un frame ed e' quindi vulnerabile a cambi di dimensione dei frame
            # TODO: si puo' sostituire usando un altra voce dumpata che da invece proprio il dt, controllo il nome e faccio prove?
            if (f >= minfs) & (f < 100000):  # there seems to be some error in the decimation if fs is too high
                channels.append(ch)
    else:  # use the provided auxiliary channel list instead
        f = open(opt.channelfile, 'r')
        L = f.readlines()
        channels = []
        for c in L:
            c = c.split()[0]
            channels.append(c)
        f.close()
    channels.sort()

    # remove the main channel from the auxiliary channels list to avoid the calculation of the coherence with itself
    channel = 'V1:' + opt.channel
    try:
        channels.remove(channel)
    except ValueError:
        pass

    channel = opt.channel
    try:
        channels.remove(channel)
    except ValueError:
        pass

    nchold = len(channels)  # number of channels before the exclusions
    print "Found %d channels\n\n" % nchold
    # make list unique, removing repeated channels, if any
    return np.unique(channels)

# -------------------------------------------------------- actual channel fetch functions -----------------------


def channel_fetch(ifo, channel, gpsb, gpse, d):
    if ifo == 'V1':
        ch, fs = getRawData(channel, gpsb, gpse - gpsb)
    else:
        buffer = d.fetch(ifo + ':' + channel, gpsb, gpse)
        ch = np.array(buffer)
        fs = len(ch) / (gpse - gpsb)
    return ch, fs