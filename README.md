Brute force coherence (Gabriele Vajente, 2015-06-11)
****************************************************

Command line arguments (with default values)
--ifo                       interferometer prefix (no default, must specify)
--channel=OAF-CAL_DARM_DQ   name of the main channel
--gpsb=1087975458           starting time
--length=180                amount of data to use (in seconds)
--outfs=8192                sampling frequency of the output results (coherence will
                            be computed up to outfs/2 if possible)
--minfs=512                 skip all channels with samplig frequency smaller than this
--naver=100                 number of averages to compute the coherence
--dir=bruco                 output directory
--top=100                   for each frequency, save to cohtab.txt and idxtab.txt this
                            maximum number of coherence channels
--webtop=20                 show this number of coherence channels per frequency, in the
                            web page summary
--plot=png                  plot format (png, pdf, none)
--nproc                     number of processes to use (if not specified, use all CPUs)
--calib                     name of a text file containing the calibration transfer 
                            function to be applied to the target channel spectrum, in 
                            a two column format (frequency, absolute value)
--xlim                      limits for the frequency axis in plots, in the format fmin:fmax
--ylim                      limits for the y axis in PSD plots, in the format ymin:ymax

Example:
./bruco.py --ifo=H1 --channel=CAL-DELTAL_EXTERNAL_DQ --calib=lho_darm_calibration.txt --gpsb=1111224016 
           --length=600 --outfs=4096 --naver=100 --dir=./bruco_1111224016 --top=100 --webtop=20 --xlim=1:1000 --ylim=1e-20:1e-14

To properly run the script, you need to setup a couple of more things. On line 52 you
must define where the lits of excluded channels is. The default is the same directory,
in a file called 'bruco_excluded_channels.txt'. On line 54 you must set up a scratch
directory where the frame library will save the temporary GWF files.


CHANGELOG:

2015-01-29 added linear detrending in PSD and coherence to improve low frequency bins
2015-02-16 - split coherence computation into PSDs and CSDs to remove redundant compuation of main channel PSD (still room to improve here)
           - removed upsampling of channels with lower sampling rate (now coherences are always computed with the minimum possible number of samples)
2015-02-18 - further splitting of coherence into primitive FFTs
2015-02-19 - allow selection of plot output format (PDF or PNG)
2015-02-24 - corrected typo in options lenght -> length
           - implemented parallel processing
2015-06-11 - using gw_data_find to locate the GWF files
2015-06-16 - added calibration transfer function option
           - added expanduser to all paths to allow use of ~
