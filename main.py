##################################################################
# coregister, crop, segment, and carve
##################################################################

#Compile the C carving program
import subprocess;
import os.path;
from sys import platform;

# only compile c programm on linux, binary for win32 included.
if platform != "win32":
    if not os.path.isfile("CarveIt.o"):
        p = subprocess.Popen("gcc -O3 CarveIt.c -lm -o CarveIt.o", shell=True);
        p.wait();

#Dummy class for storage
class Object(object):
    pass;

from HSVSegmentSeq import HSVSegmentSeq;
from TurntableCarve import TurntableCarve;
from CoregisterAndCrop import CoregisterAndCrop;

import numpy as np;

##################################################################
# initialization for 'CoregisterAndCrop'
#
# filenames
# name of original files consisting of basename, running number, and
# extension
fnin = Object();
fnin.base = 'pic/11-xl-2/';
fnin.number = list(set(range(1,37)));
fnin.extension = '.bmp';
# name of cropped and coregistered files consisting of basename,
# running number, and extension
fnroi = Object();
fnroi.base = 'pic/11-xl-2/ROI_';
fnroi.number = range(0, 360, 10);
fnroi.extension = '.png';
#
# rectangular region defining region of interest in the first input
# image, i.e. the cropping region
# crop_rect = np.array([255, 360, 215, 105]); # [xmin ymin width height].
crop_rect = np.array([245, 305, 230, 100]) # [xmin ymin width height].

# Define template and target region for coregistration
# region serving as template (region of first image)
# template_rect = np.array([220, 200, 250, 300]); # 'screw' template
template_rect = np.array([180, 250, 400, 300]) # 'tool' template
# roi in 'other' image, where template shall be sought in
# search_rect = np.array([150, 200, 600, 400]); # 'screw' region
search_rect = np.array([150, 150, 600, 400])  # 'tool' region
#
# coregister images in the input sequence, crop, and write to file
# store offsets for later usage in carving
offsets = CoregisterAndCrop(fnin,fnroi,crop_rect,search_rect,template_rect);
##################################################################

##################################################################
# initialization for 'HSVSegmentSeq'
#
# filenames
fnmask = Object();
fnmask.base = 'pic/11-xl-2/Mask_';
fnmask.number = range(0, 360, 10);
fnmask.extension = '.png';
#
# color interval of forground object in HSV space
Hint = [0, 255];
Sint = [0, 255];
Vint = [75, 255];
#
# segment seed using its HSV color value
HSVSegmentSeq(fnroi,fnmask,Hint,Sint,Vint);
##################################################################

##################################################################
# initialization for 'TurntableCarve'
#
# image and camera properties
cam = Object();
cam.orig_image_size = np.array([720, 540]); # original size of the image, needed for principal point
cam.crop_rect = crop_rect; # cropping rectangle
cam.offset = offsets; # offsets of the cropping regions
cam.alpha = range(0, -360, -10); # rotation angle
cam.PixPerMMAtZ = 145/5.74; #   #153/5.74 in home distance. 33/1.20 in the lab setup distance, 234/8.0; calibration value: pixel per mm at working depth: measure in image
cam.PixPerMMSensor = 1/0.0062; # 4.7�m pixel size (Nikon D7000, from specs)
cam.FocalLengthInMM = 12.5; # read from lens or from calibration
#
# tool in image
tool = Object();
tool.TurntableCenter = np.array([123,103]); # turntable center in roi: read from image [unit: pixel]
tool.tipY = 450; # y-position of tool tip in roi [unit: pixel]
tool.tipWidth = 20; # width of tip at its top [unit: pixel]
tool.bottomWidth = 200; # width of tip at its top [unit: pixel]
tool.margin = 10; # some extra margin to remove voxels due to machining inaccuracies of the tip etc. [unit: pixel]
#
# description of the reconstruction volume V as cuboid
V = Object();
V.VerticalOffset = 100; # Vertical offset of center of reconstruction cuboid (i.e the volume) in roi [unit: pixel]
V.VolWidth = 12.5; # width of the volume in mm (X-direction)
V.VolHeight =10.0; # height of the volume in mm (Y-direction)
V.VolDepth = 12.5; # depth of the volume in mm (Z-direction)
V.sX = 125; # number of voxels in X-direction
V.sY = 100; # number of voxels in Y-direction
V.sZ = 125; # number of voxels in Z-direction
#
# perform volume carving on mask images
volume_in_mm3 = TurntableCarve(fnmask,cam,tool,V);
##################################################################

##################################################################
# print result

print('Volume = ' + ("%0.2f" % volume_in_mm3) + 'mm^3\n');
##################################################################

