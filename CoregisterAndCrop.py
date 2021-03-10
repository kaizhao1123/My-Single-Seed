import numpy as np;
from scipy.misc import imread, imsave;
import cv2;


def CoregisterAndCrop(fnin, fnout, crop_rect, search_rect, template_rect):
    # Coregister images in a sequence using normalized cross-correlation on
    # template region as target.
    #
    # Input:
    # fnin: describes names of input files
    # fnout: describes names of output files, i.e. cropped regions
    # crop_rect: region in the first image which shall be cropped in all images
    # search_rect: search region, where template region is searched
    # template_rect: template region in first image
    #
    # Usage: Example code defining the inputs
    # # filenames
    # # name of original files consisting of basename, running number, and
    # # extension
    # fnin.base = '../images/DSC_';
    # fnin.number = [3730:3756 3758:3766];
    # fnin.extension = '.jpg';
    # # name of cropped and coregistered files consisting of basename,
    # # running number, and extension
    # fnroi.base = '../images/ROI_';
    # fnroi.number = [0:10:350];
    # fnroi.extension = '.png';
    #
    # # rectangular region defining region of interest in the first input
    # # image, i.e. the cropping region
    # crop_rect = [1955 1990 760 564]; # [xmin ymin width height].
    #
    # # Define template and target region for coregistration
    # # region serving as template (region of first image)
    # # template_rect = [3180 1200 420 540]; # 'screw' template
    # template_rect = [2197 2400 293 100]; # 'tool' template
    # # roi in 'other' image, where template shall be sought in
    # # search_rect = [3000 1070 700 800]; # 'screw' region
    # search_rect = [2168 2350 362 200]; # 'tool' region
    #
    # offsets = CoregisterAndCrop(fnin,fnroi,crop_rect,search_rect,template_rect);
    ####################################################################

    # offsets for each image
    offsets = np.zeros((2, len(fnout.number)), np.float);

    # print a point to show progress
    print('Coregister and crop\n.')

    # process first image
    im1 = ReadImage(fnin, 0);
    CropAndWrite(im1, crop_rect, fnout, 0);

    # coregister the others
    NumImgs = len(fnin.number);
    for i in range(1, NumImgs):
        # print a point to show progress
        print('.')

        # read image
        img = ReadImage(fnin, i);

        # calculate offsets by normalized cross-correlation
        offsets[:, i] = GetOffset(im1, template_rect, img, search_rect);

        # crop at offset location
        rect_crop = crop_rect + np.array([offsets[0, i], offsets[1, i], 0, 0]);
        CropAndWrite(img, rect_crop, fnout, i);

    # print end of line
    print('\n');
    return offsets;


def ReadImage(fn, idx):
    imgfilename = fn.base + ("%04d" % fn.number[idx]) + fn.extension;
    img = imread(imgfilename);
    return img;


def WriteImage(img, fn, idx):
    imgfilename = fn.base + ("%04d" % fn.number[idx]) + fn.extension;
    imsave(imgfilename, img);


def imcrop(img, rect_roi):
    return img[int(rect_roi[1]): int(rect_roi[1] + rect_roi[3]),
           int(rect_roi[0]): int(rect_roi[0] + rect_roi[2])];


def CropAndWrite(img, rect_roi, fn, idx):
    sub_roi = imcrop(img, rect_roi);
    WriteImage(sub_roi, fn, idx);


def GetOffset(im1, rect_template, im2, rect_roi):
    # closely following the Matlab example (see help on 'normxcorr2')
    sub_template = imcrop(im1, rect_template);
    sub_roi = imcrop(im2, rect_roi);

    # normalized cross correlation
    c = cv2.matchTemplate(sub_template[:, :, 0], sub_roi[:, :, 0], cv2.TM_CCORR_NORMED);

    # offset found by correlation
    corr_offset = cv2.minMaxLoc(c)[3];

    # relative offset of position of subimages
    rect_offset = np.array([(rect_roi[0] - rect_template[0]),
                            (rect_roi[1] - rect_template[1])]);

    # total offset
    offset = corr_offset + rect_offset;
    return offset;