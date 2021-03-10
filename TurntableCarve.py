import numpy as np;
from scipy.misc import imread;
import mayavi.mlab as mlab;
import matplotlib.pyplot as plt;
import subprocess;
from sys import platform;


def CarveIt(V_in, P, mask, VolWidth, VolHeight, VolDepth):
    # Write output exchange file
    out = open('carveInput.dat', 'wb');

    sY = V_in.shape[0];
    sX = V_in.shape[1];
    sZ = V_in.shape[2];

    smY = mask.shape[0];
    smX = mask.shape[1];

    width = VolWidth;
    height = VolHeight;
    depth = VolDepth;

    data = np.array([sX, sY, sZ], dtype='uint32');
    out.write(data.tobytes());

    out.write(V_in.tobytes('A'));

    out.write(P.tobytes('F'));

    data = np.array([smX, smY], dtype='uint32');
    out.write(data.tobytes());
    out.write(mask.tobytes('F'));

    data = np.array([width, height, depth], dtype='float64');
    out.write(data.tobytes());
    out.close();

    # Execute
    if platform == "win32":
        execPrefix = "CarveIt.exe";
    else:
        execPrefix = "./CarveIt.o";

    p = subprocess.Popen(execPrefix + " carveInput.dat carveResult.dat", shell=True);
    p.wait();

    # Read input exchange file
    V_in = np.fromfile('carveResult.dat', dtype='uint8');
    return V_in.reshape((sY, sX, sZ), order='F');


def TurntableCarve(fn, cam, tool, V):
    # Reconstruct volume of an object from its projection masks acquired using
    # a turntable setup
    #
    # Input:
    # fn: describes names of input mask files
    # cam: camera properties
    # tool: tool properties, describes a cone-shaped, centered at center of turntable
    # V: describes the reconstruction volume, i.e. a cuboid-shaped region
    #
    # Output:
    # vol_in_mm3: volume of the reconstructed object in mm^3
    #
    # Usage: Example code defining the inputs
    # # image and camera properties
    # cam.orig_image_size = [4928 3264]; # original size of the image, needed for principal point
    # cam.crop_rect = crop_rect; # cropping rectangle
    # cam.offset = offsets; # offsets of the cropping regions (output of 'CoregisterAndCrop')
    # cam.alpha = -[0:10:350]; # rotation angle
    # cam.PixPerMMAtZ = 524/12; # calibration value: pixel per mm at working depth: measure in image
    # cam.PixPerMMSensor = 1/0.0047; # 4.7µm pixel size (Nikon D7000, from specs)
    # cam.FocalLengthInMM = 85; # read from lens or from calibration
    # #
    # # tool in image (cone-shaped, centered at center of turntable)
    # tool.TurntableCenter = [380 485]; # turntable center in roi: read from image [unit: pixel]
    # tool.tipY = 308; # y-position of tool tip in roi [unit: pixel]
    # tool.tipWidth = 174; # width of tip at its top [unit: pixel]
    # tool.bottomWidth = 269; # width of tip at its top [unit: pixel]
    # tool.margin = 10; # some extra margin to remove voxels due to machining inaccuracies of the tip etc. [unit: pixel]
    # #
    # # description of the reconstruction volume V as cuboid
    # V.VerticalOffset = 300; # Vertical offset of center of reconstruction cuboid (i.e the volume) in roi [unit: pixel]
    # V.VolWidth = 12.5; # width of the volume in mm (X-direction)
    # V.VolHeight =10.0; # height of the volume in mm (Y-direction)
    # V.VolDepth = 12.5; # depth of the volume in mm (Z-direction)
    # V.sX = 125; # number of voxels in X-direction
    # V.sY = 100; # number of voxels in Y-direction
    # V.sZ = 125; # number of voxels in Z-direction
    # #
    # # perform volume carving on mask images
    # volume_in_mm3 = TurntableCarve(fnmask,cam,tool,V);
    ###################################################################

    # tool in image
    tool.tipX = tool.TurntableCenter[0];  # x-position of tool tip in roi [unit: pixel]
    tool.height = tool.TurntableCenter[1] - tool.tipY;  # height of the tool [unit: pixel]

    # init volume V as cuboid
    V.ImageOfOrigin = tool.TurntableCenter - [0, V.VerticalOffset];  # Center of the reconstruction cuboid
    V.dx = V.VolWidth / V.sX;  # voxelsize in X-direction
    V.dy = V.VolHeight / V.sY;  # voxelsize in Y-direction
    V.dz = V.VolDepth / V.sZ;  # voxelsize in Z-direction
    V.Voxels = [V.sY, V.sX, V.sZ];  # number of voxels in volume
    V.vol = np.ones(tuple(V.Voxels), np.uint8);  # solid filled volume

    # print info
    print('Volume carving from masks\n');

    # loop images for carving
    NumImgs = len(fn.number);  # number of images
    for i in range(NumImgs):
        # print a point to show progress
        print('.')

        # read mask image
        mask = ReadImage(fn, i);

        # calculate projection matrix for this image
        P = ProjectionMatrix(cam, V, i);

        # draw the volume outline as box in mask image
        # projectVolBox(P,mask,V,10); # uncomment to see the projected volume boxes

        # do the carving with a c implementation
        V.vol = CarveIt(V.vol, P, mask, V.VolWidth, V.VolHeight, V.VolDepth);
        # alternatively do the carving in Python -- very slow!
        # V.vol = Carve(V.vol,P,mask,V.VolWidth,V.VolHeight,V.VolDepth);

    # print end of line
    print('\n');

    # show the reconstructed object
    rotatevolume(V, 11);

    # cut the tool-region from the reconstructed object
    V.vol = CutToolFromVolume(V, tool, cam);

    # show the reconstructed object again
    # rotatevolume(V,12);

    # calculate the final volume of the object
    vol_in_mm3 = np.sum(V.vol) * V.dx * V.dy * V.dz;
    return vol_in_mm3;


##########################################################
# read an image or mask from file
def ReadImage(fn, idx):
    imgfilename = fn.base + ("%04d" % fn.number[idx]) + fn.extension;
    img = imread(imgfilename);
    return img;


##########################################################
# carving routine -- very slow implementation, better use
# the c++ implementation CarveIt
def Carve(V_in, P, mask, VolWidth, VolHeight, VolDepth):
    V = V_in;
    sz = V.shape;
    sX = sz[1];
    sY = sz[0];
    sZ = sz[2];
    dx = VolWidth / sX;
    dy = VolHeight / sY;
    dz = VolDepth / sZ;
    x0 = -VolWidth / 2 + dx / 2;
    y0 = -VolHeight / 2 + dy / 2;
    z0 = -VolDepth / 2 + dz / 2;

    sm = mask.shape;

    for iz in range(sZ):
        for ix in range(sX):
            for iy in range(sY):
                x = x0 + ix * dx;
                y = y0 + iy * dy;
                z = z0 + iz * dz;
                Q = [x, y, z, 1];
                QT = np.matrix(Q).transpose();
                q = P * QT;
                q = q / q[2];

                qx = int(np.round(q[0])[0]);
                qy = int(np.round(q[1])[0]);
                if (qx >= 1 and qx <= sm[1] and qy >= 1 and qy <= sm[0]):
                    if (mask[qy, qx] == 0):
                        V[iy, ix, iz] = 0;
                else:
                    V[iy, ix, iz] = 0;
    return V


##########################################################
# show the reconstructed volume as isosurface
def showvolume(Vin, currentfigurenum):
    mlab.figure(currentfigurenum, bgcolor=(1, 1, 1), fgcolor=(1, 1, 1));
    mlab.clf();

    p = mlab.contour3d(Vin.vol, color=(1, 0, 0));
    mlab.text(0.05, 0.95, 'Please close the window to continue calculations.', color=(0, 0, 0), width=0.9);
    mlab.text(0.3, 0.05, 'Rotate using click&drag', color=(0, 0, 0), width=0.4);

    c_scene = mlab.get_engine().current_scene;
    # c_scene.scene.light_manager.light_mode = 'vtk';
    c_scene.scene.camera.position = [0, 0, -128];
    c_scene.scene.camera.view_up = [-1, 0, 0];
    c_scene.scene.render();
    mlab.show();
    return p;


##########################################################
# show the reconstructed volume as rotating isosurface
def rotatevolume(Vin, currentfigurenum):
    p = showvolume(Vin, currentfigurenum);

    # auto rotation skipped for python version
    # rotate manually if needed


##########################################################
# draw a box in the mask image to show, where the reconstruction
# region is located
def projectVolBox(P, mask, V, currentfigurenum):
    fig = plt.figure(figsize=(6, 3.2))

    ax = fig.add_subplot(111)
    ax.set_title('projectedVolBox')
    plt.imshow(mask)
    ax.set_aspect('equal')
    # plt.colorbar(orientation='vertical')

    ## corner points of a cube, centered at the origin of the world coord
    ## system. Homogeneous coords.
    x = V.VolWidth / 2;
    y = V.VolHeight / 2;
    z = V.VolDepth / 2;
    Q1 = np.array([x, y, z, 1]).reshape((1, 4));
    q1 = P * Q1.transpose();
    q1 = q1 / q1[2];

    Q2 = np.array([x, y, -z, 1]).reshape((1, 4));
    q2 = P * Q2.transpose();
    q2 = q2 / q2[2];

    Q3 = np.array([x, -y, z, 1]).reshape((1, 4));
    q3 = P * Q3.transpose();
    q3 = q3 / q3[2];

    Q4 = np.array([x, -y, -z, 1]).reshape((1, 4));
    q4 = P * Q4.transpose();
    q4 = q4 / q4[2];

    Q5 = np.array([-x, y, z, 1]).reshape((1, 4));
    q5 = P * Q5.transpose();
    q5 = q5 / q5[2];

    Q6 = np.array([-x, y, -z, 1]).reshape((1, 4));
    q6 = P * Q6.transpose();
    q6 = q6 / q6[2];

    Q7 = np.array([-x, -y, z, 1]).reshape((1, 4));
    q7 = P * Q7.transpose();
    q7 = q7 / q7[2];

    Q8 = np.array([-x, -y, -z, 1]).reshape((1, 4));
    q8 = P * Q8.transpose();
    q8 = q8 / q8[2];

    plt.plot([float(q1[0]), float(q2[0])], [float(q1[1]), float(q2[1])]);
    plt.plot([float(q1[0]), float(q3[0])], [float(q1[1]), float(q3[1])]);
    plt.plot([float(q1[0]), float(q5[0])], [float(q1[1]), float(q5[1])]);
    plt.plot([float(q2[0]), float(q4[0])], [float(q2[1]), float(q4[1])]);
    plt.plot([float(q2[0]), float(q6[0])], [float(q2[1]), float(q6[1])]);
    plt.plot([float(q3[0]), float(q4[0])], [float(q3[1]), float(q4[1])]);
    plt.plot([float(q3[0]), float(q7[0])], [float(q3[1]), float(q7[1])]);
    plt.plot([float(q4[0]), float(q8[0])], [float(q4[1]), float(q8[1])]);
    plt.plot([float(q5[0]), float(q6[0])], [float(q5[1]), float(q6[1])]);
    plt.plot([float(q5[0]), float(q7[0])], [float(q5[1]), float(q7[1])]);
    plt.plot([float(q6[0]), float(q8[0])], [float(q6[1]), float(q8[1])]);
    plt.plot([float(q7[0]), float(q8[0])], [float(q7[1]), float(q8[1])]);
    plt.show()


##########################################################
# calculate projection matrix P from given camera and image
# information
def ProjectionMatrix(cam, V, i):
    # pricipal point in image
    camH = cam.offset[:, i];
    p = cam.orig_image_size / 2 + 0.5 - cam.crop_rect[0:2] + 1 - camH;

    # Z = f*(M/m), FocalLength is f/m, m in [mm/pix], 1/M is PixPerMMAtZ
    f = cam.FocalLengthInMM * cam.PixPerMMSensor;
    Z = f / cam.PixPerMMAtZ;  # distance of rotation axis, i.e. origin of world coords.

    X = -(V.ImageOfOrigin[0] - p[0]) / cam.PixPerMMAtZ;
    Y = -(V.ImageOfOrigin[1] - p[1]) / cam.PixPerMMAtZ;

    K = np.matrix([[f, 0, p[0]], [0, f, p[1]], [0, 0, 1]]);  # calibration matrix
    c = np.cos(cam.alpha[i] / 180.0 * np.pi);
    s = np.sin(cam.alpha[i] / 180.0 * np.pi);
    R = np.matrix([[c, 0, -s], [0, 1, 0], [s, 0, c]]);  # rotation matrix
    t = [X, Y, -Z];  # translation vector, where camera is in world coords

    tT = np.matrix(t).transpose();
    P = K * np.matrix(np.concatenate((R, -tT), axis=1));
    return P


##########################################################
# cut a rotation symmetric cone-shaped tool with flat tip
# from the reconstructed volume

# Dummy class for storage
class Object(object):
    pass;


def CutToolFromVolume(V, tool, cam):
    # convert image-based information in mm and world coordinate
    # information
    toolInMM = Object();
    toolInMM.tipX = (tool.tipX - V.ImageOfOrigin[0]) / cam.PixPerMMAtZ;
    toolInMM.tipY = (tool.tipY - V.ImageOfOrigin[1]) / cam.PixPerMMAtZ;
    toolInMM.tipRadius = tool.tipWidth / cam.PixPerMMAtZ / 2;
    toolInMM.bottomRadius = tool.bottomWidth / cam.PixPerMMAtZ / 2;
    toolInMM.height = tool.height / cam.PixPerMMAtZ;
    toolInMM.margin = tool.margin / cam.PixPerMMAtZ;

    # border line of cone as r = A*y + B + margin
    A = (toolInMM.bottomRadius - toolInMM.tipRadius) / toolInMM.height;
    B = toolInMM.tipRadius - A * toolInMM.tipY + toolInMM.margin;
    Y_tcp = toolInMM.tipY;

    # the output volume and position of its center, i.e. the origin
    vol = V.vol;
    x0 = V.VolWidth / 2 + V.dx / 2;
    y0 = V.VolHeight / 2 + V.dy / 2;
    z0 = V.VolDepth / 2 + V.dz / 2;

    # loop volume
    for iz in range(int(V.sZ)):
        for ix in range(int(V.sX)):
            for iy in range(int(V.sY)):
                x = ix * V.dx - x0;
                y = iy * V.dy - y0;
                z = iz * V.dz - z0;
                r = np.sqrt(x * x + z * z);  # radius
                if (y > Y_tcp and r < A * y + B):  # below tool tip and within tool radius?
                    vol[iy, ix, iz] = 0;  # setting vol to 0 means deleting the object at this loaction
    return vol