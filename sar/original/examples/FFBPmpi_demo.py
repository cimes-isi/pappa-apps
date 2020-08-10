##############################################################################
#                                                                            #
#  This is a demonstration of the Fast Factorized Backprojection algorithm.  #
#  Data sets can be switched in and out by commenting/uncommenting the lines #
#  of code below.                                                            #
#                                                                            #
##############################################################################

#Add include directories to default path list
import sys
sys.path.append('../')
sys.path.append('./dictionaries')
if len(sys.argv) > 1 and sys.argv[1] == '--batch':
    batch = True
else:
    batch = False

#Include Dictionaries
from SARplatform import plat_dict

#Include standard library dependencies
import matplotlib.pylab as plt
from time import time

#Include SARIT toolset
from ritsar import phsTools
from ritsar import phsRead
from ritsar import imgTools
from mpi4py import MPI

if __name__ == '__main__':
    #AFRL DSBP demo
    ###############################################################################
    #Define top level directory containing *.mat file
    #and choose polarization and starting azimuth
    pol = 'HH'
    directory = './data/AFRL/pass1'
    start_az = 1

    comm = MPI.COMM_WORLD
    # if comm.Get_size() != 4:
    #     print('Must use exactly 4 MPI ranks')
    #     exit(1)

    #Import phase history and create platform dictionary
    [phs, platform] = phsRead.AFRL(directory, pol, start_az, n_az = 4)

    #Create image plane dictionary
    img_plane = imgTools.img_plane_dict(platform, res_factor = 1.0, upsample = True, aspect = 1.0)

    #Fast-factorized backprojection with MPI
    start = time()
    img_FFBPmpi = imgTools.FFBPmpi(phs, platform, img_plane, taylor = 17, factor_max = 2)
    fbpmpi_time = time()-start

    # Only rank 0 should test other algorithms and produce plots
    if comm.Get_rank() != 0:
        exit()

    # #full backprojection
    # start = time()
    # img_bp   = imgTools.backprojection(phs, platform, img_plane, taylor = 17, upsample = 2)
    # bp_time = time()-start

    # #Fast-factorized backprojection without multi-processing
    # start = time()
    # img_FFBP = imgTools.FFBP(phs, platform, img_plane, taylor = 17, factor_max = 2)
    # fbp_time = time()-start

    # #Fast-factorized backprojection with multi-processing
    # start = time()
    # img_FFBPmp = imgTools.FFBPmp(phs, platform, img_plane, taylor = 17, factor_max = 2)
    # fbpmp_time = time()-start

    #Output image
    u = img_plane['u']; v = img_plane['v']
    extent = [u.min(), u.max(), v.min(), v.max()]

    # plt.subplot(2,2,1)
    # plt.title('Full Backprojection \n \
    # Runtime = %i s'%bp_time)
    # imgTools.imshow(img_bp, dB_scale = [-30,0], extent = extent)
    # plt.xlabel('meters'); plt.ylabel('meters')

    # plt.subplot(2,2,2)
    # plt.title('Fast Factorized Backprojection \n w/o multi-processing \n \
    # Runtime = %i s'%fbp_time)
    # imgTools.imshow(img_FFBP, dB_scale = [-30,0], extent = extent)
    # plt.xlabel('meters'); plt.ylabel('meters')

    # plt.subplot(2,2,3)
    # plt.title('Fast Factorized Backprojection \n w/ multi-processing \n \
    # Runtime = %i s'%fbpmp_time)
    # imgTools.imshow(img_FFBPmp, dB_scale = [-30,0], extent = extent)
    # plt.xlabel('meters'); plt.ylabel('meters')
    # plt.tight_layout()

    # plt.subplot(2,2,4)
    plt.title('Fast Factorized Backprojection \n w/ MPI \n \
    Runtime = %i s'%fbpmpi_time)
    imgTools.imshow(img_FFBPmpi, dB_scale = [-30,0], extent = extent)
    plt.xlabel('meters'); plt.ylabel('meters')
    plt.tight_layout()
    if batch:
        plt.savefig("FFBPmpi_demo.png")
    else:
        plt.show()

    