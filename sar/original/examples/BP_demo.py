##############################################################################
#                                                                            #
#  This is a demonstration of the Backprojection algorithm.                  #
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

#AFRL DSBP demo
###############################################################################
#Define top level directory containing *.mat file
#and choose polarization and starting azimuth
pol = 'HH'
directory = './data/AFRL/pass1'
start_az = 1

#Import phase history and create platform dictionary
[phs, platform] = phsRead.AFRL(directory, pol, start_az, n_az = 4)

#Create image plane dictionary
img_plane = imgTools.img_plane_dict(platform, res_factor = 1.0, upsample = True, aspect = 1.0)

#full backprojection
start = time()
img_bp   = imgTools.backprojection(phs, platform, img_plane, taylor = 17, upsample = 2)
bp_time = time()-start

#Output image
u = img_plane['u']; v = img_plane['v']
extent = [u.min(), u.max(), v.min(), v.max()]

plt.title('Full Backprojection \n \
Runtime = %i s'%bp_time)
imgTools.imshow(img_bp, dB_scale = [-30,0], extent = extent)
plt.xlabel('meters'); plt.ylabel('meters')

if batch:
    plt.savefig("BP_demo.png")
else:
    plt.show()
