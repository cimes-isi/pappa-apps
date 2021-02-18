##############################################################################
#                                                                            #
#  This is a demonstration of the ritsar toolset using Sandia data.          #
#  Algorithms can be switched in and out by commenting/uncommenting          #
#  the lines of code below.                                                  #
#                                                                            #
##############################################################################

#Add include directories to default path list
import sys
sys.path.append('../')
if len(sys.argv) > 1 and sys.argv[1] == '--batch':
    batch = True
else:
    batch = False

import matplotlib.pylab as plt

#Include SARIT toolset
from ritsar import phsRead
from ritsar import phsTools
from ritsar import imgTools

#Define directory containing *.au2 and *.phs files
directory = './data/Sandia/'

#Import phase history and create platform dictionary
[phs, platform] = phsRead.Sandia(directory)

#Correct for residual video phase
phs_corr = phsTools.RVP_correct(phs, platform)

#Import image plane dictionary from './parameters/img_plane'
img_plane = imgTools.img_plane_dict(platform,
                           res_factor = 1.0, n_hat = platform['n_hat'])

#Apply polar format algorithm to phase history data
#(Other options not available since platform position is unknown)
#img_pf = imgTools.polar_format(phs_corr, platform, img_plane, taylor = 30)
img_bp = imgTools.backprojection(phs_corr, platform, img_plane, taylor = 30, upsample=2)

#Output image
imgTools.imshow(img_bp, [-45,0])
if batch:
    plt.savefig("Sandia_demo.png")
else:
    plt.show()
