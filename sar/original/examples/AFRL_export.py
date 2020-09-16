#Add include directories to default path list
import os
import sys
sys.path.append('../')
sys.path.append('./dictionaries')

#Include Dictionaries
from SARplatform import plat_dict

#Include SARIT toolset
from ritsar import phsRead

import numpy as np

if __name__ == '__main__':
    #AFRL DSBP demo
    ###############################################################################
    #Define top level directory containing *.mat file
    #and choose polarization and starting azimuth
    pol = 'HH'
    directory = './data/AFRL/pass1'
    start_az = 1

    os.mkdir('HH_npy')
    os.mkdir('HH_npy_txt')

    #Import phase history and create platform dictionary
    [phs, platform] = phsRead.AFRL(directory, pol, start_az, n_az = 4)
    print('phs: ' + str(type(phs)))
    print('phs[0]: ' + str(type(phs[0])))
    print('phs[0][0]: ' + str(type(phs[0][0])))
    np.save('HH_npy/phs.npy', phs)
    np.savetxt('HH_npy_txt/phs.npy', phs, fmt='%0.5f')
    for k in platform:
        print(str(k) + ': ' + str(type(platform[k])))
        print(str(k) + '[0]: ' + str(type(np.atleast_1d(platform[k])[0])))
        if str(k) == 'pos':
            print(str(k) + '[0][0]: ' + str(type(np.atleast_1d(platform[k])[0][0])))
        np.save('HH_npy/' + str(k) + '.npy', np.atleast_1d(platform[k]))
        np.savetxt('HH_npy_txt/' + str(k) + '.npy', np.atleast_1d(platform[k]), fmt='%0.5f')
