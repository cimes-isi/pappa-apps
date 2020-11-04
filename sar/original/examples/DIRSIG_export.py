#Add include directories to default path list
import os
from sys import path
path.append('../')

#Include SARIT toolset
from ritsar import phsRead
from ritsar import phsTools
from ritsar import imgTools

import numpy as np

if __name__ == '__main__':
    #Import phase history and create platform dictionary
    #Define directory containing *.img files
    directory = './data/DIRSIG/'

    #Import phase history and create platform dictionary
    [phs, platform] = phsRead.DIRSIG(directory)

    #Correct for reisdual video phase
    phs_corr = phsTools.RVP_correct(phs, platform)

    #Demodulate phase history with constant reference, if needed
    phs_fixed = phsTools.phs_to_const_ref(phs_corr, platform, upchirp = 1)

    os.mkdir('npy')
    os.mkdir('npy_txt')

    print('--PHS--')
    print('phs: ' + str(type(phs_corr)))
    print('phs[0]: ' + str(type(phs_corr[0])))
    print('phs[0][0]: ' + str(type(phs_corr[0][0])))
    np.save('npy/phs.npy', phs_corr)
    np.savetxt('npy_txt/phs.npy', phs_corr, fmt='%0.5f')
    print('--PLATFORM--')
    for k in platform:
        print(str(k) + ': ' + str(type(platform[k])))
        print(str(k) + '[0]: ' + str(type(np.atleast_1d(platform[k])[0])))
        if str(k) == 'pos':
            print(str(k) + '[0][0]: ' + str(type(np.atleast_1d(platform[k])[0][0])))
        if str(k) == 'metadata':
            # These are plaintext values, but nothing we need to export
            for mk, mv in platform[k].items():
                print(str(mk) + ': ' + mv)
        else:
            np.save('npy/' + str(k) + '.npy', np.atleast_1d(platform[k]))
            np.savetxt('npy_txt/' + str(k) + '.npy', np.atleast_1d(platform[k]), fmt='%0.5f')

    os.mkdir('ip_npy')
    os.mkdir('ip_npy_txt')

    img_plane = imgTools.img_plane_dict(platform, res_factor = 1.0, aspect = 1.0)
    print('--IMG_PLANE--')
    for k in img_plane:
        print(str(k) + ': ' + str(type(img_plane[k])))
        print(str(k) + '[0]: ' + str(type(np.atleast_1d(img_plane[k])[0])))
        if str(k) == 'pixel_locs':
            print(str(k) + '[0][0]: ' + str(type(np.atleast_1d(img_plane[k])[0][0])))
        np.save('ip_npy/' + str(k) + '.npy', np.atleast_1d(img_plane[k]))
        np.savetxt('ip_npy_txt/' + str(k) + '.npy', np.atleast_1d(img_plane[k]), fmt='%0.5f')
