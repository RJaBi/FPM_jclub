from jclub import plaquette, plaquette_kolaorder, readopenqcd, readopenqcd_kolaorder, polyakov, polyakov_kolaorder

import numpy as np


def main():

    #gfFile = 'Gen2_8x24n7'
    gfFile = '/home/ryan/Documents/2024/conf/Gen2L/128x32/Gen2l_128x32n1'

    NS = 32
    NT = 128
    #"""
    U = readopenqcd(gfFile, NS, NS, NS, NT)  #nt,nx,ny,nz,mu,:,:
    print(f'U shape is {U.shape}')
    #sumTrP, nP, time = plaquette(U, 2, 4, 4)
    #print(f'U Plaquette for {gfFile} is {sumTrP/nP} and took {time} seconds')
    sumTrP, nP, time = polyakov(U)
    print(f'U Polyakov for {gfFile} is {sumTrP/nP} and took {time} seconds')
    sumTrP, nP, time = plaquette(U, 1, 4, 4)
    print(f'U Plaquette for {gfFile} is {sumTrP/nP} and took {time} seconds')
    #"""
    print('\n\n\n')
    UC = readopenqcd_kolaorder(gfFile, NS, NS, NS, NT)  #:,:,mu,nx,ny,nz,nt
    print(f'UC shape is {UC.shape}')
    sumTrP, nP, time = polyakov_kolaorder(UC)
    print(f'UC Polyakov for {gfFile} is {sumTrP/nP} and took {time} seconds')

    sumTrP, nP, time = plaquette_kolaorder(UC, 1, 4, 4)
    print(f'UC Plaquette for {gfFile} is {sumTrP/nP} and took {time} seconds')


if __name__ == '__main__':
    main()
