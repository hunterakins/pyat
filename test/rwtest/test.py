import numpy as np
from matplotlib import pyplot as plt
from pyat.pyat.readwrite import *

'''
Description:
Test read and write functionality of pyat
Should add more tests but things are working alright

Author: Hunter Akins

'''

if __name__ == '__main__':
    print('Running test for env files')
    [title, freq, ssp, bdry, pos, beam, cint, RMax] = read_env('../krak_test/py_env.env', 'KRAKEN')
    a = np.linspace(0, 150, 100)
    plt.plot(a,ssp.sspf(a))
    plt.show()


   # [title, freq, ssp, bdry, pos, beam, cInt, RMax] = read_env('test/bellhop_test/py_env.env', 'BELLHOP')
   # a = np.linspace(0, 5000, 100)
   # plt.plot(a,ssp.sspf(a))
   # plt.show()

    [title, freq, ssp, bdry, pos, beam, cint, RMax] = read_env('../bellhop_test/swellex.env', 'KRAKEN')
    a = np.linspace(0, 1040, 100)
    plt.plot(a, ssp.sspf(a))
    plt.show()

    write_env('test_file', 'KRAKEN', title, freq, ssp, bdry, pos, beam, cint, RMax)
    
