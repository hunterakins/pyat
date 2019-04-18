import numpy as np
import sys
sys.path.append("../")
from write_fieldflp import write_fieldflp
from equally_spaced import equally_spaced
from read_shd import read_shd
from os import system
from matplotlib import pyplot as plt

class Empty:
    def __init__(self): 
        pass
Pos = Empty()
Pos.r = Empty()
Pos.s = Empty()
Pos.r.range =  np.arange(0,10, 10/1e3)
rd = [94.125, 99.755, 105.38, 111.00, 116.62, 122.25, 127.88, 139.12, 144.74, 150.38, 155.99, 161.62, 167.26, 172.88, 178.49, 184.12, 189.76, 195.38, 200.99, 206.62, 212.25]
Pos.r.depth = np.arange(0,150,1)
Pos.s.depth = [80,100]
write_fieldflp('test_file', 'R', Pos)
system("/home/hunter/Downloads/at/bin/kraken.exe test_file")
system("/home/hunter/Downloads/at/bin/field.exe test_file")
[x,x,x,x,Pos1,pressure]= read_shd('test_file.shd',166)
pressure = np.squeeze(pressure)
pressure = abs(pressure)
pressure = 10*np.log10(pressure / np.max(pressure))
print(np.shape(pressure))
plt.contourf(Pos1.r.range, Pos1.r.depth,pressure[1,:,:])
plt.gca().invert_yaxis()
plt.show()

