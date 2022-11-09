import MDAnalysis as mda
import os
from MDAnalysis.tests.datafiles import PSF, DCD
from MDAnalysis.analysis import encore
from MDAnalysis.analysis.encore.clustering import ClusteringMethod as clm
from MDAnalysis.analysis.encore.dimensionality_reduction import DimensionalityReductionMethod as drm
from MDAnalysis.coordinates.memory import MemoryReader
import numpy as np
import matplotlib.pyplot as plt
import getopt
import sys

args=sys.argv[1:]                                                                   

ref_structure = ''

ref_traj =''

num_windows = '10'

output = 'ces_conv'

try:

    opts, args = getopt.getopt(args,"h:s:t:n:o:",["help","ref_structure =",          
                                    "ref_traj =",                                   
                                    "num_windows =",
                                    "output ="])

except getopt.GetoptError:

   print ('analysis.py -s <md.gro> -t <md.xtc> -n <number_of_windows> -o <output> ')

   sys.exit(2)                                                                      



for opt, arg in opts:                                                               

   if opt == '-h':

      print ('analysis.py -s <md.gro> -t <md.xtc> -n <number_of_windows> -o <output> ')

      sys.exit()

   elif opt in ("-s", "--ref_structure"):

      ref_structure = str(arg)

   elif opt in ("-t", "--ref_traj"):

      ref_traj = str(arg)

   elif opt in ("-n", "--num_windows"):

      num_windows = int(arg)

   elif opt in ("-o", "--output"):

      output = str(arg)

u = mda.Universe(ref_structure,ref_traj, in_memory=True)
# u = mda.Universe(ref_structure,ref_traj)

coordinates = u.trajectory.timeseries(u.atoms)
u2 = mda.Universe(ref_structure, coordinates, format=MemoryReader, order='afc')
ces_conv = encore.ces_convergence(u2,
                                  num_windows,  # window size
                                  select='name CA')  # default

count = 0
with open(output, 'w') as f:
    for i in range(len(ces_conv)):
        f.write(str(count) + "  " + str(ces_conv[i]))
        f.write('\n')
        count+=1
os.system("sed -i 's/\[/ /g' *dat")
os.system("sed -i 's/\]/ /g' *dat")
print("Cool! It's completed!")

