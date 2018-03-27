#========================================================================================
#   LIBRARIES
#========================================================================================
import numpy as np
import matplotlib.pyplot as plt
import BH_Physics as BHP
import sys
plt.close('all')


#========================================================================================
#   LIBRARIES
#========================================================================================
#Data folder
DataFolder = '/home/bustamsn/PhD/Data/cosmological_BH/Sims256/'
#DataFolder = '/home/bustamsn/bustamsn/cosmological_BH/'
DataResults = '.'
Simulation = 'cosmobh05'
#Values of anisotropy to compute
K_values = [0,1,2,10]

#========================================================================================
#   LIBRARIES
#========================================================================================
BH = BHP.black_hole_sim( simulation = Simulation, snapbase = 'snapshot', n_snap = 1000, datafolder = DataFolder, resultsfolder = DataResults )
indexes = np.loadtxt('%s%s/analysis/BH_IDs.txt'%(DataFolder,Simulation))
#Creating folder to store analytical profiles
os.system('mkdir %s%s/analysis/analytical'%(DataFolder,Simulation))

mask_chaotic = indexes[:,-1] == 1
N_files = 1.0*np.sum(mask_chaotic)
#Argument 1 in console: number of current chunk
chunk = int(sys.argv[1])
#Argument 2 in console: total number of chunks in which data is divided (parallel computing)
total_chunks = int(sys.argv[2])

N0 = int(N_files*(chunk-1)/total_chunks)
N1 = int(N_files*(chunk)/total_chunks)
if chunk == total_chunks:
    N1 = N_files

for i in indexes[mask_chaotic][N0:N1,0].astype(int):
    print chunk, i
    ndata = BH.loading_spin( i )
    for k in K_values:
        if ndata>10:
            BH.spin_evolution( mode='auto', rad_eff = -1, alpha = 0.1, kparam = k, a0 = 0.05, Nmax = 5e4 )
            np.savetxt( '%s%s/analysis/analytical/spin_%d_k%d.txt'%(DataFolder,Simulation,i,k), np.array( [BH.tms, BH.aps, BH.Mbh(BH.tms), BH.modeps] ).T, fmt='%e %e %e %d' )
