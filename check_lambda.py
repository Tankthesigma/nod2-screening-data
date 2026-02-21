import numpy as np
sched = np.load('C:/Users/vasud/nod2-screening-data/fep_pmx/wt_complex/lambda_schedule.npy')
print('Window | Elec | Sterics | Restraints')
print('-------|------|---------|------------')
for i in [15,16,17,18,19]:
    print(f'  {i:02d}   | {sched[i,0]:.2f} |  {sched[i,1]:.2f}   |   {sched[i,2]:.2f}')
