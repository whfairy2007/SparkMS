import numpy

from assignValue1 import assignValue1
from retentionCluster2 import retentionCluster2


def dataAlignmentIon4(ion_sel,nfile,mz_bound,nfeature,dt_cut,rt_cut,ionNumThresh,Dif_mz,Dif_dt,Dif_rt):
    alignData_i = []

    #  pay attention to empty arrays
    cellsz2_x = map(lambda x: len(x), ion_sel)
    # flag = sum(map(lambda x: x>0, cellsz2_x))

    if sum(cellsz2_x)>0:# there are ions in the bin
        flag1 = map(lambda x: x>=1, cellsz2_x)
        # index = [idx for idx in range(len(cellsz2_x)) if cellsz2_x[idx] > 0]
        flag2 = map(lambda x: x==1, cellsz2_x)

        if sum(flag1) == 1: # only one sample has one ion and other samples are empty
            if sum(flag2)==1:
                alignData_i = assignValue1(nfile,cellsz2_x,nfeature,ion_sel)
        else: #at least one file has more than one ion
            mz_size = mz_bound[1]-mz_bound[0]

            # if sum(cellsz2_x)<ionNumThresh% not implemented in this version
            alignData_i = retentionCluster2(ion_sel,nfeature,mz_size,dt_cut,rt_cut,cellsz2_x,Dif_mz,Dif_dt,Dif_rt)

    return alignData_i