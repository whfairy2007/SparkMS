import numpy
import math


from dataAlignmentIon4 import dataAlignmentIon4
from integrateOverlap2 import integrateOverlap2

def select_Ion(dataSelFeature, mz_bound):
    ion_sel = []
    for item in dataSelFeature:
        item_s = filter(lambda x: x[0]> mz_bound[0] and x[0]<= mz_bound[1], item)
        ion_sel.append(item_s)
    return ion_sel

def dataAlignment4(mz_bin1,dataSelFeature,nfeature,dt_cut,rt_cut,mz_bin_size,nfile,ionNumThresh,Dif_mz,Dif_dt,Dif_rt):
    nbin = len(mz_bin1)
    alignData = []
    for iter in range(nbin):
        print iter
        mz_bound = mz_bin1[iter]
        # select ions in current bin
        ion_sel = select_Ion(dataSelFeature, mz_bound)

        alignData_bin = dataAlignmentIon4(ion_sel,nfile,mz_bound,nfeature,dt_cut,rt_cut,ionNumThresh,Dif_mz,Dif_dt,Dif_rt)
        mz_bound_L = mz_bound-mz_bin_size/2# may need to round

        # select ions in the left bin
        if iter == 0:
            ion_sel_L = select_Ion(dataSelFeature, mz_bound_L)
            alignData_L = dataAlignmentIon4(ion_sel_L,nfile,mz_bound_L,nfeature,dt_cut,rt_cut,ionNumThresh,Dif_mz,Dif_dt,Dif_rt)
        else:
            alignData_L = data_R

        # select ions in the right bin
        mz_bound_R = mz_bound+mz_bin_size/2
        ion_sel_R = select_Ion(dataSelFeature,mz_bound_R)
        alignData_R = dataAlignmentIon4(ion_sel_R,nfile,mz_bound_R,nfeature,dt_cut,rt_cut,ionNumThresh,Dif_mz,Dif_dt,Dif_rt)


        # do overlapping adjustment
        if len(alignData_L)==0:# if left one is empty, don't need to do overlapping adjustment
            alignData_i = []
            data_R = alignData_R
        elif len(alignData_L)>0 and len(alignData_bin)==0:# if the left one is not empty but the middle one is empty, don't need to do adjustment
            alignData_i = alignData_L
            data_R = alignData_R
        else:
            [alignData_i,data_R] = integrateOverlap2(0,alignData_bin,alignData_L,alignData_R,mz_bound_L,nfeature)

        # save results
        if len(alignData)>0 and len(alignData_i)>0:
            alignData = numpy.vstack((alignData,alignData_i))
        elif len(alignData)==0 and len(alignData_i)>0:
            alignData = alignData_i

    alignData_o = numpy.asarray(filter(lambda x: sum(x)>0, alignData))
    if len(data_R)>0:
        alignData_o = numpy.vstack((alignData_o,data_R))

    return alignData_o
