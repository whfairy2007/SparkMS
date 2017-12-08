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


def dataAlignment4_par(dataSelFeature_T,parm1):
    dataSelFeature = dataSelFeature_T[1]
    n = dataSelFeature_T[0]

    N_dup =  parm1["N_dup"]
    nfeature = parm1["nfeature"]
    dt_cut = parm1["dt_cut"]
    rt_cut = parm1["rt_cut"]
    mz_bin_size = parm1["mz_bin_size"]
    N = parm1["N"]
    # n = parm1["n"]
    nfile = parm1["numFile"]
    ionNumThresh = parm1["ionNumThresh"]
    Dif_mz = parm1["Dif_mz"]
    Dif_dt = parm1["Dif_dt"]
    Dif_rt = parm1["Dif_rt"]

    mz_bin_split = parm1["mz_bin_split"]
    mz_bin = parm1["mz_bin"]

    mz_bin1 = filter(lambda x: x[0]>= mz_bin_split[n][0]and x[1]<= mz_bin_split[n][1], mz_bin)

    #
    nbin = len(mz_bin1)
    BLOCK_SIZE = 2000 #initial capacity (& increment size)
    listSize = BLOCK_SIZE # current list capacity

    alignData = numpy.zeros((listSize, nfile*nfeature))                  # actual list
    listPtr = 0                                # pointer to last free position

    data_R_o = []
    data_L_o = []

    data_R = []
    #
    for iter in range(nbin):
        # print (n,iter)
        # if iter==277:
        #     t=1

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
        N_row_i = len(alignData_i)
        if N_row_i>0:
            if n>0 and n<N-1:
                if iter>=N_dup and iter < nbin-N_dup:
                    plc = listPtr + N_row_i # size of non-zeros entry
                    if plc > listSize:# add new block of memory if needed
                        listSize = listSize + BLOCK_SIZE#add new BLOCK_SIZE slots
                        # alignData(listPtr+1:listSize,:) = 0
                        alignData = numpy.vstack((alignData,numpy.zeros((BLOCK_SIZE,nfeature*nfile))))

                    alignData[range(listPtr,plc),:] = alignData_i # store new item
                    listPtr = plc                    # increment position pointer
                elif iter < N_dup:
                    if (len(data_L_o)>0 and len(alignData_i)>0):
                        data_L_o = numpy.vstack((data_L_o,alignData_i))
                    elif (len(data_L_o)==0 and len(alignData_i)>0):
                        data_L_o = alignData_i
                else:
                    if (len(data_R_o)>0 and len(alignData_i)>0):
                        data_R_o = numpy.vstack((data_R_o,alignData_i))
                    elif (len(data_R_o)==0 and len(alignData_i)>0):
                        data_R_o = alignData_i

            elif n==0:
                if iter < nbin-N_dup:
                    plc = listPtr + N_row_i # size of non-zeros entry
                    if plc > listSize:# add new block of memory if needed
                        listSize = listSize + BLOCK_SIZE#add new BLOCK_SIZE slots
                        # alignData(listPtr+1:listSize,:) = 0
                        alignData = numpy.vstack((alignData,numpy.zeros((BLOCK_SIZE,nfeature*nfile))))

                    alignData[range(listPtr,plc),:] = alignData_i # store new item
                    listPtr = plc                    # increment position pointer
                else:
                    if (len(data_R_o)>0 and len(alignData_i)>0):
                        data_R_o = numpy.vstack((data_R_o,alignData_i))
                    elif (len(data_R_o)==0 and len(alignData_i)>0):
                        data_R_o = alignData_i

            elif n==N-1:
                if iter>=N_dup:
                    plc = listPtr + N_row_i # size of non-zeros entry
                    if plc > listSize:# add new block of memory if needed
                        listSize = listSize + BLOCK_SIZE#add new BLOCK_SIZE slots
                        # alignData(listPtr+1:listSize,:) = 0
                        alignData = numpy.vstack((alignData,numpy.zeros((BLOCK_SIZE,nfeature*nfile))))

                    alignData[range(listPtr,plc),:] = alignData_i # store new item
                    listPtr = plc                    # increment position pointer

                else:
                    if (len(data_L_o)>0 and len(alignData_i)>0):
                        data_L_o = numpy.vstack((data_L_o,alignData_i))
                    elif (len(data_L_o)==0 and len(alignData_i)>0):
                        data_L_o = alignData_i

    alignData_o = numpy.asarray(filter(lambda x: sum(x)>0, alignData))
    if len(data_R)>0:
        if n==N-1:
            alignData_o = numpy.vstack((alignData_o,data_R))
        else:
            data_R_o = numpy.vstack((data_R_o,data_R))

    return(alignData_o, data_R_o,data_L_o)
