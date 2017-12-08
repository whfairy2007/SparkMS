import os
import math
import numpy
import scipy.io as sio

import csv

from readInDRAMI import readInDRAMI
from splitData import splitData
from createHeaderOutput import createHeaderOutput
from dataAlignment4_par import dataAlignment4_par
from dataAlignment4 import dataAlignment4

def ion_organise(pre_slot,numFile,nfeature):
    pre_slot_cell = []
    for i in range(numFile):
        di = pre_slot[:,range(i*nfeature,(i+1)*nfeature)]
        # di_mz = di[:,1]
        pre_slot_cell.append(filter(lambda x: x[0]>0, di))#di(di_mz>0,:);

    return pre_slot_cell

def submain_matching_par(N,FileName,short_name,OutputPath,mz_bin_size1,rt_cut,dt_cut,listVarMatch,SelShowVar,Align_RT_name,selectIO,ionNumThresh):

    #attention ... window size is 2 times bin size
    mz_bin_size = mz_bin_size1*2

    #create path to save results
    path_Ord = os.path.join(OutputPath, 'MZ_ordered_file')
    if not os.path.exists(path_Ord):
        os.makedirs(path_Ord)

    path_Mat = os.path.join(OutputPath, 'Splitted_file')
    if not os.path.exists(path_Mat):
        os.makedirs(path_Mat)

    path_Out = os.path.join(OutputPath, 'Matching_result')
    if not os.path.exists(path_Out):
        os.makedirs(path_Out)

    #
    numFile = len(FileName)
    nfeature = len(listVarMatch)+len(SelShowVar)# total number of selected features

    #step1: read in dataset & select low energy part & order each file by MZ

    mz_min,mz_max, dt_min,dt_max, rt_min,rt_max = readInDRAMI(FileName,short_name,path_Ord,listVarMatch,SelShowVar,Align_RT_name,numFile)
    Dif_mz = mz_max-mz_min
    Dif_dt = dt_max-dt_min
    Dif_rt = rt_max-rt_min

    # step2: split dataset
    # find m/z bin
    mz_temp = 1/mz_bin_size
    mz_start1 = math.floor(mz_min*mz_temp)/mz_temp
    mz_end1 = math.ceil(mz_max*mz_temp)/mz_temp

    mz_range1 = numpy.arange(mz_start1, mz_end1, mz_bin_size)
    mz_bin = numpy.vstack((mz_range1[:-1],mz_range1[1:])).T
    # split data
    # dataLow_all,file_split_size,mz_bin_split = splitData(path_Ord,path_Mat,mz_bin,N,numFile)
    file_split_size,mz_bin_split = splitData(path_Ord,path_Mat,mz_bin,N,numFile)


    # create the header of output file
    listVarMatch1 = listVarMatch
    listVarMatch1[2] = Align_RT_name
    header = createHeaderOutput(FileName,short_name,numFile,nfeature ,listVarMatch1,SelShowVar)

    # step 3 matching in parallel

    pre_R = []
    alignData_R_all = []
    alignData_L_all = []
    N_dup = 10


    for n in range(N):

        n_1 = n+1
        if n_1 < 10:
            name_num = '000'+str(n_1)
        elif n_1 >=10 and n_1 < 100:
            name_num = '00'+str(n_1)
        elif n_1 >=100 and n_1 < 1000:
            name_num = '0'+str(n_1)
        else:
            name_num = str(n_1)

        file_name = path_Mat+'\\'+'dataLowNew'+name_num+'.mat'

        dataLow_dic = sio.loadmat(file_name)
        dataLow1 = dataLow_dic["dataLow"]
        dataLow = dataLow1[0]

        parm1 = dict([('N_dup',N_dup),('mz_bin_split', mz_bin_split), ('mz_bin', mz_bin), ('nfeature', nfeature),('dt_cut', dt_cut),('rt_cut', rt_cut),('mz_bin_size', mz_bin_size),\
                 ('N', N),('n', n),('numFile', numFile),('ionNumThresh', ionNumThresh),('Dif_mz', Dif_mz),('Dif_dt', Dif_dt),('Dif_rt', Dif_rt)])

        alignData,alignData_R,alignData_L = dataAlignment4_par(dataLow,parm1)

        output_filename = path_Out+'\\'+'dataLowOut'+name_num+'.csv'
        numpy.savetxt(output_filename, alignData, delimiter=",")
        alignData_R_all.append(alignData_R)
        alignData_L_all.append(alignData_L)
        # numpy.savetxt(output_filename, numpy.vstack((header, alignData)), delimiter=",")



    # import scipy.io as sio
    # sio.savemat('temp.mat', locals())

    # sio.loadmat('temp.mat')
    for n in range(N-1):
        pre_slot = alignData_R_all[n]# from the previous file
        next_slot = alignData_L_all[n+1] # from the next file
        pre_slot_cell = ion_organise(pre_slot,numFile,nfeature)
        next_slot_cell = ion_organise(next_slot,numFile,nfeature)

        dat12 = []
        for i in range(numFile):
            dat12.append(numpy.vstack((pre_slot_cell[i],next_slot_cell[i])))

        mz_low =  mz_bin_split[n,1]-mz_bin_size*(N_dup+1)
        mz_up = mz_bin_split[n,1]+mz_bin_size*(N_dup+1)
        # sel_bin = mz_bin[:,0] >= mz_low & mz_bin[:,1]<= mz_up
        mz_bin1 = filter(lambda x: x[0]>= mz_low and x[1]<= mz_up, mz_bin)

        dat12_int = dataAlignment4(mz_bin1,dat12,nfeature,dt_cut,rt_cut,mz_bin_size,numFile,ionNumThresh,Dif_mz,Dif_dt,Dif_rt)

        n_1 = n+1
        if n_1 < 10:
            name_num = '000'+str(n_1)
        elif n_1 >=10 and n_1 < 100:
            name_num = '00'+str(n_1)
        elif n_1 >=100 and n_1 < 1000:
            name_num = '0'+str(n_1)
        else:
            name_num = str(n_1)
        output_filename = path_Out+'\\'+'dataLowOut'+name_num+'.csv'
        # numpy.savetxt(output_filename, alignData, delimiter=",")
        with open(output_filename,'a') as f_handle:
            numpy.savetxt(f_handle,dat12_int, delimiter=",")

