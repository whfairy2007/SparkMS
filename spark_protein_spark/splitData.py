from fnmatch import fnmatch
from file_name_find import file_name_find
import scipy.io as sio
import numpy

def load_mat_mz(file_name,N,mz_bin):
    mat_mz = []
    for fn in file_name:

        mat_i = sio.loadmat(fn)
        dataLowSel = mat_i['dataLowSel']
        mat_mz.extend(dataLowSel[:,0])
        # print len(dataLowSel[:,0])

    mat_mz_s = sorted(mat_mz)
    len_mz = len(mat_mz_s)
    step_len = int(round(float(len_mz)/N))
    bin_end = []
    for j in range(1,N):
        mz_temp = mat_mz_s[j*step_len]
        mz_bin_c2 = mz_bin[:,1]# second column
        index = min(range(len(mz_bin_c2)), key=lambda i: abs(mz_bin_c2[i]-mz_temp))

        bin_end.append(mz_bin_c2[index])

    bin_end.append(mz_bin_c2[-1])

    bin_start = bin_end[:-1]
    bin_start.insert(0,mz_bin[0,0])
    mz_bin2 = numpy.vstack((numpy.asarray(bin_start),numpy.asarray(bin_end))).T

    return mz_bin2

def reOrgMat(dataLow):
    dataLowReOrg = []

    for i in range(len(dataLow)):
        dataLowEach = numpy.asarray(dataLow[i])
        nrow = len(dataLowEach)
        ncol = len(dataLowEach[0])
        dataLowEachVec = dataLowEach.ravel()
        dataLowEachVec1 = numpy.append([nrow,ncol],dataLowEachVec)

        dataLowReOrg.append(dataLowEachVec1)

    return dataLowReOrg




def splitData(pathSel,pathMat,mz_bin,N,numFile):
    # list all mat file
    FileList,FileList_short = file_name_find(pathSel,'*.mat')

    # use all files to split
    mz_bin2 = load_mat_mz(FileList,N,mz_bin)# find the splitting boundary

    # file_split_size = [[0 for j in range(numFile)] for i in range(N)]# initialize array
    file_split_size = numpy.zeros((numFile,N))
    dataLow_all = []
    for n in range(N):
        dataLow =[]
        mz_min = mz_bin2[n,0]
        mz_max = mz_bin2[n,1]

        for i in range(numFile):
            fn = FileList[i]
            mat_i = sio.loadmat(fn)
            dataLowSel = mat_i['dataLowSel']
            mat_mz = dataLowSel[:,0]
            dataLowSel_s = filter(lambda x: x[0] >= mz_min and x[0]<= mz_max, dataLowSel)
            dataLow.append(dataLowSel_s)
            file_split_size[i,n] = len(dataLowSel_s)

        n_1 = n+1
        if n_1 < 10:
            name_num = '000'+str(n_1)
        elif n_1 >=10 and n_1 < 100:
            name_num = '00'+str(n_1)
        elif n_1 >=100 and n_1 < 1000:
            name_num = '0'+str(n_1)
        else:
            name_num = str(n_1)

        file_name = pathMat+'\\'+'dataLowNew'+name_num+'.txt'

        # sio.savemat(file_name, {'dataLow':dataLow})
        # dataLow_all.append(dataLow)
        dataLowReOrg = reOrgMat(dataLow)
        dataLowReOrg2 = numpy.asarray(dataLowReOrg)
        with open(file_name,'a') as f_handle:
                f_handle.write("\n".join(" ".join(map(str, x)) for x in dataLowReOrg2))
        f_handle.close()

    # return(dataLow_all,file_split_size,mz_bin2)

    return(file_split_size,mz_bin2)