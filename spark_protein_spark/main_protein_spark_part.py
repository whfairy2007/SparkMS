from pyspark import SparkContext, SparkConf
import math
import numpy
from dataAlignment4_par import dataAlignment4_par
from dataAlignment4 import dataAlignment4
from pyspark.sql import SQLContext

def splitRow(row):
    newRow = row.split(' ')
    return newRow

def numberOfFile(idFile):
    nameSplit = idFile.split('New')
    name = nameSplit[1].split('.tx')
    nameInt = int(name[0])-1
    return nameInt

def floatify(array):
    newArray = [ float(x) for x in array ]
    return newArray

def reshapeFun(rowArray):
    nrow = rowArray[0]
    ncol = rowArray[1]
    vecArray = numpy.asarray(rowArray[2::])
    matArray = vecArray.reshape((nrow,ncol))

    return matArray

def buildMatrix(matrixString):
    matrixRows = matrixString.strip().split('\r\n')
    matrixArray = map(splitRow, matrixRows)
    dataArray = map(floatify, matrixArray)

    dataArrayReshape = map(reshapeFun, dataArray)

    return dataArrayReshape

def convertStr(array):
    res = []
    for i in range(len(array)):
        tmp  = [ str(x) for x in array[i]]
        res.append(','.join( tmp ))
    return res

def ion_organise(pre_slot,numFile,nfeature):
    pre_slot_cell = []
    for i in range(numFile):
        di = pre_slot[:,range(i*nfeature,(i+1)*nfeature)]
        # di_mz = di[:,1]
        pre_slot_cell.append(filter(lambda x: x[0]>0, di))#di(di_mz>0,:);

    return pre_slot_cell
#####################################
# main program                      #
#####################################

if __name__ == "__main__":
    conf = SparkConf().setAppName("PathwayEnrichment")
    sc = SparkContext(conf=conf)
    # path of input folder
    InputPath = "/user/ubuntu/"
    OutputPath1 = "/user/ubuntu/result1/"
    OutputPath2 = "/user/ubuntu/result2/"
    # MZ bin size
    mz_bin_size = 0.1*2# don'e foget to *2
    # RT bin size
    rt_cut = 2
    # DT bin size
    dt_cut = 2.5
    # the whole dataset will be splitted into N pieces
    N = 10
    # if the number of ions in a cluster is larger than ionNumThresh then further split
    ionNumThresh = 10000
    # number of repetitions on boundary
    N_dup = 10

    # results from the first part of program
    # mz_min = 248.7972
    # mz_max = 1689.3893
    # dt_min = 25.926
    # dt_max = 198.33
    # rt_min = 16.0627096276
    # rt_max = 40.070847675
    # nfeature = 21
    # numFile = 4
    # mz_bin_split = [[  248.6,   405.2], [  405.2,  462.8], [  462.8,   515.2], [  515.2,   561.2],
    #                  [  561.2,   606.2], [  606.2,   658.6], [  658.6,   722.4], [  722.4,   792. ], [  792.,    888.4], [  888.4,  1689.4]]
    sqlContext = SQLContext(sc)
    allVariable = sqlContext.read.json(InputPath+'/myVariable.json')
        # print allVariable.collect()
    allVariableCollected = allVariable.collect()
    allVariableDict = allVariableCollected[0].asDict()
    # filesName = [ str(x) for x in allVariableDict['header'][0] if x not in '']
    # headers = [ str(x) for x in allVariableDict['header'][1] if x not in '']
    mz_min = allVariableDict['mz_min']
    mz_max = allVariableDict['mz_max']
    dt_min = allVariableDict['dt_min']
    dt_max =  allVariableDict['dt_max']
    rt_min =  allVariableDict['rt_min']
    rt_max =  allVariableDict['rt_max']
    nfeature =  allVariableDict['nfeature']
    numFile =  allVariableDict['numFile']
    mz_bin_split =  allVariableDict['mz_bin_split']
    # print mz_min
    # print mz_bin_split


    ####
    Dif_mz = mz_max-mz_min
    Dif_dt = dt_max-dt_min
    Dif_rt = rt_max-rt_min

    # find m/z bin
    mz_temp = 1/mz_bin_size
    mz_start1 = math.floor(mz_min*mz_temp)/mz_temp
    mz_end1 = math.ceil(mz_max*mz_temp)/mz_temp

    mz_range1 = numpy.arange(mz_start1, mz_end1, mz_bin_size)
    mz_bin = numpy.vstack((mz_range1[:-1],mz_range1[1:])).T


    # prepare global variable
    parm1 = dict([('N_dup',N_dup),('mz_bin_split', mz_bin_split), ('mz_bin', mz_bin), ('nfeature', nfeature),('dt_cut', dt_cut),('rt_cut', rt_cut),('mz_bin_size', mz_bin_size),\
                 ('N', N),('numFile', numFile),('ionNumThresh', ionNumThresh),('Dif_mz', Dif_mz),('Dif_dt', Dif_dt),('Dif_rt', Dif_rt)])

    # read in and reorganize files
    FileName_data = sc.wholeTextFiles(InputPath+'*.txt')
    dataLow_all_T2 = FileName_data.map(lambda pair: (numberOfFile(pair[0]), buildMatrix(pair[1])))

    ######.....
    #
    # main matching part

    alignData_all3 = dataLow_all_T2.map(lambda listArray: dataAlignment4_par(listArray,parm1))
    print numpy.shape(alignData_all3.first()[0])



    # mat = alignData_all3.map(lambda tuple:[ str(x) for x in tuple[0]])
    #saveAsTextFile(OutputPath)
    mat = alignData_all3.map(lambda tuple: convertStr(tuple[0]))
    dataLow_all_T2.unpersist()
    print numpy.shape(mat.first())

    mat.saveAsTextFile(OutputPath1)


    alignData_R_all = alignData_all3.map(lambda tuple:tuple[1]).collect()
    alignData_L_all = alignData_all3.map(lambda tuple:tuple[2]).collect()
    alignData_all3.unpersist()

    dat12_int_all = []
    for n in range(N-1):
        pre_slot = alignData_R_all[n]# from the previous file
        next_slot = alignData_L_all[n+1] # from the next file
        pre_slot_cell = ion_organise(pre_slot,numFile,nfeature)
        next_slot_cell = ion_organise(next_slot,numFile,nfeature)

        dat12 = []
        for i in range(numFile):
            dat12.append(numpy.vstack((pre_slot_cell[i],next_slot_cell[i])))

        mz_low =  mz_bin_split[n][1]-mz_bin_size*(N_dup+1)
        mz_up = mz_bin_split[n][1]+mz_bin_size*(N_dup+1)
        # sel_bin = mz_bin[:,0] >= mz_low & mz_bin[:,1]<= mz_up
        mz_bin1 = filter(lambda x: x[0]>= mz_low and x[1]<= mz_up, mz_bin)

        dat12_int = dataAlignment4(mz_bin1,dat12,nfeature,dt_cut,rt_cut,mz_bin_size,numFile,ionNumThresh,Dif_mz,Dif_dt,Dif_rt)
        dat12_int_all.append(dat12_int)

    dat12_int_all_T = sc.parallelize(dat12_int_all)

    mat2 = dat12_int_all_T.map(lambda array: convertStr(array))
    mat2.saveAsTextFile(OutputPath2)