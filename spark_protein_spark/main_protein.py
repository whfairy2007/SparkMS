# from pyspark import SparkContext, SparkConf

import os, glob
import numpy
from fnmatch import fnmatch

from main_matching import main_matching
from file_name_find import file_name_find
# import scipy



#####################################
# main program                      #
#####################################

if __name__ == "__main__":
    # conf = SparkConf().setAppName("PathwayEnrichment")
    # sc = SparkContext(conf=conf)

    # path of input folder
    # InputPath = "/user/spark/Ecolidata"
    InputPath = "C:\Xian Yang\Ecoli_data_half"
    # path of output folder
    # InputPath = "/user/spark/result"
    OutputPath = "C:\Xian Yang\\result_half"

    # MZ bin size
    mz_bin_size = 0.1
    # RT bin size
    rt_cut = 2
    # DT bin size
    dt_cut = 2.5

    # the whole dataset will be splitted into N pieces
    N = 10

    # if the number of ions in a cluster is larger than ionNumThresh then further split
    ionNumThresh = 10000


    file_format_c = '*3dion_BatchLowesFitTrim.csv'

    # find file names with file_format_c
    FileName,short_name = file_name_find(InputPath,file_format_c)

    # call main matching file
    main_matching(mz_bin_size,rt_cut,dt_cut,N,ionNumThresh,FileName,short_name,OutputPath)


    t=1





