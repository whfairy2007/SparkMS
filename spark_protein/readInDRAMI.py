import numpy
import pandas as pd
# from scipy.io import savemat
import scipy.io as sio


def return_indices_of_a(a, b):
    list_index = []
    for item in a:
        # ind = numpy.where(b==item)
        ind = b.tolist().index(item)
        list_index.append(ind)
    return list_index

def order_file_mz(mat,path_Ord,filename,nfile,mz_min,mz_max, dt_min,dt_max, rt_min,rt_max):

    # mat_mz = mat.iloc[:,0]
    mat_mz = mat[:,0]
    mz_min = min(min(mat_mz),mz_min)
    mz_max = max(max(mat_mz),mz_max)
    mat_dt = mat[:,1]
    dt_min = min(min(mat_dt),dt_min)
    dt_max = max(max(mat_dt),dt_max)
    mat_rt = mat[:,2]
    rt_min = min(min(mat_rt),rt_min)
    rt_max = max(max(mat_rt),rt_max)

    dataLowSel = mat[mat[:,0].argsort()]
    index = nfile+1
    if index <10:
        name_num = '000'+str(index)
    elif index >=10 and index < 100:
        name_num = '00'+str(index)
    elif index >=100 and index < 1000:
        name_num = '0'+str(index)
    else:
        name_num = str(index)

    FileName_ord = filename.replace('.csv','')
    file_name = path_Ord+'\\'+name_num+FileName_ord+'.mat'

    # attention: dataLowSel = unique(dataLowSel,'rows');
    sio.savemat(file_name, {'dataLowSel':dataLowSel})

    return(mz_min,mz_max, dt_min,dt_max, rt_min,rt_max)


def readInDRAMI(FileName,short_name,path_Ord,listVarMatch,SelShowVar,Align_RT_name,nfile):

    mz_min = 1000
    mz_max = 0
    dt_min = 1000
    dt_max = 0
    rt_min = 1000
    rt_max = 0

    for n in range(nfile):
        s1 = short_name[n]
        FileName_ord = s1.replace('.csv','')
        # my_data = numpy.recfromcsv(FileName[n], delimiter=',', filling_values=numpy.nan, case_sensitive=True, deletechars='', replace_space=' ')
        A = pd.read_csv(FileName[n],header=0)


        # # save the low energy part
        # sel_r = A['Function'] == 1
        # A3 = A[sel_r]
        # output_filename0 = FileName[n]
        # output_filename = output_filename0.replace('_Pep3DAMRT.csv','_Pep3DAMRT_f1.csv')


        #  check whether this is aligned IO
        A_column = A.columns.values
        if sum(A_column == Align_RT_name):
            listVarMatch1 = listVarMatch
            listVarMatch1[2] = Align_RT_name
            var_match_ind = return_indices_of_a(listVarMatch1,A_column)

        var_other_ind = return_indices_of_a(SelShowVar,A_column)
        sel_ind = var_match_ind + var_other_ind

        # order file by MZ
        A4 = pd.np.array(A.iloc[:,sel_ind])#change data type
        print n
        mz_min,mz_max, dt_min,dt_max, rt_min,rt_max = order_file_mz(A4,path_Ord,s1,n,mz_min,mz_max, dt_min,dt_max, rt_min,rt_max)

    return (mz_min,mz_max, dt_min,dt_max, rt_min,rt_max)
