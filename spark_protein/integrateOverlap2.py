import numpy

def cmp_3d_dif(data1_L_mz_ij,mz_ij,data1_L_dt_ij,dt_ij,data1_L_rt_ij,rt_ij,j):# this part is an improvement on Matlab code
    dif_mz_L = max(abs(numpy.delete(data1_L_mz_ij,j)-data1_L_mz_ij[0][j]))
    dif_rt_L = max(abs(numpy.delete(data1_L_rt_ij,j)-data1_L_rt_ij[0][j]))
    dif_dt_L = max(abs(numpy.delete(data1_L_dt_ij,j)-data1_L_dt_ij[0][j]))

    dif_mz = max(abs(numpy.delete(mz_ij,j)-mz_ij[j]))
    dif_rt = max(abs(numpy.delete(rt_ij,j)-rt_ij[j]))
    dif_dt = max(abs(numpy.delete(dt_ij,j)-dt_ij[j]))

    f_mz = 1 if dif_mz_L<dif_mz else 0
    f_rt = 1 if dif_rt_L<dif_rt else 0
    f_dt = 1 if dif_dt_L<dif_dt else 0

    if (f_mz+f_dt+f_rt)>1:
        return 1
    else:
        return 0


def integrateOverlap2(ind,alignData_bin,alignData_L,alignData_R,mz_bound_L,nfeature):
    # ind,data2_sel,data1_sel_L,data1_sel_R,bin1_L,nfeature)
    nl_2d = numpy.shape(alignData_bin)# number of columns in current bin
    nl = nl_2d[1]
    # find the range of mz in each row
    data2_sel_mz = alignData_bin[:,range(ind,nl,nfeature)]
    data2_sel_mz_max = numpy.amax(data2_sel_mz,axis=1)
    data2_sel_mz[numpy.where(data2_sel_mz==0)] = 10000
    data2_sel_mz_min = numpy.amin(data2_sel_mz,axis=1)

    #the boundary between left and right bins
    mz_thresh = mz_bound_L[1]
    # find rows in current bin that cross boundary
    sel_in1 = data2_sel_mz_min<=mz_thresh
    sel_in2 = data2_sel_mz_max>=mz_thresh
    sel_in12 = numpy.where(sel_in1 & sel_in2)
    sel_in = sel_in12[0]

    if len(sel_in)==0: # empty
        alignData_i = alignData_L
        data_R = alignData_R
    else:
        data2_rep = alignData_bin[sel_in,:]

        data2_rep_mz = data2_rep[:, range(0,nl,nfeature)]
        data2_rep_dt = data2_rep[:, range(1,nl,nfeature)]
        data2_rep_rt = data2_rep[:, range(2,nl,nfeature)]

        data1_sel_L_mz = alignData_L[:, range(0,nl,nfeature)]
        data1_sel_L_dt = alignData_L[:, range(1,nl,nfeature)]
        data1_sel_L_rt = alignData_L[:, range(2,nl,nfeature)]

        data1_sel_R_mz = alignData_R[:, range(0,nl,nfeature)]
        data1_sel_R_dt = alignData_R[:, range(1,nl,nfeature)]
        data1_sel_R_rt = alignData_R[:, range(2,nl,nfeature)]

        yl_2d = numpy.shape(data2_rep_mz)
        yl = yl_2d[0]

        for i in range(yl):
            mz_ij = data2_rep_mz[i,:]
            dt_ij = data2_rep_dt[i,:]
            rt_ij = data2_rep_rt[i,:]

            sel0_0 = numpy.where(mz_ij>0)#non-zero entry
            sel0 = sel0_0[0]
            num_ion = len(sel0)#number of non-zero ions in the row
            f0 = numpy.ones(num_ion)*num_ion

            fLR = numpy.zeros(num_ion)
            fLR_3D_cmp =  numpy.zeros(num_ion)
            for k in range(num_ion): # for each non-zero entry in current row
                j = sel0[k]
                mz_k = mz_ij[j]
                dt_k = dt_ij[j]
                rt_k = rt_ij[j]
                if ind==0:
                    val_k = mz_k
                elif ind==1:
                    val_k = dt_k
                elif ind==2:
                    val_k = rt_k

                if val_k<=mz_thresh:#left
                    tf1 = data1_sel_L_mz[:,j]==mz_k
                    tf2 = data1_sel_L_dt[:,j]==dt_k
                    tf3 = data1_sel_L_rt[:,j]==rt_k

                    flag_L_0 = numpy.where(tf1&tf2&tf3)
                    flag_L = flag_L_0[0]
                    data1_sel_L_mz_s = data1_sel_L_mz[flag_L,:]>0
                    fLR[k] =  sum(data1_sel_L_mz_s[0]) #number of non-zero ions in the row containing this ion in the left bin

                    data1_L_mz_ij = data1_sel_L_mz[flag_L]
                    data1_L_dt_ij = data1_sel_L_dt[flag_L]
                    data1_L_rt_ij = data1_sel_L_rt[flag_L]

                    fLR_3D_cmp[k] = cmp_3d_dif(data1_L_mz_ij,mz_ij,data1_L_dt_ij,dt_ij,data1_L_rt_ij,rt_ij,j)

                else:#right
                    tf1 = data1_sel_R_mz[:,j]==mz_k
                    tf2 = data1_sel_R_dt[:,j]==dt_k
                    tf3 = data1_sel_R_rt[:,j]==rt_k
                    flag_R_0 = numpy.where(tf1&tf2&tf3)
                    flag_R = flag_R_0[0]
                    data1_sel_R_mz_s = data1_sel_R_mz[flag_R,:]>0
                    fLR[k] =  sum(data1_sel_R_mz_s[0])# number of non-zero ions in the row containing this ion in the RIGHT bin

            num_ion_ini = num_ion
            flag_LR_0_le = 0
            while 1:
                # flag_LR_0_0 = numpy.where(fLR>f0)
                # flag_LR_0 = flag_LR_0_0[0]
                flag_LR_0 = [i_f for i_f in range(num_ion_ini) if (fLR[i_f]>f0[i_f] or (fLR[i_f]==f0[i_f] and fLR_3D_cmp[i_f]==1))]

                if len(flag_LR_0)>0:
                    if len(flag_LR_0)>flag_LR_0_le:
                        num_ion = num_ion_ini-len(flag_LR_0)
                        flag_LR_0_le = len(flag_LR_0)
                        f0 = numpy.ones(num_ion_ini)*num_ion
                    else:
                        break

                else:
                    break

            # Delete ions
            # delete ions in data2_rep with sel0(flag_LR_0(p))
            for p in range(len(flag_LR_0)):
                delete_index = sel0[flag_LR_0[p]]
                data2_rep[i,range(delete_index*nfeature,(delete_index+1)*nfeature)] = numpy.zeros(nfeature)# attention some zero rows are created


            # delete ions in left or right bin
            flag_LR_LR_0 = numpy.where(fLR<=f0)
            flag_LR_LR = flag_LR_LR_0[0]
            if len(flag_LR_LR)>0:
                for p in range(len(flag_LR_LR)):
                    pj = sel0[flag_LR_LR[p]]
                    mz_pj = mz_ij[pj]
                    rt_pj = rt_ij[pj]
                    dt_pj = dt_ij[pj]
                    if ind==0:
                        val_pj = mz_pj
                    elif ind==1:
                        val_pj = dt_pj
                    elif ind==2:
                        val_pj = rt_pj

                    if val_pj<=mz_thresh:# left
                        tf1_2 = data1_sel_L_mz[:,pj]==mz_pj
                        tf2_2 = data1_sel_L_dt[:,pj]==dt_pj
                        tf3_2 = data1_sel_L_rt[:,pj]==rt_pj
                        flag_L_0_2 = numpy.where(tf1_2&tf2_2&tf3_2)
                        flag_L_2 = flag_L_0_2[0]
                        alignData_L[flag_L_2,range(pj*nfeature,(pj+1)*nfeature)] = numpy.zeros(nfeature)
                    else:
                        tf1_2 = data1_sel_R_mz[:,pj]==mz_pj
                        tf2_2 = data1_sel_R_dt[:,pj]==dt_pj
                        tf3_2 = data1_sel_R_rt[:,pj]==rt_pj
                        flag_R_0_2 = numpy.where(tf1_2&tf2_2&tf3_2)
                        flag_R_2 = flag_R_0_2[0]
                        alignData_R[flag_R_2,range(pj*nfeature,(pj+1)*nfeature)] = numpy.zeros(nfeature)



                t=1

        data2_rep_f = numpy.asarray(filter(lambda x: sum(x)>0, data2_rep))
        alignData_L_f = numpy.asarray(filter(lambda x: sum(x)>0, alignData_L))
        if len(data2_rep_f)>0 and len(alignData_L_f)>0:
            alignData_i = numpy.vstack((alignData_L_f,data2_rep_f))
        elif len(data2_rep_f)==0 and len(alignData_L_f)>0:
            alignData_i = alignData_L_f
        else:
            alignData_i = data2_rep_f

        data_R = numpy.asarray(filter(lambda x: sum(x)>0, alignData_R))

    return (alignData_i,data_R)