import numpy
import math
from scipy.spatial.distance import pdist
from scipy.cluster.hierarchy import linkage
from scipy.cluster.hierarchy import fcluster

def cell2mat(ion_sel):
    frame_rt = []
    for item in ion_sel:
        frame_rt.extend(item)# check extend or append

    return numpy.asarray(frame_rt)


def distfun(x,y,mz_size,dt_cut,rt_cut):
    dis_xy = (x[0]-y[0])**2+(x[1]-y[1])**2+(x[2]-y[2])**2
    flag1 = x[3]==y[3] or  math.fabs(x[4]-y[4])>mz_size/2 or math.fabs(x[5]-y[5])>dt_cut/2 or math.fabs(x[6]-y[6])>rt_cut/2
    if flag1:
        dis_xy = dis_xy+10**6

    return dis_xy


def retentionCluster2(ion_sel,nfeature,mz_size,dt_cut,rt_cut,cellsz2_x,Dif_mz,Dif_dt,Dif_rt):
    # alignData_i_j = retentionCluster2(x1,nfeature,mz_size,cutoff_th,rt_cut,cellsz,Dif_mz,Dif_dt,Dif_rt)
    # number of files
    nfile = len(ion_sel)
    frame_rt = cell2mat(ion_sel)


    # scale to RT
    frame_rt1 = frame_rt[:,range(3)]
    frame_rt2 = numpy.zeros(frame_rt1.shape)
    frame_rt2[:,0] = frame_rt1[:,0]*(Dif_rt/Dif_mz)
    frame_rt2[:,1] = frame_rt1[:,1]*(Dif_rt/Dif_dt)
    frame_rt2[:,2] = frame_rt1[:,2]

    #ions from the sample will be assigned with large distance

    frame_rt3 = []
    for k in range(nfile):
        if cellsz2_x[k]>0:
            frame_rt3.extend(k*numpy.ones(cellsz2_x[k]))
        # if cellsz2_x[k]>1:# test
        #     t=1

    frame_rt3_a = numpy.asarray(frame_rt3)
    frame_rt3_at= frame_rt3_a.T
    # frame_rt_all = numpy.concatenate((frame_rt2, frame_rt3_at,frame_rt1), axis=1)# may have error message here
    frame_rt_all = numpy.column_stack((frame_rt2,frame_rt3,frame_rt1))# may have error message here


    # define distance
    dis_rt = pdist(frame_rt_all,lambda x,y: distfun(x,y,mz_size,dt_cut,rt_cut))
    # Y = pdist(frame_rt_all,'euclidean')
    Tree_rt = linkage(dis_rt, method='average')
    clusters_rt = fcluster(Tree_rt,3*rt_cut**2,'distance')

    # save results to alignData_i_j
    n_clusters_rt = len(numpy.unique(clusters_rt))
    alignData_i_j = numpy.zeros((n_clusters_rt,nfeature*nfile))
    tol_s2 = 0
    tol_e2 = 0
    for i in range(nfile):
        tol_e2 = tol_e2 + cellsz2_x[i]
        clusters_label = clusters_rt[range(tol_s2,tol_e2)]
        tol_s2 = tol_e2
        temp = ion_sel[i]
        for j in range(n_clusters_rt):
            # print (i,j)

            if sum(clusters_label==(j+1))>0:
                idx = numpy.where(clusters_label==(j+1))
                alignData_i_j[j,range(i*nfeature,(i+1)*nfeature)]= temp[idx[0]]


    return alignData_i_j