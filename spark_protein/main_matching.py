from submain_matching_par import submain_matching_par

def main_matching(mz_bin_size,rt_cut,dt_cut,N,ionNumThresh,FileName,short_name,OutputPath):

    listVarMatch = ['ion_m_z','ion_drift','ion_rt','ion_z','ion_iso']
    SelShowVar = ['ion_area','ion_counts','clusterID','mwHPlus', 'rt_min',
        'Counts','charge','clust_drift', 'Intensity','ion_driftFWHM',
        'spectrumID','ion_msFWHM','ion_ID','atInflectUpRT','atInflectDownRT','ion_rt']

    Align_RT_name = 'aligned ion_rt'
    selectIO = [1,0,1,1] # criterion to select peptides of interest

    # call ion matching program

    submain_matching_par(N,FileName,short_name,OutputPath,mz_bin_size,rt_cut,dt_cut,listVarMatch,SelShowVar,Align_RT_name,selectIO,ionNumThresh)