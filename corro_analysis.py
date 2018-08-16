import pandas as pd
from scipy.stats import pearsonr
import os
import numpy as np



if __name__ == '__main__':
    data = pd.read_csv('LC vs BC pstable_mas5.txt', sep='\t', index_col=0)
    signalcol = [i for i in data.columns if '.CEL' in i]
    dt = data[signalcol]
    dt_T = dt.T
    probecol = dt_T.columns

    groupfile = pd.read_csv('celfile_list.txt', sep='\t', index_col=0)
    dt_T_g = pd.merge(dt_T, groupfile, how='inner', right_index=True, left_index=True)
    for i in dt_T_g.index:
        if dt_T_g.loc[i, 'group'] == 'LC':
            dt_T_g.loc[i, 'groupnum'] = 2
        else:
            dt_T_g.loc[i, 'groupnum'] = 1
    dict = {}
    for i in probecol:
        cor, pval = pearsonr(dt_T_g.loc[:,i], dt_T_g.loc[:,'groupnum'])
        cor_2 = cor**2
        dict[i] = {'corr':'%.4f'%cor, 'pvalue':'%.4f'%pval, 'corr^2':'%.4f'%cor_2}
    corr_res = pd.DataFrame.from_dict(dict, orient='index')
    corr_res.sort(['corr^2'], ascending=False, inplace=True)
    corr_res.to_csv('corr_probe_group.txt', index=True, sep='\t')
    # corrpd = dt_T_g[probecol].corr()
    # corrpd.to_csv('corr_probe_probe.txt', index=True, sep='\t')
    corr_60 = corr_res.head(60).index
    corrpd_60 = dt_T_g[corr_60].corr()
    corrpd_60.to_csv('corr_probe_60.txt', index=True, sep='\t')
    list = []
    for i in corrpd_60.index:
        aa = corrpd_60.loc[i,:]
        aa.sort(inplace=True)
        max_d = aa[-2]
        max_in = aa.index[-2]
        min_d = aa[0]
        min_in = aa.index[0]
        new = str(i) + '\t' + str(max_in) + '\t' + str('%.4f'%max_d) + '\t' + str(min_in) + '\t' + str('%.4f'%min_d) + '\n'
        list.append(new)
    with open('corr_60_max_min.txt','w') as file:
        file.write('probe1\tprobe_max\tmax\tprobe_min\tmin\n')
        file.writelines(list)





