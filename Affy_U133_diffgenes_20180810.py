import os
import subprocess
import pandas as pd
import scipy.stats as st
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import math
from glob import glob

def anno_stable(data):
    ##读取芯片的注释信息
    u133_anno = pd.read_csv(r'D:\Affy_annotation\Hsa_U133_plus_2_annotation.csv', sep=',', index_col=0)
    ##读取稳定表达的探针文件
    stable_probe = pd.read_csv(r'D:\mypackage\probeset_stable.txt', sep='\t', index_col=0)

    data_anno = pd.merge(data, u133_anno, how='left', left_index=True, right_index=True)
    data_anno_stable = pd.merge(data_anno, stable_probe, how='inner', left_index=True, right_index=True).iloc[:, :-1]
    # data_anno_filter = data_anno[data_anno.apply(lambda x: '///' not in x['Gene Symbol'] and '---' != x['Gene Symbol'], axis=1)]
    return data_anno_stable, data_anno


def TTest_diff(data, control, case, ctrname, casename, normethd):
    control = [i for i in control]
    case = [i for i in case]
    dt_gr = data.loc[:, control+case]
    dt_gr.loc[:, ctrname+'_mean'] = dt_gr.loc[:, control].mean(axis=1)
    dt_gr.loc[:, casename+'_mean'] = dt_gr.loc[:, case].mean(axis=1)
    dt_gr.loc[:, 'logFC_%s_%s'%(casename, ctrname)] = dt_gr.loc[:, casename+'_mean'] - dt_gr.loc[:, ctrname+'_mean']
    dt_gr.loc[:, 'FC_%s_%s'%(casename, ctrname)] = np.exp2(dt_gr.loc[:, 'logFC_%s_%s'%(casename, ctrname)])

    for probe in dt_gr.index:
        dt_gr.loc[probe, 'Pvalue_%s_%s'%(casename, ctrname)] = st.ttest_ind(dt_gr.loc[probe, case], dt_gr.loc[probe, control]).pvalue

    dt_gr_p = dt_gr.loc[dt_gr['Pvalue_%s_%s'%(casename, ctrname)] <= 0.05]
    dt_gr_p_fc = dt_gr_p.loc[dt_gr_p['FC_%s_%s' % (casename, ctrname)] > 1.2].append(dt_gr_p.loc[dt_gr['FC_%s_%s' % (casename, ctrname)] < 0.8])

    dt_gr_p_as, dt_gr_p_anno = anno_stable(dt_gr_p)
    dt_gr_p_as.to_csv('%s vs %s pstable_%s.txt'%(casename, ctrname, normethd), sep='\t', index=True, index_label='Probe Set ID')  #保存p<= 0.05且属于稳定表达探针的数据

    dt_gr_pfc_as, dt_gr_pfc_anno = anno_stable(dt_gr_p_fc)
    dt_gr_pfc_as.to_csv('%s vs %s pFCfilter_%s.txt'%(casename, ctrname, normethd), sep='\t', index=True, index_label='Probe Set ID') #保存p<= 0.05 & FC> 1.2或FC< 0.8且稳定表达探针的数据


if __name__ == '__main__':
    ##读取包含有分组信息的cel文件列表,分组表头为'group'
    groups = pd.read_csv('celfile_list.txt', sep='\t')
    grouped = pd.groupby(groups, by='groupID')
    groupname = ['LC VS BC'] ## case VS control

    for file in glob('*summary*.txt'):
        filename = os.path.basename(file)
        normethd = filename.split('.')[0]
        data_nor = pd.read_csv(file, sep='\t', index_col=0)
        if 'mas5' in normethd:
            arraycol = [i for i in data_nor.columns if '.mas5-Signal' in i]
            newdatanor = data_nor[arraycol]
            data_nor = np.log2(newdatanor)
            # data_nor = np.log2(data_nor)
            data_nor.columns = data_nor.columns.str.replace('.mas5-Signal', '.CEL')
        for grname in groupname:
            print(grname)
            casename = grname.split('VS')[0].strip()
            controlname = grname.split('VS')[1].strip()
            TTest_diff(data_nor, control=grouped.get_group(controlname)['cel_files'], case=grouped.get_group(casename)['cel_files'], ctrname=controlname, casename=casename, normethd=normethd)









