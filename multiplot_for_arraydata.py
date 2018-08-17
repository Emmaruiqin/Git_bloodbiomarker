import pandas as pd
import seaborn as sns
import os
import numpy as np
import scipy.stats as st

def scatterplot(x, y, xlable, ylable, FC, name, fmts=None, dpi=None):
    sns.set_context('paper')
    t = math.log2(int(FC))
    if t:
        c = []
        for i, j in zip(x, y):
            if i - j >= t:
                c.append('red')
            elif i - j <= -t:
                c.append('green')
            else:
                c.append('blue')
    else:
        c = 'blue'
    xmin,xmax,ymin,ymax = min(x),max(x),min(y),max(y)
    plt.axis([xmin, xmax, ymin, ymax])
    plt.scatter(x=x, y=y, c=c, marker='.', edgecolors=c, linewidths=1, alpha=0.8, s=1)
    plt.plot([xmin, xmax], [xmin + t, xmax + t], 'k-', [xmin, xmax], [xmin - t, xmax - t], 'k-', [xmin, xmax],
             [xmin, xmax], 'k-', lw=0.5)
    plt.ylabel(ylable)
    plt.xlabel(xlable)
    plt.title('Scatter Plot')
    for fmt in fmts:
        if dpi and fmt != 'pdf':
            plt.savefig(name + '_ScatterPlot.' + fmt, dpi=dpi)
        else:
            plt.savefig(name + '_ScatterPlot.' + fmt)
    plt.close()
#

def main(arrayfile, groupfile, case, control):
    groupdata = pd.read_csv(groupfile, sep='\t', index_col=0)
    # groupdict = groupdata.to_dict(orient='index')
    grouped = groupdata.groupby('groupID')
    casesample = grouped.get_group(case).index.tolist()
    controlsample = grouped.get_group(control).index.tolist()

    arraydata = pd.read_csv(arrayfile, sep='\t', index_col=0)

    analysisdata = arraydata_analysis(arraydata, case, casesample, control, controlsample)
    scatterplot(x=analysisdata.loc[:, '%s_mean'%case], y=analysisdata.loc[:, '%s_mean'%control])




def arraydata_analysis(arraydata, case, casesample, control, controlsample):
    cols = [i for i in arraydata.columns if '.mas5-Signal' in i]
    arraydata = arraydata[cols]
    arraydata.columns = arraydata.columns.str.replace('.mas5-Signal', '.CEL')
    arraydata = np.log2(arraydata)
    arraydata.loc[:, case + '_mean'] = arraydata.loc[:, casesample].mean(axis=1)
    arraydata.loc[:, control + '_mean'] = arraydata.loc[:, controlsample].mean(axis=1)
    arraydata.loc[:, 'logFC'] = arraydata[case + '_mean'] - arraydata[control + '_mean']
    arraydata.loc[:, 'FC'] = np.exp2(arraydata.loc[:, 'logFC'])
    arraydata.loc[:, 'pvalue'] = st.ttest_ind(arraydata.loc[:, casesample], arraydata.loc[:, controlsample],
                                              axis=1).pvalue
    return arraydata


if __name__ == '__main__':
    arrayfile = 'C:\\Users\\Administrator\\Desktop\\LC_BC_39\\LC_BC_重新分析结果\\mas5.summary.txt'
    groupfile = 'C:\\Users\\Administrator\\Desktop\\LC_BC_39\\LC_BC_重新分析结果\\celfile_list.txt'
    case = 'LC'
    control = 'BC'
    main(arrayfile=arrayfile, groupfile=groupfile,  case=case, control=control)