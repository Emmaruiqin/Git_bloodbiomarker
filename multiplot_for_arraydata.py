import pandas as pd
import seaborn as sns
import os



def main(arrayfile, groupfile, arraydatatype):
    arraydata = pd.read_csv(arrayfile, sep='\t', index_col=0)
    if arraydatatype == 'mas5':
        cols = [i for i in arraydata.columns if '.mas5-Signal']
        arraydata = arraydata[cols]
    else:
        arraydata = arraydata
    groupdata = pd.read_csv(groupfile, sep='\t', index_col=0)
    groupdict = groupdata.to_dict(orient='index')




if __name__ == '__main__':
    arrayfile = 'C:\Users\Administrator\Desktop\LC_BC_39\LC_BC_重新分析结果\mas5.summary.txt'
    groupfile = 'C:\Users\Administrator\Desktop\LC_BC_39\LC_BC_重新分析结果\celfile_list.txt'
    arraydatatype = 'mas5'
    main(arrayfile=arrayfile, groupfile=groupfile, arraydatatype)