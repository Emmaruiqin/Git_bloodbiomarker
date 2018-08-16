import pandas as pd
from sklearn import cross_validation
from sklearn import svm
import statsmodels.api as sm
from sklearn.decomposition import PCA
from sklearn.linear_model import LogisticRegression
from sklearn.linear_model import Lasso, LassoCV
from sklearn.feature_selection import RFE
from sklearn.tree import DecisionTreeRegressor
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import cross_val_score, ShuffleSplit
import numpy as np

def rfe_select(data, genes):
    lg = LogisticRegression()
    rfe = RFE(estimator=lg, step=1)
    rfe.fit(data[genes], data['groupID'])
    ranking = pd.Series(rfe.ranking_, index=genes)
    print(ranking)
    return ranking

def main(datafile, groupfile, dummy_group,):
    data = pd.read_csv(datafile, index_col=0, sep='\t')
    genes = [i for i in data.columns if 'groupID' not in i]
    group = pd.read_csv(groupfile, index_col=0, sep='\t')

    if 'groupID' not in data.columns:
        data_group = pd.merge(data, group, how='left', left_index=True, right_index=True)
    else:
        data_group = data

    rferanking = rfe_select(data=data_group, genes=genes)
    gi = [int(1) if i == 'LC' else int(0) for i in data_group['groupID']]

    lasso = LassoCV().fit(data_group[genes], gi)
    coef = pd.DataFrame(lasso.coef_, index=genes)
    coef.rename(columns={0:'lasso'}, inplace=True)
    coef['rfe'] = rferanking
    print(coef)

    ##random forest
    rfc = RandomForestRegressor()
    rfc_model = rfc.fit(data_group[genes], gi)
    coef['random'] = rfc_model.feature_importances_

    coef.to_excel('feature choose result.xlsx')


if __name__ == '__main__':
    datafile = 'genes_all_PFC_候选10.txt'
    groupfile = 'celfile_list.txt'
    dummy_group = 'LC'
    main(datafile=datafile, groupfile=groupfile, dummy_group=dummy_group)