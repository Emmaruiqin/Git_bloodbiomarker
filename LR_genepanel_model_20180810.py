import pandas as pd
import statsmodels.api as sm
import random
import numpy as np
import matplotlib.pyplot as plt
import sklearn.metrics as metrics
import itertools
import pickle
import os
import time
import shutil
from glob import glob
from sklearn import cross_validation
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestRegressor
import seaborn as sns
from sklearn.externals import joblib

def add_group(data, groupdict):
    for i in data.index:
        data.loc[i, 'response'] = groupdict[i]['groupID']
    return data

def main(datafile, groupfile):
    data = pd.read_csv(datafile, sep='\t', index_col=0)
    group = pd.read_csv(groupfile, sep='\t', index_col=0)

    genes = [i for i in data.columns if 'response' not in i]

    groupdict = group.to_dict(orient='index')
    if 'response' not in data.columns:
        data_group = add_group(data, groupdict)
    else:
        data_group = data

    group_dummy = pd.get_dummies(data_group['response'], prefix='response')
    auc_gene = []
    train_data_all, test_data_all, train_target, test_target = cross_validation.train_test_split(data_group[genes], group_dummy['response_LC'], test_size=0.3, random_state=None)
    train_data = train_data_all[genes]
    test_data = test_data_all[genes]
    lg_model = LogisticRegression()
    lg_model.fit(train_data, train_target)
    test_est = lg_model.predict(test_data)
    train_est = lg_model.predict(train_data)
    test_est_p = lg_model.predict_proba(test_data)[:,1]
    train_est_p = lg_model.predict_proba(train_data)[:,1]
    test_c = metrics.classification_report(test_target, test_est)
    train_c = metrics.classification_report(train_target, train_est)
    avg_test = '\t'.join(test_c.split('\n')[-2].split()[-4:])
    avg_train = '\t'.join(train_c.split('\n')[-2].split()[-4:])
    # loss_test = metrics.zero_one_loss(test_target, test_est)
    # loss_train = metrics.zero_one_loss(train_target,train_est)
    auc_test = metrics.roc_auc_score(test_target,test_est_p)
    auc_train = metrics.roc_auc_score(train_target, train_est_p)
    accu_train = lg_model.score(train_data, train_target)
    accu_test = lg_model.score(test_data, test_target)

    # lg_scores = cross_validation.cross_val_score(lg_model, train_data, train_target, cv=10)
    lg_scores = cross_validation.cross_val_score(lg_model, data_group[genes], group_dummy['response_LC'], cv=10)
    lg_scores_mean = lg_scores.mean()
    lg_scores_std = lg_scores.std()
    print(lg_scores_mean, lg_scores_std, accu_train, accu_test)
    # ROC_save(target_raw=train_target, target_pre=train_est_p, labname='all')
    # ROC_save(target_raw=train_target, target_pre=train_est_p, labname='all')
    # ROC_save(target_raw=train_target, target_pre=train_est_p, labname='all')

    fpr_train, tpr_train, th_train = metrics.roc_curve(train_target, train_est_p)
    fpr_test, tpr_test, th_test = metrics.roc_curve(test_target, test_est_p)
    auc_train = metrics.roc_auc_score(train_target, train_est_p)
    auc_test = metrics.roc_auc_score(test_target, test_est_p)
    plt.figure(figsize=(13, 13))
    plt.plot(fpr_train, tpr_train, color='red', lw=2, label='Training AUC(%.4f)' %auc_train)
    plt.plot(fpr_test, tpr_test, color='blue', lw=2, label='Test AUC(%.4f)' %auc_test)
    plt.plot([0, 1], [0, 1], color='black', lw=2, linestyle='--')
    plt.title('ROC curve of logit')
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.legend(loc='lower right')
    plt.savefig('ROC_best.png', dpi=800)

    print(train_c)
    print(test_c)
    with open('LR_model_classi.txt', 'w') as file:
        file.write('Train result\n')
        file.write(train_c)
        file.write('\nTest result\n')
        file.write(test_c)


    result = '\t'.join(genes) + '\t' + avg_train + '\t' +avg_test + '\t' +\
             str(auc_train) +'\t'+str(auc_test)+'\t'+str(lg_scores_mean)+'\t'+str(lg_scores_std)+'\t'+str(accu_train)+'\t'+str(accu_test)+'\n'
    auc_gene.append(result)

    with open('LR_model_assess.txt', 'w') as file:
        head1 = 'gene1\tgene2\tgene3\tgene4\tgene5\tgene6\t'
        head2 = 'precision_train\trecall_train\tf1_score_train\tsupport_train\t'
        head3 = 'precision_test\trecall_test\tf1_score_test\tsupport_test\t'
        head4 = 'auc_train\tauc_test\tscore_cv_mean\tscore_cv_std\taccuracy_train\taccuracy_test\n'
        newhead = head1+head2+head3+head4
        file.write(newhead)
        file.writelines(auc_gene)

    joblib.dump(lg_model, 'LR_model.m')


if __name__ == '__main__':
    data_all = '候选panel_xin.txt'
    groupfile = 'celfile_list.txt'
    main(datafile=data_all, groupfile=groupfile)
