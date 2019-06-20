# -*- coding: utf-8 -*-
# Created by tobias at 20.06.19

import pandas as pd
import os
import matplotlib.pyplot as plt
from exercise03 import read_data


def normalize_df(df):
    df_norm = df.copy()
    df_norm = (df_norm - df_norm.min()) / (df_norm.max() - df_norm.min())
    return df_norm


def nrmse(df):
    # normalized root mean square error
    return round((((df.iloc[:, 0] - df.iloc[:, 1]) ** 2).mean() ** .5) / df.iloc[:, 1].mean(), 4)


def nae(df):
    # normalized average error
    return round((df.iloc[:, 0].mean() - df.iloc[:, 1].mean()) / df.iloc[:, 1].mean(), 4)


def vr(df):
    # variance ratio
    return round(df.iloc[:, 0].var() / df.iloc[:, 1].var(), 4)


def ioa(df):
    # index of agreement
    nom = ((df.iloc[:, 0] - df.iloc[:, 1]) ** 2).sum()
    mean2 = df.iloc[:, 1].mean()
    denom = ((abs(df.iloc[:, 0] - mean2) + abs(df.iloc[:, 1] - mean2)) ** 2).sum()
    return round(1 - (nom / denom), 4)


def pearson_corr(df):
    # pearson correlation coefficient
    return round(df.iloc[:, 0].corr(df.iloc[:, 1]), 4)


def evaluate_model(datapath, outpath, cell):
    lpjmlpath = os.path.join(datapath, 'LPJmL')
    satpath = os.path.join(datapath, 'Satellite')

    # FAPAR
    mfapar = read_data(os.path.join(lpjmlpath, 'cell_{}'.format(cell), '{}_mfapar.txt'.format(cell)))
    mfapar = mfapar.rename(columns={'mfapar': 'LPJmL-FAPAR'}).fillna(0.0)
    fapar_sat = read_data(os.path.join(satpath, 'MOD15A2H.FPAR.forLPJcells.2000.2018.30days.txt'))
    fapar_sat_cell = fapar_sat[cell].rename('MODIS-FAPAR').fillna(0.0)
    index_min = max([mfapar.index.min(), fapar_sat_cell.index.min()])
    index_max = min([mfapar.index.max(), fapar_sat_cell.index.max()])
    fapar_comb = pd.concat([mfapar, fapar_sat_cell], axis=1, sort=True).loc[index_min:index_max]
    # fapar_comb.plot()
    # plt.show()
    fapar_var = {}
    fapar_var['NRMSE'] = nrmse(fapar_comb)
    fapar_var['NAE'] = nae(fapar_comb)
    fapar_var['VR'] = vr(fapar_comb)
    fapar_var['IoA'] = ioa(fapar_comb)
    fapar_var['Corr'] = pearson_corr(fapar_comb)
    print(fapar_var)

    # GPP/SIF
    # norm to evaluate temporal dynamic of GPP?
    mgpp = read_data(os.path.join(lpjmlpath, 'cell_{}'.format(cell), '{}_mgpp.txt'.format(cell)))
    mgpp = mgpp.rename(columns={'mgpp': 'LPJmL-GPP'}).fillna(0.0)
    mgpp_norm = normalize_df(mgpp)
    sif = read_data(os.path.join(satpath, 'GlobFluo-GOME2.SIF.forLPJcells.2007.2015.30days.txt'))
    sif_cell = sif[cell].rename('SIF').fillna(0.0)
    sif_cell_norm = normalize_df(sif_cell)
    index_min = max([mgpp.index.min(), sif_cell.index.min()])
    index_max = min([mgpp.index.max(), sif_cell.index.max()])
    gpp_sif = pd.concat([mgpp_norm, sif_cell_norm], axis=1, sort=True).loc[index_min:index_max]
    # gpp_sif.plot()
    # plt.show()

    # SWC/SSM
    mswc = read_data(os.path.join(lpjmlpath, 'cell_{}'.format(cell), '{}_mswc1.txt'.format(cell)))
    mswc = mswc.rename(columns={'mswc1': 'LPJmL-SWC'}).fillna(0.0)
    ssm = read_data(os.path.join(satpath, 'ESACCIv050.SSM.forLPJcells.1978.2017.30days.txt'))
    ssm = ssm[cell].rename('ESACCISM').fillna(0.0)
    ssm_norm = normalize_df(ssm)
    mswc_norm = normalize_df(mswc)
    index_min = max([ssm.index.min(), mswc.index.min()])
    index_max = min([ssm.index.max(), mswc.index.max()])
    ssm_comb = pd.concat([mswc_norm, ssm_norm], axis=1, sort=True).loc[index_min:index_max]
    # ssm_comb.plot()
    # plt.show()


if __name__ == '__main__':
    rootpath = os.path.dirname(os.path.realpath(__file__))
    outpath = os.path.join(rootpath, 'results', 'taskB')
    datapath = os.path.join(rootpath, 'data')
    try:
        os.makedirs(outpath)
    except:
        pass
    cells = ['6037', '7460', '25368', '32785', '53569']
    for cell in cells[:1]:
        evaluate_model(datapath, outpath, cell)
