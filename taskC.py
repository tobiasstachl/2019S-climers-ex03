# -*- coding: utf-8 -*-
# Created by tobias at 23.06.19
import os
import pandas as pd
from taskB import nrmse, pearson_corr, normalize_df
from exercise03 import read_data


def calc_metrics(df, index_name=''):
    results_dict = {}
    results_dict['NRMSE'] = nrmse(df)
    results_dict['Corr'] = pearson_corr(df)
    return pd.DataFrame(results_dict, index=[index_name])


def read_fapar(cellpath, satpath, pars, cell):
    fapar = read_data(os.path.join(cellpath, pars, '{}_{}_mfapar.txt'.format(cell, pars)))
    fapar_sat = read_data(os.path.join(satpath, 'MOD15A2H.FPAR.forLPJcells.2000.2018.30days.txt'))
    fapar_sat_cell = fapar_sat[cell].rename('MODIS-FAPAR')
    index_min = max([fapar.index.min(), fapar_sat_cell.index.min()])
    index_max = min([fapar.index.max(), fapar_sat_cell.index.max()])
    fapar_comb = pd.concat([fapar, fapar_sat_cell], axis=1, sort=True).loc[index_min:index_max]
    return fapar_comb


def read_sif(cellpath, satpath, pars, cell):
    mgpp = read_data(os.path.join(cellpath, pars, '{}_{}_mgpp.txt'.format(cell, pars)))
    mgpp = mgpp.rename(columns={'mgpp': 'LPJmL-GPP'})
    mgpp_norm = normalize_df(mgpp)
    sif = read_data(os.path.join(satpath, 'GlobFluo-GOME2.SIF.forLPJcells.2007.2015.30days.txt'))
    sif_cell = sif[cell].rename('SIF')
    sif_cell_norm = normalize_df(sif_cell)
    index_min = max([mgpp.index.min(), sif_cell.index.min()])
    index_max = min([mgpp.index.max(), sif_cell.index.max()])
    gpp_sif = pd.concat([mgpp_norm, sif_cell_norm], axis=1, sort=True).loc[index_min:index_max]
    return gpp_sif


def read_swc(cellpath, satpath, pars, cell):
    mswc = read_data(os.path.join(cellpath, pars, '{}_{}_mswc1.txt'.format(cell, pars)))
    mswc = mswc.rename(columns={'mswc1': 'LPJmL-SWC'})
    ssm = read_data(os.path.join(satpath, 'ESACCIv050.SSM.forLPJcells.1978.2017.30days.txt'))
    ssm = ssm[cell].rename('ESACCISM')
    ssm_norm = normalize_df(ssm)
    mswc_norm = normalize_df(mswc)
    index_min = max([ssm.index.min(), mswc.index.min()])
    index_max = min([ssm.index.max(), mswc.index.max()])
    ssm_comb = pd.concat([mswc_norm, ssm_norm], axis=1, sort=True).loc[index_min:index_max]
    return ssm_comb


if __name__ == '__main__':
    rootpath = os.path.dirname(os.path.realpath(__file__))
    outpath = os.path.join(rootpath, 'results', 'taskC')
    datapath = os.path.join(rootpath, 'data')
    cellpath = os.path.join(datapath, 'LPJmL', 'cell_32785')
    satpath = os.path.join(datapath, 'Satellite')
    try:
        os.makedirs(outpath)
    except:
        pass

    cell = '32785'

    pars_set = ['pars{}'.format(i) for i in range(1, 51)]
    for pars in pars_set[:1]:
        fapar_comb = read_fapar(cellpath, satpath, pars, cell)
        fapar_metrics = calc_metrics(fapar_comb, 'FAPAR')
        sif_comb = read_sif(cellpath, satpath, pars, cell)
        sif_metrics = calc_metrics(sif_comb, 'GPP-SIF')
        ssm_comb = read_swc(cellpath, satpath, pars, cell)
        ssm_metrics = calc_metrics(ssm_comb, 'SSM')
        metric_results = pd.concat([fapar_metrics, sif_metrics, ssm_metrics])
        print(metric_results)



