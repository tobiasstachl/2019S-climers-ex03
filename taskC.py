# -*- coding: utf-8 -*-
# Created by tobias at 23.06.19
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from error_metrics import KGE
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


def model_params_vs_performance_plot():
    """
    The initial plot that matthias critized.
    """
    cell = '32785'
    pars_set = ['pars{}'.format(i) for i in range(1, 51)]

    parameters = read_data(
        os.path.join(datapath, 'LPJmL', 'cell_32785_parameter-sets.txt'))
    print(parameters)

    metric_sets = {}
    param_sets = {}
    for i, pars in enumerate(pars_set):
        print()
        i += 1
        fapar_comb = read_fapar(cellpath, satpath, pars, cell)
        fapar_metrics = calc_metrics(fapar_comb, 'FAPAR')
        sif_comb = read_sif(cellpath, satpath, pars, cell)
        sif_metrics = calc_metrics(sif_comb, 'GPP-SIF')
        ssm_comb = read_swc(cellpath, satpath, pars, cell)
        ssm_metrics = calc_metrics(ssm_comb, 'SSM')
        metric_results = pd.concat([fapar_metrics, sif_metrics, ssm_metrics])

        # store metrics and params per run
        metric_sets[i] = metric_results.unstack()
        param_sets[i] = parameters.loc[pars]

    # to df
    par_set_metrics = pd.DataFrame.from_dict(metric_sets).T
    print(par_set_metrics)

    par_set_params = pd.DataFrame.from_dict(param_sets).T
    print(par_set_params)

    # plot
    fig, (ax1, ax2, ax3, ax4) = plt.subplots(nrows=4, sharex=True,
                                             figsize=(12, 8))
    par_set_metrics[['NRMSE']].plot(ax=ax1)
    par_set_metrics[['Corr']].plot(ax=ax2)
    par_set_params['WATER_BASE'].plot(ax=ax3)
    par_set_params['EMAX'].plot(ax=ax4)

    # labeling
    ax1.set_title('Parameter optimisation')
    ax1.set_ylabel('NRMSE (-)')
    ax1.legend(loc='upper right')
    ax2.set_ylabel('Pearson R (-)')
    ax2.legend(loc='upper right')
    ax3.set_ylabel('WATER_BASE (%)')
    ax4.set_ylabel('EMAX ($mm day^{-1}$)')

    plt.xlabel('Model run')
    plt.tight_layout()
    plt.savefig(os.path.join(outpath, '1_model_params_vs_performance.png'))
    plt.close(fig)


def kge_scatterplot():
    """
    Creates a scatterplot of WATER_BASE vs. EMAX parameters with the hue
    given by the Kling-Gupta efficiency (KGE).

    KGE is calculated only based on sim and obs FAPAR as we saw that the
    agreement for SSM and SIF-GPP are quite good regardless of the
    parameter settings.

    Kling-Gupta efficiencies range from -Inf to 1.
    Essentially, the closer to 1, the more accurate the model is.
    """

    cell = '32785'
    pars_set = ['pars{}'.format(i) for i in range(1, 51)]

    parameters = read_data(
        os.path.join(datapath, 'LPJmL', 'cell_32785_parameter-sets.txt'))

    metric_sets = {'kge': [],
                   'cc': [],
                   'alpha': [],
                   'beta': []}
    param_sets = {}
    for i, pars in enumerate(pars_set):
        i += 1
        fapar_comb = read_fapar(cellpath, satpath, pars, cell)
        fapar_sim = fapar_comb['mfapar'].values
        fapar_obs = fapar_comb['MODIS-FAPAR'].values

        # calc KGE for each parameter set
        kge, cc, alpha, beta = KGE(s=fapar_sim, o=fapar_obs)

        # store metrics and params per run
        metric_sets['kge'].append(kge)
        metric_sets['cc'].append(cc)
        metric_sets['alpha'].append(alpha)
        metric_sets['beta'].append(beta)
        param_sets[i] = parameters.loc[pars]

    # to df
    par_set_metrics = pd.DataFrame.from_dict(metric_sets)
    # reindex to start from 1
    par_set_metrics.index = np.arange(1, len(par_set_metrics) + 1)
    par_set_params = pd.DataFrame.from_dict(param_sets).T

    # merge -> output could be given to plot function from here on
    df_merged = pd.concat([par_set_params, par_set_metrics], axis=1)
    sorted_by_kge = df_merged.sort_values('kge', ascending=False)
    sorted_by_kge.to_csv(os.path.join(outpath,
                                      'metrics_for_param_settings.csv'))

    # plot
    fig, ax = plt.subplots(figsize=(12,5))
    x = 'WATER_BASE'
    y = 'EMAX'
    sns.scatterplot(x=x,
                         y=y,
                         size='cc',
                         hue='kge',
                         data=df_merged,
                    ax=ax,
                    legend='brief')

    # annotate setting number
    for i, txt in enumerate(df_merged.index.values):
        ax.annotate(txt, (df_merged[x].iloc[i], df_merged[y].iloc[i]),
                    color='grey')

    # Put a legend to the right of the current axis
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

    plt.title('Optimisation of model parameters via KGE')
    plt.xlabel('WATER_BASE (%)')
    plt.ylabel('EMAX ($mm day^{-1}$)')
    plt.tight_layout()
    plt.savefig(os.path.join(outpath, '2_kge_scatterplot.png'))
    plt.close(fig)

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

    kge_scatterplot()