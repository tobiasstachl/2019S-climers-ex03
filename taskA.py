import os
import matplotlib.pyplot as plt
import pandas as pd

from exercise03 import read_data


def calc_variables(lpjmlpath, outpath, cell, save=False):
    mevap = read_data(os.path.join(lpjmlpath, 'cell_{}'.format(cell), '{}_mevap.txt'.format(cell)))
    mtransp = read_data(os.path.join(lpjmlpath, 'cell_{}'.format(cell), '{}_mtransp.txt'.format(cell)))
    mnpp = read_data(os.path.join(lpjmlpath, 'cell_{}'.format(cell), '{}_mnpp.txt'.format(cell)))
    mrh = read_data(os.path.join(lpjmlpath, 'cell_{}'.format(cell), '{}_mrh.txt'.format(cell)))
    mgpp = read_data(os.path.join(lpjmlpath, 'cell_{}'.format(cell), '{}_mgpp.txt'.format(cell)))
    mfapar = read_data(os.path.join(lpjmlpath, 'cell_{}'.format(cell), '{}_mfapar.txt'.format(cell)))
    mswc = read_data(os.path.join(lpjmlpath, 'cell_{}'.format(cell), '{}_mswc1.txt'.format(cell)))
    vars_cell = pd.concat([mevap, mtransp, mnpp, mrh, mgpp, mfapar, mswc], axis=1, sort=True)
    vars_cell['met'] = vars_cell[['mevap', 'mtransp', 'minterc']].sum(axis=1)
    vars_cell['mnee'] = vars_cell['mrh'] - vars_cell['mnpp']

    if save:
        vars_cell.to_csv(os.path.join(outpath, 'cell_{}_vars.csv'.format(cell)))

    return vars_cell


def plot_vars(vars_cell, outpath, save=False):
    '''
    Plot monthly time series of FAPAR, SWC, ET, GPP and NEE for each grid cell
    '''
    fig, axes = plt.subplots(5, 1, sharex='all', figsize=(12, 12))

    # FAPAR
    vars_cell['mfapar'].plot(ax=axes[0])
    axes[0].set_title('FAPAR {}'.format(vars_cell.name))
    axes[0].set_ylabel('[-]')
    axes[0].set_ylim(0, 1)

    # SWC
    vars_cell['mswc1'].plot(ax=axes[1])
    axes[1].set_title('Soil Water Content (SWC) {}'.format(vars_cell.name))
    axes[1].set_ylabel('[-]')
    axes[1].set_ylim(0, 1)

    # ET
    vars_cell['met'].plot(ax=axes[2])
    axes[2].set_title('Evapotranspiration (ET) {}'.format(vars_cell.name))
    axes[2].set_ylabel('[mm/month]')

    # GPP
    vars_cell['mgpp'].plot(ax=axes[3])
    axes[3].set_title('Gross primary production (GPP) {}'.format(vars_cell.name))
    axes[3].set_ylabel('[gC/m2/month]')

    # NEE
    vars_cell['mnee'].plot(ax=axes[4])
    axes[4].set_title('Net ecosystem exchange (NEE) {}'.format(vars_cell.name))
    axes[4].set_ylabel('[gC/m2/month]')

    # save figure
    if save:
        plt.savefig(os.path.join(outpath, 'cell_vars_{}.png'.format(cell)), bbox_inches="tight")
    else:
        plt.show()
    plt.close()


def calc_annually(lpjmlpath, cell, title=None, out_path=None):
    '''
    Compute annual GPP and net biome productivity (NBP)
    NPB = NPP - Fire - Rh
    :return:
    '''

    mnpp = read_data(os.path.join(lpjmlpath, 'cell_{}'.format(cell), '{}_mnpp.txt'.format(cell)))
    mfirec = read_data(os.path.join(lpjmlpath, 'cell_{}'.format(cell), '{}_mfirec.txt'.format(cell)))
    mrh = read_data(os.path.join(lpjmlpath, 'cell_{}'.format(cell), '{}_mrh.txt'.format(cell)))
    mgpp = read_data(os.path.join(lpjmlpath, 'cell_{}'.format(cell), '{}_mgpp.txt'.format(cell)))
    vegc = read_data(os.path.join(lpjmlpath, 'cell_{}'.format(cell), '{}_vegc.txt'.format(cell)))
    vegc = vegc.asfreq('YS')
    fpc = read_data(os.path.join(lpjmlpath, 'cell_{}'.format(cell), '{}_fpc.txt'.format(cell)))
    vars_cell = pd.concat([mnpp, mfirec, mrh, mgpp], axis=1, sort=True)
    vars_cell['mnbp'] = vars_cell['mnpp'] - vars_cell['mfirec'] - vars_cell['mrh']
    vars_cell_yearly = vars_cell.resample('YS').sum()


    fpc_trees = fpc[["fpc TrBE", "fpc TrBR", "fpc TeNE", "fpc TeBE", "fpc TeBS", "fpc BoNE", "fpc BoBS", "fpc BoNS"]].asfreq('YS')
    dominant_tree_class = fpc_trees.mean(axis=0).sort_values(ascending=False).index[0]
    fpc_dominant_tree_ts = fpc_trees[dominant_tree_class]

    fpc_grass = fpc[["fpc TrH", "fpc TeH"]].asfreq('YS')
    dominant_grass_class = fpc_grass.mean(axis=0).sort_values(ascending=False).index[0]
    fpc_dominant_grass_ts = fpc_grass[dominant_grass_class]

    yearly_data = pd.concat([vars_cell_yearly, fpc_dominant_tree_ts, fpc_dominant_grass_ts, vegc], axis=1)
    yearly_data = yearly_data.drop(columns=['mnpp', 'mfirec', 'mrh'], axis=1)
    yearly_data = yearly_data.rename(columns={'mgpp': 'gpp', 'mnbp': 'nbp'})


    # calculate correlations
    corrs = yearly_data.corr()

    #vars_cell_yearly.plot(subplots=True)
    axes = yearly_data.plot(subplots=True, figsize=(12, 8))

    # ['mgpp', 'mnbp', 'fpc TeBE', 'fpc TeH', 'vegc']

    unit_dict = {'gpp' : '$gC m^{-2} yr^{-1}$',
                 'nbp' : '$gC m^{-2} yr^{-1}$',
                 'vegc': '$gC m^{-2} yr^{-1}$'}

    for ax in axes:
        handles, labels = ax.get_legend_handles_labels()
        ax_label = labels[0]

        if 'fpc' in ax_label:
            unit = '%'
        else:
            unit = unit_dict[ax_label]

        ax.set_ylabel('({})'.format(unit))

    if title is not None:
        axes[0].set_title(title)

    plt.tight_layout()

    if outpath is None:
        print(corrs)
        plt.show()
    else:
        corrs.to_csv(os.path.join(outpath, 'annual_corrs_cell_{}.csv'.format(cell)))
        plt.savefig(os.path.join(outpath, 'annual_cell_{}.png'.format(cell)))


if __name__ == '__main__':
    rootpath = os.path.dirname(os.path.realpath(__file__))
    outpath = os.path.join(rootpath, 'results', 'taskA')
    lpjmlpath = os.path.join(rootpath, 'data', 'LPJmL')
    try:
        os.makedirs(outpath)
    except:
        pass


    site_dict = {'32785': 'Sahel (23.75°E, 7.75°N), cropland area: 11%',
                 '6037': 'Mexico (105.25°W, 28.75°N), cropland area: 85%',
                 '7460': 'Canada (99.25°W, 52.25°N), cropland area: 0.3%',
                 '25368': 'Spain (3.25°W, 39.75°N), cropland area: 77%',
                 '53569': 'Taymyr (102.25°E, 75.75°N), cropland area: 0%'}


    cells = ['6037', '7460', '25368', '32785', '53569']
    vars_cells = {}
    for cell in cells:
        vars_cells[cell] = calc_variables(lpjmlpath, outpath, cell, False)
        vars_cells[cell].name = cell

        # get plot title for cell
        cell_title = site_dict[cell]
        plot_vars(vars_cells[cell], outpath, True)
        calc_annually(lpjmlpath, cell, cell_title, outpath)
