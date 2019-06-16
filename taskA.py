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
    vars_cell['met'] = vars_cell[['mevap', 'mtransp']].sum(axis=1)
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


def calc_annually(lpjmlpath, outpath, cell, save=False):
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
    fpc = read_data(os.path.join(lpjmlpath, 'cell_{}'.format(cell), '{}_fpc.txt'.format(cell)))
    vars_cell = pd.concat([mnpp, mfirec, mrh, mgpp], axis=1, sort=True)
    vars_cell['mnbp'] = vars_cell['mnpp'] - vars_cell['mfirec'] - vars_cell['mrh']
    vars_cell_yearly = vars_cell.resample('YS').sum()
    max_tree = fpc[["fpc TrBE", "fpc TrBR", "fpc TeNE", "fpc TeBE", "fpc TeBS", "fpc BoNE", "fpc BoBS", "fpc BoNS"]].idxmax(axis=1)
    max_grass = fpc[["fpc TrH", "fpc TeH"]].idxmax(axis=1)
    print(max_tree, max_grass)


if __name__ == '__main__':
    rootpath = os.path.dirname(os.path.realpath(__file__))
    outpath = os.path.join(rootpath, 'results', 'taskA')
    lpjmlpath = os.path.join(rootpath, 'data', 'LPJmL')
    try:
        os.makedirs(outpath)
    except:
        pass
    cells = ['6037', '7460', '25368', '32785', '53569']
    vars_cells = {}
    for cell in cells[:1]:
        vars_cells[cell] = calc_variables(lpjmlpath, outpath, cell, False)
        vars_cells[cell].name = cell
        # plot_vars(vars_cells[cell], outpath, True)
        calc_annually(lpjmlpath, outpath, cell, False)
