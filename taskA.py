import os
import numpy as np
import pandas as pd

from exercise03 import read_data


def calc_variables(lpjmlpath, outpath, cell, save=False):
    mevap = read_data(os.path.join(lpjmlpath, 'cell_{}'.format(cell), '{}_mevap.txt'.format(cell)))
    mtransp = read_data(os.path.join(lpjmlpath, 'cell_{}'.format(cell), '{}_mtransp.txt'.format(cell)))
    npp = read_data(os.path.join(lpjmlpath, 'cell_{}'.format(cell), '{}_mnpp.txt'.format(cell)))
    rh = read_data(os.path.join(lpjmlpath, 'cell_{}'.format(cell), '{}_mrh.txt'.format(cell)))
    gpp = read_data(os.path.join(lpjmlpath, 'cell_{}'.format(cell), '{}_mgpp.txt'.format(cell)))
    fapar = read_data(os.path.join(lpjmlpath, 'cell_{}'.format(cell), '{}_mfapar.txt'.format(cell)))
    swc = read_data(os.path.join(lpjmlpath, 'cell_{}'.format(cell), '{}_mswc1.txt'.format(cell)))
    vars_cell = pd.concat([mevap, mtransp, npp, rh, gpp, fapar, swc], axis=1, sort=True)
    vars_cell['met'] = vars_cell[['mevap', 'mtransp']].sum(axis=1)
    vars_cell['mnee'] = vars_cell['mrh'] - vars_cell['mnpp']

    if save:
        vars_cell.to_csv(os.path.join(outpath, 'cell_{}_vars.csv'.format(cell)))

    return vars_cell


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
    for cell in cells:
        vars_cells[cell] = calc_variables(lpjmlpath, outpath, cell, True)
    print(vars_cells)
