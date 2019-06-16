import os
import pandas as pd


def read_data(fname):
    """
    Simple reader for all .csv
    :param fname: string
        name of file to read
    :return: pandas.Dataframe
        data in file
    """
    
    df = pd.read_csv(fname, ' ')
    try:
        df.index = pd.to_datetime(df.iloc[:, 0])
    except ValueError:
        df.index = df.iloc[:, 0]
    df.drop(df.columns[0], axis=1, inplace=True)
    return df


if __name__ == '__main__':
    # get the path where this script is located
    inpath = os.path.dirname(os.path.realpath(__file__))
    # read file
    data = read_data(os.path.join(inpath, 'data', 'LPJmL', 'cell_6037', '6037_mburnt_area.txt'))




