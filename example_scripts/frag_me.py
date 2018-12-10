#!/usr/bin/env python3
import peptide_fragmentor
import pandas as pd
import pyqms
import click



if __name__ == '__main__':
    pd.set_option('display.max_columns', 500)
    pep = 'ACDEFGHI' #Phospho:2' #Acetyl:0'
    fragger = peptide_fragmentor.PeptideFragment0r(pep)
    df = fragger.df
    print(df.head(10))
    print(df.describe())
    df_by = df[df['series'].isin(['b','y'])]
    print(df_by[['name', 'modstring', 'mz']].sort_values('mz'))
