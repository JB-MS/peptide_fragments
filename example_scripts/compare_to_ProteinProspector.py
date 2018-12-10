#!/usr/bin/env python3
import peptide_fragmentor
import pandas as pd
import pyqms
import click
import sys
import os

if __name__ == '__main__':
    pep = 'KLEINER'
    fragger = peptide_fragmentor.PeptideFragment0r(pep)
    df = fragger.df
    df_by = df[df['series'].isin(['b','y'])]
    check_summary = {}
    pp_file_name = os.path.join('../data/{0}_ProteinProspector.tsv'.format(pep))
    if os.path.exists(pp_file_name) is True:
        # compare the shiznit
        pp_df = pd.read_csv(
            pp_file_name,
            sep='\t'
        )

        for index, row in df_by.iterrows():
            for sign in ['-']:
                ion = '{0}{1}{2}'.format(
                    row['name'],
                    sign,
                    row['modstring']
                )
                pf_mz = round(float(row['mz']),4)
                if ion in set(pp_df['ion']):
                    pp_mz = round(float(pp_df['mz'][pp_df['ion']==ion].values[0]),4)
                    if  ion not in check_summary.keys():
                        check_summary[ion] = {
                            'ion' : ion,
                            'pf_mz': pf_mz,
                            'pp_mz': pp_mz
                        }
                    if pp_mz == pf_mz:
                        check_summary[ion]['equal_mz']  = True

    verified_ions = []
    for ion, ion_dict in check_summary.items():
        if ion_dict['equal_mz'] is True:
            verified_ions.append(ion)
        else:
            print(
                '{0}: m/z not consistent (PP m/z: {0}; PF m/z: {1})'.format(
                    ion,
                    ion_dict['pp_mz'],
                    ion_dict['pf_mz'],
                )
            )
    if len(verified_ions) == len(check_summary.keys()):
        print('All ion m/z could be verified; {0} ions'.format(len(verified_ions)))
    else:
        print('Ion m/z verification failed')
