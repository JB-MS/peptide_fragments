#!/usr/bin/env python3
from itertools import combinations
from collections import defaultdict as ddict
import pandas as pd
import numpy as np
import pyqms
from pyqms.chemical_composition import ChemicalComposition
import copy
import pprint

import peptide_fragmentor


class PeptideFragment0r:
    def __init__(self, upep, charges=None, neutral_losses=None):
        """
        Initialize framentOr with peptide `upep`.

        Args:
            upep (str): Peptide with optional Unimod modification string in the
                format PEPTIDE#<UNIMOD_NAME>:<POS>;<UNIMOD_NAME>:<POS> ...
            charges (list, optional): Charges for frag ion creation, default
                is 1, 2, 3
            neutral_losses (list, optional): Description
        """
        if charges is None:
            self.charges = [1, 2, 3]
        else:
            self.charges = charges
        if neutral_losses is None:
            neutral_losses = peptide_fragmentor.neutral_losses
        else:
            neutral_losses += peptide_fragmentor.neutral_losses
        self.neutral_losses = neutral_losses

        self.upep_cc = ChemicalComposition(upep)
        self.upep = upep
        split = self.upep.split('#')
        self.peptide = split[0]

        self.mods = []
        if len(split) == 2:
            self.mods = split[1].split(';')

        self.internal_fragments = {
            'internal b-y': {'C': -1, 'O': -1},
            'internal a-y': {},  # TODO
        }
        self.fragment_starts_forward = {
                'a': {'C': -1, 'O': -1},
                'b': {},
                'c': {'N': +1, 'H': +3},
                # 'c(-1)': {'N': +1, 'H': +2},
                # 'c(+1)': {'N': +1, 'H': +4},
                # 'c(+2)': {'N': +1, 'H': +5},
        }
        self.fragment_starts_reverse = {
                'x': {'O': 2, 'C': 1},
                'y': {'H': 2, 'O': 1},
                'z': {'O': 1, 'N': -1, 'H': 0},
                # 'z(+1)': {'O': 1, 'N': -1, 'H': 1},
                # 'z(+2)': {'O': 1, 'N': -1, 'H': 2},
                # 'z(+3)': {'O': 1, 'N': -1, 'H': 3},
        }
        self.forward = {}
        self.reverse = {}
        for start, direction in [
                (self.fragment_starts_forward, self.forward),
                (self.fragment_starts_reverse, self.reverse)]:
            direction['pos0'] = {}
            for ion_type in start.keys():
                direction['pos0'][ion_type] = \
                    [{ 'pos': 0, 'cc': ChemicalComposition(), 'mods': []}]
                direction['pos0'][ion_type][0]['cc'] += start[ion_type]

        self._induce_fragmentation()



    def _induce_fragmentation(self):
        alread_seen_frags = set()
        for i in range(len(self.peptide)):
            groups = [
                {
                    'translated_peptide_pos': i + 1,
                    'target': self.forward,
                    'aa': self.peptide[i]
                },
                {
                    'translated_peptide_pos': len(self.peptide) - i,
                    'target': self.reverse,
                    'aa': self.peptide[::-1][i]
                }
            ]

            for grp in groups:
                cc = self.upep_cc.composition_at_pos[grp['translated_peptide_pos']]
                grp['target']['pos{0}'.format(i+1)] = ddict(list)
                # extending the dicts
                # print('self.upep_cc.unimod_at_pos', self.upep_cc.unimod_at_pos)
                for neutral_loss_dict in self.neutral_losses.get(grp['aa'], [{}]):
                    neutral_loss_can_occure = False
                    required_unimods = neutral_loss_dict.get('requires_unimod', None)
                    if required_unimods is None:
                        neutral_loss_can_occure = True
                    else:
                        for required_unimod in required_unimods:
                            if required_unimod == self.upep_cc.unimod_at_pos.get(
                                    grp['translated_peptide_pos'], ''):
                                neutral_loss_can_occure = True

                    if neutral_loss_can_occure is False:
                        continue

                    for ion_type, ion_fragments in grp['target']['pos{0}'.format(i)].items():
                        for ion_frag in ion_fragments:
                            new_ion_frag = copy.deepcopy(ion_frag)
                            new_ion_frag['pos'] += 1
                            new_ion_frag['cc'] += cc
                            new_ion_frag['cc'] += neutral_loss_dict.get('cc', {})
                            mod = neutral_loss_dict.get('name',None)
                            if mod is not None:
                                new_ion_frag['mods'].append(mod)
                            new_ion_frag['name'] = '{0}{1}'.format(ion_type,new_ion_frag['pos'])
                            new_ion_frag['hill'] = new_ion_frag['cc'].hill_notation_unimod()
                            new_ion_frag['charge'] = 1
                            new_ion_frag['perdicted intensity'] = np.NAN
                            new_ion_frag['mass'] = new_ion_frag['cc']._mass()
                            new_ion_frag['mz'] = new_ion_frag['mass'] + peptide_fragmentor.PROTON
                            new_ion_frag['series'] = ion_type
                            new_ion_frag['modstring'] = ','.join(sorted(new_ion_frag['mods']))
                            _id = '{name}{0}'.format(
                                sorted(new_ion_frag['mods']),
                                **new_ion_frag
                            )
                            if _id not in alread_seen_frags:
                                grp['target']['pos{0}'.format(i+1)][ion_type].append(new_ion_frag)
                            alread_seen_frags.add(_id)

        del self.forward['pos0']
        del self.reverse['pos0']

        all_rows = []
        for direction in [self.forward, self.reverse]:
            for pos in direction.keys():
                for ion_type in direction[pos].keys():
                    all_rows += direction[pos][ion_type]

        self.df = pd.DataFrame(all_rows)

    def _clean_up_mod_string(self, mod_string=None):
        pass

    def _create_peptide_variants(self, neutral_losses=None):
        if neutral_losses is None:
            #for testing purpose
            neutral_losses = self.neutral_losses
        variations = [""]
        for pos, aa in enumerate(self.upep):
            pass





    def _expand_charges(self, row, charges):
        rows = []
        for c in charges:
            if c == 1:
                continue
            new_row = {
                'name': row['name'],
                'cc': row['cc'],
                'charge': c,
                'predicted intensity': np.NAN,
                'series': row['series'],
                'pos': row['pos'],
                'mods': row['mods'],
                'neutral loss': row['neutral loss']
            }
            mz = (row['mz'] / c) + peptide_fragmentor.PROTON
            new_row['mz'] = mz
            rows.append(new_row)
        return rows

    def _expand_neulos(self, row, neulos):
        cc = ChemicalComposition()
        all_rows = []
        # breakpoint()
        combs = [
            (x[0][1], x[1][1]) for x in combinations(
                neulos, int(len(neulos) / 2)
            ) if (len(x) > 1) and (x[0][0] != x[1][0])
        ]
        for combi in combs:
            old_cc = ChemicalComposition(f'+{row["cc"]}')
            for cc in combi:
                old_cc.add_chemical_formula(cc)
            new_row = {
                'name': row['name'],
                'cc': old_cc.hill_notation_unimod(),
                'charge': row['charge'],
                'mz': old_cc._mass() + peptide_fragmentor.PROTON,
                'predicted intensity': row['predicted intensity'],
                'series': row['series'],
                'pos': row['pos'],
                'mods': row['mods'],
                'neutral_loss': combi
            }
            all_rows.append(new_row)
        return all_rows

    def fragment_peptide( self,ion_series = None, use_neutral_loss = True):
        """Fragment `upep` and return specified ion series.

        Args:
            ion_series (tuple, optional): Ion series to create.
                should be one of a, b, c, x, y, z or internal
            use_neutral_loss (bool, optional): Use neutral losses specified
                during init
        """
        if ion_series is None:
            ion_series = ('y', 'b')

        cc = ChemicalComposition()
        full_pep = ChemicalComposition(self.upep)
        all_rows = []
        # Start with precursor mass
        for series in ion_series:
            series_correction_cc = self.fragment_starts[series]
            cc += series_correction_cc
            mods = []
            neulos = []

            for i, aa in enumerate(self.peptide):
                pos = i + 1
                if series in 'abc':
                    i += 1
                elif series in 'xyz':
                    i = len(self.peptide) - i
                elif 'internal' in series:
                    pass
                cc += full_pep.composition_at_pos[i]
                try:
                    mod = full_pep.unimod_at_pos.get(i, '')
                    mods.append(mod)
                    neulos.append(
                        (i, self.neutral_losses[aa][mod])
                    )
                    neulos.append(
                        (i, {})
                    )  # second possibility: no nl
                except KeyError:
                    pass
                row = {
                    'name': f'{series}{pos}',
                    'cc': cc.hill_notation_unimod(),
                    'charge': 1,
                    'mz': cc._mass() + peptide_fragmentor.PROTON,
                    'predicted intensity': np.NAN,
                    'series': series,
                    'pos': i,
                    'mods': ';'.join(mods),
                    'neutral loss': [],
                }
                all_rows.append(row)
                mul_charge_rows = self._expand_charges(row, self.charges)
                all_rows += mul_charge_rows
            mul_neulo_rows = self._expand_neulos(row, neulos)
            all_rows += mul_neulo_rows
        # Precursor ion
        all_rows.append(
            {
                'name': '[MH]',
                'cc': full_pep.hill_notation_unimod(),
                'charge': 1,
                'mz': full_pep._mass() + peptide_fragmentor.PROTON,
                'predicted intensity': np.NAN,
                'series': '',
                'pos': '',
                'mods': self.mods,
                'neutral loss': '',
            }
        )
        df = pd.DataFrame(all_rows)
        return df


if __name__ == '__main__':
    mains()
