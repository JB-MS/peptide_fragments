#!/usr/bin/env python3
from itertools import combinations

import pandas as pd
import numpy as np
from pyqms.chemical_composition import ChemicalComposition

import knowledge_base as kb
from knowledge_base import neutral_losses as kb_neutral_losses


class PeptideFragment0r:
    def __init__(
        self,
        upep: str,
        charges: list = None,
        neutral_losses: list = None
    ) -> None:
        """Initialze fragmentor with peptide `upep`.

        Args:
            upep (str): Description
            charges (list, optional): Description
            neutral_losses (list, optional): Description
        """
        if charges is None:
            self.charges = [1, 2, 3]
        else:
            self.charges = charges
        if neutral_losses is None:
            neutral_losses = kb_neutral_losses
        else:
            neutral_losses += kb_neutral_losses
        self.upep = upep
        split = self.upep.split('#')
        self.peptide = split[0]
        self.mods = []
        if len(split) == 2:
            self.mods = split[1].split(';')

        self.neutral_losses = neutral_losses
        self.fragment_starts = {
            'a': {'C': -1, 'O': -1},
            'b': {},
            'c': {'N': +1, 'H': +3},
            'x': {'O': 2, 'C': 1},
            'y': {'H': 2, 'O': 1},
            'z': {'O': 1, 'N': -1, 'H': -1},
            'internal b-y': {'C': -1, 'O': -1},
            'internal a-y': {},  # TODO
        }

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
                'neutral_loss': row['neutral_loss']
            }
            mz = (row['mz'] / c) + kb.PROTON
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
                'mz': old_cc._mass() + kb.PROTON,
                'predicted intensity': row['predicted intensity'],
                'series': row['series'],
                'pos': row['pos'],
                'mods': row['mods'],
                'neutral_loss': combi
            }
            all_rows.append(new_row)
        return all_rows

    def fragment_peptide(
        self,
        ion_series: tuple = None,
        use_neutral_loss: bool = True
    ) -> pd.DataFrame:
        """Fragment `upep` and return specified ion series.

        Args:
            ion_series (tuple, optional): Ion series to create.
                should be one of a, b, c, x, y, z or interal
            use_neutral_loss (bool, optional): Use neutral losses speficied
                during init
        """
        if ion_series is None:
            ion_series = ('y', 'b')

        cc = ChemicalComposition()
        full_pep = ChemicalComposition(self.upep)
        all_rows = []
        for series in ion_series:
            series_correction_cc = self.fragment_starts[series]
            cc += series_correction_cc
            mods = []
            neulos = []

            # calc neutral_loss combinations

            for i, aa in enumerate(self.peptide):
                pos = i + 1
                if series in 'abc':
                    i += 1
                elif series in 'xyz':
                    i = len(self.peptide) - i
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
                    'mz': cc._mass() + kb.PROTON,
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

        df = pd.DataFrame(all_rows)
        return df


if __name__ == '__main__':
    pep = 'ELVISLIVES#Acetyl:0'
    fragger = PeptideFragment0r(pep)
    df = fragger.fragment_peptide(
        ion_series=('a', 'b', 'y')
    )
    print(df.set_index(['series', 'pos']))
