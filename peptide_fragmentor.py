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
            'c': {'N': 1, 'H': 3},
            'x': {'O': 2, 'C': 1},
            'y': {'H': 2, 'O': 1},
            'z': {'O': 1, 'N': 1, 'H': 1},
            'internal b-y': {'C': -1, 'O': -1},
            'internal a-y': {},
        }

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
            for i, aa in enumerate(self.peptide):
                pos = i + 1
                if series == 'b' or series == 'a':
                    i += 1
                elif series == 'y':
                    i = len(self.peptide) - i
                cc += full_pep.composition_at_pos[i]
                try:
                    mod = full_pep.unimod_at_pos[i]
                    cc += self.neutral_losses[aa][mod]
                except KeyError:
                    pass
                row = {
                    'name': f'{series}{pos}',
                    'cc': cc.hill_notation_unimod(),
                    'charge': 1,
                    'mz': cc._mass() + kb.PROTON,
                    'predicted intensity': np.NAN
                }
                all_rows.append(row)
        df = pd.DataFrame(all_rows)
        return df
