#!/usr/bin/env python3
import peptide_fragmentor




if __name__ == '__main__':
    pep = 'ELVISLIVES#Acetyl:0'
    print(dir(peptide_fragmentor))
    fragger = peptide_fragmentor.PeptideFragment0r(pep)
    df = fragger.fragment_peptide(
        ion_series=('a', 'b', 'y')
    )
    print(df.set_index(['series', 'pos']))
