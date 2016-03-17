# encoding: utf-8
"""
    pyQms_partial
        -----

        Python module for fast and accurate mass spectrometry data quantification

        This only a small part of pyQms that allows peptide_fragmentor.py to run

        :license: MIT, see LICENSE.txt for more details

        Authors:

            * Leufken, J.
            * Niehues, A.
            * Wessel, F.
            * Sarin, L.P.
            * Hippler, M.
            * Leidel, S.A.
            * Fufezan, C.

"""
from __future__ import absolute_import
import sys
import re
import pyqms
from collections import defaultdict as ddict
import copy


class ChemicalComposition(dict):
    '''
    Chemical composition class. The actual sequence or formula can be reset
    using the `add` function.

    Keyword Arguments:
        sequence (str): Peptide or chemical formula sequence
        aa_compositions (Optional[dict]): amino acid compositions
        isotopic_distributions (Optional[dict]): isotopic distributions

    Keyword argument examples:

        * **sequence** - Currently this can for example
          be::
            >>>[
                    '+H2O2H2-OH',
                    '+{0}'.format('H2O'),
                    '{peptide}'.format(pepitde='ELVISLIVES'),
                    '{peptide}+{0}'.format('PO3', peptide='ELVISLIVES'),
                    '{peptide}#{unimod}:{pos}'.format(
                        peptide = 'ELVISLIVES',
                        unimod = 'Oxidation',
                        pos = 1
                    )
                ]

    Examples:
        >>> c = pyqms.ChemicalComposition()
        >>> c.use("ELVISLIVES#Acetyl:1")
        >>> c.hill_notation()
        'C52H90N10O18'
        >>> c.hill_notation_unimod()
        'C(52)H(90)N(10)O(18)'
        >>> c
        {'O': 18, 'H': 90, 'C': 52, 'N': 10}
        >>> c.composition_of_mod_at_pos[1]
        defaultdict(<class 'int'>, {'O': 1, 'H': 2, 'C': 2})
        >>> c.composition_of_aa_at_pos[1]
        {'O': 3, 'H': 7, 'C': 5, 'N': 1}
        >>> c.composition_at_pos[1]
        defaultdict(<class 'int'>, {'O': 4, 'H': 9, 'C': 7, 'N': 1})

        >>> c = pyqms.ChemicalComposition('+H2O2H2')
        >>> c
        {'O': 2, 'H': 4}
        >>> c.subtract_chemical_formula('H3')
        >>> c
        {'O': 2, 'H': 1}

    Note:
        We did not include mass calculation, since pyQms will do it much
        more accurately using unimod and other element enrichments.

    '''
    def __init__(self, sequence=None, aa_compositions=None,
                 isotopic_distributions=None):

        self._unimod_parser = None
        self.composition_of_mod_at_pos = {}
        """dict: chemical composition of unimod modifications at given position
        (if peptide sequence was used as input or using the `use` function)

        Note:
            Numbering starts at position 1, since all PSM search engines
            use this nomenclature.
        """
        self.composition_of_aa_at_pos = {}
        """dict: chemical composition of amino acid at given peptide position
        (if peptide sequence was used as input or using the `use` function)

        Note:
            Numbering starts at position 1, since all PSM search engines
            use this nomenclature.

            Examples::

                c.composition_of_mod_at_pos[1] = {
                    '15N': 2, '13C': 6, 'N': -2, 'C': -6
                }

        """
        self.composition_at_pos = {}
        """dict: chemical composition at given peptide position incl modifications
        (if peptide sequence was used as input or using the `use` function)

        Note:
            Numbering starts at position 1, since all PSM search engines
            use this nomenclature.
        """
        self.peptide = None
        self.addon = None
        self.unimod_at_pos = {}
        # self.regex_patterns = {
        #     ':pos' : re.compile( r''':(?P<pos>[0-9]*)''' ),
        #     'aaN'  : re.compile( r'''(?P<aa>[A-Z]{1})(?P<N>[0-9]*)''' ),
        # }
        if aa_compositions is None:
            self.aa_compositions = pyqms.knowledge_base.aa_compositions
        else:
            self.aa_compositions = aa_compositions
        if isotopic_distributions is None:
            self.isotopic_distributions = \
                pyqms.knowledge_base.isotopic_distributions
        else:
            self.isotopic_distributions = isotopic_distributions
        if sequence is not None:
            self.use( sequence )

    def __add__( self, other_cc ):
        '''Experimental'''
        tmp = copy.deepcopy( self )
        for other_key, other_value in other_cc.items():
            tmp[ other_key ] += other_value
        return tmp

    def __missing__( self, key ):
        if key not in self.keys():
            self[key] = 0
        return self[key]

    def __repr__(self):
        return self.hill_notation()

    def clear( self ):
        '''
        Resets all lookup dictionaries and self

        One class instance can be used analysing a series of sequences, thereby
        avoiding class instantiation overhead.

        Warning:

            Make sure to reset when looping over sequences and use the class.
            Chemical formulas (elemental compositions) will accumulate if not
            resetted.
        '''

        self.composition_of_mod_at_pos.clear()
        self.composition_of_aa_at_pos.clear()
        self.composition_at_pos.clear()
        self.unimod_at_pos.clear()

        self.peptide = None
        self.addon = None
        for k in list(self.keys()):
            del self[ k ]

    def _parse_sequence_old_style( self, sequence ):
        '''
        Adaptor for obsolete piqDB format.
        '''
        positions = [ len(sequence) ]
        for sign in ['+', '-']:
            if sign in sequence:
                positions.append(sequence.index(sign))
        minPos       = min(positions)
        peptide      = sequence[:minPos]
        addon        = sequence[minPos:]
        self.peptide = peptide
        self.addon   = addon
        if peptide != '':
            self.add_peptide(peptide)
            self['O'] += 1
            self['H'] += 2

        chemical_formula_blocks = re.compile(r'''
                        [+|-]{1}
                        [^-+]*
                        ''', re.VERBOSE).findall(addon)
        for cb in chemical_formula_blocks:
            if cb[0] == '+':
                self.add_chemical_formula(cb[1:])
            else:
                self.subtract_chemical_formula(cb[1:])
        return

    def _parse_sequence_unimod_style( self, sequence ):
        '''
        Sequence and modification parser in current format. Can hold also the
        modification information i.e. fixed or variable mods.

        Note:

            Sequences must not have two modifications at the same position!

        Example:

            '{peptide}#{unimod}:{pos}'.format(
                peptide = 'ELVISLIVES',
                unimod = 'Oxidation',
                pos = 1
            )
        '''
        minPos = sequence.index("#")
        peptide      = sequence[:minPos]
        addon        = sequence[minPos + 1:]
        self.peptide = peptide
        if peptide != '':
            self.add_peptide(peptide)
            self['O'] += 1
            self['H'] += 2
        self.addon = addon
        unimods    = self.addon.split(';')
        # pattern = self.regex_patterns[':pos']
        pattern    = re.compile( r''':(?P<pos>[0-9]*$)''' )
        for unimod in unimods:
            if unimod == '':
                continue
            unimod = unimod.strip()
            if ':' not in unimod:
                print(
                    'This unimod: {0} requires positional information'.format(
                        unimod
                    )
                )
                exit(1)
            for occ, match in enumerate( pattern.finditer( unimod )):
                try:
                    unimodcomposition = self._unimod_parser.name2composition(
                        unimod[ :match.start() ]
                    )
                except:
                    print(
                        'Can not map unimod {0}. extracted position argument {1}'.format(
                            unimod,
                            match.start()
                        )
                    )
                    sys.exit(1)
                # if occ >= 1:
                position = int(match.group('pos'))
                if position in self.unimod_at_pos.keys():
                    sys.exit('{0} <<- Two unimods at the same position ? '.format(
                        sequence
                    ))
                self.unimod_at_pos[ position ] = unimod[ :match.start() ]
            # match = re.search( position_re_pattern, unimod)
            # if match is not None:
            #     end = match.start()
            #     print( '>>>>', match)
            # else:
            #     end = len( unimod )
            # try:
            #     unimodcomposition = self._unimod_parser.name2composition(
            #         unimod[:end ]
            #     )
            # except:
            #     print(
            #         'Unimod error:', unimod,'>>', unimod[:end],
            #         re.search( position_re_pattern , unimod),
            #         re.search( position_re_pattern , unimod).start()
            #     )
            #     exit(1)
            # print( self , 'peptide only')
            # print( 'Unimod:', unimod, unimod[:end] , )
            # Full addition
            # print( unimodcomposition , '<<<<<<')
            for k,v in unimodcomposition.items():
                self[ k ] += v
            # storage position related modifications
            position = int(match.group('pos'))
            if position == 0:
                # E.g. Acetylation at pos 0 indicates N-Term
                # but has to be counted for position 1 in this class
                position = 1

            if position not in self.composition_of_mod_at_pos.keys():
                self.composition_of_mod_at_pos[ position ] = ddict(int)
            if position not in self.composition_at_pos.keys():
                self.composition_at_pos[ position ] = ddict(int)
            for k, v in unimodcomposition.items():
                self.composition_of_mod_at_pos[ position ][ k ] += v
                self.composition_at_pos[ position ][ k ] += v

        return

    def use(self, sequence):
        '''
        Re-initialize the class with a new sequence

        This is helpful if one ones to use the same class instance
        for multiple sequence since it remove class instantiation overhead.

        Args:
            sequence (str) - See top for possible input formats.

        Note:

            Will clear the current chemical composition dict!
        '''

        self.clear()
        # reset the shiznit
        if '#' in sequence:
            # Unimod Style format
            if self._unimod_parser is None:
                self._unimod_parser = pyqms.UnimodMapper()
            self._parse_sequence_unimod_style( sequence )
        else:
            self._parse_sequence_old_style( sequence )

    def add_chemical_formula(self, chemical_formula):
        """
        Adds chemical formula to the instance

        Chemical formula can be a string or a dictionary with the element
        count.


        For example::

            chemical_formula = 'C18H36N9O18'
            chemical_formula = {
                'C' : 18,
                'H' : 36,
                'N' : 9,
                'O' : 18
            }

        """
        self._merge(chemical_formula, mode='addition')
        return

    def add_peptide(self, peptide):
        """
            Adds peptide sequence to the instance

            Note:

                Only standard amino acids can be processed. If one uses special
                amino acids like (U of F) they have to be added to the
                knowledge_base.

        """
        # pattern = self.regex_patterns['aaN']
        pattern = re.compile( r'''(?P<aa>[A-Z]{1})(?P<N>[0-9]*)''' )
        # pattern = re.compile(r'[A-Z]{1}[0-9]*')
        number_offset = 0
        # this are the count for e.g. SILAC aa, i.e. R0 R1 or C0 and so on ...
        # print( peptide, type( peptide ))
        for aaN_match in pattern.finditer( peptide ):
            aa = aaN_match.group('aa')
            N = aaN_match.group('N')
            pos = int(aaN_match.start()) - number_offset + 1
            if N != '':
                number_offset += len( N )
            try:
                aa_compo = self.aa_compositions[ aa+N ]
            except:
                sys.exit('Do not know aa composition for {0}'.format( aa+N ))
            self.add_chemical_formula( aa_compo )

            composition = self._chemical_formula_to_dict( aa_compo )
            self.composition_of_aa_at_pos[ pos ] = composition
            if pos not in self.composition_at_pos.keys():
                self.composition_at_pos[ pos ] = ddict(int)
            for k, v in composition.items():
                self.composition_at_pos[ pos ][k] += v

    def _chemical_formula_to_dict(self, chemical_formula):
        '''
        Internal function to convert a chemical formula as string to a
        dictionary.
        '''
        if isinstance( chemical_formula, pyqms.ChemicalComposition):
            chem_dict = chemical_formula
        else:
            chem_dict = {}
            # print( chemical_formula , type( chemical_formula ))
            pattern = re.compile(r'(?P<element>[A-Z][a-z]*)(?P<count>[0-9]*)')
            for match in pattern.finditer(chemical_formula):
                if match.group('count') == '':
                    count = 1
                else:
                    count = int(match.group('count'))
                if match.group('element') not in chem_dict.keys():
                    chem_dict[match.group('element')] = 0
                chem_dict[match.group('element')] += count
        return chem_dict

    def hill_notation(self, include_ones=False, cc=None):
        '''
        Formats chemical composition into `Hill notation`_ string.

        .. _Hill Notation:
            https://en.wikipedia.org/wiki/Hill_system

        Args:
            cc (dict, optional): can format other element dicts as well.

        Returns:
            str: Hill notation format of self.
                For examples::
                    C50H88N10O17
        '''
        MAJORS = ['C', 'H']
        s = ''
        if cc is None:
            cc_dict = self
        else:
            cc_dict = cc

        for major in MAJORS:
            if major in cc_dict.keys():
                if cc_dict[major] == 0:
                    continue
                s += major
                if include_ones or cc_dict[major] > 1:
                    s += str(cc_dict[major])
        for k in sorted(cc_dict.keys()):
            if k not in MAJORS:
                if cc_dict[k] == 0:
                    continue
                s += k
                if include_ones or cc_dict[k] > 1:
                    s += str(cc_dict[k])
        return s

    # def hill_notation(self, include_ones=False):
    def hill_notation_unimod( self, cc=None ):
        '''
        Formats chemical composition into `Hill notation`_ string
        adding `unimod`_ features.

        .. _Hill Notation:
            https://en.wikipedia.org/wiki/Hill_system

        .. _unimod:
            http://www.unimod.org/fields.html

        Args:
            cc (dict, optional): can format other element dicts as well.

        Returns:
            str: Hill notation format including unimod format rules of self.
                For example::
                    C(50)H(88)N(10)O(17)
                    C(50)H(88)14N(1)N(9)(17)

        '''
        MAJORS = ['C', 'H']
        s = ''
        if cc is None:
            cc_dict = self
        else:
            cc_dict = cc

        for major in MAJORS:
            if major in cc_dict.keys():
                if cc_dict[major] == 0:
                    continue
                s += '{0}({1})'.format(
                    major.replace('(','').replace(')',''),
                    cc_dict[major]
                )
        for k in sorted(cc_dict.keys()):
            if k not in MAJORS:
                if cc_dict[k] == 0:
                    continue
                s += '{0}({1})'.format(
                    k.replace('(','').replace(')',''),
                    cc_dict[k]
                )
        return s

    def _mass(self, cc=None):
        '''
        Calculate the mass of the chemical composition.
        Optional cc can be specified, i.e. a cc dict in the style of
        { element : count , ... }

        This does however not work with enriched elements, e.g. 15N
        Rather use pyqms.isotopologue_library for more accurate mass
        calculations.
        '''
        mass = 0
        if cc is None:
            cc_mass_dict = self
        else:
            cc_mass_dict = cc
        for element, count in cc_mass_dict.items():
            mass += count * self.isotopic_distributions[element][0][0]

        return mass

    def _merge(self, chemical_formula, mode='addition'):
        '''
            Generalized function that allows addition and subtraction
        '''
        if mode == 'addition':
            sign = +1
        else:
            sign = -1
        if isinstance(chemical_formula, str):
            chemical_formula = self._chemical_formula_to_dict(chemical_formula)
        for element, count in chemical_formula.items():
            self[element] = self[element] + sign * count
        # else:
        #     print(chemical_formula, type(chemical_formula))
        #     sys.exit('Do not know the format of the chemical formula')
        return

    def subtract_peptide(self, peptide):
        '''
            Subtract peptide from instance
        '''
        for aa in peptide:
            self.subtract_chemical_formula(self.aa_compositions[aa])

    def subtract_chemical_formula(self, chemical_formula):
        '''
            Subtract chemical formula from instance
        '''
        self._merge(chemical_formula, mode='subtraction')
        return

    def generate_cc_dict( self ):
        tmp = {}
        tmp.update( self )
        return tmp


