#!/usr/bin/env python3.4
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
version_info  = (0, 0, 0, 'alpha')
version = '0.0.0-alpha'

from . import knowledge_base
from .chemical_composition import ChemicalComposition
from .unimod_mapper import UnimodMapper
