"""Knowledge base for peptide fragmentor.

Attributes:
    neutral_losses (dict): Description
    PROTON (float): Mass of a proton in dalton

"""
PROTON = 1.007276466583

# keep it alphabetically sorted
neutral_losses = {
    'A' : [{}],
    'C' : [{}],
    'D': [
        {
            'name': '-H2O',
            'cc' : {'H': -2, 'O': -1}
        },
        {}
    ],
    'E': [
        {
            'name': '-H2O',
            'cc' : {'H': -2, 'O': -1}
        },
        {}
    ],
    'F' : [{}],
    'G' : [{}],
    'H':[
        {
        'name': '+H2O',
        'cc': {'H': +2, 'O': +1},
        'available_in_series': ['b']
        },
        {}
    ],
    'I': [{}],
    'K': [
        {
            'name': '-NH3',
            'cc': {'N': -1, 'H': -3},
        },
        {
            'name': '+H2O',
            'cc': {'H': +2, 'O': +1},
            'available_in_series': ['b']
        },
        {}
    ],
    'L': [{}],
    'M':[
        {
            'name': '-SOCH4',
            'requires_unimod': ['Oxidation'],
            'cc': {'H': -4, 'C': -1, 'O': -1, 'S': -1},
        },
        {} # Additional empty dict {} indicates loss is optional
    ],
    'N': [
        {
            'name': '-NH3',
            'cc': {'N': -1, 'H': -3},
        },
        {}
    ],
    'P': [{}],
    'Q': [
        {
            'name': '-NH3',
            'cc': {'N': -1, 'H': -3},
        },
        {}
    ],
    'R': [
        {
            'name': '-NH3',
            'cc': {'N': -1, 'H': -3},
        },
        {
            'name': '+H2O',
            'cc': {'H': +2, 'O': +1},
            'available_in_series': ['b']
        },
        {}
    ],
    'S': [
        {
            'name': '-P',
            'requires_unimod': ['Phospho'],
            'cc': {'H': -3, 'O': -4, 'P': -1},
        },
        {
            'name': '-H2O',
            'requires_unimod': [''],
            'cc' : {'H': -2, 'O': -1}
        },
        {} # Additional empty dict {} indicates loss is optional
    ],
    'T': [
        {
            'name': '-P',
            'requires_unimod': ['Phospho'],
            'cc': {'H': -3, 'O': -4, 'P': -1},
        },
        {
            'name': '-H2O',
            'cc' : {'H': -2, 'O': -1}
        },
        {} # Additional empty dict {} indicates loss is optional
    ],
    'V': [{}],
    'W': [{}],
    'Y': [
        {
            'name': '-P',
            'requires_unimod': ['Phospho'],
            'cc': {'H': -3, 'O': -4, 'P': -1},
        },
        {} # Additional empty dict {} indicates loss is optional
    ],
}

"""
A 71.037113785
C 103.00918495900001
D 115.026943025
E 129.0425930894
F 147.06841391380001
G 57.021463720599996
H 137.0589118574
I 113.0840639782
K 128.09496301439998
L 113.0840639782
M 131.0404850878
N 114.04292744119999
P 97.0527638494
Q 128.0585775056
R 156.1011110224
S 87.032028405
T 101.04767846940001
V 99.0684139138
W 186.07931295
Y 163.06332853380002
"""