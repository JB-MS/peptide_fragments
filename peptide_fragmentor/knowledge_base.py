"""Knowledge base for peptide fragmentor.

Attributes:
    neutral_losses (dict): Description
    PROTON (float): Mass of a proton in dalton

"""
PROTON = 1.007276466583

# first level == aa, second_level == PSI-MS mod name
neutral_losses = {
    'S': [
        {
            'name': 'P',
            'requires_unimod': ['Phospho'],
            'cc': {'H': -3, 'O': -4, 'P': -1},
        },
        {
            'name': 'H2O',
            'requires_unimod': [''],
            'cc' : {'H': -2, 'O': -1}
        },
        {} # Additional empty dict {} indicates loss is optional
    ],
    'R': [
        {
            'name': 'NH3',
            'cc': {'N': -1, 'H': -3},
        },
        {}
    ],
    'K': [
        {
            'name': 'NH3',
            'cc': {'N': -1, 'H': -3},
        },
        {}
    ],
    'Q': [
        {
            'name': 'NH3',
            'cc': {'N': -1, 'H': -3},
        },
        {}
    ],
    'N': [
        {
            'name': 'NH3',
            'cc': {'N': -1, 'H': -3},
        },
        {}
    ],
    'T': [
        {
            'name': 'H2O',
            'cc' : {'H': -2, 'O': -1}
        },
        {}
    ],
    'D': [
        {
            'name': 'H2O',
            'cc' : {'H': -2, 'O': -1}
        },
        {}
    ],
    'E': [
        {
            'name': 'H2O',
            'cc' : {'H': -2, 'O': -1}
        },
        {}
    ],

    # 'M': {
    #     'Oxidation': {
    #         'H': -4,
    #         'C': -1,
    #         'O': -1,
    #         'S': -1
    #     }
    # },
}
