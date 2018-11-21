"""Knowledge base for peptide fragmentor.

Attributes:
    neutral_losses (dict): Description
    PROTON (float): Mass of a proton in dalton

"""
PROTON = 1.007276466583

# first level == aa, second_level == PSI-MS mod name
neutral_losses = {
    'S': {
        'Phospho': {'H': -3, 'O': -4, 'P': -1},
        '': {},
    },
    'R': {
        '': {'N': -1, 'H': -3},
    },
    'K': {
        '': {'N': -1, 'H': -3},
    },
    'Q': {
        '': {'N': -1, 'H': -3},
    },
    'N': {
        '': {'N': -1, 'H': -3},
    },

    # 'M': {
    #     'Oxidation': {
    #         'H': -4,
    #         'C': -1,
    #         'O': -1,
    #         'S': -1
    #     }
    # },
}
