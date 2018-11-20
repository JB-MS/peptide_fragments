"""Knowledge base for peptide fragmentor.

Attributes:
    neutral_losses (dict): Description
    PROTON (float): Mass of a proton in dalton

"""
PROTON = 1.007276466583


neutral_losses = [
    {
        'pre': 'S#Phospho',
        'post': 'S-H',
        'name': 'Phospor neutral loss',
        'mode': 'opt',
    },
    {
        'pre': 'K',
        'post': 'K+H2O',
        'name': 'Water addition Lysin',
        'mode': 'fix',
    }
]

neutral_losses = {
    'S': {
        'Phospho': {'H': -1}
    },
    'M': {
        'Oxidation': {
            'H': 4,
            'C': 1,
            'O': 1,
            'S': 1
        }
    },
}
