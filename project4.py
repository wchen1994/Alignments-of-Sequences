"""
Algorithmic Thinking Project 4
Matching RNA or DNA sequences
"""

def build_scoring_matrix(alphabet, diag_score, off_diag_score, dash_score):
    """
    [Input]
        alphabet: set of characters
        diag_score: score for diagonal
        off_diag_score: score for not in diagonal
        dash_score: score for match with a dash
    [Output]
        dictionary of dictionaies
    """
    alphabet = alphabet.copy()
    alphabet.add('-')
    dict1 = {}
    for letter1 in alphabet:
        dict1[letter1] = {}
        for letter2 in alphabet:
            if letter1 == '-' or letter2 == '-':
                dict1[letter1][letter2] = dash_score
            elif letter1 == letter2:
                dict1[letter1][letter2] = diag_score
            else:
                dict1[letter1][letter2] = off_diag_score
    return dict1
