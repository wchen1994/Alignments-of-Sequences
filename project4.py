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

def compute_alignment_matrix(seq_x, seq_y, scoring_matrix, global_flag):
    """
    compute the score for every alignment of seq_x and seq_y according to
    the scoring matrix.
    if global_flag == True:
        can score negative number. It is useful for global alignment
    if global_flag == False:
        every socre less than 0 will become 0. It is useful for local
        alignment
    """
    len_seq_x = len(seq_x)
    len_seq_y = len(seq_y)
    lowest_score = float('-inf') if global_flag else 0

    # initialize matrix[0][0] = 0
    align_matrix = [[float('nan') for dummy_col in range(len_seq_y + 1)] 
                            for dummy_row in range(len_seq_x +1)]
    align_matrix[0][0] = 0

    # initialize alignment of '-' and seq_y, or matrix first row
    for col in range(1, len_seq_y + 1):
        score = align_matrix[0][col-1] + scoring_matrix[seq_y[col-1]]['-']
        align_matrix[0][col] = max(score, lowest_score)

    # initialize alignment of seq_x and '-', or matrix first column
    for row in range(1, len_seq_x + 1):
        score = align_matrix[row-1][0] +  scoring_matrix['-'][seq_x[row-1]]
        align_matrix[row][0] = max(score, lowest_score)

    # fill in the rest according to the previous calculation
    for row in range(1, len_seq_x + 1):
        for col in range(1, len_seq_y + 1):
            score_list = [lowest_score]
            score_list.append( align_matrix[row-1][col-1] + scoring_matrix[seq_x[row-1]][seq_y[col-1]] )
            score_list.append( align_matrix[row][col-1] + scoring_matrix['-'][seq_y[col-1]] )
            score_list.append( align_matrix[row-1][col] + scoring_matrix[seq_x[row-1]]['-'] )
            align_matrix[row][col] = max(score_list)

    return align_matrix

def compute_global_alignment(seq_x, seq_y, scoring_matrix, alignment_matrix):
    """
    [Input]
        seq_x, seq_y: the two different sequence to be aligned
        scoring _matrix: the dictionary of dictionary of the score
        alignment_matrix: the matrix of score of alignment
    [Output]
        tuple(score, align_x, align_y)
        score: the score of the alignment
        align_x, align_y: the aligned sequence
    """
    idx_x = len(seq_x)
    idx_y = len(seq_y)
    align_x = ""
    align_y = ""
    score = 0

    while idx_x !=0 and idx_y != 0:
        # match is alphabet and alphabet
        if alignment_matrix[idx_x][idx_y] == alignment_matrix[idx_x-1][idx_y-1] + scoring_matrix[seq_x[idx_x-1]][seq_y[idx_y-1]]:
            align_x = seq_x[idx_x-1] + align_x
            align_y = seq_y[idx_y-1] + align_y
            score += scoring_matrix[seq_x[idx_x-1]][seq_y[idx_y-1]]
            idx_x -= 1
            idx_y -= 1

        # match is alphabet and '-'
        elif alignment_matrix[idx_x][idx_y] == alignment_matrix[idx_x-1][idx_y] + scoring_matrix[seq_x[idx_x-1]]['-']:
            align_x = seq_x[idx_x-1] + align_x
            align_y = '-' + align_y
            score += scoring_matrix[seq_x[idx_x-1]]['-']
            idx_x -= 1

        # match is '-' and alphabet
        elif alignment_matrix[idx_x][idx_y] == alignment_matrix[idx_x][idx_y-1] + scoring_matrix['-'][seq_y[idx_y-1]]:
            align_x = '-' + align_x
            align_y = seq_y[idx_y-1] + align_y
            score += scoring_matrix['-'][seq_y[idx_y-1]]
            idx_y -= 1

    # append '-' on seq_y until same length
    while idx_x != 0:
        align_x = seq_x[idx_x-1] + align_x
        align_y = '-' + align_y
        score += scoring_matrix[seq_x[idx_x-1]]['-']
        idx_x -= 1

    # append '-' on seq_x until same length
    while idx_y != 0:
        align_x = '-' + align_x
        align_y = seq_y[idx_y-1] + align_y
        score += scoring_matrix['-'][seq_y[idx_y-1]]
        idx_y -= 1

    return (score, align_x, align_y)

