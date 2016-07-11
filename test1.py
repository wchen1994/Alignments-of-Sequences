import project4

def test_build_scoring_matrix():
    alphabet = {'A', 'U', 'C', 'G'}
    matrix = project4.build_scoring_matrix(alphabet, 10, 4, -6)
    assert matrix['A']['A'] == 10
    assert matrix['A']['U'] == 4
    assert matrix['A']['-'] == -6
    assert matrix['-']['-'] == -6

def test_compute_alignment_matrix():
    alphabet = {'A', 'T', 'C', 'G'}
    score_matrix = project4.build_scoring_matrix(alphabet, 6, 2, -4)
    align_matrix = project4.compute_alignment_matrix(['A'], ['A'], score_matrix, False)
    assert align_matrix == [[0,0],[0,6]]

def print_matrix(matrix):
    for row in matrix:
        print row

test_build_scoring_matrix()
test_compute_alignment_matrix()
