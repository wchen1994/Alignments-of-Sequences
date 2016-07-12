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
    align_matrix = project4.compute_alignment_matrix('A', 'A', score_matrix, False)
    assert align_matrix == [[0,0],[0,6]]

def test_compute_global_alignment():
    alphabet = {'A', 'T', 'C', 'G'}
    score_matrix = project4.build_scoring_matrix(alphabet, 6, 2, -4)

    result = project4.compute_global_alignment('', '', score_matrix, [[0]])
    assert result == (0,'','')

    result = project4.compute_global_alignment('A', 'A', score_matrix, [[0, -4], [-4, 6]]) 
    assert result == (6,'A','A')

    seq_x = 'ATG'
    seq_y = 'ACG'
    align_matrix = project4.compute_alignment_matrix(seq_x,seq_y,score_matrix, True)
    result = project4.compute_global_alignment(seq_x,seq_y,score_matrix,align_matrix)
    assert result == (14,'ATG','ACG')
    

def print_matrix(matrix):
    for row in matrix:
        print row

test_build_scoring_matrix()
test_compute_alignment_matrix()
test_compute_global_alignment()
