import project4

def test_build_scoring_matrix():
    alphabet = {'A', 'U', 'C', 'G'}
    matrix = project4.build_scoring_matrix(alphabet, 10, 4, -6)
    assert matrix['A']['A'] == 10
    assert matrix['A']['U'] == 4
    assert matrix['A']['-'] == -6
    assert matrix['-']['-'] == -6

test_build_scoring_matrix()
