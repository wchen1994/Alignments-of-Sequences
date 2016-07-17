import url_file_read as provided
import math
import random
import project4
import matplotlib.pyplot as plt

protein_human = provided.read_protein(provided.HUMAN_EYELESS_URL)
protein_fly = provided.read_protein(provided.FRUITFLY_EYELESS_URL)
scoring_matrix = provided.read_scoring_matrix(provided.PAM50_URL)

"""
caluclate the score of alignment of human and fruit fly eyeless protein
"""
def align_human_fly_protein():
    alignment_matrix = project4.compute_alignment_matrix(protein_human, protein_fly, scoring_matrix, False)
    result = project4.compute_local_alignment(protein_human, protein_fly, scoring_matrix, alignment_matrix)

    return result

#score, human_protein_aligned, fly_protein_aligned = align_human_fly_protein()
#print score

"""
calculate the ratio of similar the human of fly compare to consensus protein
"""
def calculate_similar_ratio():
    result = align_human_fly_protein()
    sequence_human = result[1].replace('-', '')
    sequence_fly = result[2].replace('-', '')
    
    protein_consensus = provided.read_protein(provided.CONSENSUS_PAX_URL)
    alignment_matrix = project4.compute_alignment_matrix(sequence_human, protein_consensus, scoring_matrix, True)
    result = project4.compute_global_alignment(sequence_human, protein_consensus, scoring_matrix, alignment_matrix)
    
    mark = 0
    for idx in range(len(result[1])):
        if result[1][idx] == result[2][idx]:
            mark += 1
    print mark / float(len(result[1]))
    
    protein_consensus = provided.read_protein(provided.CONSENSUS_PAX_URL)
    alignment_matrix = project4.compute_alignment_matrix(sequence_fly, protein_consensus, scoring_matrix, True)
    result = project4.compute_global_alignment(sequence_fly, protein_consensus, scoring_matrix, alignment_matrix)
    mark = 0
    for idx in range(len(result[1])):
        if result[1][idx] == result[2][idx]:
            mark += 1
    print mark / float(len(result[1]))

#calculate_diff_ratio()

"""
ploting histogram about score of alignment to the disorder sequence
"""
def save_dict(dictionary):
    dict_file = open("distribution.csv", "w")
    for key, value in dictionary.items():
        dict_file.write(str(key)+","+str(value)+"\n")
    dict_file.close()

def read_dict(fname):
    dict_file = open(fname, "r")
    dictionary = {}
    for line in dict_file:
        line = line.strip()
        key, value = line.split(",")
        dictionary[int(key)] = int(value)
    return dictionary

def generate_null_distribution(seq_x, seq_y, scoring_matrix, num_trials):
    distribution = {}
    bar = progressbar.ProgressBar(max_value=1000)
    for progress in range(num_trials):
        bar.update(progress)
        rand_y = list(seq_y)
        random.shuffle(rand_y)
        alignment_matrix = project4.compute_alignment_matrix(seq_x, rand_y, scoring_matrix, False)
        score = project4.compute_local_alignment(seq_x, rand_y, scoring_matrix, alignment_matrix)[0]
        distribution[score] = distribution.get(score,0) + 1
    save_dict(distribution)
    return distribution

def plot_histogram():
    READ = True
    if READ:
        dist =read_dict("distribution.csv")
    else:
        dist = generate_null_distribution(protein_human, protein_fly, scoring_matrix, 1000)

    x = dist.keys()
    y = dist.values()

    y_normal = [idx/1000.0 for idx in y]
    plt.bar(x, y_normal)
    plt.title("Null distribution using 1000 trials")
    plt.xlabel("Scores")
    plt.ylabel("Fraction of trials")
    plt.show()

#plot_histogram()

"""
calculate the mean and stdard deviation of the distribution of score over disorder sequence
"""
def cal_mean_stdv():
    score_list = []
    dist =read_dict("distribution.csv")
    for score, appearance in dist.items():
        score_list += [score] * appearance
    mean = sum(score_list)/float(len(score_list))
    stdv = math.sqrt(sum([(value - mean) ** 2 for value in score_list])/float(len(score_list)))
    return mean, stdv

#print cal_mean_stdv()

"""
Spelling Checking
"""
word_list = provided.read_words(provided.WORD_LIST_URL)

def check_spelling(checked_word, dist, word_list):
    # scoring matrix for edit distaion
    # edit distance = |x| + |y| - score(X,Y)
    # diag_socre = 2, off_diag_score = 1, dash_score = 0
    alphabets = set("abcdefghijklmnopqrstuvwxyz")
    scoring_matrix = project4.build_scoring_matrix(alphabets,2,1,0)
    string_set = set([])
    for word in word_list:
        alignment_matrix = project4.compute_alignment_matrix(checked_word ,word, scoring_matrix, True)
        score, _, _ = project4.compute_global_alignment(checked_word, word, scoring_matrix, alignment_matrix)
        score = len(checked_word) + len(word) - score
        if score <= dist:
            string_set.add(word)
    return string_set

def fast_check_spelling(checked_word, dist, word_list):
    word_set = set(word_list)
    align_set = set([])
    for word in word_set:
        if check_valid(checked_word, word, dist):
            align_set.add(word)
    return align_set

def check_valid(checked_word, word, dist):
    if dist < 0:
        return False
    elif checked_word and word:
        if checked_word[0] == word[0]:
            return check_valid(checked_word[1:], word[1:], dist)
        else:
            return check_valid(checked_word, word[1:], dist - 1) or check_valid(checked_word[1:], word, dist - 1) or check_valid(checked_word[1:], word[1:], dist - 1)
    elif checked_word:
        if dist - len(checked_word) < 0:
            return False
        else:
            return True
    elif word:
        if dist - len(word) < 0:
            return False
        else:
            return True
    else:
        return True

#print fast_check_spelling("humble", 1, word_list)
#print fast_check_spelling("firefly", 2, word_list)

