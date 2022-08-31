from utils import *


def calculate_fv_sequence(list_subsequence, list_profiles):
    """
    This function calculates the feature values of a protein's subsequences

    Args:
        list_subsequence (list): the list of subsequences of a protein sequence
        list_profiles (list): list of the profile dictionaries

    Return:
        list_fv_sequence (list): the list of the calculated features for a protein sequence
    """
    list_fv_sequence = list()
    for profile in list_profiles:
        max_value = -99999999
        for subsequence in list_subsequence:
            if "X" in subsequence:
                continue
            sum_value = 0
            for i in range(len(subsequence)):
                sum_value += profile[subsequence[i]][i]
            if sum_value > max_value:
                max_value = sum_value
        if max_value < -15:
            max_value = 0.0
        else:
            max_value = math.exp(max_value)
        list_fv_sequence.append(max_value)
    return list_fv_sequence


def form_dict_fv(list_profiles, fasta_dict, sub_seq_len):
    """
    This function is to construct a dictionary whose keys are protein ids and values are SPMap feature vectors.

    Args:
        list_profiles (list): list of the profile dictionaries
        fasta_dict (dict): is dictionary of fasta file that we are to construct SPMap feature vectors.
                    (keys: protein ids, values: sequences)
        sub_seq_len (int): is an integer that indicates the length of subsequences

    Return:
        dict_fv (dict): is a dictionary whose keys are protein ids and values are SPMap feature vectors.
    """

    dict_fv = {}
    for prot_id in fasta_dict:
        list_subsequence = extractSubsequences(fasta_dict[prot_id], sub_seq_len)
        list_fv_sequence = calculate_fv_sequence(list_subsequence, list_profiles)
        dict_fv[prot_id] = list_fv_sequence
    return dict_fv
