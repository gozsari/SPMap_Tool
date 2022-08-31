import argparse

from generateProfiles import constructCluster, constructProfiles
from generateFeatureVectors import form_dict_fv
from utils import readFasta, extractAllSubsequences, blosum62_reader, writeProfiles2File, read_profiles, \
    write_feature_vector

if __name__ == "__main__":
    alphabet = "ARNDCQEGHILKMFPSTWYV"
    similarityThreshold = 8
    parser = argparse.ArgumentParser(description="arguments of SPMap feature extracter")
    parser.add_argument('--generateProfile', type=bool, default=True, help='if profile files are to be generated or '
                                                                           'not, True if the profile file needs to be '
                                                                           'generated')
    parser.add_argument('--path', type=str, default="input_folder", help='path to fasta file directory')
    parser.add_argument('--fastaFile_P', type=str, default="CYT_pos.fasta", help='fasta file name to construct profiles')
    parser.add_argument('--minSeqLen', type=int, default=20,
                        help="protein sequences less than this value will not be considered")
    parser.add_argument('--profileFile', type=str, default="CYT_pos_profile.txt", help='profile file name')
    parser.add_argument('--subSeqLen', type=int, default=5, help='the length of subsequences')
    parser.add_argument('--fastaFile_O', type=str, default="CYT_golden_positive.fasta", help='fasta file name whose features will be '
                                                                          'extracted')

    args = parser.parse_args()

    bool_profile = args.generateProfile
    path_input = args.path
    min_seq_len = args.minSeqLen
    sub_seq_len = args.subSeqLen

    if bool_profile:
        path_input = args.path
        fasta_file_profile = args.fastaFile_P
        fastaDict = readFasta(path_input + "/" + fasta_file_profile, min_seq_len)
        NUMBER_SEQS = fastaDict.__len__()
        listSubsequence = extractAllSubsequences(fastaDict, sub_seq_len)
        blo62Dict = blosum62_reader('blo62.csv', alphabet)
        clustersDict = constructCluster(alphabet, sub_seq_len, blo62Dict, listSubsequence, similarityThreshold)
        profileDict = constructProfiles(clustersDict, NUMBER_SEQS, sub_seq_len, alphabet)
        profileFile = "profiles/" + fasta_file_profile.split('.')[0] + "_profile.txt"
        writeProfiles2File(profileDict, profileFile, alphabet)
    else:
        profileFile = "profiles/" + args.profileFile

    fasta_file_out = args.fastaFile_O
    fasta_dict = readFasta(path_input + "/" + fasta_file_out, min_seq_len)
    list_profiles = read_profiles(profileFile, alphabet)
    fv_dict = form_dict_fv(list_profiles, fasta_dict, sub_seq_len)
    filename_fv = "output_folder/" + fasta_file_out.split('.')[0] + "_SPMap.csv"
    write_feature_vector(filename_fv, fv_dict)
