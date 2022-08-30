import math


def readFasta(fileName, min_seq_len):
    """
    This function is to read fasta files and returns a dict of fasta file

    Args:
        fileName (string): is the path to the fasta file
        min_seq_len (int): is a cut value for the sequences

    Return:
        fastaDict (dict): It is a dictionary whose keys are UniProt protein ids and values are protein sequences
    """
    fastaDict = {}
    with open(fileName) as fp:
        protId = ''
        sequence = ''
        for line in fp:
            if line[0] == '>':
                if len(sequence) > min_seq_len:
                    if sequence.find('U') != -1:
                        sequence = sequence.replace("U", "C")
                    fastaDict[protId] = sequence
                sequence = ''
                firstLine = line.split("|")
                protId = firstLine[1]
                continue
            line = line.strip()
            sequence = sequence + line
        if len(sequence) > min_seq_len:
            fastaDict[protId] = sequence
    fp.close()
    return fastaDict


def extractAllSubsequences(fastaDict, sub_seq_len):
    """
    This function extracts subsequences whose length is sub_seq_len

    Args:
        fastaDict (dict): It is a dictionary whose keys are UniProt protein ids and values are protein sequences
        sub_seq_len (int): is an integer that indicates the length of subsequences

    Return:
        listSubsequence (list): is a list of subsequences
    """
    listSubsequence = [fastaDict[protId][i:i + sub_seq_len]
                       for protId in fastaDict for i in range(len(fastaDict[protId]) - (sub_seq_len - 1))]
    return listSubsequence


def extractSubsequences(sequence, sub_seq_len):
    listSubsequence = [sequence[i:i + sub_seq_len]
                       for i in range(len(sequence) - (sub_seq_len - 1))]
    return listSubsequence


def blosum62_reader(fileName, alphabet):
    """
    This function is to read blosum62 matrix from a file and to store as a 20 by 20 matrix

    Args:
        fileName (string): is the path to the csv file of blosum62 matrix
        alphabet (string): is the string of combined letters of valid amino acids (20 letters: ARNDCQEGHILKMFPSTWYV)

    Return:
        blosum62Dict (dict): is the dictionary whose keys are amino acid pairs and values are blosum62 values.
    """
    lineNo = 0
    blosum62Dict = {}
    with open(fileName) as fp:
        for line in fp:
            line = line.strip()
            arrayBlo62 = line.split(",")
            pos = 0
            for letter in alphabet:
                blosum62Dict[(alphabet[lineNo], letter)] = int(arrayBlo62[pos])
                pos += 1
            lineNo = lineNo + 1
    return blosum62Dict


def PSSM_initializer(alphabet, sub_seq_len):
    """
    This function is to initialize Position Specific Scoring Matrix(PSMM) matrix of each cluster of subsequences.

    Args:
        alphabet (string): is the string of combined letters of valid amino acids (20 letters: ARNDCQEGHILKMFPSTWYV).
        sub_seq_len (int): is an integer that indicates the length of subsequences.

    Return:
        PSSM (dict): is a dictionary whose keys are the pairs of amino acid letters and their position in the
        subsequence, values are the counts of amino acids in that position.
    """
    PSSM = dict()
    for letter in alphabet:
        for i in range(sub_seq_len):
            PSSM[(letter, i)] = 0
    return PSSM


def PSSM_updater(subsequence, PSSM):
    """
    This function is to update PSSM of a cluster if subsequence is added to that cluster

    Args:
        subsequence (string): a subsequence of the protein sequence
        PSSM (dict): is a dictionary whose keys are the pairs of amino acid letters and their position in the
        subsequence, values are the counts of amino acids in that position.

    Return:
        PSSM (dict): is a dictionary whose keys are the pairs of amino acid letters and their position in the
        subsequence, values are the counts of amino acids in that position.
    """
    for i, letter in enumerate(subsequence):
        PSSM[(letter, i)] += 1
    return PSSM


def calculateMaxSimilarity(blosum62Dict, listCenters, subsequence):
    """
    This function is to find the most similar cluster to the subsequence by calculating the similarity between the
    center subsequence of the cluster and the subsequence according to blosum62 matrix.

    Args:
        blosum62Dict (dict): is the dictionary whose keys are amino acid pairs and values are blosum62 values.
        listCenters (list): is the list of the center subsequences of all clusters.
        subsequence (string): is a subsequence of protein sequence.

    Return:
        maxSimilarity  (int): maxSimilarity is the maximum similarity value between the center subsequence of the cluster
        possibleCluster (string): possibleCluster is the center subsequence that used to indicate the cluster.
    """
    maxSimilarity = -9999999
    possibleCluster = listCenters[0]
    for center in listCenters:
        similarityScore = 0
        for i in range(len(subsequence)):
            similarityScore += blosum62Dict[(center[i], subsequence[i])]
        if similarityScore > maxSimilarity:
            maxSimilarity = similarityScore
            possibleCluster = center
    return maxSimilarity, possibleCluster


def PSSM2profile(PSSM, si, sub_seq_len, alphabet):
    """
    This function is to convert PSSM to profile for each PSSM of clusters

    Args:
        PSSM:
        si:
        sub_seq_len:
        alphabet:

    Return:
    """
    log = math.log
    profile = [[log((PSSM[(letter, j)] + 0.01) / (si + 0.2))
                for j in range(sub_seq_len)]
               for letter in alphabet]
    return profile


def writeProfiles2File(profileDict, fileName, alphabet):
    """
    This function is to write all profiles to a file

    Args:
        profileDict:
        fileName:
        alphabet:

    Return:

    """
    with open(fileName, 'w') as fp:
        write = fp.write
        number = 0
        for key in profileDict:
            number += 1
            write('NUMBER OF CLUSTER : %s\n' % number)
            write('CLUSTER CENTER : %s\n' % key)
            aa_index = -1
            for item in profileDict[key]:
                aa_index = aa_index + 1
                write("%s : " % alphabet[aa_index])
                write(" %s\n" % item)
            write('\n')
    fp.close()


def read_profiles(profile_file, alphabet):
    """

    :param profile_file:
    :param alphabet:
    :return:
    """
    file_profiles = open(profile_file, "r")
    list_profiles = []
    profile_each_cluster = {}
    for line in file_profiles:
        if line == "\n":
            continue
        parts = line.split(":")
        part0 = parts[0].strip()
        part1 = parts[1].strip()
        if part0 == "NUMBER OF CLUSTER":
            continue
        if part0 in alphabet:
            list_values = []
            values = part1.split(",")
            for value in values:
                value = value.strip()
                value = value.strip("[")
                value = value.strip("]")
                list_values.append(float(value))
            profile_each_cluster[part0] = list_values
        if part0 == "V":
            list_profiles.append(profile_each_cluster)
            profile_each_cluster = {}
    return list_profiles


def write_feature_vector(filename_fv, fv_dict):
    file = open(filename_fv, "w")
    for prot_id in fv_dict:
        file.write(f"{prot_id}: ")
        for value in fv_dict[prot_id]:
            file.writelines([str(value), ","])
        file.write("\n")
    file.close()
