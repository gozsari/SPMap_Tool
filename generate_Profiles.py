from utils import *



def constructCluster(alphabet, sub_seq_len, blo62Dict, listSubsequence, similarityThreshold):
    """
    This function is to construct clusters

    Args:
        alphabet:
        sub_seq_len:
        blo62Dict:
        listSubsequence:

    Return:
    """

    similarityThreshold = 8
    PSSM = PSSM_initializer(alphabet, sub_seq_len)
    clustersDict = dict()
    subsequenceCount = 1
    PSSM = PSSM_updater(listSubsequence[0], PSSM)
    clustersDict[listSubsequence[0]] = (subsequenceCount, PSSM)
    for i in range(1, listSubsequence.__len__()):
        if "X" in listSubsequence[i] or "B" in listSubsequence[i] or "Z" in listSubsequence[i]:
            continue
        listCenters = list(clustersDict.keys())
        maxSimilarityScore, possibleCluster = calculateMaxSimilarity(blo62Dict, listCenters, listSubsequence[i])

        if maxSimilarityScore >= similarityThreshold:
            clustersDict[possibleCluster] = (clustersDict[possibleCluster][0] + 1,
                                             PSSM_updater(listSubsequence[i],
                                             clustersDict[possibleCluster][1]))
        else:
            PSSM = PSSM_initializer(alphabet, sub_seq_len)
            clustersDict[listSubsequence[i]] = (subsequenceCount, PSSM_updater(listSubsequence[i], PSSM))
    return clustersDict

def constructProfiles(clustersDict, NUMBER_SEQS, sub_seq_len, alphabet):
    """
    This function is to create profiles and takes cluster dictionary and number of sequences as parameter

    Args:
        clustersDict:
        NUMBER_SEQS:
        sub_seq_len:
        alphabet:

    Return:

    """
    profileDict = {}
    for key in clustersDict:
        si = clustersDict[key][0]
        if si <= (NUMBER_SEQS*0.1):
            continue
        PSSM = clustersDict[key][1]
        profileDict[key] = PSSM2profile(PSSM, si, sub_seq_len, alphabet)
    return profileDict



