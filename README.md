# SPMap: Subsequence-based feature map for protein function classification
SPMap takes into account the information coming from the subsequences of a protein. A group of protein sequences that belong to the same level of classification is decomposed into fixed-length subsequences and they are clustered to obtain a representative feature space mapping. Mapping is defined as the distribution of the subsequences of a protein sequence over these clusters. The resulting feature space representation is used to train discriminative classifiers for functional families. The aim of this approach is to incorporate information coming from important subregions that are conserved over a family of proteins while avoiding the difficult task of explicit motif identification. 

## SPMap tool:
It is sequence-based feature extraction tool based on the subsequences profiles obtained from trainin data (in fasta format).

1. Users should first form a profile from the respective training dataset.
2. Users can then extract protein features using the profile(s).

### How to use:
```
 python runSPMap.py --generateProfile True --path 'input_folder' --fastaFile_P CYT_pos.fasta --minSeqLen 20 --subSeqLen 5 --fastaFile_O CYT_golden_positive.fasta
```
Table 1: SPMap tool's arguments

| Arguments                               | Description                                                                                 | Values
-----------------------------------------|---------------------------------------------------------------------------------------------|---------------------
 generateProfile                    | If profile files are to be generated or not, True if the profile file needs to be generated |True or False
 path                  | path to fasta file directory                                                                | default: "input_folder"
fastaFile_P | fasta file name to construct profiles (fasta file of training data)                         | default: "CYT_pos.fasta"
minSeqLen      | protein sequences shorter than this value will not be considered                            | default: 20
profileFile |profile file name to be generated |default: "CYT_pos_profile.txt"
subSeqLen | the length of subsequences | default: 5
fastaFile_O | fasta file name whose features will be extracted | default: "CYT_golden_positive.fasta"
<br/>


## References
Sarac, O. S., Gürsoy-Yüzügüllü, Ö., Cetin-Atalay, R., & Atalay, V. (2008). Subsequence-based feature map for protein function classification. Computational biology and chemistry, 32(2), 122-130.
