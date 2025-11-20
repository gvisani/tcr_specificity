# tcr_specificity

Model: ProteinMPNN

Input: TCR-pMHC structure, locations of CDR3s to design, locations of peptide

Output: tuples of CDR3s and peptide PWMs


Two steps:
1) Sample CDR3 sequences. Condition either on the native peptide sequence of the structure, or mask the peptide
2) For each sampled CDR3 sequence, compute the peptide PWMs

Notes:
- The second step should probably be parallelized across different GPUs
- The original peptide sequence can be scores in each PWM, and PWMs can be ordered based on this as well as entropy or off-target propensity
- For example, use the MageA3 structure as input and Titin as an off-target sequence we want to diverge from

Research questions:
- What are the characteristics of CDR3s that bind on-target (e.g. MAGEA3) but not off-target (e.g. Titin)?
- What are the characteristics of CDR3s that have high vs. low predicted specificity? Az measured by low vs. high entropy of the peptide PWM

