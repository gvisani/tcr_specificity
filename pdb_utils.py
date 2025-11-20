
import os
import sys

# PATH_TO_TCRDOCK = "/gscratch/spe/gvisan01/TCRdock-copy"
# sys.path.append(PATH_TO_TCRDOCK)

import tcrdock


def get_tcr_pmhc_info(pdbfile, mhc_class = 1, organism = 'human'):

    pose = tcrdock.pdblite.pose_from_pdb(pdbfile)

    num_chains = len(pose['chains'])
    if mhc_class==1:
        if num_chains == 5:
            # remove B2M
            print(f'removing chain 2 from a 5-chain MHC class I pose; residue numbers '
                  'in parsing output will not include this chain')
            pose = tcrdock.pdblite.delete_chains(pose, [1]) # 0-indexed chain number
            num_chains = len(pose['chains'])
        else:
            assert num_chains==4, \
                f'MHC-I pdbfile {pdbfile} should have 4 or 5 chains, see --help message'
        cs = pose['chainseq'].split('/')
        mhc_aseq, pep_seq, tcr_aseq, tcr_bseq = cs
        mhc_bseq = None
    else:
        assert num_chains==5, \
            f'MHC-II pdbfile {pdbfile} should have 5 chains, see --help message'
        cs = pose['chainseq'].split('/')
        mhc_aseq, mhc_bseq, pep_seq, tcr_aseq, tcr_bseq = cs

    tdinfo = tcrdock.tcrdock_info.TCRdockInfo().from_sequences(
        organism, mhc_class, mhc_aseq, mhc_bseq, pep_seq, tcr_aseq, tcr_bseq).to_dict()
    
    return pose, tdinfo


def get_cdr3_resids(pdbfile, mhc_class = 1, organism = 'human'):
    '''
    note: resids contain insertion codes as well
    '''

    pose, tdinfo = get_tcr_pmhc_info(pdbfile, mhc_class=mhc_class, organism=organism)

    cdr3a_seq = tdinfo["tcr"][0][2]
    cdr3b_seq = tdinfo["tcr"][1][2]

    num_cdrs = len(tdinfo["tcr_cdrs"]) # always a multiple of 2!
    cdr3a_loc = tdinfo["tcr_cdrs"][num_cdrs // 2 - 1]
    cdr3b_loc = tdinfo["tcr_cdrs"][num_cdrs - 1]

    ## sanity check
    assert pose["sequence"][cdr3a_loc[0] : cdr3a_loc[1] + 1] == cdr3a_seq
    assert pose["sequence"][cdr3b_loc[0] : cdr3b_loc[1] + 1] == cdr3b_seq

    cdr3a_resids = pose["resids"][cdr3a_loc[0] : cdr3a_loc[1] + 1]
    cdr3b_resids = pose["resids"][cdr3b_loc[0] : cdr3b_loc[1] + 1]

    return (cdr3a_resids, cdr3a_seq), (cdr3b_resids, cdr3b_seq)




if __name__ == '__main__':

    get_cdr3_resids("wt_pdbs/5brz.pdb")
    get_cdr3_resids("/gscratch/spe/gvisan01/tcr_pmhc/all_pdbs/1ao7.pdb")

