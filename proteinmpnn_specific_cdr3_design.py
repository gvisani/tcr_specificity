

import os
import argparse

from pdb_utils import get_cdr3_resids
from proteinmpnn_class import ProteinMPNN

'''

class ProteinMPNN(object):

    def __init__(
            self,
            model_name: str,
            proteinmpnn_path: str = "/gscratch/spe/gvisan01/ProteinMPNN-copy/"
        ) -> None:
        self.model_name = model_name
        self.proteinmpnn_path = proteinmpnn_path
    
    def run(
            self,
            output_dir: str,
            pdbfile: Optional[str] = None, # either this or pdbdir can be specified
            pdbdir: Optional[str] = None, # either this or pdbfile can be specified
            chains_or_positions_to_design_or_score: Optional[Union[List[str], Dict[str, int]]] = None, # specifies positions to design, assuming only non-empty icodes are to be designed! if None, everything is designed
            fix_pdbfile_offset: bool = True,
            conditional_probs_only: bool = False,
            num_seq_per_target: int = 10,
            sampling_temp: float = 0.1,
            batch_size: int = 1,
        ) -> None:

'''


def sample_cdr3s_with_proteinmpnn(tcr_pmhc_pdbfile: str, mhc_class: int = 1, organism: str = 'human'):
    
    ## 1) get resids of the cdr3s, form them into the dictionary of design specifics that will be used by the ProteinMPNN class
    (cdr3a_resids, cdr3a_seq), (cdr3b_resids, cdr3b_seq) = get_cdr3_resids(tcr_pmhc_pdbfile, mhc_class, organism)

    # check that there aren't any insertion codes since we can't design those yet (just a limit of the implementation)
    positions_to_design = {}
    for resids in [cdr3a_resids, cdr3b_resids]:
        chains, resnums_and_icodes_str = zip(*resids)
        unique_chains = list(set(chains))
        assert len(unique_chains) == 1
        chain = unique_chains[0]
        resnums = []
        for resnum_and_icode_str in resnums_and_icodes_str:
            assert resnum_and_icode_str[-1] == " ", f"non-empty icode found in CDR3 - we can't design this yet unfortunately"
            resnum = int(resnum_and_icode_str.strip(" "))
            resnums.append(resnum)
        positions_to_design[chain] = sorted(resnums)

    ## 2) run proteinmpnn in sampling mode!
    mpnn = ProteinMPNN('v_48_002')

    mpnn.run(
            'test_output_5brz',
            pdbfile = tcr_pmhc_pdbfile,
            chains_or_positions_to_design_or_score = positions_to_design,
            fix_pdbfile_offset = True,
            conditional_probs_only = False,
            num_seq_per_target = 10,
            sampling_temp = 0.1,
            batch_size = 10,
    )


if __name__ == '__main__':
    sample_cdr3s_with_proteinmpnn("./wt_pdbs/5brz.pdb")



