

import os
import numpy as np
from copy import deepcopy
import argparse
from typing import Dict, List, Tuple, Optional, Union
import time

from scipy.special import softmax

from pdb_utils import get_cdr3_resids
from seq_utils import read_fasta, get_designed_chains, get_fixed_chains, get_field
from utils import read_jsonl, write_jsonl
from proteinmpnn_class import ProteinMPNN

## as used by ProteinMPNN
IDX_TO_AA = 'ACDEFGHIKLMNPQRSTVWYX'
AA_TO_IDX = dict(zip(IDX_TO_AA, range(21)))


def sample_cdr3s_with_proteinmpnn(
    output_dir: str,
    num_cdr3s_to_sample: int,
    model_name: str, # e.g. v_48_002, v_48_002, v_48_030
    tcr_pmhc_pdbfile: str,
    mhc_class: int = 1,
    organism: str = 'human',
    sampling_temp: float = 0.1,
    batch_size: int = 10,
) -> None:
    
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
    mpnn = ProteinMPNN(model_name)

    mpnn.run(
            output_dir,
            pdbfile = tcr_pmhc_pdbfile,
            chains_or_positions_to_design_or_score = positions_to_design,
            fix_pdbfile_offset = True,
            conditional_probs_only = False, # sampling mode
            num_seq_per_target = num_cdr3s_to_sample,
            sampling_temp = sampling_temp,
            batch_size = batch_size,
    )

    ## TODO: save the cdr3 locations somewhere maybe? just need a way to extract them easily



def build_peptide_proteinmpnn_inputs_with_different_cdr3s_bypassing_parsing(
        pdbid: str,
        cdr3_output_dir: str,
        peptide_output_dir: str,
        peptide_chain: str
    ) -> None:
    '''
    Instead of changing the sequence in the structure, and then parsing it with the proteinmpnn parser,
    we're going to simply create new jsonl files with multiple versions of the TCR-pMHC structure,
    each version containing a different cdr3 sequence.
    Note: this only works if we assume no changes in the structure.
    '''

    os.makedirs(os.path.join(peptide_output_dir, pdbid), exist_ok=True)

    ## 1) read the designs --> full chains but with CDR3s differing, makes it easy actually

    samples_fastapath = os.path.join(cdr3_output_dir, pdbid, "seqs", f"{pdbid}.fa")

    header_and_seq_list = read_fasta(samples_fastapath)

    # first item is the original sequence, use the header to extract names of the designed chains
    first_header = header_and_seq_list[0][0]
    fixed_chains = get_fixed_chains(first_header)
    designed_chains = get_designed_chains(first_header)
    all_chains = sorted(fixed_chains + designed_chains) # sorting is key I think
    
    headers = [tup[0] for tup in header_and_seq_list[1:]]
    seq_samples = [tup[1] for tup in header_and_seq_list[1:]]

    # create sample_ids
    sample_id_list = [f"cdr3_sample_{get_field(header, 'sample')}" for header in headers]

    # this is to effectively remove duplicates
    seq_to_sample_id = {seq: sample_id for sample_id, seq in zip(sample_id_list, seq_samples)}

    # make dict[sample_id -> dict[chain -> seq]], for ease of creating parsed_pdbs.jsonl
    sample_id_to_chain_to_seq = {sample_id: dict(zip(designed_chains, seq.split('/'))) for seq, sample_id in seq_to_sample_id.items()}


    ## 2) create parsed_pdbs.jsonl
    
    # extract base data
    parsed_pdbs_base_list = read_jsonl(os.path.join(os.path.join(cdr3_output_dir, pdbid, "parsed_pdbs.jsonl")))
    assert len(parsed_pdbs_base_list) == 1

    parsed_pdb_base_dict = parsed_pdbs_base_list[0]

    parsed_pdbs_list = []
    for sample_id, chain_to_seq in sample_id_to_chain_to_seq.items():
        parsed_pdb_dict = {}

        # set the name with sample_id
        parsed_pdb_dict["name"] = f"{pdbid}_{sample_id}"

        # copy over num_of_chains
        parsed_pdb_dict["num_of_chains"] = deepcopy(parsed_pdb_base_dict["num_of_chains"])

        # copy over all the coords
        for chain in all_chains:
            parsed_pdb_dict[f"coords_chain_{chain}"] = deepcopy(parsed_pdb_base_dict[f"coords_chain_{chain}"])
        
        # copy over the seq_chain of fixed_chains
        for chain in fixed_chains:
            parsed_pdb_dict[f"seq_chain_{chain}"] = deepcopy(parsed_pdb_base_dict[f"seq_chain_{chain}"])
        
        # put in the designed_chains
        for chain in designed_chains:
            parsed_pdb_dict[f"seq_chain_{chain}"] = chain_to_seq[chain]

        # construct the new sequence
        # NOTE: assuming this is the order in which sequences are concatenated in the "seq" entry of each dictionary
        parsed_pdb_dict["seq"] = "".join([parsed_pdb_dict[f"seq_chain_{chain}"] for chain in all_chains])

        parsed_pdbs_list.append(parsed_pdb_dict)
    
    write_jsonl(parsed_pdbs_list, os.path.join(os.path.join(peptide_output_dir, pdbid, "parsed_pdbs.jsonl")))


    ## 3) create assigned_pdbs.jsonl --> designing the peptide chain
    pep_designed_chains = [peptide_chain]
    pep_fixed_chains = [chain for chain in all_chains if chain not in pep_designed_chains]
    assigned_pdbs_dict = {f"{pdbid}_{sample_id}": [pep_designed_chains, pep_fixed_chains] for sample_id in sample_id_to_chain_to_seq.keys()}
    write_jsonl([assigned_pdbs_dict], os.path.join(os.path.join(peptide_output_dir, pdbid, "assigned_pdbs.jsonl")))

    # NOTE: no need to create fixed_pdbs.jsonl since we're going to be designin/scoring a whole chain! (i.e. the whole peptide)



def generate_proteinmpnn_pwms(
    output_dir: str,
    model_name: str, # e.g. v_48_002, v_48_002, v_48_030
    num_seq: int = 100, # I think ProteinMPNN encodes at random, geting slightly different PWMs as a result when multiple amino-acids are being decoded. Thus num_seq makes the model average over multiple decoding patterns.
    sampling_temp: float = 0.1, # shouln't matter... keeping it just to check
    batch_size: int = 10,
) -> None:
    
    mpnn = ProteinMPNN(model_name)

    mpnn.run_on_prepared_inputs(
        output_dir,
        specified_chains_to_design = True,
        conditional_probs_only = True,
        num_seq_per_target = num_seq,
        sampling_temp = sampling_temp,
        batch_size = batch_size,
    )


def get_log_p(cond_proba_npz_path: str) -> np.ndarray:
    arr = np.load(cond_proba_npz_path)
    log_p = np.mean(arr['log_p'], axis=0)[arr['design_mask'].astype(np.bool_)]
    return log_p

def per_site_entropy(pwm_N20: np.ndarray) -> np.ndarray:
    return - np.sum(pwm_N20 * np.log2(pwm_N20), axis=1) # entropy_N

def parse_peptide_pwms(output_dir: str):

    # conditional_probs_only

    cond_probs_dir = os.path.join(output_dir, "conditional_probs_only")

    name_to_pwm = {}
    for npzfile in os.listdir(cond_probs_dir):
        name = npzfile[:-4] # exclude ".npz"
        log_p = get_log_p(os.path.join(cond_probs_dir, npzfile))
        name_to_pwm[name] = softmax(log_p, axis=1)
    
    



if __name__ == '__main__':

    cdr3_output_dir = "cdr3_samples"
    os.makedirs(cdr3_output_dir, exist_ok=True)

    pdbid = "5brz"

    pdbfile = f"./wt_pdbs/{pdbid}.pdb"
    model_name = "v_48_030"
    num_cdr3s_to_sample = 100

    peptide_output_dir = "peptide_pwms"
    peptide_chain = "B"

    # print('\nSampling multiple CDR3s...', flush=True)

    # start = time.time()
    # sample_cdr3s_with_proteinmpnn(
    #     os.path.join(cdr3_output_dir, f"{pdbid}"),
    #     num_cdr3s_to_sample,
    #     model_name,
    #     pdbfile,
    #     mhc_class = 1,
    #     organism = 'human',
    #     sampling_temp = 0.2,
    #     batch_size = 32,
    # )
    # print(f'\n\tElapsed time: {time.time() - start:.3f}', flush=True)

    # print('\nBuilding inputs for peptides..', flush=True)

    # start = time.time()
    # build_peptide_proteinmpnn_inputs_with_different_cdr3s_bypassing_parsing(
    #     pdbid,
    #     cdr3_output_dir,
    #     peptide_output_dir,
    #     peptide_chain
    # )
    # print(f'\n\tElapsed time: {time.time() - start:.3f}', flush=True)

    # print('\nGenerating peptide pwms...', flush=True)

    # start = time.time()
    # generate_proteinmpnn_pwms(
    #     os.path.join(peptide_output_dir, f"{pdbid}"),
    #     model_name,
    #     num_seq = 64,
    #     sampling_temp = 0.1,
    #     batch_size = 64,
    # )
    # print(f'\n\tElapsed time: {time.time() - start:.3f}', flush=True)


    parse_peptide_pwms(os.path.join(peptide_output_dir, f"{pdbid}"))


