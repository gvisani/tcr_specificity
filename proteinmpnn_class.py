
'''

The following class can be used to run proteinmpnn to do the following tasks:
- Sample sequence designs for a provided set of residues
- Output the PWM for a provided set of residues

'''

import os
from typing import Optional, Union, Tuple, Dict, List
from Bio.PDB import PDBParser

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
        
        Args:
        ----
        output_dir: str
            where to store the outputs of proteinmpnn

        pdbfile: Optional[str]
            path to a pdb file to run proteinmpnn on, either this or `pdbdir` must be specified, and only one can be specified
        
        pdbdir: Optional[str]
            path to a directory containing pdb files to run proteinmpnn on, either this or `pdbfile` must be specified, and only one can be specified
        
        chains_or_positions_to_design_or_score: Optional[Union[List[str], Dict[str, int]]]
            if None, the whole protein is designed/scored
            if List[str], then the items represent chainids to be designed/scored
            if Dict[str, int], then keys represent chainids and values represent resnums to be designed/scored (resnums as shown in the pdb file(s), or positions as used internally by proteinmpnn).
            only empty icodes can be requested. see `fix_pdbfile_offset` if providing resnums as appearing in the pdb file(s)
        
        fix_pdbfile_offset: bool (default True)
            only applies if `chains_or_positions_to_design_or_score` is a Dict[str, int]
            essentially it converts the resnums in the pdb file that you (i.e. the user) want to design into the numbering used internally by proteinmpnn
            it involves parsing the structure with biopython, which can slow things down a bit, so only do this if you know for sure that you need it
            also, alternatively you might want to provide the positions as used by proteinmpnn to `chains_or_positions_to_design_or_score` and set this argument to False
        
        conditional_probs_only: 
        
        num_seq_per_target: int (default 10)
            proteinmpnn parameter. how many sequences to generate per pdb file
        
        sampling_temp: float (default 0.1)
            proteinmpnn parameter. inconsequential if `conditional_probs_only` is set to True

        batch_size: int (default 1)
            proteinmpnn parameter
        
        
        Returns:
        ----
        None
        
        '''

        ## initialize directory with pdbs
        ## either a pdbdir is specified, in which case ProteinMPNN is run on all pdbfiles in it,
        ##      or a single pdbfile is specified, in which case it is moved to a temporary directory within output_dir
        assert pdbfile is not None or pdbdir is not None, "One of `pdbfile` or `pdbdir` must be specified."
        assert not ((pdbfile is not None) and (pdbdir is not None)), "ONLY one of `pdbfile` and `pdbdir` can be specified, not both."

        if pdbfile is not None:
            # pdbfile case
            # make temporary directory with the desired pdb
            pdbdir_in_use = os.path.join(output_dir, "temp_pdbs")
            if os.path.exists(pdbdir_in_use):
                print(f"Warning: {pdbdir_in_use} already exists, deleting it...")
                os.system(f'rm -r {pdbdir_in_use}')
        
            os.makedirs(pdbdir_in_use)
            os.system(f'cp {pdbfile} {pdbdir_in_use}')

            pdbfile_for_fix = pdbfile
        
        elif pdbdir is not None:
            # pdbdir case
            # assuming that the same offset fix applies to all pdb files in pdbdir (true for sure if structurres are the same just with different sequences)
            #       just grabbing one representative pdbfile for the fix 
            pdbdir_in_use = pdbdir
            pdbfile_for_fix = os.path.join(pdbdir, os.listdir(pdbdir)[0])

        path_for_parsed_chains = os.path.join(output_dir, 'parsed_pdbs.jsonl')
        os.system(f'python {os.path.join(self.proteinmpnn_path, f"helper_scripts/parse_multiple_chains.py")} --input_path={pdbdir_in_use} --output_path={path_for_parsed_chains}')

        proteinmpnn_arguments = [
            f"--suppress_print 1",
            f"--model_name {self.model_name}",
            f"--jsonl_path {path_for_parsed_chains}",
            f"--out_folder {output_dir}",
            f"--num_seq_per_target {num_seq_per_target}",
            f"--batch_size {batch_size}",
            f"--sampling_temp {sampling_temp}",
        ]

        if conditional_probs_only:
            proteinmpnn_arguments.append(f"--conditional_probs_only 1")

        if chains_or_positions_to_design_or_score is not None:

            if isinstance(chains_or_positions_to_design_or_score, list):
                ## design specific chains

                chains_str = " ".join(chains_or_positions_to_design_or_score)
                path_for_assigned_chains = os.path.join(output_dir, 'assigned_pdbs.jsonl')
            
                os.system(f'python {os.path.join(self.proteinmpnn_path, f"helper_scripts/assign_fixed_chains.py")} --input_path={path_for_parsed_chains} --output_path={path_for_assigned_chains} --chain_list "{chains_str}"')

                proteinmpnn_arguments.append(f"--chain_id_jsonl {path_for_assigned_chains}")

            elif isinstance(chains_or_positions_to_design_or_score, dict):
                ## design specific positions on chains, assumig only empty icodes are being designed for now

                chains_str, resnums_str = self._format_chain_and_resnums(chains_or_positions_to_design_or_score, fix_pdbfile_offset=fix_pdbfile_offset, pdbfile=pdbfile_for_fix)
                print()
                print(chains_str)
                print(resnums_str)
                print()
                path_for_assigned_chains = os.path.join(output_dir, 'assigned_pdbs.jsonl')
                path_for_fixed_positions = os.path.join(output_dir, 'fixed_pdbs.jsonl')
            
                os.system(f'python {os.path.join(self.proteinmpnn_path, f"helper_scripts/assign_fixed_chains.py")} --input_path={path_for_parsed_chains} --output_path={path_for_assigned_chains} --chain_list "{chains_str}"')
                os.system(f'python {os.path.join(self.proteinmpnn_path, f"helper_scripts/make_fixed_positions_dict.py")} --input_path={path_for_parsed_chains} --output_path={path_for_fixed_positions} --chain_list "{chains_str}" --position_list "{resnums_str}" --specify_non_fixed')

                proteinmpnn_arguments.append(f"--chain_id_jsonl {path_for_assigned_chains}")
                proteinmpnn_arguments.append(f"--fixed_positions_jsonl {path_for_fixed_positions}")
            
            else:
                raise ValueError(f"Invalid type for `chains_or_positions_to_design_or_score`, must be None, list or dict, not {type(chains_or_positions_to_design_or_score)}.")
        
        final_command = f"python {os.path.join(self.proteinmpnn_path, f'protein_mpnn_run.py')} {' '.join(proteinmpnn_arguments)} "
        print()
        print(final_command)
        print()
        os.system(final_command)


    def _format_chain_and_resnums(
            self,
            chain_to_resnums_dict: Dict[str, int], # what's called `chains_or_positions_to_design_or_score` in self.run(), assumed to be a dictionary here
            fix_pdbfile_offset: bool = True,
            pdbfile: Optional[str] = None,
        ) -> Tuple[str, str]:

        if fix_pdbfile_offset:
            assert pdbfile is not None, f"Must specify pdbfile if fix_pdbfile_offset is set to True"
            chain_to_offset_dict = self._get_offset_between_proteinmpnn_residue_index_and_PDB_residue_index(pdbfile)
            chain_to_resnums_dict ={
                chain: [resnum + chain_to_offset_dict[chain][resnum] for resnum in resnum_list]
                    for chain, resnum_list in chain_to_resnums_dict.items()
            }

        chains_str = ' '.join(chain_to_resnums_dict.keys())
        resnums_str = ','.join([' '.join(map(str, resnums)) for resnums in chain_to_resnums_dict.values()])
        return chains_str, resnums_str

    def _get_offset_between_proteinmpnn_residue_index_and_PDB_residue_index(
            self,
            pdbfile: str,
    ) -> Dict[str, Dict[int, int]]:

        parser = PDBParser()
        structure = parser.get_structure(os.path.basename(pdbfile).strip('.pdb'), pdbfile)

        chain_to_start_offset_dict = {}
        for chain_obj in structure.get_chains():
            chain = chain_obj.get_id()
            for residue_obj in chain_obj:
                offset = residue_obj.get_id()[1] - 1
                chain_to_start_offset_dict[chain] = - offset
                break
        
        chain_to_offset_dict = {}
        for chain_obj in structure.get_chains():
            chain = chain_obj.get_id()
            chain_to_offset_dict[chain] = {}
            num_nonempty_icodes_so_far = 0
            for residue_obj in chain_obj:
                resnum, icode = residue_obj.get_id()[1], residue_obj.get_id()[2]
                if icode != ' ':
                    num_nonempty_icodes_so_far += 1
                else:
                    chain_to_offset_dict[chain][resnum] = chain_to_start_offset_dict[chain] + num_nonempty_icodes_so_far
                
        return chain_to_offset_dict


'''

python /gscratch/spe/gvisan01/ProteinMPNN-copy/protein_mpnn_run.py --suppress_print 1 --model_name v_48_002 --jsonl_path test_output_5brz/parsed_pdbs.jsonl  --chain_id_jsonl test_output_5brz/assigned_pdbs.jsonl  --fixed_positions_jsonl test_output_5brz/fixed_pdbs.jsonl  --out_folder test_output_5brz  --num_seq_per_target 10  --batch_size 10  --sampling_temp 0.1

'''
