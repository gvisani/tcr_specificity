

def extract_v_and_j_genes_and_cdr3s_from_tcrs(organism, alpha_chain, beta_chain):

    '''
    
    Returns list of two tuples. Each tuple element contains the V gene, J gene, and CDR3 sequence, respectively for the alpha and beta chains of a TCR.

    Example usage:
    (a_v_gene, a_j_gene, a_cdr3), (b_v_gene, b_j_gene, b_cdr3) = extract_v_and_j_genes_and_cdr3s_from_tcrs('human', 'MKLVTSITVLLSLGIMGDAKTTQPNSMESNEEEPVHLPCNHSTISGTDYIHWYRQLPSQGPEYVIHGLTSNVNNRMASLAIAEDRKSSTLILHRATLRDAAVYYCILDNNNDMRFGAGTRLTVKP', 'MGTSLLCWVVLGFLGTDHTGAGVSQSPRYKVTKRGQDVALRCDPISGHVSLYWYRQALGQGPEFLTYFNYEAQQDKSGLPNDRFSAERPEGSISTLTIQRTEQRDSAMYRCASSLAPGATNEKLFFGSGTQLSVL')

    '''

    assert organism in ['human', 'mouse']

    from .tcrdock_info import TCRdockInfo

    info = TCRdockInfo().from_sequences(organism, None, '', '', '', alpha_chain, beta_chain)

    return info.tcr




