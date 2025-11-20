
import ast
import re
from typing import List, Tuple

def read_fasta(filepath: str) -> List[Tuple[str, str]]:
    sequences = []
    header = None
    seq_chunks = []

    with open(filepath, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                # save previous record
                if header is not None:
                    sequences.append((header, "".join(seq_chunks)))
                header = line[1:]          # remove ">"
                seq_chunks = []
            else:
                seq_chunks.append(line)

        # save the last record
        if header is not None:
            sequences.append((header, "".join(seq_chunks)))

    return sequences

def get_fixed_chains(s: str):
    # Match fixed_chains=[...anything...]
    m = re.search(r"fixed_chains=\[([^\]]*)\]", s)
    if not m:
        return None
    # Safely parse the content inside the brackets as Python literal
    return ast.literal_eval("[" + m.group(1) + "]")

def get_designed_chains(s: str):
    # Match designed_chains=[...anything...]
    m = re.search(r"designed_chains=\[([^\]]*)\]", s)
    if not m:
        return None
    # Safely parse the content inside the brackets as Python literal
    return ast.literal_eval("[" + m.group(1) + "]")


def get_field(s: str, key: str):
    """
    Extract the value for `key` using regex and parse it with ast.literal_eval.
    Falls back to raw string if literal_eval fails.
    """
    # Match key=VALUE up to comma or end of string
    pattern = rf"{re.escape(key)}=([^,]+)"
    m = re.search(pattern, s)
    if not m:
        return None

    value_str = m.group(1).strip()

    # Try safe Python literal parsing
    try:
        return ast.literal_eval(value_str)
    except Exception:
        return value_str
    