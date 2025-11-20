
import json
from typing import List, Dict

def read_jsonl(path: str) -> List[Dict]:
    """Read a JSONL file and return a list of Python objects."""
    data = []
    with open(path, "r") as f:
        for line in f:
            data.append(json.loads(line))
    return data

def write_jsonl(data: List[Dict], path: str):
    """Write a list of Python objects to a JSONL file."""
    with open(path, "w") as f:
        for obj in data:
            f.write(json.dumps(obj) + "\n")
