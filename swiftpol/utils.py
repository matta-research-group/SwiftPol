### Utilities for Swiftpol/Polyply MARTINI builder
### Under development

# Perceive sequences
def perceive_sequences(sys):
    import warnings
    sequence_dict = {}
    for polymer in sys.chain_rdkit:
        chain_id = sys.chain_rdkit.index(polymer)
        residue_ids = []
        residue_names = []
        for atom in polymer.GetAtoms():
            info = atom.GetPDBResidueInfo()
            try:
                residue_names.append(info.GetResidueName())
                residue_ids.append(info.GetResidueNumber())
            except:
                warnings.warn(f"Atom {atom.GetIdx()} in chain {chain_id} does not have valid residue information. Skipping this atom.")
                continue
        residue_ids_dict = {residue_ids[i]: residue_names[i] for i in range(len(residue_ids))} 
        for i in residue_ids_dict:
            if residue_ids_dict[i].find('A') > 0:
                residue_ids_dict[i] = 'A'
            elif residue_ids_dict[i].find('B') > 0:
                residue_ids_dict[i] = 'B'
            elif residue_ids_dict[i].find('C') > 0:
                residue_ids_dict[i] = 'C'       
            elif residue_ids_dict[i].find('D') > 0:
                residue_ids_dict[i] = 'D' 
            elif residue_ids_dict[i].find('S') > 0:
                residue_ids_dict[i] = 'S'                               
        sorted_residue_ids_dict = dict(sorted(residue_ids_dict.items(), key=lambda item: int(item[0])))
        sequence = ''.join(sorted_residue_ids_dict.values())
        # Terminal info for ester and carboxyl groups
        if sys.terminals:
            if sys.terminals == 'ester':
                sequence += 'Z'
            elif sys.terminals == 'carboxyl':
                sequence += 'X'
            else:
                warnings.warn(f"Unsupported terminal type '{sys.terminals}'. No terminal information will be appended to the sequence.")
        
        sequence_dict[chain_id] = sequence

    return sequence_dict


def export_seq_to_json(dict_seq, resname_map, chirality_map={"A": "S", "S": "R", "B": "S"}, terminals_map={"Z": "E", "X": "C"}):
    """
    Export sequences to JSON format with support for terminals.

    Parameters
    ----------
    dict_seq : dict
        Dictionary of sequences where keys are sequence IDs and values are sequences.
    resname_map : dict
        Mapping of sequence letters to residue names.
    chirality_map : dict, optional
        Mapping of sequence letters to chirality (default: {"A": "S", "S": "R", "B": "S"}).
    terminals_map : dict, optional
        Mapping of terminal names to their corresponding sequence letters (default: {"ESTER": "E", "ACID": "C"}).
    """
    import json

    for s in dict_seq:
        sequence = dict_seq[s]
        nodes = []
        edges = []

        # Check for terminals at the start and end of the sequence
        start_terminal = None
        end_terminal = None

        if sequence.startswith(tuple(terminals_map.values())):
            start_terminal = sequence[0]
            sequence = sequence[1:]  # Remove the start terminal from the sequence

        if sequence.endswith(tuple(terminals_map.values())):
            end_terminal = sequence[-1]
            sequence = sequence[:-1]  # Remove the end terminal from the sequence

        # Build nodes
        for i, letter in enumerate(sequence):
            nodes.append({
                "resname": resname_map[letter],
                "seqid": 0,
                "chiral": chirality_map.get(letter, None),  # Default to None if chirality is not defined
                "id": i
            })

        # Add start terminal as the first node, if present
        if start_terminal:
            nodes.insert(0, {
                "resname": resname_map[start_terminal],
                "seqid": -1,  # Use -1 for terminal residues
                "chiral": None,  # Terminals typically don't have chirality
                "id": 0
            })

        # Add end terminal as the last node, if present
        if end_terminal:
            nodes.append({
                "resname": resname_map[end_terminal],
                "seqid": len(sequence),  # Use the sequence length as the ID for the end terminal
                "chiral": None,
                "id": len(nodes)  # ID is the current length of the nodes list
            })

        # Build linear edges
        for i in range(len(nodes) - 1):
            edges.append({
                "source": i,
                "target": i + 1
            })

        # Full graph structure
        graph = {
            "directed": False,
            "multigraph": False,
            "graph": {},
            "nodes": nodes,
            "edges": edges
        }

        # Write to JSON file
        with open(f"sequence_{s}.json", "w") as f:
            json.dump(graph, f, indent=2)