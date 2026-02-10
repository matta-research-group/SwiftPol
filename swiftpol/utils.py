### Utilities for Swiftpol/Polyply MARTINI builder
### Under development

# Perceive sequences
def perceive_sequences(sys):
    sequence_dict = {}
    for polymer in sys.chain_rdkit:
        chain_id = sys.chain_rdkit.index(polymer)
        residue_ids = []
        residue_names = []
        for atom in polymer.GetAtoms():
            info = atom.GetPDBResidueInfo()
            residue_names.append(info.GetResidueName())
            residue_ids.append(info.GetResidueNumber())
        residue_ids_dict = {residue_ids[i]:residue_names[i] for i in range(len(residue_ids))} 
        for i in residue_ids_dict:
            if residue_ids_dict[i].find('A')>0:
                residue_ids_dict[i] = 'A'
            elif residue_ids_dict[i].find('B')>0:
                residue_ids_dict[i] = 'B'
            elif residue_ids_dict[i].find('C')>0:
                residue_ids_dict[i] = 'C'       
            elif residue_ids_dict[i].find('D')>0:
                residue_ids_dict[i] = 'D' 
            elif residue_ids_dict[i].find('S')>0:
                residue_ids_dict[i] = 'S'                               
        sorted_residue_ids_dict = dict(sorted(residue_ids_dict.items(), key=lambda item: int(item[0])))
        sequence = ''.join(sorted_residue_ids_dict.values())
        sequence_dict[chain_id] = sequence
    return sequence_dict


def export_seq_to_json(dict_seq, resname_map, chirality_map = {"A": "S", "S": "R", "B":"S"}):
    import json
    for s in dict_seq:
        sequence = dict_seq[s]
        nodes = []
        edges = []
    
        # build nodes
        for i, letter in enumerate(sequence):
            nodes.append({
                "resname": resname_map[letter],
                "seqid": 0,
                "chiral": chirality_map[letter],
                "id": i
            })
        # build linear edges
        for i in range(len(sequence) - 1):
            edges.append({
                "source": i,
                "target": i + 1
            })
        # full graph structure
        graph = {
            "directed": False,
            "multigraph": False,
            "graph": {},
            "nodes": nodes,
            "edges": edges
        }

        # write to json file
        with open(f"sequence_{s}.json", "w") as f:
            json.dump(graph, f, indent=2)
    