import networkx as nx
import pandas as pd
import numpy as np
import sys
import getopt

def assign_proteins(data: pd.DataFrame, pept_dict: dict) -> (pd.DataFrame, dict):
    """
    Assign psms to proteins.
    This function appends the dataframe with a column 'n_possible_proteins' which indicates how many proteins a psm could be matched to.
    It returns the appended dataframe and a dictionary `found_proteins` where each protein is mapped to the psms indices.

    Args:
        data (pd.DataFrame): psms table of scored and filtered search results from alphapept.
        pept_dict (dict): dictionary that matches peptide sequences to proteins

    Returns:
        pd.DataFrame: psms table of search results from alphapept appended with the number of matched proteins.
        dict: dictionary mapping psms indices to proteins.

    """

    data = data.reset_index(drop=True)

    data['n_possible_proteins'] = data['sequence'].apply(lambda x: len(pept_dict[x]))
    unique_peptides = (data['n_possible_proteins'] == 1).sum()
    shared_peptides = (data['n_possible_proteins'] > 1).sum()

    logging.info(f'A total of {unique_peptides:,} unique and {shared_peptides:,} shared peptides.')

    sub = data[data['n_possible_proteins'] == 1]
    psms_to_protein = sub['sequence'].apply(lambda x: pept_dict[x])

    found_proteins = {}
    for idx, _ in enumerate(psms_to_protein):
        idx_ = psms_to_protein.index[idx]
        p_str = 'p' + str(_[0])
        if p_str in found_proteins:
            found_proteins[p_str] = found_proteins[p_str] + [str(idx_)]
        else:
            found_proteins[p_str] = [str(idx_)]

    return data, found_proteins


def get_shared_proteins(data: pd.DataFrame, found_proteins: dict, pept_dict: dict) -> dict:
    """
    Assign peptides to razor proteins.

    Args:
        data (pd.DataFrame): psms table of scored and filtered search results from alphapept, appended with `n_possible_proteins`.
        found_proteins (dict): dictionary mapping psms indices to proteins
        pept_dict (dict): dictionary mapping peptide indices to the originating proteins as a list

    Returns:
        dict: dictionary mapping peptides to razor proteins

    """

    G = nx.Graph()

    sub = data[data['n_possible_proteins'] > 1]

    for i in range(len(sub)):
        seq, score = sub.iloc[i][['sequence', 'score']]
        idx = sub.index[i]
        possible_proteins = pept_dict[seq]

        for p in possible_proteins:
            G.add_edge(str(idx), 'p' + str(p), score=score)

    connected_groups = np.array([list(c) for c in sorted(nx.connected_components(G), key=len, reverse=True)],
                                dtype=object)
    n_groups = len(connected_groups)


    # Solving with razor:
    found_proteins_razor = {}
    for a in connected_groups[::-1]:
        H = G.subgraph(a).copy()
        shared_proteins = list(np.array(a)[np.array(list(i[0] == 'p' for i in a))])

        while len(shared_proteins) > 0:
            neighbors_list = []

            for node in shared_proteins:
                shared_peptides = list(H.neighbors(node))

                if node in G:
                    if node in found_proteins.keys():
                        shared_peptides += found_proteins[node]

                n_neigbhors = len(shared_peptides)

                neighbors_list.append((n_neigbhors, node, shared_peptides))

            # Check if we have a protein_group (e.g. they share the same everythin)
            neighbors_list.sort()

            # Check for protein group
            node_ = [neighbors_list[-1][1]]
            idx = 1
            while idx < len(neighbors_list):  # Check for protein groups
                if neighbors_list[-idx][0] == neighbors_list[-idx - 1][0]:  # lenght check
                    if set(neighbors_list[-idx][2]) == set(neighbors_list[-idx - 1][2]):  # identical peptides
                        node_.append(neighbors_list[-idx - 1][1])
                        idx += 1
                    else:
                        break
                else:
                    break

            # Remove the last entry:
            shared_peptides = neighbors_list[-1][2]
            for node in node_:
                shared_proteins.remove(node)

            for _ in shared_peptides:
                if _ in H:
                    H.remove_node(_)

            if len(shared_peptides) > 0:
                if len(node_) > 1:
                    node_ = tuple(node_)
                else:
                    node_ = node_[0]

                found_proteins_razor[node_] = shared_peptides

    return found_proteins_razor


def get_protein_groups(data: pd.DataFrame, pept_dict: dict, fasta_dict: dict, decoy=False, callback=None,
                       **kwargs) -> pd.DataFrame:
    """
    Function to perform protein grouping by razor approach.
    This function calls `assign_proteins` and `get_shared_proteins`.
    ToDo: implement callback for solving
    Each protein is indicated with a p -> protein index

    Args:
        data (pd.DataFrame): psms table of scored and filtered search results from alphapept.
        pept_dict (dict): A dictionary mapping peptide indices to the originating proteins as a list.
        fasta_dict (dict): A dictionary with fasta sequences.
        decoy (bool, optional): Defaults to False.
        callback (bool, optional): Defaults to None.

    Returns:
        pd.DataFrame: alphapept results table now including protein level information.
    """
    data, found_proteins = assign_proteins(data, pept_dict)
    found_proteins_razor = get_shared_proteins(data, found_proteins, pept_dict)

    report = data.copy()

    assignment = np.zeros(len(report), dtype=object)
    assignment[:] = ''
    assignment_pg = assignment.copy()

    assignment_idx = assignment.copy()
    assignment_idx[:] = ''

    razor = assignment.copy()
    razor[:] = False

    if decoy:
        add = 'REV__'
    else:
        add = ''

    for protein_str in found_proteins.keys():
        protein = int(protein_str[1:])
        protein_name = add + fasta_dict[protein]['name']
        indexes = [int(_) for _ in found_proteins[protein_str]]
        assignment[indexes] = protein_name
        assignment_pg[indexes] = protein_name
        assignment_idx[indexes] = str(protein)

    for protein_str in found_proteins_razor.keys():
        indexes = [int(_) for _ in found_proteins_razor[protein_str]]

        if isinstance(protein_str, tuple):
            proteins = [int(_[1:]) for _ in protein_str]
            protein_name = ','.join([add + fasta_dict[_]['name'] for _ in proteins])
            protein = ','.join([str(_) for _ in proteins])

        else:
            protein = int(protein_str[1:])
            protein_name = add + fasta_dict[protein]['name']

        assignment[indexes] = protein_name
        assignment_pg[indexes] = protein_name
        assignment_idx[indexes] = str(protein)
        razor[indexes] = True

    report['protein'] = assignment
    report['protein_group'] = assignment_pg
    report['razor'] = razor
    report['protein_idx'] = assignment_idx

    return report


def perform_protein_grouping(data: pd.DataFrame) -> pd.DataFrame:
    """
    Wrapper function to perform protein grouping by razor approach

    Args:
        data (pd.DataFrame): psms table of scored and filtered search results from alphapept.
        pept_dict (dict): A dictionary mapping peptide indices to the originating proteins as a list.
        fasta_dict (dict): A dictionary with fasta sequences.

    Returns:
        pd.DataFrame: alphapept results table now including protein level information.
    """
    data_sub = data[['db_idx','sequence', 'score', 'decoy']]
    data_sub = data[data['score']>0.99]
    data_sub_unique = data_sub.groupby(['db_idx','sequence', 'decoy'], as_index=False).agg({"score": "max"})
    targets = data_sub_unique[data_sub_unique.decoy == False]
    targets = targets.reset_index(drop=True)
    #protein_targets = get_protein_groups(targets, pept_dict, fasta_dict, **kwargs)

    #protein_targets['decoy_protein'] = False

    decoys = data_sub_unique[data_sub_unique.decoy == True]
    decoys = decoys.reset_index(drop=True)
    #protein_decoys = get_protein_groups(decoys, pept_dict, fasta_dict, decoy=True, **kwargs)

    #protein_decoys['decoy_protein'] = True

    #protein_groups = pd.concat([protein_targets, protein_decoys])
    #protein_groups_app = protein_groups[
    #    ['sequence', 'decoy', 'protein', 'protein_group', 'razor', 'protein_idx', 'decoy_protein',
    #     'n_possible_proteins']]
    #protein_report = pd.merge(data,
    #                          protein_groups_app,
    #                          how='inner',
    #                          on=['sequence', 'decoy'],
    #                          validate="many_to_one")

    return targets

if __name__ == '__main__':
    argv=sys.argv[1:]
    try:
        opts, args = getopt.getopt(argv, "hi:m:o:")
    except:
        print("Error Option, using -h for help information.")
        sys.exit(1)
    if len(opts)==0:
        print("\n\nUsage:\n")
        print("-i\t xml format output file for conversion")
        print("-o\t Output for WinnowNet's input\n")
        sys.exit(1)
    input_file=""
    output_file=""
    for opt, arg in opts:
        if opt in ("-h"):
            print("\n\nUsage:\n")
            print("-i\t xml format output file for conversion")
            print("-o\t Output for WinnowNet's input\n")
            sys.exit(1)
        elif opt in ("-i"):
            input_file=arg
        elif opt in ("-o"):
            output_file=arg

    df = pd.read_csv(input_file)
    df=perform_protein_grouping(df)
    df.to_csv(output_file)

    print('All done.')