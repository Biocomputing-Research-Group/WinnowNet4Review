import pandas as pd
import sys
import getopt
import numpy as np


def get_q_values(fdr_values: np.ndarray) -> np.ndarray:
    """
    Calculate q-values from fdr_values.

    Args:
        fdr_values (np.ndarray): np.ndarray of fdr values.

    Returns:
        np.ndarray: np.ndarray of q-values.
    """
    q_values = np.zeros_like(fdr_values)
    min_q_value = np.max(fdr_values)
    for i in range(len(fdr_values) - 1, -1, -1):
        fdr = fdr_values[i]
        if fdr < min_q_value:
            min_q_value = fdr
        q_values[i] = min_q_value

    return q_values


def cut_fdr(df: pd.DataFrame, fdr_level: float = 0.01, cut: bool = True) -> (float, pd.DataFrame):
    """
    Cuts a dataframe with a given fdr level

    Args:
        df (pd.DataFrame): psms table of search results from alphapept.
        fdr_level (float, optional): fdr level that should be used for filtering. The value should lie between 0 and 1. Defaults to 0.01.
        plot (bool, optional): flag to enable plot. Defaults to 'True'.
        cut (bool, optional): flag to cut above fdr threshold. Defaults to 'True'.

    Returns:
        float: numerical value of the applied score cutoff
        pd.DataFrame: df with psms within fdr

    """

    df["target"] = ~df["decoy"]

    df = df.sort_values(by=["score", "decoy"], ascending=False)
    df = df.reset_index()

    df["target_cum"] = np.cumsum(df["target"])
    df["decoys_cum"] = np.cumsum(df["decoy"])

    df["fdr"] = df["decoys_cum"] / df["target_cum"]
    df["q_value"] = get_q_values(df["fdr"].values)

    last_q_value = df["q_value"].iloc[-1]
    first_q_value = df["q_value"].iloc[0]

    if last_q_value <= fdr_level:
        cutoff_index = len(df) - 1

    elif first_q_value >= fdr_level:
        cutoff_index = 0

    else:
        cutoff_index = df[df["q_value"].gt(fdr_level)].index[0] - 1

    cutoff_value = df.loc[cutoff_index]["score"]

    if cut:
        cutoff = df[df["score"] >= cutoff_value]
    else:
        cutoff = df

    targets = df.loc[cutoff_index, "target_cum"]
    decoy = df.loc[cutoff_index, "decoys_cum"]

    fdr = df.loc[cutoff_index, "fdr"]

if __name__ == '__main__':
    argv=sys.argv[1:]
    try:
        opts, args = getopt.getopt(argv, "hi:r:o:")
    except:
        print("Error Option, using -h for help information.")
        sys.exit(1)
    if len(opts)==0:
        print("\n\nUsage:\n")
        print("-i\t xml format output file for conversion")
        print("-r\t re-scoring file for conversion")
        print("-o\t Output for WinnowNet's input\n")
        sys.exit(1)
    input_file=""
    rescore_file=""
    output_file=""
    for opt, arg in opts:
        if opt in ("-h"):
            print("\n\nUsage:\n")
            print("-i\t xml format output file for conversion")
            print("-r\t re-scoring file for conversion")
            print("-o\t Output for WinnowNet's input\n")
            sys.exit(1)
        elif opt in ("-i"):
            input_file=arg
        elif opt in ("-r"):
            rescore_file=arg
        elif opt in ("-o"):
            output_file=arg

    df=pd.read_csv(input_file)
    rescore_list=[]
    with open(rescore_file) as f:
        for line in f:
            rescore_list.append(float(line.strip()))

    df['score']=rescore_list
    cut_fdr(df,fdr_level=0.01,cut=True)
    df.to_csv(output_file)



    print('All done.')