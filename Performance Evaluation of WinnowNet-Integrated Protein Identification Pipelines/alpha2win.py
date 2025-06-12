import pandas as pd
import sys
import getopt

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
    psm_list = []
    prefix=input_file.split('.')[0]
    df=pd.read_csv(input_file)
    psm_list=[]
    for line_id, line in enumerate(df['sequence_naked']):
        scan=str(df['scan_no'][line_id])
        charge=str(df['charge'][line_id])
        rank='1'
        psmid='_'.join([prefix, scan, charge, rank])
        peptide='-.'+line.replace('ox','~')+'.-'
        protein = str(df['db_idx'][line_id])
        psm_list.append([psmid, peptide, protein])

    with open(output_file,'w') as f:
        for line in psm_list:
            f.write('\t'.join(line))


    print('All done.')
