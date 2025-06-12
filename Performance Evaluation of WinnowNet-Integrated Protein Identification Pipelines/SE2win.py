import re
import sys
import getopt

class peptide:
    def __init__(self):
        self.identified_pep = ""
        self.protein_list = []


class scan:
    def __init__(self):
        self.fidx = ""               # Filename
        self.measured_mass = None    # MeasuredParentMass
        self.scan_number = 0         # ScanNumber
        self.charge = 0              # ParentCharge
        self.pep_list = []           # list of peptide objects

    def add_pep(self, pep):
        self.pep_list.append(pep)

    def get_winnownet_output(self):
        common_prefix = f"{self.fidx}_{self.scan_number}_{self.charge}"
        results = []
        for pepid, pep in enumerate(self.pep_list):
            row = [
                f"{common_prefix}_{pepid+1}",
                pep.identified_pep.replace("+16", "~"),
                ",".join(pep.protein_list)
            ]
            results.append(row)
        return results



def parse_psm_file(path_to_file):
    scans = []
    current_scan = None

    with open(path_to_file, "r") as fin:
        for line_id, line in enumerate(fin):
            if line_id<161:
                continue
            line = line.strip()
            if not line:
                continue
            if line.startswith("+\tFilename"):
                continue
            if line.startswith("*\tIdentifiedPeptide"):
                continue

            if line.startswith("+\t"):
                parts = line.split("\t")
                fname = parts[1].split('.')[0]
                scan_num = int(parts[2])
                parent_charge = int(parts[3])
                measured_mass = float(parts[4])
                current_scan = scan()
                current_scan.fidx = fname
                current_scan.scan_number = scan_num
                current_scan.charge = parent_charge
                current_scan.measured_mass = measured_mass

                scans.append(current_scan)
                continue

            if line.startswith("*\t"):

                parts = line.split("\t")

                identified = parts[1].replace('[','.').replace(']','.')
                proteins   = parts[7].strip()
                proteins = proteins.strip()
                if proteins.startswith("{") and proteins.endswith("}"):
                    proteins = proteins[1:-1]
                protein_list = proteins.split(",") if proteins else []
                pep = peptide()
                pep.identified_pep = identified

                pep.charge = current_scan.charge
                pep.protein_list = protein_list

                if current_scan is None:
                    raise RuntimeError("Found a '*' line before any '+' line.")
                current_scan.add_pep(pep)
                continue

            continue

    return scans

def write_output(psm_list, output_file_str):
    with open(output_file_str, 'w') as fw:
        for v in psm_list:
            if v.pep_list:
                ll = v.get_winnownet_output()
                for l in ll:
                    fw.write('\t'.join(l))
                    fw.write('\n')

    print('write_output is done.')


if __name__ == "__main__":
    argv=sys.argv[1:]
    try:
        opts, args = getopt.getopt(argv, "hi:o:")
    except:
        print("Error Option, using -h for help information.")
        sys.exit(1)
    if len(opts)==0:
        print("\n\nUsage:\n")
        print("-i\t MS-GF+'s output file for conversion")
        print("-o\t Output for WinnowNet's input\n")
        sys.exit(1)
    input_file=""
    output_file=""
    for opt, arg in opts:
        if opt in ("-h"):
            print("\n\nUsage:\n")
            print("-i\t MS-GF+'s format output file for conversion")
            print("-o\t Output for WinnowNet's input\n")
            sys.exit(1)
        elif opt in ("-i"):
            input_file=arg
        elif opt in ("-o"):
            output_file=arg

    all_scans = parse_psm_file(input_file)
    write_output(all_scans,output_file)