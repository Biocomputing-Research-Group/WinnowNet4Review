import xml.etree.ElementTree as ET
import sys
import getopt
W_hydrogen=1.00784
#5_S_57.021464,7_V_15.994900
#GMGRCPM[147]M
#<mod_aminoacid_mass position="5" mass="160.030649" static="57.021464"/>
#<mod_aminoacid_mass position="7" mass="147.035385" variable="15.994900" source="param"/>
class peptide:

    def __init__(self):
        self.identified_pep = ""
        self.original_pep = ""
        self.calculated_mass = None
        self.mz_score = None
        self.mvh = None
        self.mvh_diff=0
        self.protein_list = []
        self.modifications=None
        self.prev_aa=''
        self.next_aa=''

class scan:

    def __init__(self):
        self.fidx = 0
        self.measured_mass = None
        self.charge = None
        self.scan_number = 0
        self.pep_list = []

    def add_pep(self, pep):
        self.pep_list.append(pep)

    def find_diff(self):
        pep_sorted_list = sorted(self.pep_list, key=lambda pep: -pep.mvh)
        for i in range(len(pep_sorted_list) - 1):
            pep_sorted_list[i].mvh_diff = pep_sorted_list[i].mvh - pep_sorted_list[i + 1].mvh

    def get_ms2rescore_output(self):
        # peptidoform\tspectrum_id\trun\tcollection\tis_decoy\tscore\n
        results = []
        for pep in self.pep_list:
            l = [pep.identified_pep.replace('~', '[U:Oxidation]')+'/'+self.charge]
            l.append(str(self.scan_number))
            l.append(self.fidx)
            l.append('')
            Label=False
            for p in pep.protein_list:
                if 'Rev_' not in p:
                    Label=True
                    break
            if Label:
                l.append('False')
            else:
                l.append('True')
            l.append(str(pep.mvh))
            l.append(str(pep.protein_list))
            results.append(l)

        return results




def read_one_pep_xml_myrimatch(input_file_str, psm_list, fidx):
    global xcorr_minimum
    tree = ET.parse(input_file_str)
    root = tree.getroot()
    msms_run_summary = root.find('{http://regis-web.systemsbiology.net/pepXML}msms_run_summary')
    for spectrum_query in msms_run_summary.iter('{http://regis-web.systemsbiology.net/pepXML}spectrum_query'):
        # scan information
        split_list = spectrum_query.attrib['spectrum'].split('.')
        one_scan = scan()
        one_scan.fidx = fidx
        one_scan.scan_number = int(split_list[-2]) 
        one_scan.measured_mass = spectrum_query.attrib['precursor_neutral_mass']
        one_scan.charge = spectrum_query.attrib['assumed_charge']
        for search_result in spectrum_query.iter('{http://regis-web.systemsbiology.net/pepXML}search_result'):
            for search_hit in search_result.iter('{http://regis-web.systemsbiology.net/pepXML}search_hit'):
                pep = peptide()
                # set the original peptide
                pep.original_pep = search_hit.attrib['peptide']
                # get the prev aa
                pep.prev_aa = search_hit.attrib['peptide_prev_aa']
                # get the next aa
                pep.next_aa = search_hit.attrib['peptide_next_aa']
                # get the calculated mass
                pep.calculated_mass = search_hit.attrib['calc_neutral_pep_mass']
                # get the identified peptide
                pos_list = []
                modification_list=[]
                for modification_info in search_hit.findall('{http://regis-web.systemsbiology.net/pepXML}modification_info'):
                    for mod_aminoacid_mass in modification_info.iter('{http://regis-web.systemsbiology.net/pepXML}mod_aminoacid_mass'):
                        pos_int = int(mod_aminoacid_mass.attrib['position'])
                        mass = mod_aminoacid_mass.attrib['mass']
                        if mass == '160.030649':
                            modification_list.append(str(pos_int)+'_'+'S_57.021464')
                        elif mass == '147.035385':
                            modification_list.append(str(pos_int)+'_'+'V_15.994900')
                        if pep.original_pep[pos_int-1] == 'M':
                            pos_list.append(pos_int)
                pep.modifications=','.join(modification_list)
                pep.identified_pep = pep.original_pep
                if len(pos_list) > 0:
                    pos_list = sorted(pos_list)
                    for idx, v in enumerate(pos_list):
                        pep.identified_pep = pep.identified_pep[:(v+idx)] + '~' + pep.identified_pep[(v+idx):]
                # get the negative log e valeu score
                for search_score in search_hit.iter('{http://regis-web.systemsbiology.net/pepXML}search_score'):
                    if search_score.attrib['name']=='mvh':
                        pep.mvh = float(search_score.attrib['value'])
                    if search_score.attrib['name']=='mzFidelity':
                        pep.mz_score= search_score.attrib['value']
                # get proteins
                pep.protein_list.append(search_hit.attrib['protein'])
                for alternative_protein in search_hit.iter('{http://regis-web.systemsbiology.net/pepXML}alternative_protein'):
                    pep.protein_list.append(alternative_protein.attrib['protein'])
                # add to this scan
                one_scan.add_pep(pep)
        if len(one_scan.pep_list) > 0:
            psm_list.append(one_scan)
    

def write_output(psm_list, output_file_str):
    with open(output_file_str, 'w') as fw:
        header_str = 'peptidoform\tspectrum_id\trun\tcollection\tis_decoy\tscore\tprotein_list'
        fw.write(header_str)
        fw.write('\n')
        for v in psm_list:
            if v.pep_list:
                ll = v.get_ms2rescore_output()
                for l in ll:
                    fw.write('\t'.join(l))
                    fw.write('\n')

    print('write_output is done.')

def main():
    argv=sys.argv[1:]
    try:
        opts, args = getopt.getopt(argv, "hi:m:o:")
    except:
        print("Error Option, using -h for help information.")
        sys.exit(1)
    if len(opts)==0:
        print("\n\nUsage:\n")
        print("-i\t MyriMatch's pepXML format output file for conversion")
        print("-o\t Output for Q-ranker's input\n")
        sys.exit(1)
    input_file=""
    output_file=""
    for opt, arg in opts:
        if opt in ("-h"):
            print("\n\nUsage:\n")
            print("-i\t MyriMatch's pepXML format output file for conversion")
            print("-o\t Output for Q-ranker's input\n")
            sys.exit(1)
        elif opt in ("-i"):
            input_file=arg
        elif opt in ("-o"):
            output_file=arg
    psm_list = []
    prefix=input_file.split('.')[-2]
    read_one_pep_xml_myrimatch(input_file,psm_list,prefix)
    write_output(psm_list, output_file)

    print('All done.')

if __name__ == '__main__':
    sys.exit(main())

