import xml.etree.ElementTree as ET
import sys
import getopt
import math
W_proton=1.00724
W_hydrogen=1.00784

class peptide:

    def __init__(self):
        self.identified_pep = ""
        self.original_pep = ""
        self.calculated_mass = None
        self.charge = None
        self.EValue=None
        self.score_diff = 0
        self.protein_list = []
        self.modifications = None
        self.prev_aa = ''
        self.next_aa = ''
        self.modifications=None
        self.prev_aa=''
        self.next_aa=''

class scan:

    def __init__(self):
        self.fidx = 0
        self.measured_mass = None
        self.scan_number = 0
        self.charge=0
        self.pep_list = []

    def add_pep(self, pep):
        self.pep_list.append(pep)

    def get_winnownet_output(self):
        common_prefix=self.fidx+'_'+str(self.scan_number)+'_'+self.charge
        results = []
        for pepid, pep in enumerate(self.pep_list):
            l=[common_prefix+'_'+str(pepid+1)]
            l.append(pep.identified_pep.replace('+16','~'))
            l.append(','.join(pep.protein_list))
            results.append(l)

        return results



def pep2proMap(input_file_str,Map):
    tree = ET.parse(input_file_str)
    root = tree.getroot()
    SequenceCollection=root.find('{http://psidev.info/psi/pi/mzIdentML/1.1}SequenceCollection')
    for DBSequence in SequenceCollection.iter('{http://psidev.info/psi/pi/mzIdentML/1.1}DBSequence'):
        Map[DBSequence.attrib['id']]=DBSequence.attrib['accession']

def prepostDBref(input_file_str,prepostMap):
    tree = ET.parse(input_file_str)
    root = tree.getroot()
    SequenceCollection=root.find('{http://psidev.info/psi/pi/mzIdentML/1.1}SequenceCollection')
    for PeptideEvidence in SequenceCollection.iter('{http://psidev.info/psi/pi/mzIdentML/1.1}PeptideEvidence'):
        prepostMap[PeptideEvidence.attrib['id']]=[PeptideEvidence.attrib['pre'],PeptideEvidence.attrib['post'],PeptideEvidence.attrib['dBSequence_ref']]

def read_one_pep_xml_msgf(input_file_str, psm_list,fidx):
    Map=dict()
    prepostMap=dict()
    pep2proMap(input_file_str,Map)
    prepostDBref(input_file_str,prepostMap)
    tree = ET.parse(input_file_str)
    root = tree.getroot()
    DataCollection=root.find('{http://psidev.info/psi/pi/mzIdentML/1.1}DataCollection')
    AnalysisData=DataCollection.find('{http://psidev.info/psi/pi/mzIdentML/1.1}AnalysisData')
    SpectrumIdentificationList=AnalysisData.find('{http://psidev.info/psi/pi/mzIdentML/1.1}SpectrumIdentificationList')
    for SpectrumIdentificationResult in SpectrumIdentificationList.iter('{http://psidev.info/psi/pi/mzIdentML/1.1}SpectrumIdentificationResult'):
        scan_charge=set()
        first=True
        one_scan=scan()
        one_scan.fidx = fidx
        for cvParam in SpectrumIdentificationResult.iter('{http://psidev.info/psi/pi/mzIdentML/1.1}cvParam'):
            one_scan.scan_number=cvParam.attrib['value']
        for SpectrumIdentificationItem in SpectrumIdentificationResult.iter('{http://psidev.info/psi/pi/mzIdentML/1.1}SpectrumIdentificationItem'):
            if cvParam.attrib['value']+'_'+SpectrumIdentificationItem.attrib['chargeState'] in scan_charge:
                scan_charge.add(cvParam.attrib['value']+'_'+SpectrumIdentificationItem.attrib['chargeState'])
            else:
                if first:
                    first=False
                    scan_charge.add(cvParam.attrib['value'] + '_' + SpectrumIdentificationItem.attrib['chargeState'])
                else:
                    psm_list.append(one_scan)
                    one_scan = scan()
                    one_scan.fidx = fidx
                    one_scan.scan_number = cvParam.attrib['value']
                    one_scan.charge=SpectrumIdentificationItem.attrib['chargeState']
            pep=peptide()
            pep.charge=SpectrumIdentificationItem.attrib['chargeState']
            one_scan.charge=SpectrumIdentificationItem.attrib['chargeState']
            one_scan.measured_mass=(float(SpectrumIdentificationItem.attrib['experimentalMassToCharge'])-W_proton)*int(pep.charge)
            pep.calculated_mass=(float(SpectrumIdentificationItem.attrib['calculatedMassToCharge'])-W_proton)*int(pep.charge)
            p='.'+SpectrumIdentificationItem.attrib['peptide_ref'].replace('Pep_','')+'.'
            for PeptideEvidenceRef in SpectrumIdentificationItem.iter('{http://psidev.info/psi/pi/mzIdentML/1.1}PeptideEvidenceRef'):
                pep.identified_pep=prepostMap[PeptideEvidenceRef.attrib['peptideEvidence_ref']][0]+p+prepostMap[PeptideEvidenceRef.attrib['peptideEvidence_ref']][1]
                pep.original_pep=SpectrumIdentificationItem.attrib['peptide_ref'].replace('Pep_','')
                seq=SpectrumIdentificationItem.attrib['peptide_ref'].replace('Pep_','').replace('M+16','*')
                mods=[]
                for i, c in enumerate(seq):
                    if c=='C':
                        mods.append(str(i+1)+'_S_57.021464')
                    elif c=='*':
                        mods.append(str(i+1)+'_V_15.994900')
                pep.modifications=','.join(mods)
                pep.prev_aa=prepostMap[PeptideEvidenceRef.attrib['peptideEvidence_ref']][0]
                pep.next_aa=prepostMap[PeptideEvidenceRef.attrib['peptideEvidence_ref']][1]
                pep.protein_list.append(Map[prepostMap[PeptideEvidenceRef.attrib['peptideEvidence_ref']][2]])
            for cvP in SpectrumIdentificationItem.iter('{http://psidev.info/psi/pi/mzIdentML/1.1}cvParam'):
                if cvP.attrib['name']=="MS-GF:EValue":
                    pep.EValue=-math.log(float(cvP.attrib['value']))
            one_scan.add_pep(pep)
        if len(one_scan.pep_list)>0:
            psm_list.append(one_scan)

def write_output(psm_list, output_file_str):
    with open(output_file_str, 'w') as fw:
        for v in psm_list:
            if v.pep_list:
                ll = v.get_winnownet_output()
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
    psm_list = []
    prefix=input_file.split('.')[0]
    read_one_pep_xml_msgf(input_file,psm_list,prefix)
    write_output(psm_list, output_file)

    print('All done.')

if __name__ == '__main__':
    sys.exit(main())
