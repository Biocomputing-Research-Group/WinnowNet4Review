import xml.etree.ElementTree as ET
from xml.dom import minidom
import sys
import getopt

def pretty_xml(elem):
    rough = ET.tostring(elem, 'utf-8')
    dom   = minidom.parseString(rough)
    pretty = dom.documentElement.toprettyxml(indent="  ")
    lines = [line for line in pretty.splitlines() if line.strip()]
    return "\n".join(lines)

def reconstruct_xml(input_file_str,rescore_file,output_file):
    rescore_list=[]
    with open(rescore_file) as f:
        for line in f:
            rescore_list.append(line.strip())
    ns = {'pep': 'http://regis-web.systemsbiology.net/pepXML'}
    ET.register_namespace('', 'http://regis-web.systemsbiology.net/pepXML')
    tree = ET.parse(input_file_str)
    root = tree.getroot()
    for sq_id, sq in enumerate(root.findall('.//pep:spectrum_query', ns)):
        hit = sq.find('.//pep:search_hit', ns)
        analysis = ET.Element('analysis_result', analysis='peptideprophet')
        pp=ET.SubElement(analysis, 'peptideprophet_result', probability=rescore_list[sq_id], all_ntt_prob="(0,0,0)")
        score_sum = ET.SubElement(pp, 'search_score_summary')
        params = {
            'fval': '0.0',
            'ntt': '2',
            'nmc': '0',
            'massd': '0',
            'isomassd': '0'
        }
        for name, val in params.items():
            ET.SubElement(score_sum, 'parameter',
                          name=name, value=val)
        hit.append(analysis)

    with open(output_file, 'w') as f:
        f.write(open('ProphetHeader.txt').read())
        for idx, sq in enumerate(root.findall('.//pep:spectrum_query', ns), start=1):
            xml_str = pretty_xml(sq)
            f.write(xml_str)
        f.write('\n</msms_run_summary>\n</msms_pipeline_analysis>\n')


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
    reconstruct_xml(input_file,rescore_file,output_file)