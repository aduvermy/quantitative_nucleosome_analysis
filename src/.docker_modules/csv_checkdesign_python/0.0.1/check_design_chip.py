#!/usr/bin/env python3

import sys


design=sys.argv[1] ## input csv  file
design_checked=sys.argv[2] ## output csv file: same as design input but make sure there is no space in the row
                            ## avoid error in nextflow process


def csv_parser(csv, csv_out):
    ERROR_STR = 'ERROR: Please check design file\n'
    def check_line_format(lspl):
        ## CHECK VALID NUMBER OF COLUMNS PER SAMPLE
        numCols = len(lspl)
        if numCols not in [2]:
            print (f"{ERROR_STR}Invalid number of columns (should be 2)!\nLine: {lspl}\nWARN: Colomn separator must be \t")
            sys.exit(1)

        IP,INPUT = lspl[0],lspl[1]
        if str(IP)==str(INPUT):
            print (f"{ERROR_STR}Same file specified as IP and INPUT!\nLine: {lspl}\n")
            sys.exit(1)
        return IP,INPUT

    HEADER = ['ip', 'input']

    header = csv.readline().strip().split('\t')
    if header != HEADER:
        print (f"{ERROR_STR} header:{header}!= {HEADER}\nWARN: Colomn separator must be \t")
        sys.exit(1)
    csv_out.write(f"ip\tinput\n")
    sampleMappingDict = {}
    for line in csv:
        line_sple = [elt.strip() for elt in line.strip().split("\t")]
        IP,INPUT=check_line_format(line_sple)
        IP=str(IP)
        INPUT=str(INPUT)
        ## CHECK UNICITY COUPLE [IP,INPUT]
        couple=f"{IP};{INPUT}"
        if couple not in sampleMappingDict:
            sampleMappingDict[couple]=""
        else:
            print (f"WARN: couple IP vs INPUT specified multiple times, migth be an error!\nLine: {line_sple}")
            sys.exit(1)
        csv_out.write(f"{IP}\t{INPUT}\n")
    return None



if __name__ == "__main__":
    # execute only if run as a script
    try:
        design_f=open(design,'r')
    except IOError:
        print('Cannot open', design)
        exit()
    try:
        output_f=open(design_checked,'w')
    except IOError:
        print('Cannot open', design_checked)
        exit()

    csv_parser(design_f,output_f)
