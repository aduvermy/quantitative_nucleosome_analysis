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
        if numCols not in [4]:
            print (f"{ERROR_STR}Invalid number of columns (should be 4)!\nLine: {lspl}\nWARN: Colomn separator must be \t")
            sys.exit(1)

        IP_w,WCE_w,IP_m,WCE_m = lspl[0],lspl[1],lspl[2],lspl[3]
        if str(IP_w)==str(IP_m) or str(WCE_w)==str(WCE_m) or str(IP_w)==str(WCE_w) or str(IP_m)==str(WCE_m) or str(IP_w)==str(WCE_m) or str(WCE_w)==str(IP_m) :
            print (f"{ERROR_STR}Same file specified multiple times on line:\n {lspl}\n")
            sys.exit(1)
        return IP_w,WCE_w,IP_m,WCE_m

    HEADER = ['IP_w', 'WCE_w','IP_m','WCE_m']

    header = csv.readline().strip().split('\t')
    if header != HEADER:
        print (f"{ERROR_STR} header:{header}!= {HEADER}\nWARN: Colomn separator must be \t")
        sys.exit(1)
    csv_out.write(f"IP_w\tWCE_w\tIP_m\tWCE_m\n")
    sampleMappingDict = {}
    for line in csv:
        line_sple = [elt.strip() for elt in line.strip().split("\t")]
        IP_w,WCE_w,IP_m,WCE_m=check_line_format(line_sple)
        IP_w=str(IP_w)
        WCE_w=str(WCE_w)
        IP_m=str(IP_m)
        WCE_=str(WCE_m)
        ## CHECK UNICITY EXPERIMENT [IP_w,WCE_w,IP_m,WCE_m]
        exp=f"{IP_w};{WCE_w};{IP_m};{WCE_m}"
        #exp=str(str(IP_w)+";"+str(WCE_w)+";"+str(IP_m)+";"+str(WCE_m))
        if exp not in sampleMappingDict:
            sampleMappingDict[exp]=""
        else:
            print (f"WARN: experiment specified multiple times, migth be an error!\nLine: {line_sple}")
            sys.exit(1)
        csv_out.write(f"{IP_w}\t{WCE_w}\t{IP_m}\t{WCE_m}\n")
        #csv_out.write(str(IP_w)+"\t"+str(WCE_w)+"\t"+str(IP_m)+"\t"+str(WCE_m)+"\n")
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
