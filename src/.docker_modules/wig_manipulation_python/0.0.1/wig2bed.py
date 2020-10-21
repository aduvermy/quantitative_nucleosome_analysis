#!/usr/bin/env python3

import sys


wig_file_path=sys.argv[1] ## input wig file
length_min_NDR=int(sys.argv[2]) ## min length NDR
cutoff=float(sys.argv[3]) ##cutoff use to define ND position
bed_file_path=sys.argv[4] ## output bed file


def mean_occupancy(wig_file):
    occupancy_per_position=()
    sum_occupancy=0
    nb_position=0
    for line in wig_file:
        if 'chrom' in line : ## new chromosome
            continue
        occupancy=float(line.strip()) ##occupancy value
        sum_occupancy+=occupancy
        nb_position+=1
    return float(sum_occupancy/nb_position)



def wig_parser(wig_file, L_min, c_max):
    NDR=[] ##list [chromosome, position debut, position fin] NDR
    NDR_bool="FALSE" ## NDR boolean init to false
    length_NDR=0 ##init length NDR to 0

    def upload_NDR_list(NDR, chrom, start, length): ##add new NDR to NDR_list
        end=start+length-1 #end position
        NDR.append([chrom,start,end]) ##NDR chromosome and position list

    for line in wig_file:
        if 'chrom' in line : ## new chromosome
            if NDR_bool=="TRUE" and length_NDR >= L_min: ##about previous chromosome check if we were in NDR
                upload_NDR_list(NDR,chromosome,pos_start,length_NDR)
            try:
                chromosome=line.strip().split(' ')[1].split('=')[1] ##chromosome name
            except IndexError:
                print("File format might be not suitable")
                print("Incorrect parsing line:", line.strip())
            NDR_bool="FALSE" ##init NDR boolean to false
            length_NDR=0 ##init NDR length to 0
            k=1 ##chromosome start position init to 1
            continue

        occupancy=float(line.strip()) ##occupancy value
        if float(occupancy) < c_max: ##occupancy test
            length_NDR+=1 ## NDR length increase
            if NDR_bool=="FALSE": ##check if we start NDR
                pos_start=k ##NDR start position init to 1
                NDR_bool="TRUE" ##we are in NDR

        else: ##position k doesn't satisfy occupancy test

            if NDR_bool=="TRUE" and length_NDR >= L_min: ##check if we were in NDR and NDR long enough
                upload_NDR_list(NDR,chromosome,pos_start, length_NDR)

            length_NDR=0 ##init NDR length to 0
            NDR_bool="FALSE" ##init NDR boolean to false

        k+=1 ##next position
    return NDR






def build_bed_file(NDR_l, bed_file):
    for elt in NDR_l:
        chromosome=elt[0]
        start=elt[1]
        end=elt[2]
        bed_file.write(f"{chromosome}\t{start}\t{end}\n")
    return None




if __name__ == "__main__":
    # execute only if run as a script
    try:
        wig_f=open(wig_file_path,'r')
    except IOError:
        print('Cannot open', wig_file_path)

    try:
        bed_f=open(bed_file_path,'w')
    except IOError:
        print('Cannot open', bed_file_path)


    NDR_list=wig_parser(wig_f, length_min_NDR, cutoff)
    build_bed_file(NDR_list, bed_f)
    wig_f.close()
    bed_f.close()
