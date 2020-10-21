#!/usr/bin/env python3

import sys


wig_file_path=sys.argv[1] ## input wig file
wig_norm_file_path=sys.argv[2] ## output bed file




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


def build_wig_average_normalized(norm_factor, wig_file, wig_norm_file):
    for line in wig_file:
        if 'chrom' in line : ## new chromosome
            wig_norm_file.write(line)
            continue
        occupancy=float(line.strip()) ##occupancy value
        occupancy_norm=float(occupancy/norm_factor)
        wig_norm_file.write(f"{occupancy_norm}\n")
    return None




if __name__ == "__main__":
    # execute only if run as a script
    try:
        wig_f=open(wig_file_path,'r')
    except IOError:
        print('Cannot open', wig_file_path)
        exit()

    try:
        wig_norm_f=open(wig_norm_file_path,'w')
    except IOError:
        print('Cannot open', wig_norm_file_path)
        exit()

    mean_occ=mean_occupancy(wig_f)
    wig_f.seek(0)
    build_wig_average_normalized(mean_occ, wig_f, wig_norm_f)
    wig_f.close()
    wig_norm_f.close()
