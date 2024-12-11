import argparse
import os
import subprocess
import pandas as pd 
from Bio import SeqIO
from utils import *

def parameter_input():
    parser = argparse.ArgumentParser(description='EnzMicroMap (Enzyme and Microbe Mapping in Metagenomics)')
    parser.add_argument('-1', '--input', help='Path to the fastq(.gz) file of read1')
    parser.add_argument('-2', '--Input', help='Path to the fastq(.gz) file of read2')
    parser.add_argument('-o', '--output', help='Save file path', default=f'{os.getcwd()}')
    parser.add_argument('-w', '--threads', type=int, help='Threads used to run this pipeline (default:1)', default=1)
    parser.add_argument('-x', '--host', help='Name(s) of Bowtie2 index(es) built before, please separate different names with "," when removing more than one hosts (for example: "human,mouse") (default: None, which means this step will be skipped)', default=None) 
    parser.add_argument('-s', '--specific', type=int, help='If you want to run a specific step in the pipeline, you can specify this parameter. 1: run_fastp, 2: run_bowtie2, 3: acquire_contigs, 4: run_bowtie2_build_index, 5: run_bins, 6: run_checkm, 7: select_bins, info_prodigal, 8: run_prodigal, select_prodigal_partial00', default=0)   
    parser.add_argument('-f', '--fast', type=int, help='Choose the pipeline\'s running speed: for a faster option (parameter: 1), select MEGAHIT to assemble contigs [with good quality]; for a slower option (parameter: 0, default), select SPAdes to assemble contigs [with better quality].', default=0)
    args = parser.parse_args()
    return args 

def get_sample_name(input1):
    name = os.path.basename(os.path.abspath(input1))
    name = name.split('_')[0]
    return name

def read_checkpoint(output):
    if os.path.exists(f"{output}/checkpoint.txt"):
        with open(f"{output}/checkpoint.txt", "r") as f:
            return f.read().strip()
    return None

def write_checkpoint(step,output):
    with open(f"{output}/checkpoint.txt", "w") as f:
        f.write(step)

def main():
    args = parameter_input()
    threads = args.threads
    input1 = args.input
    input2 = args.Input
    output = args.output
    host = args.host
    name = get_sample_name(input1)

    if args.specific == 0:
        last_checkpoint = read_checkpoint(output)
        if last_checkpoint is None:
            run_fastp(input1, input2, threads, output, name)
            write_checkpoint("1.fastp",output)

        last_checkpoint = read_checkpoint(output)
        if last_checkpoint == "1.fastp":
            index_path = [f"{sys.path[0]}/bowtie2_index/{host}/{host}"]
            run_bowtie2(output, threads, name, index_path)
            write_checkpoint("2.bowtie2",output)

        last_checkpoint = read_checkpoint(output)
        if last_checkpoint == "2.bowtie2":
            if args.fast == 0:
                run_spades(output, threads, name)
                write_checkpoint("3.contigs",output)
            else:
                run_megahit(output, threads, name)
                write_checkpoint("3.contigs",output)

        last_checkpoint = read_checkpoint(output)
        if last_checkpoint == "3.contigs":
            run_bowtie2_build_index(output, threads, name)
            write_checkpoint("4.index",output)  

        last_checkpoint = read_checkpoint(output)
        if last_checkpoint == "4.index":
            run_bins(output, threads, name)
            write_checkpoint("5.bins",output) 

        last_checkpoint = read_checkpoint(output)
        if last_checkpoint == "5.bins":
            run_checkm(output, threads, name)
            write_checkpoint("6.checkm",output) 

        last_checkpoint = read_checkpoint(output)
        if last_checkpoint == "6.checkm":
            high_id = select_bins(output, name)
            if not high_id:  # 检查列表是否为空
                print(f"{name} high_id is empty")
                exit()  # 终止程序
            info_prodigal(output,name)
            write_checkpoint("7.high_quality_bins",output)

        last_checkpoint = read_checkpoint(output)
        if last_checkpoint == "7.high_quality_bins":
            run_prodigal(output,name)
            select_prodigal_partial00(output,name)
            write_checkpoint("9.prodigal",output)
    
    elif args.specific == 1:
        run_fastp(input1, input2, threads, output, name)

    elif args.specific == 2:
        index_path = [f"{sys.path[0]}/bowtie2_index/{host}/{host}"]
        run_bowtie2(output, threads, name, index_path)        
    
    elif args.specific == 3:
        if args.fast == 0:
            run_spades(output, threads, name)
        else:
            run_megahit(output, threads, name)

    elif args.specific == 4:
        run_bowtie2_build_index(output, threads, name)
    
    elif args.specific == 5:
        run_bins(output, threads, name)
    
    elif args.specific == 6:  
        run_checkm(output, threads, name)

    elif args.specific == 7:
        high_id = select_bins(output, name)
        if not high_id:  # 检查列表是否为空
            print(f"{name} high_id is empty")
            exit()  # 终止程序
        info_prodigal(output,name)
    
    elif args.specific == 8:
        run_prodigal(output,name)
        select_prodigal_partial00(output,name)

if __name__ =="__main__":
    main()
