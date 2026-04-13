import os
import subprocess
import sys
import pandas as pd
from Bio import SeqIO

def run_in_conda_env(cmd_list, env_name, error_msg="Error"):
    """在指定 conda 环境中运行命令"""
    full_cmd = ["conda", "run", "-n", env_name] + cmd_list
    print(f"Running in {env_name} environment: {' '.join(cmd_list)}")
    ret = subprocess.call(full_cmd)
    if ret != 0:
        sys.exit(f"{error_msg} in {env_name} environment")
    return ret


def run_fastp(input1, input2, threads, output, name):
    """Perform quality control and adapter trimming on raw sequencing data.
    
    Args:
        input1 (str): Path to read1 FASTQ file
        input2 (str): Path to read2 FASTQ file
        threads (int): Number of CPU threads to use
        output (str): Output directory path
        name (str): Sample name for output files
    """
    print("Run fastp")
    if os.path.exists(f"{output}/1.fastp") is True:
        subprocess.call([f"rm -rf {output}/1.fastp"], shell=True)
    subprocess.call([f"mkdir {output}/1.fastp"], shell=True)
    ret = subprocess.call([f"fastp -i {input1} -I {input2} -o {output}/1.fastp/{name}_trim_1.fastq.gz -O {output}/1.fastp/{name}_trim_2.fastq.gz --detect_adapter_for_pe -h {output}/1.fastp/fastp_report.html -w {threads} -j {output}/1.fastp/fastp_report.json"], shell=True)
    if ret != 0:
        print("Warning: fastqc error")

def run_bowtie2(output, threads, name, index_path):
    """Perform host DNA removal using Bowtie2 alignment.
    
    Args:
        output (str): Output directory path
        threads (int): CPU threads for alignment
        name (str): Sample identifier
        index_path (list): Bowtie2 index paths for host genomes
    """
    if len(index_path) == 0:
        print("No need for running bowtie2")
    else:
        print("Run bowtie2")
    if os.path.exists(f"{output}/2.bowtie2") is True:
        subprocess.call([f"rm -rf {output}/2.bowtie2"], shell=True)
    subprocess.call([f"mkdir {output}/2.bowtie2"], shell=True)
    subprocess.call([f"mv {output}/1.fastp/{name}_trim_1.fastq.gz {output}/2.bowtie2/{name}_trim_1.fastq.gz"], shell=True)
    subprocess.call([f"mv {output}/1.fastp/{name}_trim_2.fastq.gz {output}/2.bowtie2/{name}_trim_2.fastq.gz"], shell=True)
    for na in index_path:
        ret = subprocess.call([f"bowtie2 -p {threads} -x {na} -1 {output}/2.bowtie2/{name}_trim_1.fastq.gz -2 {output}/2.bowtie2/{name}_trim_2.fastq.gz --un-conc-gz {output}/2.bowtie2/tmp > {output}/2.bowtie2/tmp.sam"], shell=True)
        if ret != 0:
            sys.exit("Error: bowtie2 error")
        subprocess.call([f"mv {output}/2.bowtie2/tmp.1 {output}/2.bowtie2/{name}_1.fastq.gz"], shell=True)
        subprocess.call([f"mv {output}/2.bowtie2/tmp.2 {output}/2.bowtie2/{name}_2.fastq.gz"], shell=True)
    if len(index_path) != 0:
        subprocess.call([f"rm {output}/2.bowtie2/tmp.sam"], shell=True)

def run_spades(output, threads, name):
    """Execute metagenomic assembly using SPAdes with quality-controlled reads.
    
    Args:
        output (str): Output directory path
        threads (int): CPU threads for assembly
        name (str): Sample identifier
    """
    print("Run spades")
    if os.path.exists(f"{output}/3.contigs") is True:
        subprocess.call([f"rm -rf {output}/3.contigs"], shell=True)
    subprocess.call([f"mkdir {output}/3.contigs"], shell=True)
    ret = subprocess.call([f"spades.py -1 {output}/2.bowtie2/{name}_1.fastq.gz -2 {output}/2.bowtie2/{name}_2.fastq.gz -t {threads} -o {output}/3.contigs --meta -m 500"], shell=True)
    if ret != 0:
        sys.exit("Error: spades error")
    subprocess.call([f"find {output}/3.contigs -mindepth 1 ! -name 'contigs.fasta' -delete"], shell=True)

def run_megahit(output, threads, name):
    """Execute metagenomic assembly using MEGAHIT for rapid contig generation.
    
    Args:
        output (str): Output directory path
        threads (int): CPU threads for assembly
        name (str): Sample identifier
    """
    print("Run megahit")
    if os.path.exists(f"{output}/3.contigs") is True:
        subprocess.call([f"rm -rf {output}/3.contigs"], shell=True)
    # subprocess.call([f"mkdir {output}/3.contigs"], shell=True)
    ret = subprocess.call([f"megahit -1 {output}/2.bowtie2/{name}_1.fastq.gz -2 {output}/2.bowtie2/{name}_2.fastq.gz -t {threads} -o {output}/3.contigs"], shell=True)
    if ret != 0:
        sys.exit("Error: megahit error")
    subprocess.call([f"find {output}/3.contigs -mindepth 1 ! -name 'final.contigs.fa' -delete"], shell=True)
    subprocess.call([f"mv {output}/3.contigs/final.contigs.fa {output}/3.contigs/contigs.fasta"], shell=True)

def run_bowtie2_build_index(output, threads, name):
    """Build Bowtie2 index from assembled contigs for read alignment.
    
    Args:
        output (str): Output directory path
        threads (int): CPU threads for index construction
        name (str): Sample identifier (unused in current implementation)
    """
    print("Run bowtie2_build index")
    if os.path.exists(f"{output}/4.index") is True:
        subprocess.call([f"rm -rf {output}/4.index"], shell=True)
    subprocess.call([f"mkdir {output}/4.index"], shell=True)
    ret = subprocess.call([f"bowtie2-build -f {output}/3.contigs/contigs.fasta {output}/4.index/final --threads {threads}"], shell=True)
    if ret != 0:
        sys.exit("Error: bowtie2-build error") 

def run_bins_metabat2(output, threads, name):
    """Execute metagenomic binning process using Bowtie2 and MetaBAT2.
    
    Args:
        output (str): Output directory path
        threads (int): CPU threads for parallel processing
        name (str): Sample identifier
    """
    print("Run bins")
    if os.path.exists(f"{output}/5.bins") is True:
        subprocess.call([f"rm -rf {output}/5.bins"], shell=True)
    subprocess.call([f"mkdir {output}/5.bins"], shell=True)
    ret = subprocess.call([f"bowtie2 -1 {output}/2.bowtie2/{name}_1.fastq.gz -2 {output}/2.bowtie2/{name}_2.fastq.gz -p {threads} -x {output}/4.index/final -S {output}/5.bins/{name}.sam"], shell=True)
    if ret != 0:
        sys.exit("Error: bowtie2 error")
    ret = subprocess.call([f"samtools view -@ {threads} -b -S {output}/5.bins/{name}.sam -o {output}/5.bins/{name}.bam"], shell=True)
    if ret != 0:
        sys.exit("Error: samtools view  error")
    ret = subprocess.call([f"samtools sort -@ {threads} -l 9 -O BAM {output}/5.bins/{name}.bam -o {output}/5.bins/{name}.sorted.bam"], shell=True)
    if ret != 0:
        sys.exit("Error: samtools sort error")
    ret = subprocess.call([f"samtools index -@ {threads} {output}/5.bins/{name}.sorted.bam {output}/5.bins/{name}.sorted.bam.bai"], shell=True)
    if ret != 0:
        sys.exit("Error: samtools index error")
    ret = subprocess.call([f"samtools idxstats -@ {threads} {output}/5.bins/{name}.sorted.bam > {output}/5.bins/{name}.read"], shell=True)
    if ret != 0:
        sys.exit("Error: samtools idxstats error")
    ret = subprocess.call([f"jgi_summarize_bam_contig_depths --outputDepth {output}/5.bins/{name}.depth.txt {output}/5.bins/{name}.sorted.bam"], shell=True)
    if ret != 0:
        sys.exit("Error: jgi_summarize_bam_contig_depths error")
    ret = subprocess.call([f"metabat2 -m 1500 -t {threads} -i {output}/3.contigs/contigs.fasta -a {output}/5.bins/{name}.depth.txt -o {output}/5.bins/{name}_bin"], shell=True)
    if ret != 0:
        sys.exit("Error: metabat2 error")

def run_samtools(output, threads, name):
    """Process SAM/BAM files for metagenomic binning workflow.
    
    Args:
        output (str): Output directory path
        threads (int): CPU threads for processing
        name (str): Sample identifier
    """
    print("Run samtools")
    if os.path.exists(f"{output}/5.bins") is True:
        subprocess.call([f"rm -rf {output}/5.bins"], shell=True)
    subprocess.call([f"mkdir {output}/5.bins"], shell=True)
    if os.path.exists(f"{output}/bam_folder") is True:
        subprocess.call([f"rm -rf {output}/bam_folder"], shell=True)
    subprocess.call([f"mkdir {output}/bam_folder"], shell=True)
    ret = subprocess.call([f"bowtie2 -1 {output}/2.bowtie2/{name}_1.fastq.gz -2 {output}/2.bowtie2/{name}_2.fastq.gz -p {threads} -x {output}/4.index/final -S {output}/5.bins/{name}.sam"], shell=True)
    if ret != 0:
        sys.exit("Error: bowtie2 error")
    ret = subprocess.call([f"samtools view -@ {threads} -b -S {output}/5.bins/{name}.sam -o {output}/5.bins/{name}.bam"], shell=True)
    if ret != 0:
        sys.exit("Error: samtools view  error")
    ret = subprocess.call([f"samtools sort -@ {threads} -l 9 -O BAM {output}/5.bins/{name}.bam -o {output}/5.bins/{name}.sorted.bam"], shell=True)
    if ret != 0:
        sys.exit("Error: samtools sort error")
    ret = subprocess.call([f"samtools index -@ {threads} {output}/5.bins/{name}.sorted.bam {output}/5.bins/{name}.sorted.bam.bai"], shell=True)
    if ret != 0:
        sys.exit("Error: samtools index error")
    ret = subprocess.call([f"samtools idxstats -@ {threads} {output}/5.bins/{name}.sorted.bam > {output}/5.bins/{name}.read"], shell=True)
    if ret != 0:
        sys.exit("Error: samtools idxstats error")
    ret = subprocess.call([f"cp {output}/5.bins/{name}.sorted.bam {output}/bam_folder/"], shell=True)
    if ret != 0:
        sys.exit("Error: cp bam file error")
    ret = subprocess.call([f"cp {output}/5.bins/{name}.sorted.bam.bai {output}/bam_folder/"], shell=True)
    if ret != 0:
        sys.exit("Error: cp bam.bai file error")

def run_bins_comebin(output, threads, name):
    """Execute metagenomic binning process using ComeBin (在 comebin_env 中运行)"""
    print("Run bins with ComeBin")
    if os.path.exists(f"{output}/5.tmp") is True:
        subprocess.call([f"rm -rf {output}/5.tmp"], shell=True)

    # 使用 conda run 在 comebin_env 中执行 run_comebin.sh
    comebin_cmd = [
        "run_comebin.sh", 
        "-a", f"{output}/3.contigs/contigs.fasta", 
        "-o", f"{output}/5.tmp", 
        "-p", f"{output}/bam_folder/", 
        "-t", str(threads)
    ]
    run_in_conda_env(comebin_cmd, "comebin_env", "run_comebin error")

    subprocess.call([f"rm -rf {output}/bam_folder"], shell=True)

    # 重命名（在主环境执行即可）
    bin_dir = f"{output}/5.tmp/comebin_res/comebin_res_bins"
    files = sorted(
        [f for f in os.listdir(bin_dir) if f.endswith('.fa')],
        key=lambda x: int(x.split('.')[0]) if x.split('.')[0].isdigit() else 999999
    )
    
    for i, old_name in enumerate(files, start=1):
        old_path = os.path.join(bin_dir, old_name)
        new_path = os.path.join(bin_dir, f"test_bin.{i}.fa")
        print(f"renamed '{old_path}' -> '{new_path}'")
        os.rename(old_path, new_path)

    subprocess.call([f"mv {output}/5.tmp/comebin_res/comebin_res_bins/*.fa {output}/5.bins/"], shell=True)
    subprocess.call([f"rm -rf {output}/5.tmp"], shell=True)

# def run_bins_comebin(output, threads, name):
#     """Execute metagenomic binning process using ComeBin.
    
#     Args:
#         output (str): Output directory path
#         threads (int): CPU threads for parallel processing
#         name (str): Sample identifier
#     """
#     print("Run bins")
#     if os.path.exists(f"{output}/5.tmp") is True:
#         subprocess.call([f"rm -rf {output}/5.tmp"], shell=True)
#     ret = subprocess.call([f"run_comebin.sh -a {output}/3.contigs/contigs.fasta -o {output}/5.tmp -p {output}/bam_folder/ -t {threads}"], shell=True)
#     if ret != 0:
#         sys.exit("Error: run_comebin error")
#     ret = subprocess.call([f"rm -rf {output}/bam_folder"], shell=True)
#     rename_cmd = f'''
#         i=1
#         for f in $(ls -v "{output}/5.tmp/comebin_res_bins/"*.fa); do
#             mv -v "$f" "{output}/5.tmp/comebin_res_bins/test_bin.$i.fa"
#             ((i++))
#         done
#     '''
#     ret = subprocess.call(rename_cmd, shell=True)
#     if ret != 0:
#         sys.exit("Error: rename bins failed")
#     subprocess.call([f"mv {output}/5.tmp/comebin_res_bins/ {output}/5.bins/"], shell=True)
#     subprocess.call([f"rm -rf {output}/5.tmp"], shell=True)
    
def run_checkm(output, threads, name):
    """Execute CheckM quality assessment workflow for metagenomic bins.
    
    Args:
        output (str): Output directory path
        threads (int): CPU threads for parallel processing
        name (str): Sample identifier
    
    Workflow Steps:
    1. Lineage-specific marker analysis
    2. Phylogenetic tree construction
    3. Quality assessment reporting
    4. Coverage profiling
    """
    print("Run checkm")
    if os.path.exists(f"{output}/6.tmp") is True:
        subprocess.call([f"rm -rf {output}/6.tmp"], shell=True)
    subprocess.call([f"mkdir {output}/6.tmp"], shell=True)
    if os.path.exists(f"{output}/6.checkm") is True:
        subprocess.call([f"rm -rf {output}/6.checkm"], shell=True)
    subprocess.call([f"mkdir {output}/6.checkm"], shell=True)
    ret = subprocess.call([f"checkm lineage_wf -t {threads} -x fa --nt --tab_table -f {output}/6.tmp/{name}_bins_qa.txt {output}/5.bins/ {output}/6.tmp/"], shell=True)
    if ret != 0:
        sys.exit("Error: checkm error") 
    subprocess.call([f"mv {output}/6.tmp/lineage.ms {output}/6.checkm/"], shell=True)
    subprocess.call([f"mv {output}/6.tmp/{name}_bins_qa.txt {output}/6.checkm/"], shell=True)
    subprocess.call([f"rm -rf {output}/6.tmp"], shell=True)
    subprocess.call([f"mkdir {output}/6.tmp"], shell=True)
    ret = subprocess.call([f"checkm tree -t {threads} -x fa --nt {output}/5.bins/ {output}/6.tmp/"], shell=True)
    if ret != 0:
        sys.exit("Error: checkm error") 
    ret = subprocess.call([f"checkm tree_qa {output}/6.tmp/ -o 1 -f {output}/6.checkm/{name}_qa_result.tsv --tab_table"], shell=True)
    if ret != 0:
        sys.exit("Error: checkm error") 
    subprocess.call([f"rm -rf {output}/6.tmp"], shell=True)
    subprocess.call([f"mkdir {output}/6.tmp"], shell=True)
    ret = subprocess.call([f"checkm coverage -t {threads} -x fa {output}/5.bins/ {output}/6.tmp/coverage.tsv {output}/5.bins/{name}.sorted.bam"], shell=True)
    if ret != 0:
        sys.exit("Error: checkm error") 
    ret = subprocess.call([f"checkm profile {output}/6.tmp/coverage.tsv -f {output}/6.checkm/{name}_profile_result.tsv --tab_table"], shell=True)
    if ret != 0:
        sys.exit("Error: checkm error") 
    subprocess.call([f"rm -rf {output}/6.tmp"], shell=True)

# def run_checkm2(output, threads, name):
#     """Execute CheckM2 quality assessment workflow for metagenomic bins.
    
#     Args:
#         output (str): Output directory path
#         threads (int): CPU threads for parallel processing
#         name (str): Sample identifier
#     """
#     print("Run checkm2")
#     if os.path.exists(f"{output}/6.tmp") is True:
#         subprocess.call([f"rm -rf {output}/6.tmp"], shell=True)
#     subprocess.call([f"mkdir {output}/6.tmp"], shell=True)
#     ret = subprocess.call([f"checkm2 predict  -x .fa --threads {threads} --input {output}/5.bins/ --output-directory{output}/6.tmp/"], shell=True)
#     if ret != 0:
#         sys.exit("Error: checkm2 error") 
#     data = pd.read_table(f'{output}/6.tmp/quality_report.tsv', sep='\t')
#     data = data.rename(columns={"Name": "Bin Id"})
#     subprocess.call([f"rm -rf {output}/6.checkm/{name}_bins_qa.txt"], shell=True)
#     data.to_csv(f'{output}/6.checkm/{name}_bins_qa.txt', sep='\t', index=False)

def run_checkm2(output, threads, name):
    """Execute CheckM2 quality assessment (在 checkm2_env 中运行)"""
    print("Run checkm2")
    if os.path.exists(f"{output}/6.tmp") is True:
        subprocess.call([f"rm -rf {output}/6.tmp"], shell=True)
    subprocess.call([f"mkdir {output}/6.tmp"], shell=True)

    # 使用 conda run 在 checkm2 中执行 checkm2
    checkm2_cmd = [
        "checkm2", "predict",
        "-x", ".fa",
        "--threads", str(threads),
        "--input", f"{output}/5.bins/",
        "--output-directory", f"{output}/6.tmp/",
        "--database_path", f"{sys.path[0]}/CheckM2_database/uniref100.KO.1.dmnd"
    ]
    run_in_conda_env(checkm2_cmd, "checkm2", "checkm2 error")

    # 后续处理仍在主环境（pandas 操作）
    data = pd.read_table(f'{output}/6.tmp/quality_report.tsv', sep='\t')
    data = data.rename(columns={"Name": "Bin Id"})
    subprocess.call([f"rm -rf {output}/6.checkm/{name}_bins_qa.txt"], shell=True)
    data.to_csv(f'{output}/6.checkm/{name}_bins_qa.txt', sep='\t', index=False)
    subprocess.call([f"rm -rf {output}/6.tmp"], shell=True)

def select_bins(output,name):
    """Filter and organize metagenomic bins based on CheckM quality metrics.
    
    Args:
        output (str): Output directory path
        name (str): Sample identifier
        
    Returns:
        list: IDs of high-quality bins
        
    Quality Thresholds:
        High: ≥90% completeness, ≤5% contamination
        Medium: ≥50% completeness, ≤10% contamination
    """
    print("Select high quality and medium quality bins")
    if os.path.exists(f"{output}/7.high_quality_bins") is True:
        subprocess.call([f"rm -rf {output}/7.high_quality_bins"], shell=True)
    subprocess.call([f"mkdir {output}/7.high_quality_bins"], shell=True)
    if os.path.exists(f"{output}/8.medium_quality_bins") is True:
        subprocess.call([f"rm -rf {output}/8.medium_quality_bins"], shell=True)
    subprocess.call([f"mkdir {output}/8.medium_quality_bins"], shell=True)
    evaluation = pd.read_table(f'{output}/6.checkm/{name}_bins_qa.txt')
    high_id = list(evaluation[(evaluation['Completeness'] >= 90.00) & (evaluation['Contamination'] <= 5.00)]['Bin Id'])
    medium_id = list(evaluation[(evaluation['Completeness'] >= 50.00) & (evaluation['Contamination'] <= 10.00)]['Bin Id'])
    for i in high_id:
        ret = subprocess.call([f"cp {output}/5.bins/{i}.fa {output}/7.high_quality_bins"], shell=True)
        if ret != 0:
            sys.exit("Error: cp high quality bins error") 
    for i in medium_id:
        ret = subprocess.call([f"cp {output}/5.bins/{i}.fa {output}/8.medium_quality_bins"], shell=True)
        if ret != 0:
            sys.exit("Error: cp medium quality bins error") 
    return  high_id

def info_prodigal(output,name):
    """Generate metadata mapping between predicted genes and their source contigs/bins.
    
    Args:
        output (str): Output directory path
        name (str): Sample identifier
    
    Output File:
        Creates 'high_bin_prodigal.txt' with columns:
        - sample_name: Identifier from sequencing run
        - file_name: Source bin filename
        - contig_name: Original contig identifier
    """
    print("information of prodigal contigs and bins")
    with open(f'{output}/7.high_quality_bins/high_bin_prodigal.txt', 'w') as out_f:
        out_f.write('sample_name\tfile_name\tcontig_name\n')          
        for filename in os.listdir(f'{output}/7.high_quality_bins'):
            if filename.endswith('.fasta') or filename.endswith('.fa'):  
                filepath = os.path.join(f'{output}/7.high_quality_bins', filename)
                with open(filepath, 'r') as fasta_file:
                    for record in SeqIO.parse(fasta_file, 'fasta'):
                        out_f.write(f'{name}\t{filename}\t{record.id}\n')

def run_prodigal(output,name):
    """Execute Prodigal gene prediction on high-quality metagenomic bins.
    
    Args:
        output (str): Output directory path
        name (str): Sample identifier
        
    Output Files:
        - GFF3 file: Structural annotations
        - FNA file: Predicted gene nucleotide sequences
        - FAA file: Predicted protein translations
    """
    print("Run prodigal")
    if os.path.exists(f"{output}/9.prodigal") is True:
        subprocess.call([f"rm -rf {output}/9.prodigal"], shell=True)
    subprocess.call([f"mkdir {output}/9.prodigal"], shell=True)  
    subprocess.call([f"cat {output}/7.high_quality_bins/*.fa > {output}/7.high_quality_bins/{name}_high_bin.fasta "], shell=True)  
    ret = subprocess.call([f"prodigal -i {output}/7.high_quality_bins/{name}_high_bin.fasta -f gff -o {output}/9.prodigal/{name}_high_bin.gff3 -d {output}/9.prodigal/{name}_high_bin_gene.fna -a {output}/9.prodigal/{name}_high_bin_protein.faa -p meta"], shell=True)
    if ret != 0:
        sys.exit("Error: prodigal error") 

def select_prodigal_partial00(output,name):
    """Filter and save complete gene predictions from Prodigal results.
    
    Args:
        output (str): Output directory path
        name (str): Sample identifier
    
    Biological Significance:
        Selects only complete genes (no partial fragments) based on Prodigal's 
        'partial=00' flag in sequence descriptions, indicating full CDS calls
        from start to stop codon without truncation.
    """
    with open(f'{output}/9.prodigal/{name}_high_bin_protein_partial00.fasta', 'w') as out_f: 
        with open(f'{output}/9.prodigal/{name}_high_bin_protein.faa', 'r') as fasta_file: 
            for record in SeqIO.parse(fasta_file, 'fasta'):
                if 'partial=00' in record.description:
                    SeqIO.write(record, out_f, 'fasta')
    with open(f'{output}/9.prodigal/{name}_high_bin_gene_partial00.fasta', 'w') as out_f: 
        with open(f'{output}/9.prodigal/{name}_high_bin_gene.fna', 'r') as fasta_file: 
            for record in SeqIO.parse(fasta_file, 'fasta'):
                if 'partial=00' in record.description:
                    SeqIO.write(record, out_f, 'fasta')    
    
    subprocess.call([f"rm -rf {output}/5.bins/{name}.bam"], shell=True)
    subprocess.call([f"rm -rf {output}/5.bins/{name}.sorted.bam"], shell=True)
    subprocess.call([f"rm -rf {output}/5.bins/{name}.sorted.bam.bai"], shell=True)
    subprocess.call([f"rm -rf {output}/5.bins/{name}.sam"], shell=True)


# def run_fedkea(output,name):
#     """Execute FEDKEA enzyme annotation pipeline on complete gene predictions.
    
#     Args:
#         output (str): Output directory path
#         name (str): Sample identifier
        
#     Output:
#         Creates directory structure:
#         - 10.fedkea/data/: Intermediate data storage
#         - 10.fedkea/result/: Final annotation results
        
#     Parameters:
#         -b 32: Batch size for parallel processing
#         -a 1: Algorithm selection for enzyme prediction
#     """
#     print("Run fedkea")
#     if os.path.exists(f"{output}/10.fedkea") is True:
#         subprocess.call([f"rm -rf {output}/10.fedkea"], shell=True)
#     subprocess.call([f"mkdir {output}/10.fedkea"], shell=True)  
#     ret = subprocess.call([f"python {sys.path[0]}/FEDKEA/main.py -i {output}/9.prodigal/{name}_high_bin_protein_partial00.fasta -b 32 -d {output}/10.fedkea/data/ -o {output}/10.fedkea/result/ "], shell=True)
#     if ret != 0:
#         sys.exit("Error: fedkea error") 

