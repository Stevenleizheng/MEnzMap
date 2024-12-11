import argparse
import os
import subprocess
import pandas as pd 
from Bio import SeqIO

def parameter_input():
    parser = argparse.ArgumentParser(description='iVirP (Integrative virome pipeline)')
    parser.add_argument('-1', '--input', help='Path to the fastq(.gz) file of read1')
    parser.add_argument('-o', '--output', help='Save file path', default=f'{os.getcwd()}')
    args = parser.parse_args()
    return args 

def get_sample_name(input1):
    name = os.path.basename(os.path.abspath(input1))
    name = name.split('_')[0]
    return name

def calculate_rpkm(output,name):
    print("Calculate rpkm")
    if os.path.exists(f"{output}/11.enzyme_rpkm") is True:
        subprocess.call([f"rm -rf {output}/11.enzyme_rpkm"], shell=True)
    subprocess.call([f"mkdir {output}/11.enzyme_rpkm"], shell=True) 
    map_data = pd.read_table(f'{output}/5.bins/{name}.read',header=None)
    total_reads = sum(map_data.iloc[:,2]) + sum(map_data.iloc[:,3])
    rpkm = (map_data.iloc[:,2] + map_data.iloc[:,3])*(1e9)/(map_data.iloc[:,1] * total_reads)
    rpkm = pd.concat([map_data.iloc[:,0], rpkm], axis=1)
    rpkm = rpkm.drop(rpkm.index[-1])
    new_column_names = rpkm.columns.tolist()
    new_column_names[0] = 'Accession'
    new_column_names[1] = name
    rpkm.columns = new_column_names
    high_bin_prodigal = pd.read_table(f'{output}/7.high_quality_bins/high_bin_prodigal.txt')
    high_bin_prodigal = high_bin_prodigal.rename(columns={'contig_name': 'Accession'})
    rpkm = pd.merge(rpkm, high_bin_prodigal, on='Accession', how='left')
    rpkm = rpkm.drop(columns=['sample_name'])
    rpkm = rpkm.rename(columns={'file_name': 'Bin Id'})
    rpkm = rpkm.dropna(subset=['Bin Id'])
    rpkm['Bin Id'] = rpkm['Bin Id'].str.replace(r'\.fa$', '', regex=True)
    taxo_result = pd.read_csv(f'{output}/6.checkm/{name}_qa_result.tsv',sep='\t')
    rpkm = pd.merge(rpkm, taxo_result, on='Bin Id', how='left')
    rpkm =rpkm.drop(columns=['# unique markers (of 43)'])
    rpkm =rpkm.drop(columns=['# multi-copy'])
    rpkm['Taxonomy'] = rpkm['Taxonomy'].str.split(';').str[-1].str.split('__').str[-1]
    profile_result = pd.read_csv(f'{output}/6.checkm/{name}profile_result.tsv',sep='\t')
    rpkm = pd.merge(rpkm, profile_result, on='Bin Id', how='left')
    rpkm =rpkm.drop(columns=['Bin size (Mbp)'])
    rpkm =rpkm.drop(columns=[f'{name}.sorted: mapped reads'])
    rpkm =rpkm.drop(columns=[f'{name}.sorted: % mapped reads'])
    rpkm =rpkm.drop(columns=[f'{name}.sorted: % binned populations'])
    # construct taxonomy: percentage
    tax_per_df = rpkm[['Taxonomy',f'{name}.sorted: % community']].drop_duplicates(subset='Taxonomy')
    taxo_list = list(tax_per_df['Taxonomy'])
    tax_per_dict = dict(zip(tax_per_df['Taxonomy'], tax_per_df[f'{name}.sorted: % community']))
    enzyme_data = pd.read_csv(f'{output}/10.fedkea/enzyme_result.csv')
    data = pd.DataFrame() 
    data['Accession'] = enzyme_data['Accession'].str.replace(r'_(\d+)$', '', regex=True)
    data['EC'] = enzyme_data['FinalResult']
    data = pd.merge(data, rpkm[['Accession', f'{name}','Taxonomy']], on='Accession', how='left')
    grouped_data = data[['EC',f'{name}']].groupby('EC', as_index=False).sum()
    store_list = []
    EC_list = list(grouped_data['EC'])
    for i in EC_list:
        length = len(tax_per_dict)
        a = dict(data[data['EC']==i]['Taxonomy'].value_counts())
        tmp = []
        for i in range(length):
            taxo_i = str(list(tax_per_dict.keys())[i])
            if taxo_i in a.keys():
                tmp.append(f"{a[taxo_i]}*{tax_per_dict[taxo_i]/100:.4f}")
            else:
                tmp.append(0)
        store_list.append(tmp)
    df = pd.DataFrame(store_list, columns=list(tax_per_dict.keys()))
    grouped_data = pd.concat([grouped_data,df],axis=1)
    print(grouped_data)
    grouped_data.to_csv(f'{output}/11.enzyme_rpkm/enzyme_contig_abundance_rpkm.csv',index=None)
        
def main():
    args = parameter_input()
    input1 = args.input
    output = args.output
    name = get_sample_name(input1)

    calculate_rpkm(output,name)

if __name__ =="__main__":
    main()