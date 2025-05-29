import pandas as pd
import os
import glob

def get_mutation_type(info_column):
    mutation_type = "SNV"
    csq_field = [field.split('|') for field in info_column.split(';') if field.startswith('CSQ=')]
    if csq_field:
        for csq in csq_field[0]:
            if 'insertion' in csq or 'deletion' in csq:
                mutation_type = "Indel"
    return mutation_type


def process_quantumclone_results(output_dir, vcf, sample_name):
    """
    Processes QuantumClone results to extract clonal mutations for a specific sample.

    Args:
        output_dir (str): Path to directory with sample subdirectories.
        vcf (str): Path to root directory with original VCFs, organized by tumor type or flat.
        sample_name (str): Name of the sample to process.
    
    Saves:
        truncalmutations.txt in the sample's directory.
    """
    sample_path = os.path.join(output_dir, sample_name)

    try:
        filtered_path = os.path.join(sample_path, "filtered.csv")
        cluster_path = os.path.join(sample_path, "clustering.csv")
        centers_path = os.path.join(sample_path, "centers.csv")

        filtered = pd.read_csv(filtered_path, index_col=0)
        clusters = pd.read_csv(cluster_path, index_col=0)
        centers = pd.read_csv(centers_path, index_col=0)

        clonal_cluster = int(centers["X..i.."].idxmax())
        filtered = filtered.sort_values(by="id")
        filtered['Cluster'] = clusters['Number']
        quantum_clonal = filtered[filtered['Cluster'] == clonal_cluster]
        quantum_clonal = quantum_clonal.rename(columns={'Chr': '#CHROM', 'Start': 'POS'})

        vcf_files = []
        if os.path.isfile(vcf):
            # If a specific VCF file was provided
            vcf_files = [vcf]
        elif os.path.isdir(vcf):
            # If a directory of VCF files was provided
            vcf_files = glob.glob(os.path.join(vcf, f"*{sample_name}*.vcf"))
        
        if not vcf_files:
            print(f"No VCF file found for {sample_name}, skipping.")
            return

        vcf_data = {}
        with open(vcf_files[0], 'r') as vcf_file:
            for line in vcf_file:
                if line.startswith('#'):
                    continue
                fields = line.strip().split('\t')
                chrom, pos = fields[0].replace("chr", ""), fields[1]
                info_column = fields[7]
                mutation_type = get_mutation_type(info_column)
                vcf_data[(chrom, pos)] = mutation_type

        mutation_types = []
        for _, row in quantum_clonal.iterrows():
            chrom = str(row['#CHROM']).replace("chr", "")
            pos = str(row['POS'])
            mutation_types.append(vcf_data.get((chrom, pos), "SNV"))

        quantum_clonal['MutationType'] = mutation_types
        out_file = os.path.join(sample_path, "truncalmutations.txt")
        quantum_clonal.to_csv(out_file, sep='\t', index=False)

        print(f"{sample_name}: {len(quantum_clonal)} truncal mutations saved")

    except Exception as e:
        print(f"Error processing {sample_name}: {e}")

