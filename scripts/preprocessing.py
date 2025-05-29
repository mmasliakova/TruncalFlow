import os
import pandas as pd
import csv
import argparse
import glob
import numpy as np

from vcf_filtering import filter_variants, remove_zeroCN_regions

def extract_sample_id(filename):
    return os.path.basename(filename).split('.')[0]


def parse_format_field(row, field_name, sample_col='TUMOR'):
    sample_col_upper = sample_col.upper()
    sample_values = row[sample_col_upper].split(':')
    format_fields = row['FORMAT'].split(':')
    field_index = format_fields.index(field_name)
    return sample_values[field_index]

def parse_vcf(vcf_file):
    with open(vcf_file, 'r') as file:
        lines = []
        for line in file:
            if line.startswith('##'):
                continue
            elif line.startswith('#'):
                header = line.strip().split('\t')
            else:
                lines.append(line.strip().split('\t'))
    return pd.DataFrame(lines, columns=header)


def preprocess_vcf(vcf, sample_id, sample_col='TUMOR', default_genotype=False):
    vcf.columns = vcf.columns.str.upper()
    sample_col_upper = sample_col.upper()
    
    vcf = vcf[vcf['FILTER'] == 'PASS']
    
    vcf = vcf[vcf['FILTER'] == 'PASS'].copy()
    vcf['POS'] = pd.to_numeric(vcf['POS'], errors='coerce')
    vcf = vcf.dropna(subset=['POS'])

    if vcf.empty:
        return None

    df = pd.DataFrame({
        'Sample': sample_id,
        'Chr': vcf['#CHROM'].str.replace('chr', '', regex=True),
        'Start': vcf['POS']
    })

    df['Depth'] = vcf.apply(lambda row: parse_format_field(row, 'DP', sample_col_upper), axis=1)
    df['Alt'] = vcf.apply(lambda row: parse_format_field(row, 'AD', sample_col_upper).split(',')[1], axis=1)

    df['Depth'] = df['Depth'].astype(int)
    df['Alt'] = df['Alt'].astype(int)
    
    if default_genotype:
        df['Genotype'] = 'AB'

    return df


def convert_major_minor_to_freec(cnv_file, exclude_cn0=True):
    cnv_df = pd.read_csv(cnv_file, sep='\t')
    if exclude_cn0:
        cn0 = remove_zeroCN_regions(cnv_file)
        cnv_df = cnv_df[~cnv_df.apply(
            lambda row: any(
                (row["Chromosome"] == cn_row["Chromosome"]) and
                (cn_row["Start"] <= row["Start"] <= cn_row["End"])
                for _, cn_row in cn0.iterrows()
            ),
            axis=1
        )]

    header = ["Chromosome", "Start", "Ratio", "MedianRatio", "CopyNumber", "BAF", "EstimatedBAF", "Genotype", "UncertaintyOfGT"]
    freec = []

    for _, row in cnv_df.iterrows():
        if row['Copy_Number'] == 0:
            continue
        ratio = row['Copy_Number'] / 2
        BAF = row['Minor_Copy_Number'] / row['Copy_Number']
        corr_BAF = abs(round(BAF, 2) - 0.5)

        genotype = (
            "0" if (row['Major_Copy_Number'] == 0 and row['Minor_Copy_Number'] == 0)
            else "A" if row['Major_Copy_Number'] == 1 and row['Minor_Copy_Number'] == 0
            else "A" * int(row['Major_Copy_Number']) + "B" * int(row['Minor_Copy_Number'])
        )

        freec.append([row['Chromosome'].replace('chr', ''), row['Start'], ratio, ratio,
                      row['Copy_Number'], corr_BAF, round(BAF, 2), genotype, 0])
    return pd.DataFrame(freec, columns=header)


def convert_battenberg_to_freec(battenberg_file):
    header = ["Chromosome", "Start", "Ratio", "MedianRatio", "CopyNumber", "BAF", "EstimatedBAF", "Genotype", "UncertaintyOfGT"]
    freec = []

    with open(battenberg_file, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            ratio = 2 ** float(row['LogR'])
            corr_baf = abs(round(float(row['BAF']), 2) - 0.5)
            copy_number = (int(row['nMaj1_A']) + int(row['nMin1_A'])) * float(row['frac1_A'])
            if row['frac2_A'] != 'NA':
                copy_number += (int(row['nMaj2_A']) + int(row['nMin2_A'])) * float(row['frac2_A'])

            copy_number = round(copy_number)

            if row['frac2_A'] != 'NA' and float(row['frac2_A']) > 0:
                genotype = "AB" if round(copy_number) >= 2 else "A"
            else:
                genotype = (
                    "0" if (int(row['nMaj1_A']) == 0 and int(row['nMin1_A']) == 0)
                    else "A" if (int(row['nMaj1_A']) == 1 and int(row['nMin1_A']) == 0)
                    else "A" * int(row['nMaj1_A']) + "B" * int(row['nMin1_A'])
                )

            uncertainty = round(100 * float(row['SDfrac_A'])) if row['SDfrac_A'] != "NA" else 0

            freec.append([row['chr'], row['startpos'], ratio, ratio, copy_number,
                          corr_baf, round(float(row['BAF']), 2), genotype, uncertainty])

    return pd.DataFrame(freec, columns=header)
    

def limit_mutations(vcf_df, max_mutations=65536, key_mutations_path=None):
    if len(vcf_df) <= max_mutations:
        return vcf_df

    print(f"[INFO] Limiting to {max_mutations} mutations while preserving chromosome distribution...")

    #chrom_counts = vcf_df['#CHROM'].value_counts(normalize=True)
    #sampled = []

    #for chrom, frac in chrom_counts.items():
        #chrom_df = vcf_df[vcf_df['#CHROM'] == chrom]
        #n_samples = int(max_mutations * frac)
        #sampled.append(chrom_df.sample(n=min(n_samples, len(chrom_df)), random_state=42))

    #final_df = pd.concat(sampled)

    # If total is slightly over/under due to rounding, resample exact amount
    #if len(final_df) > max_mutations:
        #final_df = final_df.sample(n=max_mutations, random_state=42)
    #elif len(final_df) < max_mutations:
        #print(f"[WARN] Only {len(final_df)} variants retained after chrom-balanced sampling (target was {max_mutations})")

    #return final_df
    
    if keep_mutations_path:
        keep_df = pd.read_csv(keep_mutations_path, sep=None, engine='python')
        keep_df.columns = [c.strip().upper() for c in keep_df.columns]
        keep_df.rename(columns={'CHROM': '#CHROM'}, inplace=True)

        key_cols = ['#CHROM', 'POS', 'REF', 'ALT']
        merged = pd.merge(vcf_df, keep_df[key_cols], on=key_cols, how='inner')
        print(f"[INFO] Found {len(merged)} priority mutations to keep.")
        vcf_df = vcf_df.drop(merged.index)  # remove them from original
    else:
        merged = pd.DataFrame(columns=vcf_df.columns)

    remaining_quota = max_mutations - len(merged)

    if remaining_quota <= 0:
        print("[WARN] Number of priority mutations exceeds max limit, returning only those.")
        return merged.iloc[:max_mutations]

    # Continue with chrom-balanced sampling for the rest
    chrom_counts = vcf_df['#CHROM'].value_counts(normalize=True)
    sampled = []

    for chrom, frac in chrom_counts.items():
        chrom_df = vcf_df[vcf_df['#CHROM'] == chrom]
        n_samples = int(remaining_quota * frac)
        sampled.append(chrom_df.sample(n=min(n_samples, len(chrom_df)), random_state=42))

    final_sample = pd.concat(sampled)

    if len(final_sample) > remaining_quota:
        final_sample = final_sample.sample(n=remaining_quota, random_state=42)

    final_df = pd.concat([merged, final_sample])
    return final_df


def process_sample(vcf_path, cnv_path, sample_id, cnv_format, output_dir, mutation_types=None, cn0_dict=None, functional_filter=None, filtering_mode=None, key_mutations=None):
    vcf = parse_vcf(vcf_path)

    # === FILTERING applied to raw VCF ===
    if mutation_types or functional_filter:
        if cn0_dict is None:
            cn0_dict = {}
        vcf = filter_variants(
            vcf, 
            mutation_types=mutation_types, 
            functional_filter=functional_filter,
            cn0_regions=cn0_dict.get(sample_id, None),
            filtering_mode=filtering_mode
        )
    
    vcf = limit_mutations(vcf, key_mutations_path=key_mutations)
    
    default_genotype = cnv_path is None
    snvs = preprocess_vcf(vcf, sample_id, default_genotype=default_genotype)

    # === OUTPUT ===
    sample_output_dir = os.path.join(output_dir, sample_id)
    os.makedirs(sample_output_dir, exist_ok=True)

    if snvs is not None:
        snvs.to_csv(os.path.join(sample_output_dir, f"{sample_id}_SNVlist.txt"), sep='\t', index=False)

    # === CNV PROCESSING ===
    if cnv_path:
        if cnv_format == 'battenberg':
            freec = convert_battenberg_to_freec(cnv_path)
        else:
            freec = convert_major_minor_to_freec(cnv_path, exclude_cn0=True)

        freec.to_csv(os.path.join(sample_output_dir, f"{sample_id}_freec.txt"), sep='\t', index=False)

    #if not default_genotype:
     #   if cnv_format == 'battenberg':
      #      freec = convert_battenberg_to_freec(cnv_path)
       # else:
        #    freec = convert_major_minor_to_freec(cnv_path, exclude_cn0=True)
        
        #freec.to_csv(os.path.join(sample_output_dir, f"{sample_id}_freec.txt"), sep='\t', index=False)
