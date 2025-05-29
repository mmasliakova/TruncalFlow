import pandas as pd

def has_relevant_consequence(info_field, functional_filter):
    csq_prefix = [field for field in info_field.split(';') if field.startswith('CSQ=')]
    if not csq_prefix:
        return False

    csq_entries = csq_prefix[0].replace('CSQ=', '').split(',')
    for entry in csq_entries:
        fields = entry.split('|')
        
        if len(fields) < 8:  # Ensure there are enough fields to check (0 to 7)
            continue

        for field in fields[:8]:  # Checking fields from 0 to 7 (inclusive)
            if any(f in functional_filter for f in field.split('&')):  # Check each part of the field, as some fields may have multiple values
                return True
    return False



def filter_variants(vcf_df, mutation_types=None, cn0_regions=None, functional_filter=None, filtering_mode=None):
    if vcf_df.empty:
        return vcf_df
        
        
    info_fields = vcf_df['INFO'].dropna().astype(str)

    # Check for mutation type filtering capability
    if mutation_types:
        has_mutation_info = (
            info_fields.str.contains('VT=', case=False).any() or
            info_fields.str.contains('CSQ=|ANN=', case=False).any()
        )
        if not has_mutation_info:
            raise ValueError(
                "VCF file lacks information required for mutation type filtering."
            )
            
    if functional_filter:
        has_functional_info = info_fields.str.contains('CSQ=|ANN=', case=False).any()
        if not has_functional_info:
            raise ValueError(
                "VCF file lacks information required for functional consequence filtering."
            )
            
    if mutation_types and functional_filter:
        mode = filtering_mode or 'single'
    else:
        # If only one filter is used, use that one only
        mode = 'single'  # but will be ignored effectively
        filtering_mode = None

    filtered_rows = []

    for idx, row in vcf_df.iterrows():
        info = row['INFO']
        chrom = str(row['#CHROM']).replace('chr', '')  # Normalize chromosome names
        pos = int(row['POS'])

        # --- CN=0 filtering ---
        if cn0_regions is not None:
            overlaps_cn0 = cn0_regions[
                (cn0_regions['Chromosome'].astype(str).str.replace('chr', '') == chrom) &
                (cn0_regions['Start'] <= pos) &
                (cn0_regions['End'] >= pos)
            ]
            if not overlaps_cn0.empty:
                continue  # skip if in CN=0 region
        
        mutation_passes = True
        # --- Mutation type filtering ---
        if mutation_types:
            info_lower = info.lower()
            # Skip if 'snv' is required but not present
            if 'snv' in mutation_types and not any(term in info_lower for term in ['snv', 'snp']):
                mutation_passes = False
            
            # Skip if 'indel' is required but neither 'insertion' nor 'deletion' is present
            if 'indel' in mutation_types and not any(term in info_lower for term in ['insertion', 'deletion']):
                mutation_passes = False
        
        consequence_passes = True
        # --- Functional consequence filtering ---
        if functional_filter:
            consequence_found = False

            if 'CSQ=' in info:
                if has_relevant_consequence(info, functional_filter):
                    consequence_found = True

            elif 'ANN=' in info:
                ann_data = [entry for entry in info.split(';') if entry.startswith('ANN=')]
                if ann_data:
                    ann_entries = ann_data[0].split('=')[1].split(',')
                    for ann_entry in ann_entries:
                        fields = ann_entry.split('|')
                        if len(fields) > 1:
                            consequences = fields[1].split('&')
                            if any(c in functional_filter for c in consequences):
                                consequence_found = True
                                break

            if not consequence_found:
                consequence_passes = False

        if filtering_mode == 'both':
            if mutation_types and functional_filter:
                if mutation_passes and consequence_passes:
                    filtered_rows.append(row)
            elif mutation_types and mutation_passes:
                filtered_rows.append(row)
            elif functional_filter and consequence_passes:
                filtered_rows.append(row)
        else:
            if (mutation_types and mutation_passes) or (functional_filter and consequence_passes):
                filtered_rows.append(row)

    return pd.DataFrame(filtered_rows)



def remove_zeroCN_regions(cnv_file):
    cnv_df = pd.read_csv(cnv_file, sep='\t')
    cn0 = cnv_df[cnv_df["Copy_Number"] == 0]
    return cn0
