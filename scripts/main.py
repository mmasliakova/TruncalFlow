import os
import glob
import argparse
import time
import subprocess
import sys

from preprocessing import process_sample, extract_sample_id
from run_clustering import submit_clustering_jobs, check_clustering_errors
from results_processing import process_quantumclone_results

def wait_and_postprocess(output_dir, vcf_input, postprocess_fn, timeout_hours=10, poll_interval=60):
    """
    Waits for clustering outputs and starts postprocessing as soon as each sample finishes.

    Args:
        output_dir (str): Base output directory.
        vcf_input (str): VCF input path (file or folder) for process_quantumclone_results().
        postprocess_fn (function): The function to call for postprocessing.
        timeout_hours (int): Max time to wait before raising TimeoutError.
        poll_interval (int): Seconds to wait between checks.
    """
    
    sample_dirs = [d for d in os.listdir(output_dir) if os.path.isdir(os.path.join(output_dir, d))]
    completed_samples = set()
    errored_samples = set()
    #expected_files = [os.path.join(output_dir, s, "clustering.csv") for s in sample_dirs]
    start_time = time.time()

    while True:
        for sample in sample_dirs:
            sample_path = os.path.join(output_dir, sample)
            cluster_file = os.path.join(sample_path, "clustering.csv")
            err_file = os.path.join(sample_path, "pbs_jobs", f"RunClustering_{sample}.err")
    
            if sample in completed_samples or sample in errored_samples:
                continue
    
            if os.path.exists(cluster_file):
                print(f"[INFO] Clustering complete for {sample}. Starting postprocessing...")
                try:
                    postprocess_fn(output_dir, vcf_input, sample)
                    completed_samples.add(sample)
                except Exception as e:
                    print(f"[ERROR] Postprocessing failed for {sample}: {e}")
                    traceback.print_exc()
                    errored_samples.add(sample)
            elif os.path.exists(err_file):
                with open(err_file) as f:
                    content = f.read()
                    if "OOM" in content or "Killed" in content or "error" in content.lower():
                        print(f"[ERROR] Clustering job failed for {sample}. Error:\n{content}")
                        errored_samples.add(sample)
    
        print(f"[INFO] Waiting for clustering... {len(completed_samples)}/{len(sample_dirs)} done.")
    
        if len(completed_samples) + len(errored_samples) == len(sample_dirs):
            break
    
        if time.time() - start_time > timeout_hours * 3600:
            print("[ERROR] Timeout exceeded while waiting for clustering jobs.")
            break
    
        time.sleep(poll_interval)

    if len(completed_samples) == 0:
        raise RuntimeError("[ERROR] No successful clustering jobs.")
    else:
        print(f"[INFO] Finished postprocessing {len(completed_samples)} sample(s).")


def main():
    parser = argparse.ArgumentParser(description="QuantumClone pipeline runner")
    parser.add_argument('--vcf', required=True, help='VCF file or folder')
    parser.add_argument('--cnv', required=False, help='Optional CNV file or folder')
    parser.add_argument('--cnv_format', required=False, choices=['battenberg', 'major_minor_format'])
    parser.add_argument('--output_dir', required=True)
    #parser.add_argument('--vcf_root_dir', required=False, help='Directory where original VCFs are stored (for result annotation)')
    parser.add_argument('--mutation_types', nargs='+', choices=['snv', 'indel'], help='Filter for mutation types: snv, indel')
    parser.add_argument('--functional_filter', nargs='+', choices=['protein_coding', 'missense', 'nonsense', 'transcript_ablation',
    'splice_acceptor_variant', 'deletion', 'insertion', 'splice_donor_variant', 'stop_gained', 'frameshift_variant',
    'stop_lost', 'start_lost', 'transcript_amplification', 'feature_elongation', 'feature_truncation'], help='Filter for functional consequence')
    parser.add_argument('--filtering_mode', choices=['both', 'single'], default='single', required=False,
    help="Combine mutation_types and functional_filter using 'both' (both must match) or 'single' (either can match). Default is 'single'.")
    parser.add_argument('--key_mutations', required=False, help='Path to a TXT file listing mutations to retain (CHROM, POS, REF, ALT)')
    parser.add_argument('--time', required=False, default='10:00:00',
    help='Walltime required to complete clustering. Increase if a clustering job failed due to runtime issue (format hours:minutes:seconds), default 10:00:00.')


    args = parser.parse_args()
    os.makedirs(args.output_dir, exist_ok=True)
    

    # === STEP 1: Preprocessing ===
    if os.path.isfile(args.vcf):
    
        sample_id = extract_sample_id(args.vcf)
        cnv_file = args.cnv if args.cnv and os.path.isfile(args.cnv) else None
        
        if cnv_file:
            cnv_sample_id = extract_sample_id(cnv_file)
            if cnv_sample_id != sample_id:
                print(f"[ERROR] Sample ID mismatch: VCF is '{sample_id}' but CNV is '{cnv_sample_id}'")
                return
        
        print(f"[INFO] Preprocessing single sample: {sample_id}")
        process_sample(args.vcf, cnv_file, sample_id, args.cnv_format, args.output_dir, mutation_types=args.mutation_types, functional_filter=args.functional_filter, filtering_mode=args.filtering_mode, key_mutations=args.key_mutations)

    elif os.path.isdir(args.vcf):
        vcf_files = glob.glob(os.path.join(args.vcf, '*.vcf'))
        cnv_dict = {}
        
        if args.cnv and os.path.isdir(args.cnv):
            cnv_files = glob.glob(os.path.join(args.cnv, '*'))
            cnv_dict = {extract_sample_id(f): f for f in cnv_files}
    
            vcf_sample_ids = {extract_sample_id(f) for f in vcf_files}
            cnv_sample_ids = set(cnv_dict.keys())
            unmatched_cnvs = cnv_sample_ids - vcf_sample_ids
    
            if unmatched_cnvs:
                print(f"[ERROR] CNV files found with no matching VCFs: {', '.join(unmatched_cnvs)}")
                return
        
        if not cnv_dict:
            print(f"[WARN] No CNV files found in the specified CNV folder: {args.cnv}. Continuing without CNV data.")


        for vcf_file in vcf_files:
            sample_id = extract_sample_id(vcf_file)
            cnv_file = cnv_dict.get(sample_id)
            print(f"[INFO] Preprocessing sample: {sample_id}")
            process_sample(vcf_file, cnv_file, sample_id, args.cnv_format, args.output_dir, mutation_types=args.mutation_types, functional_filter=args.functional_filter, filtering_mode=args.filtering_mode, key_mutations=args.key_mutations)

    else:
        print("[ERROR] Invalid input. Please provide VCF file/folder and optionally CNV file/folder.")
        return
    
    job_ids = []
    # === STEP 2: Submit Clustering Jobs ===
    try:
        job_ids = submit_clustering_jobs(args.output_dir, args.time)
        
        # === STEP 3: Wait for Clustering Completion ===
        wait_and_postprocess(
            args.output_dir,
            args.vcf,
            postprocess_fn=process_quantumclone_results,
            timeout_hours=10
        )
    
        # === STEP 4: Process Results ===
        #sample_dirs = [d for d in os.listdir(args.output_dir) if os.path.isdir(os.path.join(args.output_dir, d))]
        #successful_samples = [s for s in sample_dirs if s not in failed_samples]
        
        #if successful_samples:
         #   for sample in successful_samples:
          #      print(f"[INFO] Postprocessing results for: {sample}")
           #     process_quantumclone_results(args.output_dir, args.vcf, sample)
        #else:
         #   print("[ERROR] All clustering jobs failed. Exiting.")
    except KeyboardInterrupt:
        print("\n[WARN] KeyboardInterrupt received during job submission. Canceling jobs...")
        for job_id in job_ids:
            try:
                subprocess.run(f"qdel {job_id}", shell=True)
                print(f"[INFO] Canceled job {job_id}")
            except Exception as e:
                print(f"[ERROR] Failed to cancel job {job_id}: {e}")
        sys.exit(1)


if __name__ == '__main__':
    main()