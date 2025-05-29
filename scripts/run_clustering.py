import os
import glob
import re
import sys
import subprocess

def submit_clustering_jobs(output_dir, time):
    """
    Submits QuantumClone clustering jobs for each sample in output_dir.
    
    Args:
        output_dir (str): Path to directory containing sample subdirectories.
        time (str): Walltime required to complete clustering.
    """
    output_dir = os.path.abspath(output_dir)
    walltime = time

    script_dir = os.path.dirname(os.path.abspath(__file__))
    clustering_script = os.path.join(script_dir, "clustering.R")

    sample_folders = [d for d in glob.glob(os.path.join(output_dir, "*")) if os.path.isdir(d)]
    submitted_job_ids = []

    for sample in sample_folders:
        sample_name = os.path.basename(sample)
        
        jobs_dir = os.path.join(sample, "pbs_jobs")
        os.makedirs(jobs_dir, exist_ok=True)

        snv_file = glob.glob(os.path.join(sample, "*_SNVlist.txt"))
        if not snv_file:
            print(f"Skipping {sample_name} - No SNVlist file found.")
            continue

        expected_freec_file = glob.glob(os.path.join(sample, "*_freec.txt"))
        if not expected_freec_file:
            print(f"Note: CNV file not found for {sample_name}, proceeding without FREEC.")
            # Still submit, R will handle it
            expected_freec_file = None

        job_script_path = os.path.join(jobs_dir, f"job_{sample_name}.sh")
        
        script_dir = os.path.dirname(os.path.abspath(__file__))
        r_script = os.path.join(script_dir, "clustering.R")
        sif_path = os.path.join(script_dir, "thesis_2025_marina_latest.sif")

        with open(job_script_path, "w") as fh:
            fh.write(f"""#!/bin/bash
#PBS -N QC_{sample_name}
#PBS -l nodes=1:ppn=12
#PBS -l walltime={walltime}
#PBS -l mem=200gb
#PBS -o {jobs_dir}/RunClustering_{sample_name}.out
#PBS -e {jobs_dir}/RunClustering_{sample_name}.err
#PBS -d {sample}


apptainer exec -e --no-home --env R_LIBS_USER=/usr/local/lib/R/library {sif_path} Rscript {r_script} "{output_dir}" "{sample_name}"
""")

        submit_cmd = f"qsub {job_script_path}"
        result = subprocess.run(submit_cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        if result.returncode != 0:
            print(f"[ERROR] Failed to submit clustering job for {sample_name}")
            print(result.stderr)
            sys.exit(1)
        else:
            job_id = result.stdout.strip().split('.')[0]
            print(f"[INFO] Submitted clustering job for {sample_name}")
            submitted_job_ids.append(job_id)
            
    return submitted_job_ids
            

def check_clustering_errors(output_dir):
    sample_dirs = [d for d in os.listdir(output_dir) if os.path.isdir(os.path.join(output_dir, d))]
    failed_samples = []

    for sample in sample_dirs:
        err_file = os.path.join(output_dir, sample, "pbs_jobs", f"RunClustering_{sample}.err")
        if os.path.exists(err_file):
            with open(err_file, 'r') as f:
                content = f.read().lower()
                # Match keywords that indicate serious failure
                if re.search(r"\b(error|oom|killed)\b", content):
                    print(f"[ERROR] Clustering failed for sample {sample}. Check: {err_file}")
                    failed_samples.append(sample)

    return failed_samples


