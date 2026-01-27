####bin/referenceGenerate.py
import os
import subprocess
import sys
import glob
import shlex

def generateReference(gtf,fasta,protein=False,thread=1):
    print(f"Starting generateReference: gtf={gtf}, fasta={fasta}, thread={thread}")
    
    scriptPath=sys.path[0]+'/referenceGenerate/'
    gtf=os.path.abspath(gtf)
    fasta=os.path.abspath(fasta)
    
    if protein:
        # ISSUE 1: Path with special chars in Rscript command
        process_cmd=f'Rscript {shlex.quote(scriptPath+"0.gtf_process.R")} -g {shlex.quote(gtf)}'
        print(f"Running: {process_cmd}")
        subprocess.run(process_cmd,shell=True,check=True)
        gtf=gtf.replace('.gtf','_proteinCoding.gtf')
    
    outfa=gtf.replace('.gtf','.fa')
    
    # ISSUE 2: Path with special chars in gffread command  
    process_cmd=f'gffread -w {shlex.quote(outfa)} -g {shlex.quote(fasta)} {shlex.quote(gtf)}'
    print(f"Running: {process_cmd}")
    subprocess.run(process_cmd,shell=True,check=True)
    
    outpath=os.path.dirname(outfa)
    print(f"Output path: {outpath}")
    
    if thread == 1:
        # ISSUE 3: Path with special chars in single-thread ORF generation
        orf_script = scriptPath+'1.ORF_generate.py'
        process_cmd=f'python {shlex.quote(orf_script)} -i {shlex.quote(outfa)} -o {shlex.quote(outpath)}'
        print(f"Running: {process_cmd}")
        subprocess.run(process_cmd,shell=True,check=True)
        
    elif thread > 1:
        # Multi-thread processing
        split_transcript_dir = os.path.join(outpath, 'split', 'trancript')
        split_orf_dir = os.path.join(outpath, 'split', 'ORF')
        
        os.makedirs(split_transcript_dir, exist_ok=True)
        os.makedirs(split_orf_dir, exist_ok=True)
        
        # ISSUE 4: Path with special chars in R split script
        split_script = scriptPath+'2.split_fasta.R'
        process_cmd=f'Rscript {shlex.quote(split_script)} -i {shlex.quote(outfa)} -s {str(thread)} -o {shlex.quote(split_transcript_dir)}'
        print(f"Running: {process_cmd}")
        subprocess.run(process_cmd,shell=True,check=True)
        
        # ISSUE 5: Path with special chars in batch script preparation
        original_script = scriptPath+'ORF_generate.sh'
        temp_script = scriptPath+'tmp_ORF_generate.sh'
        
        if not os.path.exists(original_script):
            print(f"ERROR: Batch script not found: {original_script}")
            return
            
        # Read and modify the batch script
        with open(original_script, 'r') as f:
            script_content = f.read()
        
        # Replace placeholders - be careful with paths containing special characters
        script_content = script_content.replace('1.ORF_generate.py', shlex.quote(scriptPath+'1.ORF_generate.py'))
        
        fa_split_base = os.path.join(split_transcript_dir, os.path.basename(outfa).replace('.fa',''))
        script_content = script_content.replace('transcript', shlex.quote(fa_split_base))
        script_content = script_content.replace('path', shlex.quote(split_orf_dir+'/'))
        script_content = script_content.replace('number', str(thread))
        
        with open(temp_script, 'w') as f:
            f.write(script_content)
        
        # ISSUE 6: Path with special chars in batch script execution
        process_cmd=f'bash {shlex.quote(temp_script)}'
        print(f"Running: {process_cmd}")
        subprocess.run(process_cmd,shell=True,check=True)
        
        # ISSUE 7: The main problem - cat command with special characters
        # Replace shell cat command with Python file merging
        process_cmd=f'cat {shlex.quote(outpath)}/split/ORF/* > {shlex.quote(outpath)}/ORF.fa'
        subprocess.run(process_cmd,shell=True,check=True)
    
    # ISSUE 8: Path with special chars in redundancy removal script
    orf_file = os.path.join(outpath, 'ORF.fa')
    redundancy_script = scriptPath+'4.remove_ORF_redanduncy.py'
    
    if os.path.exists(orf_file) and os.path.exists(redundancy_script):
        if os.path.getsize(orf_file) > 0:
            process_cmd=f'python {shlex.quote(redundancy_script)} -i {shlex.quote(orf_file)} -o {shlex.quote(outpath)}'
            print(f"Running: {process_cmd}")
            subprocess.run(process_cmd,shell=True,check=True)
        else:
            print("ORF file is empty, creating empty ORF_no_redundancy.fa")
            with open(os.path.join(outpath, 'ORF_no_redundancy.fa'), 'w') as f:
                pass
    
    # ISSUE 9: Path with special chars in cleanup
    split_dir = os.path.join(outpath, 'split')
    if os.path.exists(split_dir):
        # Use Python's shutil instead of shell rm command
        import shutil
        shutil.rmtree(split_dir)
        print("Cleaned up split directory")


def removeHomo(ORF, db, refFa, thread=1):
    import shutil
    scriptPath = os.path.join(sys.path[0], 'referenceGenerate')
    ORF = os.path.abspath(ORF)
    orfPath = os.path.dirname(ORF)

    # Prepare BLAST database if not provided
    if not db:
        db_dir = os.path.join(orfPath, 'db')
        os.makedirs(db_dir, exist_ok=True)
        db_path = os.path.join(db_dir, 'ref.db')
        process_cmd = f'makeblastdb -in {shlex.quote(refFa)} -dbtype prot -parse_seqids -out {shlex.quote(db_path)}'
        print(f"Running: {process_cmd}")
        subprocess.run(process_cmd, shell=True, check=True)
        db = db_path
    else:
        db = os.path.abspath(db)

    if thread == 1:
        # Single-thread BLAST
        blast_out = ORF.replace('.fa', '.blast')
        process_cmd = (
            f'blastp -query {shlex.quote(ORF)} '
            f'-out {shlex.quote(blast_out)} '
            f'-db {shlex.quote(db)} -outfmt 6 -evalue 1e-5 -num_threads 1'
        )
        print(f"Running: {process_cmd}")
        subprocess.run(process_cmd, shell=True, check=True)

    elif thread > 1:
        # Multi-thread BLAST
        split_orf_dir = os.path.join(orfPath, 'split', 'ORF')
        split_blast_dir = os.path.join(orfPath, 'split', 'blast')
        os.makedirs(split_orf_dir, exist_ok=True)
        os.makedirs(split_blast_dir, exist_ok=True)

        # Split ORF fasta
        split_script = os.path.join(scriptPath, '2.split_fasta.R')
        process_cmd = (
            f'Rscript {shlex.quote(split_script)} '
            f'-i {shlex.quote(ORF)} -s {str(thread)} '
            f'-o {shlex.quote(split_orf_dir)}'
        )
        print(f"Running: {process_cmd}")
        subprocess.run(process_cmd, shell=True, check=True)

        # Prepare batch BLAST script from template
        original_script = os.path.join(scriptPath, 'removeORFhomo.sh')
        temp_script = os.path.join(scriptPath, 'tmp_removeORFhomo.sh')

        if not os.path.exists(original_script):
            print(f"ERROR: BLAST batch script not found: {original_script}")
            return

        with open(original_script, 'r') as f:
            script_content = f.read()

        # Expect placeholders in the original script: {orf}, {path}, {db}, {n}
        fa_split_base = os.path.join(split_orf_dir, os.path.basename(ORF).replace('.fa', ''))
        blast_output_base = os.path.join(split_blast_dir, os.path.basename(ORF).replace('.fa', ''))

        script_content = script_content.format(
            orf=shlex.quote(fa_split_base),
            path=shlex.quote(blast_output_base),
            db=shlex.quote(db),
            n=str(thread)
        )

        with open(temp_script, 'w') as f:
            f.write(script_content)

        # Run batch BLAST script
        process_cmd = f'bash {shlex.quote(temp_script)}'
        print(f"Running: {process_cmd}")
        subprocess.run(process_cmd, shell=True, check=True)

        # Combine BLAST outputs
        process_cmd=f"cat {shlex.quote(split_blast_dir)}/* > {shlex.quote(orfPath)}/{os.path.basename(ORF).replace('.fa', '.blast')}"
        subprocess.run(process_cmd,shell=True,check=True)

    # Remove homologous ORFs
    blast_file = ORF.replace('.fa', '.blast')
    r_script = os.path.join(scriptPath, '3.remove_homo_ORF.R')
    process_cmd = (
        f'Rscript {shlex.quote(r_script)} '
        f'-f {shlex.quote(ORF)} '
        f'-b {shlex.quote(blast_file)} '
        f'-n ORF_unmatch_ref.fa'
    )
    print(f"Running: {process_cmd}")
    subprocess.run(process_cmd, shell=True, check=True)

    # Cleanup
    split_dir = os.path.join(orfPath, 'split')
    if os.path.exists(split_dir):
        shutil.rmtree(split_dir)
        print("Cleaned up split directory")
