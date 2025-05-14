# Snakefile

import sys
import pandas as pd

PATIENT = config.get('patient', None)
REGION = config.get('region', None)
if PATIENT is None or REGION is None:
    print("Error: Please provide a patient name and region using --config patient=PATIENT_ID region=REGION_ID")
    sys.exit(1)

#Load accession numbers for tumor samples for the patient
info_table = pd.read_csv('../../lab/data/nih_brain_tumor_data/GBM/SraRunTable.txt', sep=',')
pattern = f"{PATIENT}.{REGION}.Tumor.DNA"
ACCESSIONS = info_table[info_table['Sample Name'].str.contains(pattern, na=False)]['Run'].tolist()

if len(ACCESSIONS) == 0:
    print(f"Error: No accessions found matching pattern: {pattern}")
    sys.exit(1)
else:
    print(f"Found {len(ACCESSIONS)} accession(s) for pattern '{pattern}': {ACCESSIONS}")

rule all:
    input:
        expand("results/{patient}.{region}/{accession}.tumor.aln.srt.bam", accession=ACCESSIONS, patient=PATIENT, region=REGION),
        expand("results/{patient}.{region}/{accession}.tumor.aln.srt.bam.bai", accession=ACCESSIONS, patient=PATIENT, region=REGION)

rule bwa_map:
    input:
        ref = 'hs38DH.fa',
        r1='../../lab/data/nih_brain_tumor_data/brmet/{accession}/{accession}_1.fastq.gz',
        r2='../../lab/data/nih_brain_tumor_data/brmet/{accession}/{accession}_2.fastq.gz'
    output:
        bam='results/{patient}.{region}/{accession}.tumor.aln.bam'
    params:
        rg=r"@RG\tID:{accession}\tSM:{patient}.{region}\tPL:ILLUMINA"
    threads: workflow.cores
    shell:
        """
        module load GCC/7.3.0-2.30  OpenMPI/3.1.1
        module load bwa-mem2/2.2.1
        module load SAMtools/1.15

        echo "Running alignment for {wildcards.accession}..."
        mkdir -p results/{wildcards.patient}.{wildcards.region}/

        # Remove old output if it exists
        if [ -f {output.bam} ]; then
        echo "Warning: Overwriting existing file {output.bam}"
        rm -f {output.bam}
        fi

        bwa-mem2 mem -t{threads} -R '{params.rg}' {input.ref} {input.r1} {input.r2}\
        | samtools view -1 -Sb - > {output.bam};
        
        # Verify output was created
        if [ ! -f {output.bam} ]; then
            echo "Error: Alignment failed - output file not created"
            exit 1
        fi
        """

rule samtools_sort:
    input:
        'results/{patient}.{region}/{accession}.tumor.aln.bam'
    output:
        'results/{patient}.{region}/{accession}.tumor.aln.srt.bam'
    threads: workflow.cores
    shell:
        """
        samtools sort -T -@ {threads} 'results/{patient}.{region}/{accession} -O BAM {input} > {output}'
        """

rule samtools_index:
    input:
        'results/{patient}.{region}/{accession}.tumor.aln.srt.bam'
    output:
        'results/{patient}.{region}/{accession}.tumor.aln.srt.bam.bai'
    threads: workflow.cores
    shell:
        "samtools -@ {threads} index {input}"