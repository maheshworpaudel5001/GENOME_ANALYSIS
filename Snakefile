# Snakefile

import sys
import pandas as pd

PATIENT = config.get('patient', None)
REGION = config.get('region', None)

if PATIENT is None or REGION is None:
    print("Error: Please provide a patient name and region using --config patient=PATIENT_ID region=REGION_ID")
    sys.exit(1)

# Determine patient class
PATIENT_CLASS = 'brmet' if 'BrMET' in PATIENT else 'gbm'

# Load accession numbers for tumor samples for the patient
sra_table = 'SraRunTable.txt'
info_table = pd.read_csv(sra_table, sep=',')
tumor_pattern = f"{PATIENT}.{REGION}.Tumor.DNA"
normal_pattern = f"{PATIENT}.Normal.DNA"
TUMOR_ACCESSIONS = info_table[info_table['Sample Name'].str.contains(tumor_pattern, na=False)]['Run'].tolist()
NORMAL_ACCESSIONS = info_table[info_table['Sample Name'] == normal_pattern]['Run'].tolist()

if len(TUMOR_ACCESSIONS) == 0:
    print(f"Error: No accessions found matching pattern: {tumor_pattern}")
    sys.exit(1)
else:
    print(f"Found {TUMOR_ACCESSIONS} accessions for pattern {tumor_pattern}")

if len(NORMAL_ACCESSIONS) == 0:
    print(f"Error: No normal sample found for pattern '{normal_pattern}'")
    sys.exit(1)
else:
    print(f"Found {NORMAL_ACCESSIONS} accession(s) for pattern {normal_pattern}")

TUMOR_ACCESSION = TUMOR_ACCESSIONS[0]
NORMAL_ACCESSION = NORMAL_ACCESSIONS[0]

rule all:
    input:
        f"results/{PATIENT}.{REGION}/{TUMOR_ACCESSION}.tumor.aln.srt.bam",
        f"results/{PATIENT}.{REGION}/{TUMOR_ACCESSION}.tumor.aln.srt.bam.bai",
        f"results/{PATIENT}.normal.aln.srt.bam",
        f"results/{PATIENT}.normal.aln.srt.bam.bai",
        f"results/{PATIENT}.{REGION}/strelka/results/variants/somatic.snvs.vcf.gz",
        f"results/{PATIENT}.{REGION}/strelka/results/variants/somatic.indels.vcf.gz"

rule map_tumor:
    input:
        ref = 'hs38DH.fa'
    params:
        r1 = f"../../lab/data/nih_brain_tumor_data/{PATIENT_CLASS}/{TUMOR_ACCESSION}/{TUMOR_ACCESSION}_1.fastq.gz",
        r2 = f"../../lab/data/nih_brain_tumor_data/{PATIENT_CLASS}/{TUMOR_ACCESSION}/{TUMOR_ACCESSION}_2.fastq.gz",
        rg = f"@RG\tID:{TUMOR_ACCESSION}\tSM:{PATIENT}.{REGION}\tPL:ILLUMINA"
    output:
        bam = f'results/{PATIENT}.{REGION}/{TUMOR_ACCESSION}.tumor.aln.srt.bam',
        bai = f'results/{PATIENT}.{REGION}/{TUMOR_ACCESSION}.tumor.aln.srt.bam.bai'
    threads: workflow.cores
    shell:
        """
        module load bwa-mem2/2.2.1
        module load SAMtools/1.15

        bwa-mem2 mem -t {threads} -R '{params.rg}' {input.ref} {params.r1} {params.r2} \
        | samtools view -@ {threads} -Sb - \
        | samtools sort -@ {threads} -o {output.bam} -
        samtools index -@ {threads} {output.bam}
        """

rule map_normal:
    input:
        ref = 'hs38DH.fa'
    params:
        r1 = f"../../lab/data/nih_brain_tumor_data/{PATIENT_CLASS}/{NORMAL_ACCESSION}/{NORMAL_ACCESSION}_1.fastq.gz",
        r2 = f"../../lab/data/nih_brain_tumor_data/{PATIENT_CLASS}/{NORMAL_ACCESSION}/{NORMAL_ACCESSION}_2.fastq.gz",
        rg = f"@RG\tID:{NORMAL_ACCESSION}\tSM:{PATIENT}.normal\tPL:ILLUMINA"
    output:
        bam = f'results/{PATIENT}.normal.aln.srt.bam',
        bai = f'results/{PATIENT}.normal.aln.srt.bam.bai'
    threads: workflow.cores
    shell:
        """
        module load bwa-mem2/2.2.1
        module load SAMtools/1.15

        bwa-mem2 mem -t {threads} -R '{params.rg}' {input.ref} {params.r1} {params.r2} \
        | samtools view -@ {threads} -Sb - \
        | samtools sort -@ {threads} -o {output.bam} -
        samtools index -@ {threads} {output.bam}
        """

rule somatic_variant_call_strelka:
    input:
        ref = 'hs38DH.fa'
    output:
        snvs = f"results/{PATIENT}.{REGION}/strelka/results/variants/somatic.snvs.vcf.gz",
        indels = f"results/{PATIENT}.{REGION}/strelka/results/variants/somatic.indels.vcf.gz"
    params:
        runDir = f'results/{PATIENT}.{REGION}/strelka',
        normal = f'results/{PATIENT}.normal.aln.srt.bam',
        tumor = f'results/{PATIENT}.{REGION}/{TUMOR_ACCESSION}.tumor.aln.srt.bam'
    threads: workflow.cores
    shell:
        """
        module load GCCcore/10.2.0
        module load Python/2.7.18
        
        #configure strelka
        ../strelka-2.9.2.centos6_x86_64/bin/configureStrelkaSomaticWorkflow.py \
        --normalBam {params.normal} \
        --tumorBam {params.tumor} \
        --referenceFasta {input.ref} \
        --runDir {params.runDir}

        #run strelka
        {params.runDir}/runWorkflow.py -m local -j {threads}
        """