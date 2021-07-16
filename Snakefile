
READS_DIR = '/vortexfs1/omics/alexander/data/TARA/PRJEB4352-snakmake-output/trimmed/PRJEB4352'
READS, = glob_wildcards(os.path.join(READS_DIR, '{samples}_1.trimmed.fastq.gz'))
MAG_DIR = 'genomes'
OUTPUTDIR = 'output'

rule all:
    input: expand(os.path.join(OUTPUTDIR ,"{sample}.coverm.abundance.tab"), sample=READS), expand('/vortexfs1/scratch/halexander/coverm-tmp-202012/{sample}.mapped.bam', sample=READS)

rule coverm_genome:
    input:
        r1 = os.path.join(READS_DIR, "{sample}" + "_1.trimmed.fastq.gz"),
        r2 = os.path.join(READS_DIR, "{sample}" + "_2.trimmed.fastq.gz"),
    params:
        genome_dir = MAG_DIR,  
        tmpdir = '/vortexfs1/scratch/halexander/coverm-tmp-202012/{sample}'

    output:
        os.path.join(OUTPUTDIR ,"{sample}.coverm.abundance.tab")
    conda:
        "coverm.yaml"
    shell:
        """
        mkdir -p {params.tmpdir}
        export TMPDIR={params.tmpdir} 
        coverm genome --coupled {input.r1} {input.r2} --genome-fasta-directory {params.genome_dir} --genome-fasta-extension fna --min-read-percent-identity 0.95 --min-read-aligned-percent 0.75 --output-format dense     --min-covered-fraction 0     --contig-end-exclusion 75     --trim-min 0.05     --trim-max 0.95 --proper-pairs-only --methods count length covered_bases covered_fraction reads_per_base mean variance trimmed_mean rpkm relative_abundance --threads 16 --bam-file-cache-directory {params.tmpdir} > {output} 
        """

rule mapped_reads: 
    input: os.path.join(OUTPUTDIR ,"{sample}.coverm.abundance.tab")
    output: '/vortexfs1/scratch/halexander/coverm-tmp-202012/{sample}.mapped.bam'
    params:  folder = '/vortexfs1/scratch/halexander/coverm-tmp-202012/{sample}', bam=  '/vortexfs1/scratch/halexander/coverm-tmp-202012/{sample}/coverm-genome.{sample}_1.trimmed.fastq.gz.bam'

    conda: 
        "coverm.yaml"
    shell:
        """
        samtools view -b -F 4 {params.bam} > {output}
        rm -rf {params.folder}
        """


