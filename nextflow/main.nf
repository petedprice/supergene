reads_ch = Channel.fromFilePairs(params.reads + '/*{1,2}.fastq.gz')
gtf=file(params.gtf)
ref=file(params.ref)
ref_index=file(params.ref_index)
metadata=file(params.metadata)
adapter=file(params.adapter)
//targetsnps=params.targetsnps

ref_ch=channel
    .fromPath(metadata)
    .splitCsv()
    .map {row ->tuple(row[0],row[1])}


process trim {


    queue = "ressexcon.q"
    clusterOptions = { '-P ressexcon' }
    cpus = 4
    memory = '8 GB'
    time = '2h'
    
    tag {'trim_' + species + '_' + sid }


    //publishDir 'paired', mode: 'copy', overwrite: true, pattern: '*paired*'

    input:
    tuple val(sid), file(reads), val(species) from reads_ch.combine(ref_ch, by:0)

    output:
    tuple val(species), val(sid), file("${species}_${sid}_forward_paired.fastq.gz"), file("${species}_${sid}_reverse_paired.fastq.gz") into trimmed1
    tuple val(species), val(sid), file("${species}_${sid}_forward_paired.fastq.gz"), file("${species}_${sid}_reverse_paired.fastq.gz") into trimmed2

    script:

    """
    #!/bin/bash
    source /usr/local/extras/Genomics/.bashrc
    trimmomatic PE -phred33 $reads ${species}_${sid}_forward_paired.fastq.gz ${sid}_forward_unpaired.fastq.gz ${species}_${sid}_reverse_paired.fastq.gz ${sid}_reverse_unpaired.fastq.gz ILLUMINACLIP:$adapter:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:95
    """
}


process prep {


    output:
    tuple file('genome.ss'), file('genome.exon') into prepped
    tuple file('genome.ss'), file('genome.exon') into prepped2


    script:
    """
    #!/bin/bash
    source /usr/local/extras/Genomics/.bashrc
    hisat2_extract_splice_sites.py $gtf > genome.ss
    hisat2_extract_exons.py $gtf > genome.exon
    """

}




process index_hisat2 {

    queue = "ressexcon.q"
    cpus = 16
    memory = '192 GB'
    time = '1h'
    clusterOptions = { '-P ressexcon' }

    
    tag {'hisat index'}


    publishDir 'index', mode: 'copy', overwrite: true, pattern: '*'

    input:
    tuple file('genome.ss'), file('genome.exon') from prepped

    output:
    file('hisat_index') into hisat_indexed


    script:
    """
    #!/bin/bash
    source /usr/local/extras/Genomics/.bashrc
    ##hisat2-build --exon genome.exon --ss genome.ss $ref index -p 16
    hisat2-build $ref index -p 16
    mkdir hisat_index
    mv index*ht2 hisat_index
    """



}


process allignment_hisat2 {
    
    tag {'allign_' + '_' + sid }

    queue = "ressexcon.q"
    cpus = 4
    memory = '8 GB'
    time = '4h'
    clusterOptions = { '-P ressexcon' }
    

    publishDir 'bam', mode: 'copy', overwrite: true, pattern: '*bam'

    input:
    tuple val(species), val(sid), file("${species}_${sid}_forward_paired.fastq.gz"), file("${species}_${sid}_reverse_paired.fastq.gz") from trimmed2
    file('hisat_index') from hisat_indexed
    tuple file('genome.ss'), file('genome.exon') from prepped2

    output:
    tuple val(sid), file("${sid}.bam") into alligned1
    tuple val(sid), file("${sid}.bam") into alligned2     

    script: 
    """
    #!/bin/bash
    source /usr/local/extras/Genomics/.bashrc
    hisat2 \
	-x hisat_index/index \
	-1 ${species}_${sid}_forward_paired.fastq.gz \
	-2 ${species}_${sid}_reverse_paired.fastq.gz \
        --known-splicesite-infile genome.ss \
        --summary-file ${sid}_hisat2.summary.log \
        --threads 16 \
            | samtools view -bS -F 4 -F 8 -F 256 - > ${sid}.bam
    """

}

process RG_add {
    errorStrategy 'ignore'

    tag {'RG_' + '_' + sid }

    queue = "ressexcon.q"
    clusterOptions = { '-P ressexcon' }

    cpus = 4
    memory = '8 GB'
    time = '2h'

    input:
    tuple val(sid), file("${sid}.bam") from alligned1

    output:
    tuple val(sid), file("${sid}.RG.bam") into readgrouped1
    tuple val(sid), file("${sid}.RG.bam") into readgrouped2

    script:
    """
    #!/bin/bash
    source /usr/local/extras/Genomics/.bashrc
    source activate gatk
    picard AddOrReplaceReadGroups I=${sid}.bam O=${sid}.RG.bam SORT_ORDER=coordinate RGID=zf RGLB=zf RGPL=illumina RGSM=zf RGPU=zf CREATE_INDEX=True
    """

}



process snp_calling {
    errorStrategy 'ignore'

    tag {'haplotypecaller_' + '_' + sid }


    queue = "ressexcon.q"
    clusterOptions = { '-P ressexcon' }
    cpus = 2
    memory = '4 GB'
    time = '2h'

    input:
    tuple val(sid), file("${sid}.RG.bam") from readgrouped1
    file('ref.fasta') from ref
    file('ref.fasta.fai') from ref_index

    output:
   

    script:
    """
    #!/bin/bash
    source /usr/local/extras/Genomics/.bashrc

    source activate gatk 

    samtools faidx ref.fasta

    picard CreateSequenceDictionary R=ref.fasta O=ref.dict

    gatk --java-options "-Xmx4g" HaplotypeCaller \
        -R ref.fasta \
        -I ${sid}.RG.bam \
        -O ${sid}.vcf.gz \
       	-ERC GVCF
    """

}

process featureCounts { 
    errorStrategy 'ignore'
    tag {'featurecount_' + '_' + sid }

    input:
    tuple val(sid), file("${sid}.RG.bam") from readgrouped2

    output:


    script:
    """
    #!/bin/bash
    source /usr/local/extras/Genomics/.bashrc

    source activate subread



    """

}

