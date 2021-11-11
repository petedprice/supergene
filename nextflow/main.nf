reads_ch = Channel.fromFilePairs(params.reads + '/*{1,2}.fastq.gz')
gtf=file(params.gtf)
ref=file(params.ref)
metadata=file(params.metadata)
adapter=file(params.adapter)
//targetsnps=params.targetsnps

ref_ch=channel
    .fromPath(metadata)
    .splitCsv()
    .map {row ->tuple(row[0],row[1])}


process trim {

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

process index_star {

    cpus = 8
    memory = '40 GB'
    time = '2h'

    tag {'star index'}


    publishDir 'index', mode: 'copy', overwrite: true, pattern: '*'

    output:
    file('star_index') into star_indexed


    script:
    """
    #!/bin/bash
    source /usr/local/extras/Genomics/.bashrc
    mkdir star_index
    source activate star
    STAR --runThreadN 32 \
	--runMode genomeGenerate \
	--genomeDir star_index \
	--genomeFastaFiles $ref \
	--sjdbGTFfile $gtf \
	--sjdbOverhang 99 \
	--genomeSAindexNbases 13

    """



}


process allignment_star {
    

    cpus = 4
    memory = '16 GB'
    time = '2h'

    input:
    tuple val(species), val(sid), file("${species}_${sid}_forward_paired.fastq.gz"), file("${species}_${sid}_reverse_paired.fastq.gz") from trimmed2
    file('star_index') from star_indexed

    output:
    tuple val(sid), file("star_alligned_${sid}") into alligned     

    script: 
    """
    #!/bin/bash
    source /usr/local/extras/Genomics/.bashrc
    source activate star
    STAR --genomeDir star_index/ \
	--runThreadN 8 \
	--readFilesIn ${species}_${sid}_forward_paired.fastq.gz ${species}_${sid}_reverse_paired.fastq.gz \
	--outFileNamePrefix star_alligned_${sid} \
	--outSAMtype BAM SortedByCoordinate \
	--outSAMunmapped Within \
	--outSAMattributes Standard \
	--quantMode GeneCounts
    """

}

process snp_calling {


    input:
    tuple val(sid), file("star_alligned_${sid}") from alligned


    output:


    script:
    """
    source /usr/local/extras/Genomics/.bashrc
    gatk --java-options "-Xmx4g" HaplotypeCaller \
	-R $ref \
	-I input.bam \
	-O ${sid}.vcf.gz \
	-ERC GVCF
    """

}

/*

process allele_quant {


    input:
    tuple val(sid), file("star_alligned_${sid}") from alligned

    output:


    script:
    """
    source /usr/local/extras/Genomics/.bashrc
    gatk ASEReadCounter \ 
	-R $ref \ 
	-I ${sid}.bam \ 
	-V targetsites.vcf.gz \ 
	-O ${sid}_allele_counts.table
    """

}


*/
