reads_ch = Channel.fromFilePairs(params.reads + '/**/*{1,2}.fastq.gz')
reads_ch2 = Channel.fromPath(params.reads + '/**/*.fastq.gz')
meta_path = Channel.fromPath(params.metadata)
adapter=file(params.adapter)
metadata=file(params.metadata)


cdna_ch=Channel.fromPath(params.cdna + '/*.gz')
cds_ch=Channel.fromPath(params.cds + '/*.gz')



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
    input: 
    file(ref) from reference
    file(gtf) from annotation

    output:



    script:
    """
    #!/bin/bash
    source /usr/local/extras/Genomics/.bashrc
    source activate star
    STAR --runThreadN 6 \
	--runMode genomeGenerate \
	--genomeDir star_index \
	--genomeFastaFiles ref \
	--sjdbGTFfile gtf \
	--sjdbOverhang 99    

    """



}


process allignment_star {
    

    input:
    tuple val(species), val(sid), file("${species}_${sid}_forward_paired.fastq.gz"), file("${species}_${sid}_reverse_paired.fastq.gz") from trimmed2


    output:
    

    script: 
    """
    #!/bin/bash
    source /usr/local/extras/Genomics/.bashrc
    source activate star
    STAR --genomeDir star_index/ \
	--runThreadN 6 \
	--readFilesIn xxxxxx \
	--outFileNamePrefix star_alligned_${sid} \
	--outSAMtype BAM SortedByCoordinate \
	--outSAMunmapped Within \
	--outSAMattributes Standard
    """

}
