
// Downloading annotation file to associate the location to genes 
process downloadAnnotation {
	input : 
	val link  // Link to NCBI

	output :
	file "annotation.gff" 

	script :  // WARNING this function may need to be changed if it is used on linux (wget might not exist on some distros)
	"""
	curl -o annotation.gff "$link"
	"""

}

// Downloading the reference genome to create the index files
process downloadRefGenome {
	input : 
	val link // Link to NCBI

	output : 
	file "reference.fasta"  // I chose the name but it can be changed

	script : 
	""" 
	curl -o reference.fasta "$link"
	"""
}

process creatingGenomeIndex { // Creating the genome index that is required to run the mapping command
	input : 
	file ref_genome // fasta file

	output : 
        tuple val(params.RefName), file("*.ebwt")

	script:
	"""
	bowtie-build $ref_genome $params.RefName
	"""
}

process downloadFastq {    // Downloading fastq files from the NCBI database 

	input :
	val sraid // One SRAID

	output : 
	file "*.fastq.gz"

	shell:
	"""
	fasterq-dump !{sraid}
	gzip !{sraid}.fastq
	"""
}

// Trimming the files
process trimmingFastQ {

	input :
	tuple val(sraid),file(sample) 
	each error
	each adapter
	each quality
	each overlap_min

	output : 
	file("*.fastq.gz")
	tuple val(error), val(adapter), val(quality), val(overlap_min)

	shell:
	"""
	cutadapt -e !{error} -q !{quality} -O !{overlap_min} -a !{params.adapter_seq} -o !{sraid}_trimmed.fastq.gz !{sample}
	"""

}

// Mapping the sequences
process MappingFQ_samtools{
        input :
        tuple val(sraid),file(sample),val(RefName),file(index_files)
	tuple val(error), val(adapter), val(quality), val(overlap_min)


        output :
        file "*.bam"
	tuple val(error), val(adapter), val(quality), val(overlap_min)

        shell:
        """
        bowtie -S !{RefName} <(gunzip -c !{sample}) | samtools sort > !{sraid}.bam
        """

}


// Counting 
process FeatureCount{

        input :
        val list_bam    // All the bam files 
        file annotation // Annotation file
	tuple val(error), val(adapter), val(quality), val(overlap_min)

        output :
        file "counts.txt" // Count files 
	tuple val(error), val(adapter), val(quality), val(overlap_min)

        shell:
        """
        featureCounts -t gene -g ID -F GTF -s 1 -a !{annotation} -o counts.txt !{list_bam.join(' ')}
        """
}


// Downloading the AureoWiki file
process downloadAureoWiki {

	input :
	val data // Data for the php request
	val link_php //Link to the php script

	output :
	file "NCTC8325.csv"

	script :
	"""
	curl -X POST $link_php --data $data -o NCTC8325.csv
	"""

}


// Processing
process processR{
        input :
        file count // Count file 
        file script_r // R script 
	file NCTC8325 
	tuple val(error), val(adapter), val(quality), val(overlap_min)

        output :
        tuple file("MA_plot_complet.pdf"), file("MA_plot_enhanced.pdf"), val("${adapter}_${quality}_${overlap_min}_${error}") 

        shell:
        """
        R -f !{script_r} 
        """
}


