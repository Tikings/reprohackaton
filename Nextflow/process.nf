SRAIDs = [
"SRR10379721",
"SRR10379722",
"SRR10379723",
"SRR10379724",
"SRR10379725",
"SRR10379726",
]

linkAnnotation = "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?db=nuccore&report=gff3&id=CP000253.1"
linkRefGenome = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=CP000253.1&rettype=fasta"

params.thread_nb = 2 // To change later
params.RefName = "INDEX_S_aureus" 			// Name of the reference built by the bowtie-build command 
params.adapter_seq = "AGATCGGAAGAGC"


// Downloading annotation file to associate the location to genes 
process downloadAnnotation {
	input : 
	val link  // Link to NCBI

	output :
	file "annotation.gff" // I chose the name but it can be changed

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

	script :  // WARNING this function may need to be changed if it is used on linux (wget might not exist on some distros)
	""" 
	curl -o reference.fasta "$link"
	"""
}

process creatingGenomeIndex { // Creating the genome index that is required to run the mapping command
	input : 
	file ref_genome // fasta file

	output : 
	file "*.ebwt" 

	script:
	"""
	bowtie-build $ref_genome $params.RefName
	"""
}

process downloadFastq {    // Downloading fastq files from the NCBI database 

	input :
	val sraid // One SRAID

	output : 
	file "*.fastq"

	script:
	"""
	fasterq-dump --threads $params.thread_nb $sraid
	"""
}

process trimmingFastQ {
	input :
	tuple val(sraid),file(sample) 

	output : 
	file "*.fq" 

	shell:
	"""
	cutadapt -e 0.1 -q 20 -O 1 -a !{params.adapter_seq} !{sample} > !{sraid}_trimmed.fq  
	"""
}




workflow {

	// Downloading Reference genome --------------------------------

	RefGenomelink = channel.of(linkRefGenome)
	// RefGenomelink.view()

	// Creating the index 
	ref_genome = downloadRefGenome(RefGenomelink)
	index_files = creatingGenomeIndex(ref_genome)
	// index_files.view()

	// Downloading annotations --------------------------------

	Annotationlink = channel.of(linkAnnotation)
	// Annotationlink.view() 
	annotations = downloadAnnotation(Annotationlink)
 
	// Running locally --------------------------------

	// Retrieving the files path
	fastq_files = channel.fromPath("../data/fastq_files/*.fastq")
	// Getting the name of the fastq files
	fastq_names = fastq_files.map{v -> v.getSimpleName()}
	// Creating a tuple with the name
	tuple_fastq= fastq_names.merge(fastq_files)
	tuple_fastq.view()


	// Downloading from data_bases --------------------------------

// 	sraids = channel.fromList(SRAIDs)
// 	// Download files from database
// 	fastq_files = downloadFastq(sraids) // ça serait cool de pouvoir output les outputs de la commande dans le stdout pendant que ça fonctionne
// 	fastq_files.view()
// 	// Getting the name of the fastq files
// 	fastq_names = fastq_files.map{v -> v.getSimpleName()}
// 	// Creating a tuple with the name
// 	tuple_fastq= fastq_names.merge(fastq_files)

	// Trimming files --------------------------------

	fastq_trimmed = trimmingFastQ(tuple_fastq)
	fastq_trimmed.view()


}
