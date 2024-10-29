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

process downloadFastq {

	input :
	val sraid

	output : 
	file "*.fastq"

	script:
	"""
	fasterq-dump --threads $params.thread_nb --progress $sraid
	"""
}

process trimmingFastQ {
	input :
	file sample

	output : 
	file "*.txt"
	file "*.fq"

	script:
	"""
	trim_galore -q 20 --phred33 --length 25 $sample
	"""
}




workflow {
	RefGenomelink = channel.of(linkRefGenome)
	// RefGenomelink.view()
	ref_genome = downloadRefGenome(RefGenomelink)
	index_files = creatingGenomeIndex(ref_genome)
	// index_files.view()

	Annotationlink = channel.of(linkAnnotation)
	// Annotationlink.view() 
	annotations = downloadAnnotation(Annotationlink)

	sraids = channel.from(SRAIDs)
//	fastq_files = downloadFastq(sraids) // ça serait cool de pouvoir output les outputs de la commande dans le stdout pendant que ça fonctionne
//	fastq_files.view()

	fastq_files = channel.fromPath("../data/*.fastq")
	trimmed_fastq_files = trimmingFastQ(fastq_files)

}
