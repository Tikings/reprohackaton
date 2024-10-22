SRAIDs = [
"SRR10379721",
"SRR10379722",
"SRR10379723",
"SRR10379724",
"SRR10379725",
"SRR10379726",
]


params.thread_nb = 2
params.linkAnnotation = "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?db=nuccore&report=gff3&id=CP000253.1"
params.linkRefGenome = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=CP000253.1&rettype=fasta"


process downloadAnnotation {
	output :
	file "annotation.gff"

	script : 
	"""
	curl -o annotation.gff $params.linkAnnotation
	"""

}


process downloadRefGenome {
	output : 

	script : 
	""" 
	curl -o reference.fasta "$params.linkRefGenome"
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
	ref_genome = downloadRefGenome()
	annotations = downloadAnnotation()
	sraids = channel.from(SRAIDs)
//	fastq_files = downloadFastq(sraids) // ça serait cool de pouvoir output les outputs de la commande dans le stdout pendant que ça fonctionne
//	fastq_files.view()

	fastq_files = channel.fromPath("../data/*.fastq")
	trimmed_fastq_files = trimmingFastQ(fastq_files)


}
