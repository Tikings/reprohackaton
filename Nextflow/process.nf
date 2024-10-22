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


workflow {
	sraids = channel.from(SRAIDs)
	fastq_files = downloadFastq(sraids) // ça serait cool de pouvoir output les outputs de la commande dans le stdout pendant que ça fonctionne

}
