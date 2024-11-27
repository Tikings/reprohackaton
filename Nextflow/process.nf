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

process trimmingFastQ {
	input :
	tuple val(sraid),file(sample) 

	output : 
	file "*.fastq.gz" 

	shell:
	"""
	cutadapt -e 0.1 -q 20 -O 1 -a !{params.adapter_seq} -o !{sraid}_trimmed.fastq.gz !{sample}
	"""
}


process MappingFQ_samtools{
        input :

        tuple val(sraid),file(sample),val(RefName),file(index_files)


        output :
        file "*.bam"

        shell:
        """
        bowtie -S !{RefName} <(gunzip -c !{sample}) | samtools sort > !{sraid}.bam
        """

}


process FeatureCount{
        input :

        val list_bam
        file annotation

        output :
        file "counts.txt"

        shell:
        """
        featureCounts -t gene -g ID -F GTF -s 1 -a !{annotation} -o counts.txt !{list_bam.join(' ')}
        """
}


process processR{
        input :
        file count
        file script_r
        output :
        file "output.log"


        shell:
        """
        R -f !{script_r}  > output.log
        """
}





workflow {

        // Downloading Reference genome --------------------------------
        RefGenomelink = channel.of(linkRefGenome)

        // Creating the index --------------------------------
        ref_genome = downloadRefGenome(RefGenomelink)
        index_files = creatingGenomeIndex(ref_genome)


        // Downloading annotations --------------------------------
        Annotationlink = channel.of(linkAnnotation)
        annotations = downloadAnnotation(Annotationlink)


        // Downloading from data_bases --------------------------------

        sraids = channel.fromList(SRAIDs)
        // Download files from database
        fastq_files = downloadFastq(sraids) 
        // Creating a tuple with the name
        tuple_fastq= fastq_files.map{v -> [v.getSimpleName(), v]}

        // Trimming files --------------------------------

        fastq_trimmed = trimmingFastQ(tuple_fastq)

        // Mapping ————————————————————

        tuple_fastq_trimmed= fastq_trimmed.map{v -> [v.getSimpleName() , v]}
        combined_tuple=tuple_fastq_trimmed.combine(index_files) // pour pouvoir utiliser l'index plusieurs fois

        bam_files=MappingFQ_samtools(combined_tuple)
        list_bam=bam_files.collect()

        count_file=FeatureCount(list_bam,annotations)

        // R Analysis ————————————————————

        //output_file=processR(count_file)

        script_r=channel.fromPath("./script.R")
        script_r.view()
        output_2=processR(count_file,script_r)
}
