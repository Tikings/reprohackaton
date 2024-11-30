// Dossier pour l'enregistrement des données finales 
params.pub_dir = "../results/"

// Tous les SRAIDs des données disponibles sur NCBI
SRAIDs = [ "SRR10379721", "SRR10379722", "SRR10379723", "SRR10379724", "SRR10379725", "SRR10379726" ]

// ______ Liens utiles ______

// Fichier d'annotation du génome
linkAnnotation = "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?db=nuccore&report=gff3&id=CP000253.1"

// Génome de référence
linkRefGenome = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=CP000253.1&rettype=fasta"

// Méthode pour récupérer les données AureoWiki
link_Aureo_wiki_php = "'https://aureowiki.med.uni-greifswald.de/extensions/AureoDownload/download.php'"
data_aureo_request = "'gsiStrainSelection%5B%5D=NCTC8325%3A0&gsiColumnSelection%5B%5D=panlocus&gsiColumnSelection%5B%5D=pansymbol&gsiColumnSelection%5B%5D=symbol&gsiColumnSelection%5B%5D=synonym&gsiColumnSelection%5B%5D=geneGI&includeLinesWithoutData=1&download_token_hidden_field=1732555537548&button=Download+as+Text+%28.tsv%29'"

// Nom que l'on va donner à l'index bowtie 
params.RefName = "INDEX_S_aureus"

// Séquence d'adapter
params.adapter_seq = "AGATCGGAAGAGC"


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

	output : 
	file "*.fastq.gz" 

	shell:
	"""
	cutadapt -e 0.1 -q 20 -O 1 -a !{params.adapter_seq} -o !{sraid}_trimmed.fastq.gz !{sample}
	"""
}

// Mapping the sequences
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


// Counting 
process FeatureCount{
	publishDir params.pub_dir

        input :

        val list_bam    // All the bam files 
        file annotation // Annotation file

        output :
        file "counts.txt" // Count files 

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
	publichDir params.pub_dir

        input :
        file count // Count file 
        file script_r // R script 
	file NCTC8325 

        output :
        file "*.pdf" 

        shell:
        """
        R -f !{script_r} 
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

	link_php = channel.of(link_Aureo_wiki_php)
	data_request = channel.of(data_aureo_request)
	NCTC8325 = downloadAureoWiki(data_request, link_php)

        script_r=channel.fromPath("../Plotting_script.R")


	// Plotting
	MAplot = processR(count_file,script_r, NCTC8325)

}
