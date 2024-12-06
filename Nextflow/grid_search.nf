nextflow.preview.output = true

// File to run the workflow over multiple

include { 
	downloadAnnotation;
	downloadRefGenome;
	creatingGenomeIndex;
	downloadFastq;
	trimmingFastQ;
	MappingFQ_samtools;
	FeatureCount;
	processR;
	downloadAureoWiki
	} from "./process_gridsearch.nf"

// ______________________________ Paramètres ______________________________

// Tous les SRAIDs des données disponibles sur NCBI
SRAIDs = [ "SRR10379721", "SRR10379722", "SRR10379723", "SRR10379724", "SRR10379725", "SRR10379726" ]

// ______ Liens utiles ______

// Fichier d'annotation du génome
linkAnnotation = "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?db=nuccore&report=gff3&id=CP000253.1"

// Génome de référence
linkRefGenome = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=CP000253.1&rettype=fasta"

// Méthode pour récupérer les données AureoWiki
link_Aureo_wiki_php = "'https://aureowiki.med.uni-greifswald.de/extensions/AureoDownload/download.php'"
data_aureo_request = "'gsiStrainSelection%5B%5D=NCTC8325%3A0&gsiColumnSelection%5B%5D=pansymbol&includeLinesWithoutData=1&download_token_hidden_field=1733490093889&button=Download+as+Text+%28.tsv%29'"

workflow download_data {

	take:
	Annotationlink
	RefGenomelink
	link_php
	data_request
	sraids

	main:
	annotations = downloadAnnotation(Annotationlink)
	ref_genome = downloadRefGenome(RefGenomelink)
	index_files = creatingGenomeIndex(ref_genome)
	NCTC8325 = downloadAureoWiki(data_request, link_php)
	fastq_files = downloadFastq(sraids) 

	emit:
	annotations
	index_files
	NCTC8325
	fastq_files

}

workflow {
	main : 
		Annotationlink = channel.of(linkAnnotation)
		RefGenomelink = channel.of(linkRefGenome)
		link_php = channel.of(link_Aureo_wiki_php)
		data_request = channel.of(data_aureo_request)
		sraids = channel.fromList(SRAIDs)

		(annotations, index_files, NCTC8325, fastq_files) =  download_data(Annotationlink, RefGenomelink, link_php, data_request, sraids)


		adapter_list = channel.of("AGATCGGAAGAGCACACGTCTGAACTCCAGTCA", "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT" ,"AGATCGGAAGAGC")
		error_values = [ 0.1, 0.15, 0.2 ]
		overlap_min_list = [ 1,2,3 ]
		quality_list = [ 20,25,30 ]

		tuple_fastq = fastq_files.map{v -> [v.getSimpleName(), v]}


 		// Trimming files --------------------------------
 
 		( fastq_trimmed, tuple_settings ) = trimmingFastQ(tuple_fastq, error_values, adapter_list, quality_list, overlap_min_list)
 
 		// Mapping ————————————————————

		tuple_fastq_trimmed = fastq_trimmed.map{v -> [v.getSimpleName() , v]}
		combined_tuple = tuple_fastq_trimmed.combine(index_files) // pour pouvoir utiliser l'index plusieurs fois

		(bam_files, tuple_settings) = MappingFQ_samtools(combined_tuple, tuple_settings)
		list_bam = bam_files.collect()

		(count_file, tuple_settings) = FeatureCount(list_bam,annotations,tuple_settings)


		// R Analysis ————————————————————

		script_r=channel.fromPath("./Plotting_grid_search.R")

		// Plotting
		MAplot_and_settings = processR(count_file,script_r, NCTC8325, tuple_settings)

	publish :
		MAplot_and_settings >> "figures"
}

output {
	"figures" { 
		path { ma_plot_complet, ma_plot_enhanced, settings -> "figures/${settings}" }
	}
}
