# Le workflow et les fichiers manipulés

## Process qui peut être fait en amont

Ce sont les process qui doivent être fait en amont avant d'analyser les séquences séparement. 

- ***Téléchargement du génome de référence*** :  
```sh
curl -o reference.fasta "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=CP000253.1&rettype=fasta"
```
C'est le fichier qui contiendra tout le génome de référence sous le format `fasta`


- ***Téléchargement de l'annotation de génome de référence*** :  
```sh
curl -o reference.gff "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?db=nuccore&report=gff3&id=CP000253.1"
```
C'est le fichier qui contiendra tout le génome de référence sous le format `gff`

<span style="color : red"> NB</span> : La commande `curl` est utilisé ici car j'ai plus l'habitude de l'utiliser que wget

- ***Création de l'index du génome*** : 
La commande est la suivante :
```sh
bowtie-build reference.fasta <nom_ref>
```

Dans ce cas là, on utilise les noms de fichier donnés dans le `curl`

à l'issue de ce traitement on se retrouve avec 6 fichiers :
- `<nom_ref>.1.ebwt`
- `<nom_ref>.2.ebwt`
- `<nom_ref>.3.ebwt`
- `<nom_ref>.4.ebwt`
- `<nom_ref>.1.rev.ebwt`
- `<nom_ref>.2.rev.ebwt`

Pour comprendre ce qu'il se passe : [Documentation bowtie](https://rnnh.github.io/bioinfo-notebook/docs/bowtie.html)


---


1. Récupération des séquences grâce à SRA-toolkit :
```sh
fasterq-dump --threads <nb-CPUS> --progress <SRAID>
```

Les séquences sont les SRAID sont les suivants : 
- Persistor 1 (GSM4145661) : SRR10379721
- Persistor 2 (GSM4145662) : SRR10379722
- Persistor 3 (GSM4145663) : SRR10379723
- Control 1 (GSM4145664) : SRR10379724
- Control 2 (GSM4145665) : SRR10379725
- Control 3 (GSM4145666) : SRR10379726

--> ça pourrait être sympa de faire un fichier text pour ne pas avoir à hardcoder les séquences en dur dans le workflow. 
Comme ça on pourrait potentiellement faire tourner le workflow sur d'autres séquences. 

Après ça on se retrouve avec les différents fichiers fastq : `<SRAID>.fastq` 

2. Trimming des données reçue : 

**ça serait sympa de savoir ce que ça fait concrètement.**   
Il semble que c'est un package qui permet de supprimer les séquences qui semblent erronées (par rapport à ?).

Pour faire ça on utilise la commande : 

```sh
trim_galore -q 20 --phred33 --length 25  <FASTQ FILE>
```

À partir de ce traitement on obtient 2 fichiers par fichier `fastq` : 
- `<SRAID>.fastq_trimming_report.txt` : qui contient des informations du trimming
- `<SRAID>_trimmed.fq`

3. Mapping avec bowtie

On utilise la commande suivante : 
```sh 
bowtie -p <CPU> -S <INDEX NAME> <FASTQ FILE>
```

NB : Les fichiers FastQ doivent être unzip. Si ce n'est pas le cas on peut les unzip en utilisant la commande. 

Attention : J'ai essayé de faire tourner ça et de le mettre dans un dossier txt. 
Il me semble que la sortie ça doit être un `.bam`. 
Pour une séquence c'est un fichier de 10Go en sortie

Normalement on pipe direct le résultat dans samtools directement : 

```sh
bowtie -p <NB CPU> -S <Nom index> <Fichier fastq> | samtools sort -@ <Nb CPU> <NAME>.bam
```

Ok je ne sais pas pourquoi mais je n'arrive pas à faire fonctionner bowtie avec le pipe.
Quand je l'output dans un fichier txt ça marche très bien. Mais en dehors de ça et plus rien ne marche. 


UPDATE :
Quand je pipe ça marche, mais j'ai pas tout fait tourné et je sais pas pourquoi mais le mac me demande d'imprimer des documents qui ne sont pas lisibles par des humains. 
Mais on se retrouve aussi avec des fichiers `.bam` de ce format en sortie : 
-  `samtools.40510.7738.tmp.0000.bam`

UPDATE :

Pour l'utilisation de la commande `samtools sort` je ne suis pas passé par le pipe. 
J'ai mis la sortie de la commande `bowtie` dans un fichier txt puis je l'ai mise en input de la commande samtools et j'ai ensuite mis la sortie dans le fichier `bam`que j'ai nommé. 

Grosso modo j'ai fait : 
```sh
bowtie -p <NB_CPUs> -S <Nom index> <Fichier fastq>      > <Nom à choisir>.txt
samtools sort -@ <NB_CPUs> < <Nom à choisir>.txt > <Name>.bam
```




# References / Aide :

## Lien doc et téléchargement

- [Documentation fasterq-dump](https://rnnh.github.io/bioinfo-notebook/docs/fasterq-dump.html)
- [Documentation bowtie](https://rnnh.github.io/bioinfo-notebook/docs/bowtie.html)
- [Documentation samtools](https://rnnh.github.io/bioinfo-notebook/docs/samtools.html)
- [Lien telechargement subreads](https://sourceforge.net/projects/subread/) : Pour avoir featureCounts

## Version des outils à utiliser : 
- `fasterq-dump` : Peu importe la version comme c'est pour le téléchargement des données
- `trim_galore` : Dans le papier je ne l'ai pas trouvé, mais il me semble que c'est construit par dessus cutadapt qui a pour version 1.11
- `bowtie` : 0.12.7
- `subread` (Pour utiliser featurecounts) : 1.4.6-p3
- `R package` : 3.4.1
- `DESeq2` : 1.16


