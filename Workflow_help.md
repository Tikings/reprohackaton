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

# Description des étapes du process


## 1. Récupération des séquences grâce à SRA-toolkit :
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

Pour travailler plus facilement avec les fichiers on peut gzipper les fichiers -> La majorité des fonctions peuvent travailler avec ces fichiers. \
Pour cela on peut mettre dans le process : 
```sh
fasterq-dump --threads <nb-CPUs> --progress <SRAID>
gzip *.fastq
```

La compression des fichiers est efficace, par exemple pour SRR10379721 qui fait 10Gb. Sa compression donne un fichier de 2Gb. Cependant c'est une étape qui prend un temps comsidérable. 
Voir ce qui est le mieux.

## 2. Trimming des données reçue : 

C'est un module qui permet de supprimer les read de mauvaise qualité (c'est l'option : `--phred33`) et supprimer les adapters utilisés pour le séquençage. 
On supprime également les reads qui sont trop petits (`--length 25` qui est la longueur minimale).

Pour faire ça on utilise la commande : 

```sh
trim_galore -q 20 --phred33 --length 25  <FASTQ FILE>
```

À partir de ce traitement on obtient 2 fichiers par fichier `fastq` : 
- `<SRAID>.fastq_trimming_report.txt` : qui contient des informations du trimming
- `<SRAID>_trimmed.fq`

Ici on utilise la fonction `trim_galore` mais dans le papier ils utilisent la fonction `cutadapt` qui en fait est un module sur lequel se base `trim_galore`.
Pour se rapprocher des résultats du papier on va donc se focaliser sur utiliser cette fonction ci. 

En utilisant `cutadapt` la commande change, pour choisir les paramètres on utilise les mêmes settings que pour `trim_galore` (on les retrouve dans le report du trimming par trim_galore)

Exemple de settings : 

- Quality Phred score cutoff: 20
- Quality encoding type selected: ASCII+33
- Using Nextera adapter for trimming (count: 1). Second best hit was smallRNA (count: 0)
- Adapter sequence: 'CTGTCTCTTATA' (Nextera Transposase sequence; auto-detected) <span style="color : red"> Attention les adapters ne sont pas les mêmes parce que visiblement les machines de séquençage ne sont pas pareil</span>
- Maximum trimming error rate: 0.1 (default)
- Minimum required adapter overlap (stringency): 1 bp
- Minimum required sequence length before a sequence gets removed: 25 bp

Dans le trimming report fournit par trim_galore après le traitement du fichier, ils mentionnent la commande utilisée :
```sh
cutadapt -j 1 -e 0.1 -q 20 -O 1 -a AGATCGGAAGAGC SRR10379721.fastq
```
avec : 
- `-j` : Le nombre de coeurs à utiliser pour l'analyse (c'est la seule option qui n'est pas dans la version du papier)
- `-e` : Error-rate 
- `-q` : Quality score cutoff
- `-O`: Minimum required adapter overlap 1
- `-a` : La sequence d'adapter qui est utilisée (Attention pour SRR10379722 `trim_galore` n'a pas détecté un adapter de séquenceur illumina)

On va donc utiliser la commande suivante qui fonctionne avec la version 1.11 (celle du papier) de `cutadapt`: 

```sh
cutadapt -e 0.1 -q 20 -O 1 -a AGATCGGAAGAGC SRR10379721.fastq
```

On peut utiliser cette commande avec les fichiers compressés mais il faut préciser l'ouput du fichier : 

```sh
cutadapt -e 0.1 -q 20 -O 1 -a AGATCGGAAGAGC <SRAID>.fastq.gz > <SRAID>_trimmed.fq 
```

L'output de cette commande reste également assez volumineuse, est ce que ça vaut le coup de la compresser aussi ? \ 
Ça voudrait dire une baisse de performance (en temps) pour le workflow mais on aurait un meilleur management des données disponibles.

J'ai essayé de faire le différence des fichiers en sortie en utilisant trim_galore et cutadapt mais ça m'a pris trop de temps et j'ai du arreter... 

Avec la version 1.11 (celle du papier), la commande fonctionne bien (aussi bien avec les fichiers commpressés qu'avec ceux qui ne le sont pas). 


## 3. Mapping avec bowtie

On utilise la commande suivante : 
```sh 
bowtie -p <CPU> -S <INDEX NAME> <FASTQ FILE>
```

NB : Les fichiers FastQ doivent être unzip. Si ce n'est pas le cas on peut les unzip en utilisant la commande. 

On peut directement utiliser la commande suivante (dans le cas où les fichiers n'ont pas été zippé) en pipant le résultat de `bowtie` dans `samtools sort`

```sh
bowtie -p <NB CPU> -S <Nom index> <Fichier fastq> | samtools sort -@ <Nb CPU> > <NAME>.bam
```
Et si jamais on veut utiliser les fichiers `fastq` zippé, on peut utiliser :

```sh
bowtie -p <NB CPU> -S <Nom index> <(gunzip -c <Fichier gz>) | samtools sort -@ <Nb CPU> > <NAME>.bam
```
NB : La partie `<NAME>` correspond à un nom choisit par l'utilisateur --> Il faudrait qu'on la définissent plus tard. Je pensais à utiliser le SRAID

L'output final de ce process est un ensemble de fichiers `bam` qui contiennent les données de séquences alignées pour chaque `SRAID`. On a donc finalement: `<NAME_SRAID>.bam`

## 4. Counting des séquences : 

Une fois tous les fichiers `bam` produit on peut utiliser la dernière librairie avant de passer à l'analyse sur R. \ 
Celle ci repose sur `featureCounts` de la librairie `subread`. 
Elle permet de compter le nombre de reads alignés sur un gène donnée (obtenu grâce aux annotations de génome qu'on a téléchargé précédemment)

Le commande donnée dans les indications (cf TP_nf_sn) est la suivante : 

```sh
featureCounts --extraAttributes Name -t gene -g ID -F GTF -T <Nb-CPU> -a <annotation.gff> -o counts.txt <BAM Files>
```

Mais avec la version du papier ça ne fonctionne pas, on utilise donc la commande suivante : 

```sh
featureCounts -t gene -g ID -F GTF -T <Nb-CPU> -a <annotation.gff> -o counts.txt <BAM Files>
```

Il ajoute également une partie en plus à la commande : 

```sh
featureCounts -t gene -g ID -F GTF -s 1 -T <Nb-CPU> -a <annotation.gff> -o counts.txt <BAM Files>
```

`<BAM files>` correspond à tous les fichiers `bam` à la fois (== `*.bam`) \ 
Et l'output final est un fichier `counts.txt` qui contient les données nécessaires à l'analyse R.

# Explication de la structure des fichiers

## FastQ Files que l'on télécharge depuis NCBI

Ce sont des fichiers qui contiennent les données de séquençage RNA-seq avec les scores de qualité associés. Une séquence se présente sous la forme: 
```
@SRR10379721.1 1 length=65
CGGCAGACAAACAAAAGAATGGAATCAAAGCTAACTTCAAAATTCGCCACAACATTGAAGATGGA
+SRR10379721.1 1 length=65
FFFFFBF/FFFFFFFFFFFFFFFFFFFFFFF/FFFF/FFFFFFFFFBFFFFFFFFFFFFFBBBBB
```
La première partie représente la séquence, et la deuxième représente un chaine qui permet de qualifier la qualité des données.

## Données `<NOM>_trimmed.fastq` obtenues après le trimming

C'est l'équivalent du fichier FastQ téléchargé depuis la base donnée mais avec les séquences avec une faible qualité (On la quantifie à partir de la probabilité que la séquence contient des erreurs).
On supprime également les adapteurs utilisés pour le RNA-seq : Explication dans [ce papier](https://support.illumina.com/content/dam/illumina-support/documents/documentation/chemistry_documentation/samplepreps_truseq/truseq-stranded-mrna-workflow/truseq-stranded-mrna-workflow-reference-1000000040498-00.pdf)


## Le génome de référence

C'est uniquement un fichier texte avec la séquence du génome de référence.

## Index du génome de référence 

Ce sont les fichier `ebwt`. Ils permettent de faciliter l'alignement des séquences des fichiers fastq avec bowtie. 

## Les fichiers `bam`

C'est le format compressé des fichiers `sam` qui veut dire : "Sequence alignement map". C'est tout simplement la localisation des séquences alignées

## Le fichier d'output de `featureCounts`

C'est un fichier qui contient un matrice de l'association des alignements de chaque séquence avec le gène associé à la portion de génome sur lequel il est aligné. 
La localisation des gènes est obtenue à partir du fichier d'annotation du génome (téléchargé depuis NCBI).

<span style="color : #F97316"> TODO -> Peut être détailler les colonnes que l'on a dans l'output</span>


# References / Aide :

## Lien doc et téléchargement

- [Documentation fasterq-dump](https://rnnh.github.io/bioinfo-notebook/docs/fasterq-dump.html)
- [Documentation bowtie](https://rnnh.github.io/bioinfo-notebook/docs/bowtie.html)
- [Documentation samtools](https://rnnh.github.io/bioinfo-notebook/docs/samtools.html)
- [Lien telechargement subreads](https://sourceforge.net/projects/subread/) : Pour avoir featureCounts
- [RNA-seq workflow: gene-level exploratory analysis and differential expression](https://master.bioconductor.org/packages/release/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html#building-the-results-table) : Pour comprendre le fonctionnement du workflow
- [Trueseq Stranded mRNA Guide Reference](https://support.illumina.com/content/dam/illumina-support/documents/documentation/chemistry_documentation/samplepreps_truseq/truseq-stranded-mrna-workflow/truseq-stranded-mrna-workflow-reference-1000000040498-00.pdf) : Kit utilisé dans le papier pour qui explique le processus de séquençage


## Version des outils à utiliser : 
- `fasterq-dump` : Peu importe la version comme c'est pour le téléchargement des données
- `trim_galore` : Dans le papier je ne l'ai pas trouvé, mais il me semble que c'est construit par dessus cutadapt qui a pour version 1.11
- `bowtie` : 0.12.7
- `subread` (Pour utiliser featurecounts) : 1.4.6-p3
- `R package` : 3.4.1
- `DESeq2` : 1.16


