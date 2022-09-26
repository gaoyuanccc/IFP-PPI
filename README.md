# IFP-PPI
Interaction Frequency Predict Protein-Protein Interaction (IFP-PPI)

E-mail: 18434753515@163.com
## Introduction

IFP-PPI is a tool for building gene association networks that use the interaction frequency of chromosome segments in space to estimate the spatial interaction frequency between genes and identify chromosomal interaction regions (CIDs). With normalized or unnormalized IF matrices and bacterial genome-wide Genbank annotation files as inputs, the software will construct the gene association network, identify the CIDs, and export the network as a file.


## Installation
Requirement:
* python (3.x)
* biopython
* numpy
* pandas
* scipy

Notes:
You may alternatively install using the following `conda` commands shown below, which are derived from our built environment file `environment.yaml`.

    conda env create -f environment.yaml -n new_env_name


##  How to use IFP-PPI?

### Parameters
| Parameter |Full name| Description| Default 
|:-:|:-:|--|:-:|
|-gb|genbank| Genbank filename with annotation information for the strain. |\|
|-p|multichromosome|Used to determine if the target strain is multichromosome; if the strain is multichromosome, set to 1 and is 0 by default.|0|
|-gbc|multichromosomeGenbank|If the strain is multichromosome, enter the Genbank filenames for the chromosomes in order, separated by a '/' symbol.|\|
|-im |IFmatrix| Contact matrix filename to convert interaction frequencies between bins to interaction frequencies between genes. |\|
|-ir|IFmatrixResolution|Resolution of IFmatrix.|\|
|-i|input|Input folder name.|input||
|-b|binBeginNumber |Begin the bin serial number in the contact matrix, which is 0 by default.|0|
|-n|genepairNumber|The number of selected gene pairs and is 10000 by default.|10000|
|-d|removeDistance|Linear genomic distances to select removed gene pairs and is 300 by default.|300|
|-cm|CIDmatrix|Contact matrix filename to identify CIDs and is IFmatrix by default. If CIDmatrix = None, then it will not identify CIDs|IFmatrix|
|-cr|CIDmatrixResolution|Resolution of CIDmatrix and is IFmatrixResolution by default.|IFmatrixResolution|
|-o|output|Output folder name.|output|
|-h|help|Help documentation|\|

### Usage

    python IFP-PPI.py -i [input folder] -gb [genbank] -im [IFmatrix] -ir [IFmatrixResolution] [...]

#### For example

 - Single-chromosome
The GSM2870409 sample of *Escherichia coli* K-12 MG1655 was taken is an example to demonstrate the complete function.
````
    # Run IFP-PPI
     python IFP-PPI.py -i input -gb U00096.3.gb -im GSM2870409_1000_iced.matrix -ir 1000 -cm GSM2870409_10000_iced.matrix -cr 10000 -o output
    ------------------------------
	2022-08-13 18:32:32
	Run IFP-PPI
	
	------------------------------
	2022-08-13 18:32:32
	Convert interaction frequency
	
	------------------------------
	2022-08-13 18:33:14
	Select gene pairs

	------------------------------
	2022-08-13 18:33:14
	Remove close gene pairs

	------------------------------
	2022-08-13 18:33:14
	Identify CIDs

	Time:    42.5720648765564s
````

 - Multi-chromosome
 The SRR3180951 sample of *Vibrio cholerae*  was taken is an example to demonstrate the complete function.
 ````
 #Run IFP-PPI
python IFP-PPI.py -i input -p 1 -gbc vibrio_cholerae_1.gb/vibrio_cholerae_2.gb -im SRR3180951_1000_iced.matrix -ir 1000 -cm SRR3180951_10000_iced.matrix -cr 10000 -o output -b 1
------------------------------
2022-08-25 17:04:16
Run IFP-PPI

------------------------------
2022-08-25 17:04:16
Convert interaction frequency

------------------------------
2022-08-25 17:04:44
Select gene pairs

------------------------------
2022-08-25 17:04:44
Remove close gene pairs

------------------------------
2022-08-25 17:04:44
Identify CIDs

Time:    41.56097483634949s
 ````

#### Detailed usage

 1. -i input
Put all input files into the  'input', folder before starting the run. 'input' is the folder name where all input files are located.
 2. -o output
 'output' is the name of the folder where all output files are located. The software will output 6 files if it identifies CIDs. '*cid_range.txt*' and '*cid_boundary.txt*' are files describing CID ranges and boundaries. *’bin_t_p.txt*' is the leftward and rightward interaction preference file for bin. '*gene_association_network.txt*', '*gene_cid.txt*' and '*gene_all_cid.txt*' are files from the obtained gene association network and the CID attribute of the gene. '*gene_association_network.txt*' is visualized by Cytoscape. Its specific structure is as follows:
 ```
	+ input
		++ U00096.3.gb (genbank)
		++ GSM2870409_1000_iced.matrix (IFmatrix)
		++ GSM2870409_10000_iced.matrix (CIDmatrix)
	+ output
		+ cid
			+ cid_range
				++ cid_range.txt
				++ cid_boundary.txt
			+ identify_cid
				++ bin_t_p.txt
		+ network
			++ gene_association_network.txt
			++ gene_cid.txt
			++ gene_all_cid.txt
```
 3.  -gb genbank
 This command inputs the Genbank annotation file of the strain, which needs to contain the genome length and the position information of the gene, and its format is as follows:
 ````
 LOCUS       U00096               4641652 bp    DNA     circular BCT 23-SEP-2020
DEFINITION  Escherichia coli str. K-12 substr. MG1655, complete genome.
            ......
            ......
     gene            190..255
                     /gene="thrL"
                     /locus_tag="b0001"
                     /gene_synonym="ECK0001"
                     /db_xref="ASAP:ABE-0000006"
                     /db_xref="ECOCYC:EG11277"
     CDS             190..255
                     /gene="thrL"
                     /locus_tag="b0001"
                     /gene_synonym="ECK0001"
                     /codon_start=1
                     /transl_table=11
                     /product="thr operon leader peptide"
                     /protein_id="AAC73112.1"
                     /db_xref="UniProtKB/Swiss-Prot:P0AD86"
                     /db_xref="ASAP:ABE-0000006"
                     /db_xref="ECOCYC:EG11277"
                     /translation="MKRISTTITTTITITTGNGAG"
````
 4.  -im IFmatrix -ir IFmatrixResolution
This command inputs the contact matrix file that converts the interaction frequency, and its format includes two kinds, matrix and three column data. The resolution is this matrix resolution, 1000 is recommended.	
**Matrix format**
		`````
		65.478946	633.67299	9.645624	5.717007	5.582881	...	3337.363203
		633.67299	893.962839	1387.628772	32.182989	10.138047	...	16.641758
		9.645624	1387.628772	1368.663351	1164.003637	17.234005	...	9.224961
		5.717007	32.182989	1164.003637	988.886841	1602.016586	...	6.561214
		5.582881	10.138047	17.234005	1602.016586	503.522212	...	3.203641
		...			...			...			...			...			...		...
		3337.363203	16.641758	9.224961	6.561214	3.203641	...	532.456055
		`````
		**Three columns format**
The data consists of three columns: the first column is the serial number of the first bin, the second column is the serial number of the second bin, and the third column is the interaction frequency of the two bins.
		``````
		0	0	65.478946
		0	1	633.672990
		0	2	9.645624
		0	3	5.717007
		0	4	5.582881
		...	...	...
		4641	4641	532.456055
		``````
 5.  -cm --CIDmatrix -cr CIDmatrixResolution
 If you need to identify CIDs, provide the contact matrix in the same format as the above IFmatrix; 10000 is the suggested resolution. If you do not need to identify CIDs, enter the command `-cm None`.
 6. -p --multichromosome -gbc --multichromosomeGenbank
 If the strain is multi-chromosome, set -p 1. Genbank filenames for chromosomes are entered in order after -gbc, separated by a '/'.

 

 
 

 


