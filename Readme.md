# Pamgene phosphorylation chip analysis
Custom homebrew solution to analyse data from [Pamgene](https://pamgene.com/kinase-activity-profiling-services/) kinase activity profiling system, created due to issues with the original Bionavigator software. Workflow was reverse engineered and approximated from various presentations and posters found online and as such is less complex and complete then the one offered commercially.

## Example usage
`
pamgene.R -i raw_data_stk.csv --pamgene pamgene.csv --phosphonet phosphonet.csv --ctrl A1 --exp A2 --kinexus_score 300 --batch_correction
`

## Commandline arguments
`-i, --file_in` input file containg data from pamgene kinase activity profiling system. See preview of input formats for more information.

`-o, --dir_out` output directory

`--pamgene` csv file containing information about peptides on the pamgene chip

`--phosphonet` csv file containing information about kinases phosphorylating specific protein residues of proteins included im pamgene csv file

`-c, --ctrl` control group name as a single string or comma-separated indices of alphabetically arranged <group_name>_<chip_id> denominators

`-e, --exp` experimental group name as a single string or comma-separated indices of alphabetically arranged <group_name>_<chip_id> denominators

`--top_kinases` top kinases from phosphonet file to include in calculations, if NULL the number will be estimated from data

`--kinexus_score` threshold kinexus version 2 score that determines whether the residue is considered as being phosphorylated by a kinase. All phosphorylations below threshold are omitted.

`--no_perm` number of permutations used to calculate specificity and selectivity scores

`-b, --batch_correction` whether to remove batch effect from peptide kinetics matrix. Uses sva::ComBat for calculations.

## Example output

### Peptide statistics
| id           	| mean_ctrl 	| mean_exp 	| statistic 	| pvalue 	| A1_710250713 	| A1_710250715 	| A2_710250713 	| A2_710250715 	|
|--------------	|-----------	|----------	|-----------	|--------	|--------------	|--------------	|--------------	|--------------	|
| ACM1_421_433 	| 6.79      	| 6.38     	| 0.65      	| 0.55   	| 6.54         	| 6.32         	| 6.54         	| 6.76         	|
| ACM1_444_456 	| 6.70      	| 6.46     	| 0.37      	| 0.74   	| 6.53         	| 6.42         	| 6.85         	| 7.02         	|
| ACM4_456_468 	| 6.44      	| 6.03     	| 0.69      	| 0.52   	| 6.21         	| 5.87         	| 6.40         	| 6.31         	|
| ACM5_494_506 	| 11.02     	| 10.70    	| 0.22      	| 0.84   	| 11.39        	| 10.67        	| 11.06        	| 11.62        	|

### Kinase scores
| kinase_id 	| kinase_name          	| specificity 	| selectivity 	| total 	| sum_score 	| mean_score 	| sd_score 	|
|-----------	|----------------------	|-------------	|-------------	|-------	|-----------	|------------	|----------	|
| Q9UBS0    	| p70S6Kb (RPS6KB2)    	| 0.70        	| 1.00        	| 1.70  	| -0.28     	| -0.01      	| 0.03     	|
| P41743    	| PKCi (PRKCI)         	| 0.36        	| 1.00        	| 1.36  	| -0.66     	| -0.33      	| 0.48     	|
| P22694    	| PKACb PRKACB         	| 1.00        	| 0.34        	| 1.34  	| -0.30     	| -0.01      	| 0.02     	|
| P25098    	| BARK1 (ADRBK1; GRK2) 	| 0.20        	| 1.00        	| 1.20  	| -0.45     	| -0.23      	| 0.42     	|

## Preview of input formats

### Raw data
Csv file with one of the following structures:

| ID         	| UniprotAccession 	| \<group1\>\_\<chip_id1\>\_\<time1\> 	| \<group1\>\_\<chip_id1\>\_\<time2\>  	| etc.|
|------------	|------------	|-------------	|-------------	|-------------	|

for all groups, chips and time points, or same data transformed into long format with exactly 6 columns:

| ID         	| UniprotAccession 	| group 	| chip_id	| time 	| value 	|
|------------	|------------	|-------------	|-------------	|-------------	|-------------	|


### Pamgene

File describing the peptides included on the pamgene chip. Information can be inferred from the peptides reported in the raw data file by [aligning](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastp) the peptide protein sequences and adding proteins that pass similarity threshold and contain phosphoryltable amino acids in the corresponding frequences. Following columns are required in any order.

| id         	| uniprot_id 	| phosphosite 	|
|------------	|------------	|-------------	|
| ATF2_47_59 	| P15336     	| 69          	|
| ATF2_47_59 	| P15336     	| 71          	|
| ATF2_47_59 	| P15336     	| 73          	|
| ATF2_47_59 	| P17544     	| 51          	|


### Phosphonet
With the knowledge of peptides included in the chip, phosphonet file can be generated by [Phosphonet scraper](https://github.com/icervenka/phosphonet_scraper). Following columns are required in any order.

| substrate 	| site 	| kinase_rank 	| kinase_name   	| kinase_id 	| kinexus_score_v2 	|
|-----------	|------	|-------------	|---------------	|-----------	|------------------	|
| P15336    	| 69   	| 1           	| JNK1 (MAPK8)  	| P45983    	| 633              	|
| P15336    	| 69   	| 2           	| JNK3 (MAPK10) 	| P53779    	| 633              	|
| P15336    	| 69   	| 3           	| ERK2 (MAPK1)  	| P28482    	| 615              	|
| P15336    	| 69   	| 4           	| ERK1          	| P27361    	| 611              	|

## Known Issues and limitations
- Script is only able to calculate kinase activity comparing two groups (control vs. experimental)
- Pamgene and Phosphonet files are not supplied, users are required to generate their own
