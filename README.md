# BioE-201-230

Repo for assignments after week 3

## Assignment Week 3

### Download the Dataset Through NCBI-Dataset Package

#### Install the package

conda install -c conda-forge ncbi-datasets-cli

#### Download the dataset

datasets download genome taxon 2 --released-after 1980-01-01 --released-before 2001-12-31 --include genome

Tab-separated values (TSV) downloaded through website and uploaded through MobaXterm

| name                                                | GC Percent | Contig N50 | scaffolds count | size    |
|:---------------------------------------------------:|:----------:|:----------:|:---------------:|:-------:|
| Vibrio cholerae O1 biovar El Tor str. N16961        | 47.5       | 2961149    | 2               | 4033464 |
| Brucella melitensis bv. 1 str. 16M                  | 57         | 2117144    | 2               | 3294931 |
| Deinococcus radiodurans R1 = ATCC 13939 = DSM 20539 | 66.5       | 2648638    | 4               | 3284156 |
| Lactococcus lactis subsp. lactis Il1403             | 35.5       | 2365589    | 1               | 2365589 |
| Pasteurella multocida subsp. multocida str. Pm70    | 40.5       | 2257487    | 1               | 2257487 |
| Thermotoga maritima MSB8                            | 46         | 1860725    | 1               | 1860725 |
| Haemophilus influenzae Rd KW20                      | 38         | 1830138    | 1               | 1830138 |
| Helicobacter pylori 26695                           | 39         | 1667867    | 1               | 1667867 |
| Helicobacter pylori J99                             | 39         | 1643831    | 1               | 1643831 |
| Aquifex aeolicus VF5                                | 43.5       | 1551335    | 2               | 1590791 |
| Chlamydia pneumoniae CWL029                         | 40.5       | 1230230    | 1               | 1230230 |
| Chlamydia pneumoniae AR39                           | 40.5       | 1229853    | 1               | 1229853 |
| Treponema pallidum subsp. pallidum str. Nichols     | 53         | 1138011    | 1               | 1138011 |
| Chlamydia trachomatis D/UW-3/CX                     | 41.5       | 1042519    | 1               | 1042519 |

### Find the smallest (largest) genome

#### Biggest size ouput:

awk -F'\t' 'NR == 1 || \$5 > min {min = $5; line = $0} END {print line}' summary.tsv

**Output:** Vibrio cholerae O1 biovar El Tor str. N16961    47.5    2961149 2       4033464

#### Smallest size ouput:

awk -F'\t' 'NR == 1 || \$5 < min {min = $5; line = $0} END {print line}' summary.tsv

**Output:** Chlamydia trachomatis D/UW-3/CX 41.5    1042519 1       1042519

#### Biggest with size ouput only:

awk -F'\t' 'NR==1{next} {if(max=="" || \$5 > max){max=$5}} END{print max}' summary.tsv

**Output:** 4033464

#### Smallest with size ouput only:

awk -F'\t' 'NR==1{next} {if(max=="" || \$5 < max){max=$5}} END{print max}' summary.tsv

**Output:** 1042519

### Find the number of genomes that contain at least two “c” in the species name

awk -F'\t' 'BEGIN{count=0} {if(gsub(/[Cc]/,"c",$1) >= 2) count++} END{print count}' summary.tsv

**Output:** 7

### How many of the species names contain two or more “c” but do not contain the word “coccus”

awk -F'\t' 'BEGIN{count=0} {if(gsub(/c/,"c",\$1) >= 2 && $1 !~ /coccus/) count++} END{print count}' summary.tsv

**Output:** 3

### Use the find command to find all genome files (FASTA) larger than 3MB

awk -F'\t' '$5 > 3000000 {count++} END{print count}' summary.tsv

**Output:** 3
