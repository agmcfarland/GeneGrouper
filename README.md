
**GeneGrouper** finds gene clusters (genes that are next to each other in the genome) that contain a gene of interest and separates them into groups of gene clusters based on the similarity of their gene content. The goal of GeneGrouper is to allow for understanding how gene content can vary, if only slightly, in up to thousands of similar gene clusters, and how those gene clusters are distributed in genomes.

<<<<<<< Updated upstream
GeneGrouper is a command-line tool that finds gene clusters that contain a gene of interest in a set of genomes and bins them into groups of similar gene clusters.

=======
## Table of contents 
>>>>>>> Stashed changes

1. **Why use GeneGrouper**

2. **Example application**

3. **Installation**

4. **Usage**

5. **Example dataset and searches**

6. **Search output file structure and file descriptions**

<<<<<<< Updated upstream
1. One translated gene of interest (.faa/.fasta/.txt)

2. Two or more genomes from RefSeq (.gbff)
=======
7. **Creating a conda environment with all dependencies**

8. **FAQ**
>>>>>>> Stashed changes

9. **Citation**

<<<<<<< Updated upstream
1. For each individual search, a new folder will be outputted containing all gene clusters and their groupings
=======
10. **References**

>>>>>>> Stashed changes

# Why use GeneGrouper?

GeneGrouper searches many genomes for a query gene using BLAST. When GeneGrouper finds a hit, it uses that gene as a seed and extracts the surrounding upstream and downstream genes **(Fig. 1 A-C)**. The extracted regions are then separated into groups of regions that share similar gene content **(Fig. 1 D)**. 

In this way, GeneGrouper can show the user all regions that have different gene content (separated by group), and those that have similar gene content (found within groups) **(Fig. 1 D)**. 

This approach can be used to find whether a specific gene cluster is found in all searched genomes, and whether each gene cluster that is found has all the expected genes, or whether there is unusual gene content.

<img src="docs/overview_figure.png" alt="GeneGrouper overview figure" width=1000>

Figure 1: GeneGrouper overview (A-C). Panel D shows GeneGrouper results after searching 1,130 *Salmonella enterica* genomes for regions containing *pduA* homologs.

# Example application

**Introduction**

We wanted to know whether the catabolic Pdu gene cluster was present in 1,130 *Salmonella enterica*. This gene cluster is made up of 23 genes **(Fig. 2)**. We expected one copy of the intact Pdu gene cluster per genome.

<img src="docs/pdu_gene_cluster.png" alt="Pdu gene cluster figure" width=1000>

**Methods**

<<<<<<< Updated upstream
You only need to make a database of the genomes once.

```
GeneGrouper -g /path/to/gbff -d /path/to/output_directory \
build_database
```
=======
We used GeneGrouper  to search for all occurrences of the *pduA* gene, which is an important component of the Pdu gene cluster. We specified that if *pduA* is found, extract all genes 2,000 bp upstream and 18,000 bp downstream of it **(Fig. 2)**. This should encompass all 23 Pdu gene cluster genes. 
>>>>>>> Stashed changes

**Results**

<<<<<<< Updated upstream
Now you can search the database of genomes for gene clusters that contain your gene of interest! 

```
GeneGrouper -d /path/to/output_directory -n gene_name \
find_regions -f /path/to/seed_gene.fasta 
```
=======
We had 2,252 *pduA* hits and the region surrounding each *pduA* gene homolog was extracted as defined above. GeneGrouper separated all 2,252 regions into five different groups according to their gene content. We can see in the output (Fig. 1 D center) that group 0 is composed of 1,120 regions that have the genes we expected to find in the Pdu gene cluster. We can also see how dissimilar the gene content is for each member of each group (Fig. 1 D right). We can also see how the identity and coverage of each *pduA* gene from each member compares to our query gene **(Fig. 1 D left)**. 

Other groups with different gene architectures from the Pdu gene cluster, but containing a *pduA* gene homolog are also present **(Fig. 1D groups 1-3)**. 
>>>>>>> Stashed changes

**Conclusions**

The group 0 boxplot of region dissimilarities **(Fig. 1 D right)** indicate 99% of genomes had a region that almost exactly matched the Pdu gene cluster, with gene content differences ranging between 0-0.25, with a median of 0. We can conclude that the Pdu gene cluster is a core component of most *S. enterica* genomes, and that some of these Pdu gene clusters are undergoing gene/gain loss. [We explore the implications of these results in our publication pre-print](https://doi.org/10.1101/2021.05.27.446007). 


# Installation

Use pip to install:

```
pip install genegrouper
```

See the section 'Creating a conda environment with all dependencies' for more detailed installation instructions.

## Requirements and dependencies

[Python >= 3.6](https://www.python.org/)
[biopython](https://biopython.org/wiki/Packages)
[scikit-learn](https://scikit-learn.org/stable/)
[pandas](https://pandas.pydata.org/docs/index.html)
[matplotlib](https://matplotlib.org/)
[MMseqs2]( https://github.com/soedinglab/MMseqs2)
[MCL](https://github.com/JohannesBuchner/mcl)
[BLAST]( https://anaconda.org/bioconda/blast)
[R >= 4.0.0 (for visualizations)](https://www.r-project.org/)
[gggenes](https://cran.r-project.org/web/packages/gggenes/vignettes/introduction-to-gggenes.html)
[reshape](https://cran.r-project.org/web/packages/reshape/index.html)
[ggplot2](https://cran.r-project.org/web/packages/ggplot2/index.html)
[cowplot](https://cran.r-project.org/web/packages/cowplot/vignettes/introduction.html)
[dplyr](https://cran.r-project.org/web/packages/dplyr/vignettes/dplyr.html)
[groupdata2](https://cran.r-project.org/web/packages/groupdata2/vignettes/introduction_to_groupdata2.html)

# Usage
GeneGrouper has two required inputs:
1. A translated gene sequence in fasta format (with file extension .fasta/.txt)
2. A folder containing RefSeq GenBank-format genomes (with the file extension .gbff). 

#### Use `build_database` to make a GeneGrouper database of your RefSeq .gbff genomes
```
GeneGrouper -g /path/to/gbff -d /path/to/main_directory \
build_database
```
This database can be searched using `find_regions` as many times as you want!

#### Use `find_regions` to search for gene regions and output to a search-specific directory, 'gene_name'
```
GeneGrouper -d /path/to/main_directory -n gene_name \
find_regions -f /path/to/region_search_directory.fasta 
```
#### Use `visualize` to output visualizations of gene regions and their distribution among genomes and taxa
```
GeneGrouper -d /path/to/main_directory -n gene_name \
visualize --visual_type main
```

## Examples of different searches using `visualize`

Search for regions to 2,000 bp upstream and 18,000 bp downstream of the seed gene.
```
GeneGrouper -d /path/to/main_directory -n gene_name \
find_regions -f /path/to/region_search_directory.fasta -us 2000 -ds 18000
```
Restrict regions to those containing a seed gene with >=70% identity and >=90% coverage to the query gene.
```
GeneGrouper -d /path/to/main_directory -n gene_name \
find_regions -f /path/to/region_search_directory.fasta -i 70 -c 90
```

Allow for up to one region extracted per genome.
```
GeneGrouper -d /path/to/main_directory -n gene_name \
find_regions -f /path/to/region_search_directory.fasta -hk 1
```
Repeat the region grouping procedure 2 times.
```
GeneGrouper -d /path/to/main_directory -n gene_name \
find_regions -f /path/to/region_search_directory.fasta -re 2
```
Do all the above in one search
```
GeneGrouper -d /path/to/main_directory -n gene_name \
find_regions -f /path/to/region_search_directory.fasta -us 2000 -ds 18000 -i 70 -c 90 -hk 1 -re 2 
```

## Examples of different group visualizations using `visualize`

Visualize all members within group label 0
```
GeneGrouper -d /path/to/main_directory -n gene_name \
visualize --visual_type group --group_label 0
```
# Example dataset and searches

In this section we use GeneGrouper to search for three different genes and the gene clusters they represent.

#### Change directory to an empty folder and download the following data. Ungzip the .gbff files.
```
svn checkout https://github.com/agmcfarland/GeneGrouper/trunk/test_data/genomes
svn checkout https://github.com/agmcfarland/GeneGrouper/trunk/test_data/query_genes
gunzip ./genomes/*.gz
```

There will be 15 genomes total composed of four *Salmonella*, three *Klebsiella*, four *Citrobacter*, and four *Pseudomonas* genomes.

#### Now build the database of the genomes. 

```
GeneGrouper \
-g genomes -d example_search \
build_database
```


#### We will start for the Pdu gene cluster using the *pduA* gene in all the genomes 

A 2,000 upstream and 18,000 downstream search region will be used. The blast hit threshold will be set to 30% identity and 80% coverage relative to our *pduA* query gene.

```
# start search
GeneGrouper \
-n pdua -d example_search -g genomes -t 8 \
find_regions \
-f query_genes/pdua.txt \
-us 2000 \
-ds 18000 \
-i 30 \
-c 80
```
#### Visualize the groups that were obtained by the search.

```
# make main visualizations
GeneGrouper \
-n pdua -d example_search visualize --visual_type main

# view visualizations
open ./example_search/pdua/visualizations/*
```
Hopefully you see two distinct groups: the Pdu gene cluster (g0) and the Eut gene cluster (g1). Note that 11 genomes have the Pdu gene cluster and 4 do not. The four genomes missing a Pdu gene cluster all belong to *Pseudomonas aeruginosa*. The right hand panel of the three-part visualization shows that there is some slight variation in gene content within the members of group 0.


#### Inspect group 0 to see what kind of variation in gene content exists within the group.

```
# inspect group label 0
GeneGrouper \
-n pdua -d example_search \
visualize --visual_type group --group_label 0

# view group 0 visualiation
open ./example_search/pdua/visualizations/inspect_group_0_1.png
```

The inspect group visualization shows how many unique gene architectures are present in the group, and number of members that have that architecture. As you can see, subgroup 0 has seven members with identical architecture. Interestingly, subgroup 3 is missing the *pocR* regulator. There is instead a transposase insertion! I wonder what effects that may have on regulation of the Pdu gene cluster in this genome?

#### Now we'll explore how the five-gene Pst gene cluster, which is involved in phosphate transport, is distributed in our genomes.

We already have the database built. We just need to re-run GeneGrouper using a new query gene and whatever parameters we want.

```
# start search
GeneGrouper \
-n psts_ecoli -d example_search -g genomes -t 8 \
find_regions \
-f query_genes/psts_ecoli.txt \
-us 8000 \
-ds 8000 \
-i 10 \
-c 50 

# make main visualizations
GeneGrouper \
-n psts_ecoli -d example_search \
visualize --visual_type main

# view visualizations
open ./example_search/psts_ecoli/visualizations/*
```

Our visualizations show that four groups are created. All taxa except *Pseudomonas aeruginsa* have *stSCAB/phoU*. However, genes surrounding the Pst gene cluster differe depending on taxa. This demonstrates the selective benefit of having an intact Pst gene cluster.

#### What about genes that do not have a stable gene content surrounding them? Let's search for a horizontally-transferred carbanem-resistant beta-lactamase. 

Horizontally-transferred genes are typically in regions of the genome that have highly variable gene content, or linked to large segments of mobile DNA. We will set the minimum size of a group to just two gene regions to account for how dissimilar we expect gene regions to be.

```
# start search
GeneGrouper \
-n ndmbla -d example_search -g genomes -t 8 \
find_regions \
-f query_genes/c_bla.txt \
-us 8000 \
-ds 8000 \
-i 10 \
-c 50 \
--min_group_size 2

# make main visualizions
GeneGrouper \
-n ndmbla -d example_search \
visualize --visual_type main

# visualize group -1, which has 9 regions that could not be placed into their own group(s)
GeneGrouper \
-n ndmbla -d example_search \
visualize --visual_type group --group_label -1

# View visualizations
open ./example_search/ndmbla/visualizations/*
```
From the outputs, only two distinct groups show up. A third group, g-1, contains all gene regions that could not be placed into their own group. From looking at the gene content of all groups, we can see that there are many mobile elements, such as transposases and type-IV conjugation secretion system genes, as we expected. This supports the evidence that these beta-lactamases we've found are mobile. 

Interestingly, of the two distinct groups, we observe a taxa dependent distribution. The mobile genetic elements these beta-lactamases are associated with may be adapted to transfer between members of the same species.


# Search output file structure and file descriptions

Each GeneGrouper region search produces the following output:

```
├── main_directory
│   ├── region_search_directory
│   │   ├── group_statistics_summmary.csv
│   │   ├── representative_group_member_summary.csv
│   │   ├── group_taxa_summary.csv
│   │   ├── representative_group_member_summary.csv
│   │   ├── group_regions.csv
│   │   ├── group_region_seqs.faa
│   │   ├── visualizations
│   │   │   ├── group_summary.png
│   │   │   ├── groups_by_taxa.png
│   │   │   ├── taxa_searched.png
│   │   │   ├── inspect_group_-1.png
│   │   ├── internal_data
│   │   ├── seed_results.db
````
### Tabular files

#### group_statistics_summmary.csv

* Each row shows has the group id, number of members, and the min, mean, max, and std dev for a member's pair-wise dissimilarity to other members in the group, and the identity and coverage of a member's seed gene relative to the query gene.

#### representative_group_member_summary

* Each row shows the group id and each gene found in the reprsentative group member.

#### group_taxa_summary.csv

* Each row shows the percentage of genomes from each taxa that have a member in a group. An asterisk indicates that a genome had more than one member in that group.

#### group_regions.csv

**This is a large table that has all data generated from the run. Each row is a gene that was found in a genome and has information about:**

* Which genome it originates from, what contig it is on, and which group member it belongs to.

* Whether it is a pseudogene or not.

* The RefSeq locus tag, RefSeq gene name, and RefSeq product annotation.

* The GeneGrouper group it belongs to.

* The relative and pair-wise dissimilarity of the group member it belongs to.

* The position of the gene in the gene region.

* The start, end, and strand orientation of the gene.

### Sequence files

#### group_region_seqs.faa

* Contains all amino acid translated seqeunces for all genes from all gene regions. The header of each translated gene matches the orf_id from group_regions.csv


### Visual files

**Note: Some visual files will have a numbered suffix at the end. This is because 30 groups at a time are displayed in each visualization. So a search that produces 60 groups will have two of each visualization.**

#### group_summary.png

* Three-part visualization of the group sarch output. Left panel shows the identity (red) and coverage (blue) of each member's seed gene used to find the region relative to the query gene. Middle panel shows the gene architecture for each group's representative member. X's indicate pseudogenes. Numbers above genes are orthology identifiers. Right panel shows the relative gene content dissimilarity (blue) of each member to the representative member in the middle panel, and the mean pair-wise dissimilarity of each member to all other members in the group.

#### groups_by_taxa.png

* Heatmap showing the percentage of genomes in a taxon with at least one member per group. Asterirks indicates that at least one genome in that taxon has more than one member in that group.

#### taxa_searched.png

* The total number of genomes in a taxon searched (blue) and the number of genomes in that taxon that had at least one region extracted (red).

#### inspect_group_<group id>.png

* Three-part visualization of a group inspection output. The left panel shows the counts of each unique subgroup architecture. Middle panel shows the subgroup architecture. Right panel shows the dissimilarity of the subgroup gene content relative to subgroup 0, which is also the group representative.

# Creating a conda environment with all dependencies

Instructions for creating a self-contained conda environment for GeneGrouper with all required dependencies.


```
conda create -n GeneGrouper_env python=3.9
source activate GeneGrouper_env
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
pip install biopython scipy scikit-learn pandas matplotlib
conda install -c conda-forge -c bioconda 
conda install -c bioconda mcl blast
conda install -c bioconda blast
```


If you do not have R in your path (i.e. the command ```which R ``` does not print a /path/to/R), you can install R using conda:

```
conda install -c conda-forge r-base
```

If you already have R installed, or after installing R, install the following packages from the CRAN repository:

(This might take a while if you have a fresh installation of R!)

```
r
packages <- c("reshape", "ggplot2", "cowplot", "dplyr", "gggenes", "groupdata2")
install.packages(setdiff(packages, rownames(installed.packages()))) 
q()
```

Finally,

```
pip install GeneGrouper
```

<<<<<<< Updated upstream
 **Dependencies:**

```
conda list

# Name                    Version                   Build  Channel
_libgcc_mutex             0.1                 conda_forge    conda-forge
_openmp_mutex             4.5                       1_gnu    conda-forge
_r-mutex                  1.0.1               anacondar_1    conda-forge
binutils_impl_linux-64    2.35.1               h193b22a_2    conda-forge
binutils_linux-64         2.35                h67ddf6f_30    conda-forge
biopython                 1.78                     pypi_0    pypi
blast                     2.5.0                hc0b0e79_3    bioconda
boost                     1.76.0           py39h5472131_0    conda-forge
boost-cpp                 1.76.0               h312852a_1    conda-forge
bwidget                   1.9.14               ha770c72_0    conda-forge
bzip2                     1.0.8                h7f98852_4    conda-forge
c-ares                    1.17.1               h7f98852_1    conda-forge
ca-certificates           2020.12.5            ha878542_0    conda-forge
cairo                     1.16.0            h6cf1ce9_1008    conda-forge
certifi                   2020.12.5        py39hf3d152e_1    conda-forge
curl                      7.76.1               hea6ffbf_2    conda-forge
cycler                    0.10.0                   pypi_0    pypi
font-ttf-dejavu-sans-mono 2.37                 hab24e00_0    conda-forge
font-ttf-inconsolata      3.000                h77eed37_0    conda-forge
font-ttf-source-code-pro  2.038                h77eed37_0    conda-forge
font-ttf-ubuntu           0.83                 hab24e00_0    conda-forge
fontconfig                2.13.1            hba837de_1005    conda-forge
fonts-conda-ecosystem     1                             0    conda-forge
fonts-conda-forge         1                             0    conda-forge
freetype                  2.10.4               h0708190_1    conda-forge
fribidi                   1.0.10               h36c2ea0_0    conda-forge
gawk                      5.1.0                h7f98852_0    conda-forge
gcc_impl_linux-64         9.3.0               h70c0ae5_19    conda-forge
gcc_linux-64              9.3.0               hf25ea35_30    conda-forge
genegrouper               0.0.1                    pypi_0    pypi
gettext                   0.19.8.1          h0b5b191_1005    conda-forge
gfortran_impl_linux-64    9.3.0               hc4a2995_19    conda-forge
gfortran_linux-64         9.3.0               hdc58fab_30    conda-forge
graphite2                 1.3.13            h58526e2_1001    conda-forge
gsl                       2.6                  he838d99_2    conda-forge
gxx_impl_linux-64         9.3.0               hd87eabc_19    conda-forge
gxx_linux-64              9.3.0               h3fbe746_30    conda-forge
harfbuzz                  2.8.1                h83ec7ef_0    conda-forge
icu                       68.1                 h58526e2_0    conda-forge
jbig                      2.1               h7f98852_2003    conda-forge
joblib                    1.0.1                    pypi_0    pypi
jpeg                      9d                   h36c2ea0_0    conda-forge
kernel-headers_linux-64   2.6.32              h77966d4_13    conda-forge
kiwisolver                1.3.1                    pypi_0    pypi
krb5                      1.19.1               hcc1bbae_0    conda-forge
ld_impl_linux-64          2.35.1               hea4e1c9_2    conda-forge
lerc                      2.2.1                h9c3ff4c_0    conda-forge
libblas                   3.9.0                9_openblas    conda-forge
libcblas                  3.9.0                9_openblas    conda-forge
libcurl                   7.76.1               h2574ce0_2    conda-forge
libdeflate                1.7                  h7f98852_5    conda-forge
libedit                   3.1.20191231         he28a2e2_2    conda-forge
libev                     4.33                 h516909a_1    conda-forge
libffi                    3.3                  h58526e2_2    conda-forge
libgcc-devel_linux-64     9.3.0               h7864c58_19    conda-forge
libgcc-ng                 9.3.0               h2828fa1_19    conda-forge
libgfortran-ng            9.3.0               hff62375_19    conda-forge
libgfortran5              9.3.0               hff62375_19    conda-forge
libglib                   2.68.2               h3e27bee_0    conda-forge
libgomp                   9.3.0               h2828fa1_19    conda-forge
libiconv                  1.16                 h516909a_0    conda-forge
libidn2                   2.3.1                h7f98852_0    conda-forge
liblapack                 3.9.0                9_openblas    conda-forge
libnghttp2                1.43.0               h812cca2_0    conda-forge
libopenblas               0.3.15          pthreads_h8fe5266_1    conda-forge
libpng                    1.6.37               h21135ba_2    conda-forge
libssh2                   1.9.0                ha56f1ee_6    conda-forge
libstdcxx-devel_linux-64  9.3.0               hb016644_19    conda-forge
libstdcxx-ng              9.3.0               h6de172a_19    conda-forge
libtiff                   4.3.0                hf544144_1    conda-forge
libunistring              0.9.10               h14c3975_0    conda-forge
libuuid                   2.32.1            h7f98852_1000    conda-forge
libwebp-base              1.2.0                h7f98852_2    conda-forge
libxcb                    1.13              h7f98852_1003    conda-forge
libxml2                   2.9.12               h72842e0_0    conda-forge
lz4-c                     1.9.3                h9c3ff4c_0    conda-forge
make                      4.3                  hd18ef5c_1    conda-forge
matplotlib                3.4.2                    pypi_0    pypi
mcl                       14.137          pl5262h779adbc_6    bioconda
mmseqs2                   13.45111             h95f258a_1    bioconda
ncurses                   6.2                  h58526e2_4    conda-forge
numpy                     1.20.3           py39hdbf815f_0    conda-forge
openssl                   1.1.1k               h7f98852_0    conda-forge
pandas                    1.2.4                    pypi_0    pypi
pango                     1.48.5               hb8ff022_0    conda-forge
pcre                      8.44                 he1b5a44_0    conda-forge
pcre2                     10.36                h032f7d1_1    conda-forge
perl                      5.26.2            h36c2ea0_1008    conda-forge
pillow                    8.2.0                    pypi_0    pypi
pip                       21.1.2             pyhd8ed1ab_0    conda-forge
pixman                    0.40.0               h36c2ea0_0    conda-forge
pthread-stubs             0.4               h36c2ea0_1001    conda-forge
pyparsing                 2.4.7                    pypi_0    pypi
python                    3.9.4           hffdb5ce_0_cpython    conda-forge
python-dateutil           2.8.1                    pypi_0    pypi
python_abi                3.9                      1_cp39    conda-forge
pytz                      2021.1                   pypi_0    pypi
r-base                    4.1.0                h9e01966_1    conda-forge
readline                  8.1                  h46c0cb4_0    conda-forge
scikit-learn              0.24.2                   pypi_0    pypi
scipy                     1.6.3                    pypi_0    pypi
sed                       4.8                  he412f7d_0    conda-forge
setuptools                49.6.0           py39hf3d152e_3    conda-forge
six                       1.16.0                   pypi_0    pypi
sqlite                    3.35.5               h74cdb3f_0    conda-forge
sysroot_linux-64          2.12                h77966d4_13    conda-forge
threadpoolctl             2.1.0                    pypi_0    pypi
tk                        8.6.10               h21135ba_1    conda-forge
tktable                   2.10                 hb7b940f_3    conda-forge
tzdata                    2021a                he74cb21_0    conda-forge
wget                      1.20.1               h22169c7_0    conda-forge
wheel                     0.36.2             pyhd3deb0d_0    conda-forge
xorg-kbproto              1.0.7             h7f98852_1002    conda-forge
xorg-libice               1.0.10               h7f98852_0    conda-forge
xorg-libsm                1.2.3             hd9c2040_1000    conda-forge
xorg-libx11               1.7.1                h7f98852_0    conda-forge
xorg-libxau               1.0.9                h7f98852_0    conda-forge
xorg-libxdmcp             1.1.3                h7f98852_0    conda-forge
xorg-libxext              1.3.4                h7f98852_1    conda-forge
xorg-libxrender           0.9.10            h7f98852_1003    conda-forge
xorg-libxt                1.2.1                h7f98852_2    conda-forge
xorg-renderproto          0.11.1            h7f98852_1002    conda-forge
xorg-xextproto            7.3.0             h7f98852_1002    conda-forge
xorg-xproto               7.0.31            h7f98852_1007    conda-forge
xz                        5.2.5                h516909a_1    conda-forge
zlib                      1.2.11            h516909a_1010    conda-forge
zstd                      1.5.0                ha95c52a_0    conda-forge

 ```




## Citation

Please cite:

Density-based binning of gene clusters to infer function or evolutionary history using GeneGrouper
=======
# FAQ

#### 1. Where can I download GenBank-format RefSeq genomes with file extension .gbff?

There are acouple of simple options I use. Let's try them on downloading all clostridium genomes with complete- or chromosome-level assemblies.

**The first option** is to use [NCBI's Assembly Advanced Search Builder](https://www.ncbi.nlm.nih.gov/assembly/advanced)

* Using the search builder, select 'Organism' and input 'clostridium'. Next, select 'Assembly Level' and select both 'chromosome' and 'complete'. You should get the following code generated automatically:

`("clostridium"[Organism]) AND (("chromosome"[Assembly Level] OR "complete genome"[Assembly Level]))`

* Press 'Search' 

* Next click on the 'Download Assemblies' button.

* Make sure the source database says 'RefSeq' and file type is 'Genomic GenBank Format (.gbff)'

**The second option** is to use [ncbi-refseq-download](https://github.com/kblin/ncbi-genome-download])

Use ```pip install ncbi-refseq-download``` to download the package.

* In the command line copy the following code:

```ncbi-genome-download --section refseq --formats genbank --assembly-levels complete,chromosome --genera Clostridium bacteria```

#### 2. Can I build a database using gzipped .gbff genomes?

No, please extract all uncompressed .gbff genomes to a single folder.

#### 3. How do I choose the correct upstream/downstream coordinates ?

GeneGrouper extracts a default 10,000 bp upstream/downstream of the seed gene. However, depending on the types of adjacent genes you are trying to capture, adjust the upstream and downstream settings. A rule of thumb is that 1,000 bp is equal to 1 gene length. You can also add additional distance to accommodate for potential insertions/deletions. If you are unsure of what you expect to see, use the default values first!

Try different settings and inspect how regions have grouped. It should be fast!

#### 4. How do I choose the correct identity and coverage cutoff values for my seed genes?
This depends on your research question. If you want more distant homologs, then lower both. Generally, identity values lower than 30% suggest very low homology. More distant homologs will likely return gene regions with gene content very different from that closely related homologs to the query gene. But you never know what you might find!

#### 5. What if I want to add more .gbff genomes to my genomes database?

Go to your genomes folder and add the new .gbff genomes. Afterwards run the following command:

`GeneGrouper -g /path/to/gbff -d /path/to/output_directory \
build_database`

This will update the database and all necessary files for GeneGrouper to use.

#### 6. How many genomes can I run GeneGrouper on and how much time should it take?

GeneGrouper has been tested on 1,130 genomes using a i7 quadcore MacBook Pro. It took 6 minutes to build the database. A search takes an average of 3 minutes. 

#### 7. What operating systems does GeneGrouper work on?

GeneGrouper has been tested on Mac OS and Linux but not Windows.

#### 8. Does GeneGrouper work with non-RefSeq GenBank format genomes?

GeneGrouper has only been tested on RefSeq GenBank format (.gbff) genomes.

#### 9. Does GeneGrouper work on Eukaryotic, Archaeal, or viral genomes?
GeneGrouper has only been tested on bacterial genomes. It likely does not work at all on Eukaryotes.

# Citation

**Density-based binning of gene clusters to infer function or evolutionary history using GeneGrouper**
>>>>>>> Stashed changes

Alexander G McFarland, Nolan W Kennedy, Carolyn E Mills, Danielle Tullman-Ercek, Curtis Huttenhower, Erica M Hartmann

bioRxiv 2021.05.27.446007; doi: https://doi.org/10.1101/2021.05.27.446007

<<<<<<< Updated upstream
## Contact

Contact me at alexandermcfarland2022@u.northwestern.edu
=======
# Contact

Feel free to message me at alexandermcfarland2022@u.northwestern.edu

Follow me on twitter [@alexmcfarland_](https://twitter.com/alexmcfarland_)! 


# References


1. Cock PJA, Antao T, Chang JT, Chapman BA, Cox CJ, Dalke A, et al. Biopython: freely available Python tools for computational molecular biology and bioinformatics. Bioinformatics. 2009 Jun 1;25(11):1422–3.

2. Buitinck L, Louppe G, Blondel M, Pedregosa F, Mueller A, Grisel O, et al. API design for machine learning software: experiences from the scikit-learn project. arXiv:13090238 [cs] [Internet]. 2013 Sep 1 [cited 2021 Apr 5]; Available from: http://arxiv.org/abs/1309.0238

3. McKerns MM, Strand L, Sullivan T, Fang A, Aivazis MAG. Building a Framework for Predictive Science. arXiv:12021056 [cs] [Internet]. 2012 Feb 6 [cited 2021 Apr 5]; Available from: http://arxiv.org/abs/1202.1056

4. Camacho C, Coulouris G, Avagyan V, Ma N, Papadopoulos J, Bealer K, et al. BLAST+: architecture and applications. BMC Bioinformatics. 2009 Dec 15;10(1):421.

5. Steinegger M, Söding J. MMseqs2 enables sensitive protein sequence searching for the analysis of massive data sets. Nature Biotechnology. 2017 Nov;35(11):1026–8.

6. Enright AJ, Van Dongen S, Ouzounis CA. An efficient algorithm for large-scale detection of protein families. Nucleic Acids Res. 2002 Apr 1;30(7):1575–84.

7. gggenes @ METACRAN [Internet]. [cited 2021 Apr 5]. Available from: https://www.r-pkg.org/pkg/gggenes

8. O’Leary NA, Wright MW, Brister JR, Ciufo S, Haddad D, McVeigh R, et al. Reference sequence (RefSeq) database at NCBI: current status, taxonomic expansion, and functional annotation. Nucleic Acids Res. 2016 Jan 4;44(Database issue):D733–45.

9. Steinegger M, Söding J. Clustering huge protein sequence sets in linear time. Nature Communications. 2018 Jun 29;9(1):2542.

10. Ester M, Kriegel H-P, Xu X. A Density-Based Algorithm for Discovering Clusters in Large Spatial Databases with Noise. :6.
>>>>>>> Stashed changes

11. Caliński T, Harabasz J. A dendrite method for cluster analysis. Communications in Statistics. 1974 Jan 1;3(1):1–27.