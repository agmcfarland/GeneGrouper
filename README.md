# GeneGrouper


GeneGrouper is a command-line tool that finds gene clusters that contain a gene of interest in a set of genomes and bins them into groups of similar gene clusters.


<img src="docs/overview_figure.png" alt="GeneGrouper overview figure" width=1000>

<br><br>


---

<br><br>
## Quick Overview

#### Inputs

1. One translated gene of interest (.faa/.fasta/.txt)

2. Two or more genomes from RefSeq (.gbff)

#### Outputs

1. For each individual search, a new folder will be outputted containing all gene clusters and their groupings

#### Visualizations and data

1. GeneGrouper produces 4 different visualizations to understand the pan-genomic context of gene cluster bins.

2. Several datasets are outputted for further inspection by the user.

#### Performance

For 1,130 genomes and using a 2.2Ghz quad-core MacBook Pro, GeneGrouper:

1. Builds a database in 8 minutes (this step is only done once)

2. Neatly bins 2,300 separate gene clusters in \~4 minutes


## Quick Start


#### Use `build_database` to make a database of your RefSeq .gbff genomes

You only need to make a database of the genomes once.

```
GeneGrouper -g /path/to/gbff -d /path/to/output_directory \
build_database
```

#### Use `find_regions` to search for gene clusters and output to a search-specific directory, 'gene_name'

Now you can search the database of genomes for gene clusters that contain your gene of interest! 

```
GeneGrouper -d /path/to/output_directory -n gene_name \
find_regions -f /path/to/seed_gene.fasta 
```

#### Use `visualize` to output visualizations of gene clusters and their distribution among genomes and taxa

```
GeneGrouper -d /path/to/output_directory -n gene_name \
visualize -vt main
```


## Additional paramters and examples

#### `find_regions`


Search for gene clusters 2,000 bp upstream and 18,000 bp downstream of the seed gene.

```
GeneGrouper -d /path/to/output_directory -n gene_name \
find_regions -f /path/to/seed_gene.fasta -us 2000 -ds 18000
```

Restrict gene clusters to those containing a seed gene with >=70% identity and >=90% coverage to the query gene.

```
GeneGrouper -d /path/to/output_directory -n gene_name \
find_regions -f /path/to/seed_gene.fasta -i 70 -c 90
```

Allow for up to one gene cluster found per genome.

```
GeneGrouper -d /path/to/output_directory -n gene_name \
find_regions -f /path/to/seed_gene.fasta -hk 1
```

Repeat the gene cluster clustering procedure 2 times.

```
GeneGrouper -d /path/to/output_directory -n gene_name \
find_regions -f /path/to/seed_gene.fasta -re 2
```

Do it all at once.

```
GeneGrouper -d /path/to/output_directory -n gene_name \
find_regions -f /path/to/seed_gene.fasta -us 2000 -ds 18000 -i 70 -c 90 -hk 1 -re 2 
```

#### `visualize` 

Visualize all subclusters within cluster label 'c0'.

```
GeneGrouper -d /path/to/output_directory -n gene_name \
visualize -vt region_cluster -clab 0
```

<br><br>

---

<br><br>

## Tutorial

Coming soon! Example dataset and tutorial!

## Search dataset output descriptions

Coming soon!

<br><br>

---

<br><br>

## Commands

**Top-level commands**
```
GeneGrouper [-h] [-d PROJECT_DIRECTORY] [-n SEARCH_NAME]
                     [-g GENOMES_DIRECTORY] [-t THREADS]
                     {build_database,find_regions,visualize} ...

optional arguments:
  -h, --help            show this help message and exit
  -d PROJECT_DIRECTORY, --project_directory PROJECT_DIRECTORY
                        main directory to contain the base files used for
                        region searching and clustering
  -n SEARCH_NAME, --search_name SEARCH_NAME
                        name of the directory to contain search-specific
                        results
  -g GENOMES_DIRECTORY, --genomes_directory GENOMES_DIRECTORY
                        directory containing genbank-file format genomes with
                        the suffix .gbff
  -t THREADS, --threads THREADS
                        number of threads to use

subcommands:
  valid subcommands

  {build_database,find_regions,visualize}
                        sub-command help
    build_database      convert a set of genomes into a useuable format for
                        GeneGrouper
    find_regions        find regions given a translated gene and a set of
                        genomes
    visualize           visualize region clusters
```


**find_regions commands**
```
GeneGrouper find_regions 
  -h, --help            show this help message and exit
  -f SEED_FILE, --seed_file SEED_FILE
                        provide the absolute path to a fasta file containing a
                        translated gene sequence
  -us UPSTREAM_SEARCH, --upstream_search UPSTREAM_SEARCH
                        upstream search length in basepairs
  -ds DOWNSTREAM_SEARCH, --downstream_search DOWNSTREAM_SEARCH
                        downstream search length in basepairs
  -i SEED_IDENTITY, --seed_identity SEED_IDENTITY
                        identity cutoff for initial blast search
  -c SEED_COVERAGE, --seed_coverage SEED_COVERAGE
                        coverage cutoff for initial blast search
  -hk SEED_HITS_KEPT, --seed_hits_kept SEED_HITS_KEPT
                        number of blast hits to keep
  -re RECLUSTER_ITERATIONS, --recluster_iterations RECLUSTER_ITERATIONS
                        number of region re-clustering attempts after the
                        initial clustering
```

**visualize commands**
```
GeneGrouper visualize 
optional arguments:
  -h, --help            show this help message and exit
  -vt {main,region_cluster}, --visual_type {main,region_cluster}
  -clab CLUSTER_LABEL, --cluster_label CLUSTER_LABEL
```


## Installation

**Simple:**

Simple installation assuming you already have dependencies installed.

```pip install GeneGrouper```


**Detailed:**

Instructions for creating a self-contained conda environment for GeneGrouper with all required dependencies.



```
conda create -n GeneGrouper_env python=3.9

source activate GeneGrouper_env

conda config --add channels defaults

conda config --add channels bioconda

conda config --add channels conda-forge

pip install biopython scipy scikit-learn pandas matplotlib

conda install -c conda-forge -c bioconda mmseqs2

conda install -c bioconda mcl

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


2. Download GeneGrouper

```
pip install GeneGrouper
```

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

Alexander G McFarland, Nolan W Kennedy, Carolyn E Mills, Danielle Tullman-Ercek, Curtis Huttenhower, Erica M Hartmann

bioRxiv 2021.05.27.446007; doi: https://doi.org/10.1101/2021.05.27.446007

## Contact

Contact me at alexandermcfarland2022@u.northwestern.edu

