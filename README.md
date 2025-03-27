## Welcome
Welcome to our Arthropod page! please follow below steps to install and configure your environment variables for command line or web version. Please make sure that command line works well prior to the installation of web version.

### Requirements
#### R
- ##### R 4.1 or later (https://www.r-project.org/)
R is utilized for visualization in Arthropod software. Please install R first and make sure R and Rscript are under your environment variables.
- ##### R packages
Several R packages are needed including ggplot2, ggtree and treeio packages. Follow the installation step, or you can install the packages by yourself.
```R
install.packages("ggtree")
install.packages("ggplot2")
install.packages("treeio")
```
#### perl Modules
Several perl Modules are needed including bioperl and Log::Log4perl. 
Follow the installation step to install by different methods:
- ##### Ubuntu
```bash
sudo apt-get install liblog-log4perl-perl
sudo apt install bioperl
```
- ##### CentOS
```bash
sudo yum install perl-Log-Log4perl
```
- ##### conda
```bash
conda install bioconda::perl-log-log4perl
conda install bioconda::perl-bioperl
```
#### Necessary assembly software download
- easy353	https://github.com/plant720/Easy353
```bash
git clone https://github.com/plant720/Easy353.git
chmod +x Easy353/build_database.py
chmod +X Easy353/easy353.py
echo "export PATH=/your_path/Easy353:\$PATH" >> ~/.bashrc
source ~/.bashrc
pip install biopython psutil requests beautifulsoup4 -i https://pypi.tuna.tsinghua.edu.cn/simple
```
- captus	https://github.com/edgardomortiz/Captus
```bash
conda install bioconda::captus
```
- fastp	Quality control and data-filtering of FASTQ files
```bash
conda install bioconda::fastp
```
Necessary phylogenetic tree building software download
- mafft and muscle	Align the homologous gene sequences
```bash
conda install muscle mafft -c bioconda
```
- trimal	Trim the alignment files
```bash
conda install bioconda::trimal
```
- iqtree and raxml	Construction of phylogenetic trees for single-gene species
```bash
conda install raxml iqtree -c bioconda
```
- astral	Use coalescent-model to construct species phylogenetic tree
```bash
conda install conda-forge::astral
```

#### Installation Arthropod procedures
- Download the low-copy gene set for for the phylum Arthropod from Figshare Database (10.6084/m9.figshare.27644622)
- We could obtain the software in the Arthropod website and uncompress the Arthropod software package
```bash
wget https://github.com/Wangjiajia-1111/Arthropod-a-tool-for-phylogenomic-research-in-arthropods/blob/main/Arthropod-v1.tar
tar -xvf Arthropod-v1.tar
```
- Add lib/ to LD_LIBRARY_PATH. To do this, add the following text to ~/.bashrc:
```text
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/your_software_path:
export PERL5LIB=$PERL5LIB:/home/your_software_path:
```
- and run:
```bash
source ~/.bashrc
```
- Test if Arthropod software is installed successfully: Arthropod. If you see the following content, congratulations! Arthropod is successfully installed. If not, see if all the requirements are satisfied or contact the authors for help.
```bash
Usage: Arthropod <command> ...

Avalable commands:
- assemble    	Assemble transcriptome without reference genome and protein prediction
- build_tree	Multiple sequence alignment,sequence trim,build gene trees and the species tree
```

## Usage
We use the arthropod low-copy gene library generated in this study and the raw data of Macropsis flavida as examples
```bash
# low_copy_OG.tar (https://figshare.com/ndownloader/files/50343768)
tar -xvf low_copy_OG.tar
# raw data of Macropsis flavida
wget https://download.cncb.ac.cn/gsa2/CRA020542/CRR1382787/CRR1382787_r1.fastq.gz
wget https://download.cncb.ac.cn/gsa2/CRA020542/CRR1382787/CRR1382787_r2.fastq.gz
```
### Step 1: Assemble new sequences using the arthropod low-copy gene library as reference fasta
```bash
perl Arthropod assemble

Usage: Arthropod assemble [commands] ...

Commands:
- fastp		        Quality control and data-filtering of FASTQ files.
- easy353		Based on reference sequences filter and assemble the transcriptome, the whole genome or genome skimming sequencing data to recover - target genes.
- captus		Build phylogenomic datasets from multiple types of sequencing data based on reference sequences.
- transSeq		Convert coding sequence into amino acid sequence.
```
- fastp	Quality control and data-filtering
```bash
perl Arthropod assemble fastp -t < thread_num> <fq1> <fq2> <output_dir> <prefix:latin name>

# example
perl Arthropod assemble fastp -t 20 CRR1382787_r1.fastq.gz CRR1382787_r2.fastq.gz Macropsis_flavida_fastp Macropsis_flavida
# The data generated after quality control are Macropsis_flavida_fastp/Macropsis_flavida_R1.fq.gz and Macropsis_flavida_fastp/Macropsis_flavida_R2.fq.gz
```

- Assemble	captus or easy353

captus (Recommend)
```bash
perl Arthropod assemble captus -t <assemble threads> -T <extract threads> -c < Multi-threads> <FASTQ files directory or list> <latin name> <reference fasta directory>

# example
perl Arthropod assemble captus -t 60 -T 30 -c 4 Macropsis_flavida_fastp Macropsis_flavida low_copy_OG
# The OG sequences assembled using captus based on the low copy gene library "low_copy_OG" are stored in the Macropsis_flavida_CAPTUSmerge directory.
```

easy353
```bash
perl Arthropod assemble easy353 -t <filtering threads> -T <assembly threads> <fq1> <fq2> <reference fasta directory> <latin name>

# example
perl Arthropod assemble easy353 -t 4 -T 8  Macropsis_flavida_fastp/Macropsis_flavida_R1.fq.gz Macropsis_flavida_fastp/Macropsis_flavida_R2.fq.gz low_copy_OG Macropsis_flavida
# The OG sequences assembled using easy353 based on the low copy gene library "low_copy_OG" are stored in the Macropsis_flavida_Easy353merge directory.
```

- - Convert nucleotide sequence to amino acid sequence
```bash
perl Arthropod assemble transSeq <latin name_CAPTUSmerge>

# example
perl Arthropod assemble transSeq Macropsis_flavida_CAPTUSmerge
# The amino acid sequence of Convert is generated in the Macropsis_flavida_CAPTUSmerge_pep.
```
```bash
perl Arthropod assemble transSeq <latin name_Easy353merge>

# example
perl Arthropod assemble transSeq Macropsis_flavida_Easy353merge
# The amino acid sequence of Convert is generated in the Macropsis_flavida_Easy353Smerge_pep.
```

### Step 2: Constructing phylogenetic tree
```bash
Usage: Arthropod build_tree [commands] ...

Commands:

- - alignment_muscle		Align the homologous gene sequences with muscle.
- - alignment_mafft		Align the homologous gene sequences with mafft.
- - Trim			Trim the alignment files with trimal.
- - RAxMLtree			Use multi-species coalescent-model to build the phylogenetic trees with RAxML(gene tree) and ASTRAL(species tree).
- - iqtree		        Use multi-species coalescent-model to build the phylogenetic trees with iqtree(gene tree) and ASTRAL(species tree).
- - tree_plot		        Visualization of the phylogeny tree.
```

- Sequence alignment	muscle or mafft

muscle
```bash
perl Arthropod build_tree alignment_muscle -t <Multi-threads> -o <outdir> <OG fasta directory>
# example
perl Arthropod build_tree alignment_muscle -t 4 -o Macropsis_flavida_CAPTUStree/1_alignment_muscle_out Macropsis_flavida_CAPTUSmerge_pep
```


mafft (recommend)
```bash
perl Arthropod build_tree alignment_mafft -t <alignment threads> -c <Multi-threads> -o <outdir> <OG fasta directory>

# example
perl Arthropod build_tree alignment_mafft -t 4 -c 30 -o Macropsis_flavida_CAPTUStree/1_alignment_mafft_out Macropsis_flavida_CAPTUSmerge_pep
```

- Trim the alignment fasta
```bash
perl Arthropod build_tree trim -t <Multi-threads> -o <outdir:2_trim_out> <alignment file directory>

# example
perl Arthropod build_tree trim -t 30 -o Macropsis_flavida_CAPTUStree/2_trim_out Macropsis_flavida_CAPTUStree/1_alignment_mafft_out
```

- constructe tree

iqtree (Recommend)
```bash
perl Arthropod build_tree iqtree -t <Multi-threads> -T <iqtree_threads> -m <model> -B 1000 -o <outdir> <trimal file directory>

# example
perl Arthropod build_tree iqtree -t 15 -T 4 -m MF -B 1000 -o Macropsis_flavida_CAPTUStree/3_iqtree_out Macropsis_flavida_CAPTUStree/2_trim_out
# The final species tree is Macropsis_flavida_CAPTUStree/3_iqtree_out/astral_out/astral_coalescent.result
```

raxml
```bash
perl Arthropod build_tree RAxMLtree -t <Multi-threads> -T <raxml_threads> -m <model> -N 100 -o <outdir> <trimal file directory>

# example
perl Arthropod build_tree RAxMLtree -t 15 -T 4 -m PROTGAMMAJTT -N 100 -o Macropsis_flavida_CAPTUStree/3_raxml_out Macropsis_flavida_CAPTUStree/2_trim_out
# The final species tree is Macropsis_flavida_CAPTUStree/3_raxml_out/astral_out/astral_coalescent.result
```

- Visualization of tree
```bash
perl Arthropod build_tree tree_plot <tree file> <species group> <output prefix>

# example
Rscript tree_plot_order.r Macropsis_flavida_CAPTUStree/3_iqtree_out/astral_out/astral_coalescent.result species_color.txt Arthropod_add-Macropsis_flavida_order
# Visualization of the tree are Arthropod_add-Macropsis_flavida_order.circular_tree.pdf and Arthropod_add-Macropsis_flavida_order.rectangular_tree.pdf
```

In addition, we can upload the final species tree file to the iTOL online website (https://itol.embl.de/upload.cgi) for beautification.

