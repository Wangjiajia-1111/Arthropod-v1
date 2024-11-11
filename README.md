#### Welcome
Welcome to our Arthropod page! please follow below steps to install and configure your environment variables for command line or web version. Please make sure that command line works well prior to the installation of web version.

#### Requirements
##### R
- ###### R 4.1 or later (https://www.r-project.org/)
R is utilized for visualization in Arthropod software. Please install R first and make sure R and Rscript are under your environment variables.
- ###### R packages
Several R packages are needed including ggplot2, ggtree and treeio packages. Follow the installation step, or you can install the packages by yourself.
install.packages("ggtree")
install.packages("ggplot2")
install.packages("treeio")

#### perl Modules
Several perl Modules are needed including bioperl and Log::Log4perl. 
Follow the installation step to install by different methods:
- ###### Ubuntu
sudo apt-get install liblog-log4perl-perl
sudo apt install bioperl
- ###### CentOS
sudo yum install perl-Log-Log4perl
- ###### conda
conda install bioconda::perl-log-log4perl
conda install bioconda::perl-bioperl

#### Necessary assembly software download
- easy353	https://github.com/plant720/Easy353
git clone https://github.com/plant720/Easy353.git
chmod +x Easy353/build_database.py
chmod +X Easy353/easy353.py
echo "export PATH=/your_path/Easy353:\$PATH" >> ~/.bashrc
source ~/.bashrc
pip install biopython psutil requests beautifulsoup4 -i https://pypi.tuna.tsinghua.edu.cn/simple
- captus	https://github.com/edgardomortiz/Captus
conda install bioconda::captus
- fastp	Quality control and data-filtering of FASTQ files
conda install bioconda::fastp
Necessary phylogenetic tree building software download
- mafft and muscle	Align the homologous gene sequences
conda install muscle mafft -c bioconda
- trimal	Trim the alignment files
conda install bioconda::trimal
- iqtree and raxml	Construction of phylogenetic trees for single-gene species
conda install raxml iqtree -c bioconda
- astral	Use coalescent-model to construct species phylogenetic tree
conda install conda-forge::astral

#### Installation Arthropod procedures
- Download the low-copy gene set for for the phylum Arthropod from Figshare Database (10.6084/m9.figshare.27644622)
- We could obtain the software in the Arthropod website and uncompress the Arthropod software package
wget https://github.com/Wangjiajia-1111/Arthropod-a-tool-for-phylogenomic-research-in-arthropods/blob/main/Arthropod-v1.tar
tar -xvf Arthropod-v1.tar
- Add lib/ to LD_LIBRARY_PATH. To do this, add the following text to ~/.bashrc:
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/your_software_path:
export PERL5LIB=$PERL5LIB:/home/your_software_path:
- and run:
source ~/.bashrc
- Test if Arthropod software is installed successfully: Arthropod. If you see the following content, congratulations! Arthropod is successfully installed. If not, see if all the requirements are satisfied or contact the authors for help.
Usage: Arthropod <command> ...
Avalable commands:
	assemble    	Assemble transcriptome without reference genome and protein prediction
	build_tree	Multiple sequence alignment,sequence trim,build gene trees and the species tree

#### Usage
##### Step 1: Assemble new sequences using the arthropod low-copy gene library as reference fasta
perl Arthropod assemble
Usage: Arthropod assemble [commands] ...
Commands:
fastp		Quality control and data-filtering of FASTQ files.
easy353		Based on reference sequences filter and assemble the transcriptome, the whole genome or genome skimming sequencing data to recover target genes.
captus		Build phylogenomic datasets from multiple types of sequencing data based on reference sequences.
transSeq		Convert coding sequence into amino acid sequence.

- fastp	Quality control and data-filtering
perl Arthropod assemble fastp -h
E.g.：perl Arthropod assemble fastp -t < thread_num> <fq1> <fq2> <output_dir> <prefix:latin name>
- Assemble	captus or easy353
- - captus (Recommend)
perl Arthropod assemble captus -h
e.g.：perl Arthropod assemble captus -t <assemble threads> -T <extract threads> -c < Multi-threads> <FASTQ files directory or list> <latin name> <reference fasta directory>
- - easy353
perl Arthropod assemble easy353 -h
E.g.：perl Arthropod assemble easy353 -t <filtering threads> -T <assembly threads> <fq1> <fq2> <reference fasta directory> <latin name>
- Convert nucleotide sequence to amino acid sequence
perl Arthropod assemble transSeq -h
perl Arthropod assemble transSeq <latin name_CAPTUSmerge>
or perl Arthropod assemble transSeq <latin name_Easy353merge>

##### Step 2: Constructing phylogenetic tree
perl Arthropod build_tree
Usage: Arthropod build_tree [commands] ...
Commands:
alignment_muscle		Align the homologous gene sequences with muscle.
alignment_mafft		Align the homologous gene sequences with mafft.
Trim					Trim the alignment files with trimal.
RAxMLtree			Use multi-species coalescent-model to build the phylogenetic trees with RAxML(gene tree) and ASTRAL(species tree).
iqtree				Use multi-species coalescent-model to build the phylogenetic trees with iqtree(gene tree) and ASTRAL(species tree).
tree_plot				Visualization of the phylogeny tree.
- Sequence alignment	muscle or mafft
- - muscle
perl Arthropod build_tree alignment_muscle -h
E.g.：perl Arthropod build_tree alignment_muscle -t <Multi-threads> -o <outdir> <OG fasta directory>
- - mafft (recommend)
perl Arthropod build_tree alignment_mafft -h
E.g.：perl Arthropod build_tree alignment_mafft -t <alignment threads> -c <Multi-threads> -o <outdir> <OG fasta directory>
- Trim the alignment fasta
perl Arthropod build_tree trim -h
E.g.：perl Arthropod build_tree trim -t <Multi-threads> -o <outdir:2_trim_out> <alignment file directory>
- buildtree
- - iqtree (Recommend)
perl Arthropod build_tree iqtree -h
E.g.：perl Arthropod build_tree iqtree -t <Multi-threads> -T <iqtree_threads> -m <model> -B 1000 -o <outdir> <trimal file directory>
- - raxml
perl Arthropod build_tree RAxMLtree -h
E.g.：perl Arthropod build_tree RAxMLtree -t <Multi-threads> -T <raxml_threads> -m <model> -N 100 -o <outdir> <trimal file directory>
- Visualization of tree
perl Arthropod build_tree tree_plot -h
E.g.：perl Arthropod build_tree tree_plot <tree file> <species group> <output prefix>
In addition, we can upload the final species tree file to the iTOL online website (https://itol.embl.de/upload.cgi) for beautification.

