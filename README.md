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

Necessary assembly software download
