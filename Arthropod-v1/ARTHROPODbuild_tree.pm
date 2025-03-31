#!/usr/bin/perl
package build_tree;

# 检查是否存在程序
sub check_command{
    my ($cmd) = @_;
    my $path = `which $cmd`;
    $path=~s/\n//g;
    return $path;
}


sub build_tree{
    my $usage="
        Usage: Arthropod build_tree [commands] ...

        Commands:
        alignment_muscle        Align the homologous gene sequences with muscle.
        alignment_mafft         Align the homologous gene sequences with mafft.
        trim                    Trim the alignment files with trimal.
        RAxMLtree               Use multi-species coalescent-model to build the phylogenetic trees with 
                                RAxML(gene tree) and ASTRAL(species tree).
        iqtree                  Use multi-species coalescent-model to build the phylogenetic trees with 
                                iqtree(gene tree) and ASTRAL(species tree).
        tree_plot               Visualization of the phylogeny tree.
        ";

    die $usage if @ARGV<1;
    my $com=shift @ARGV;
    if ($com eq "alignment_muscle"){
        alignment_muscle(@ARGV);
    }
    elsif ($com eq "alignment_mafft"){
        alignment_mafft(@ARGV);
    }
    elsif ($com eq "trim"){
        trim(@ARGV);
    }
    elsif ($com eq "RAxMLtree"){
        RAxMLtree(@ARGV);
    }
    elsif ($com eq "iqtree"){
        iqtree(@ARGV);
    }
    elsif ($com eq "tree_plot"){
        tree_plot(@ARGV);
    }
    else {
        print STDERR "Unknown command: $com\n";
        die($usage);
    }
}


sub alignment_muscle{
    use strict;
    use warnings;
    use Log::Log4perl;
    use Getopt::Std;
    use vars qw($opt_h $opt_o $opt_t);
    getopts("ho:t:");

    Log::Log4perl->init("Build_tree.conf");
    my $logger = Log::Log4perl->get_logger();

    my $usage = "\nUsage: Arthropod build_tree alignment_muscle [options] <input_dir>
        
        Necessary input description:
        input_dir       <str>       This directory should contain one or more sra data.
                                    Files ended by .fa or .fasta or .fas or .faa

        Options:
        -h                          Print this usage page.
        -t              <int>       Multi-threading number(multi_threads.pl). default=10
        -o              <str>       Output directory. default=1_alignment_muscle_out
        ";

    die $usage if @ARGV!=1;
    die $usage if defined($opt_h);
    my ($input_dir)=@ARGV;

    # check the availability of muscle
    $logger->info("=========================== 1.Muscle alignment start ===========================");
    $logger->error("The muscle for alignment is not available!") unless check_command('muscle');
    my $muscle_path;
    if (check_command('muscle')){
        $muscle_path = check_command('muscle');
    } else {
        die "ERROR: The muscle for alignment is not available\n";
    }

    # create output directory
    my $out_dir = "1_alignment_mafft_out";
    $out_dir = $opt_o if defined($opt_o);
    system(qq(mkdir -p $out_dir)) unless -e $out_dir;
    # threads
    my $thread_num = 10;
    $thread_num = $opt_t if defined($opt_t);
    
    # 读取数据
    opendir(FASTA,$input_dir) || die("Error: can't open input data directory!\n");
    my @sample = readdir(FASTA);
    close FASTA;

    # 生成多线程shell脚本
    my $alignment_sh = "alignment_muscle.sh";
    $logger->info("Generate multi-thread muscle alignment script: $alignment_sh");
    open(SH,">$alignment_sh") || die("Error:can not write muscle alignment shell script:$alignment_sh\n");
    print STDOUT "The commands for alignment are:\n";
    foreach my $s(@sample){
        next if $s=~/^\./;
        my $input_file = $input_dir."/".$s;
        die("The file $input_file is not existed !") unless -e $input_file;
        my $OG_number = (split/\./,$s)[0];
        my $output_file = $out_dir."/".$OG_number.".alignment.fasta";
        print STDOUT "$muscle_path -align $input_file -output $output_file \n";
        print SH "$muscle_path -align $input_file -output $output_file \n";
    }
    close SH;

    # check multi_threads.pl
    my $multi_threads_path;
    if (check_command('multi_threads.pl')){
        $multi_threads_path=check_command('multi_threads.pl');
    } else {
        die "ERROR: multi_threads.pl was not found! \n";
    }
    $logger->info("The command of alignment multi-threads script:[ perl $multi_threads_path -c $thread_num $alignment_sh ]");
    system(qq(perl $multi_threads_path -c $thread_num $alignment_sh));
    $logger->info("=========================== Muscle alignment finished ===========================");
}


sub alignment_mafft{
    use strict;
    use warnings;
    use Log::Log4perl;
    use Getopt::Std;
    use vars qw($opt_h $opt_o $opt_t $opt_c);
    getopts("ho:t:c:");

    Log::Log4perl->init("Build_tree.conf");
    my $logger = Log::Log4perl->get_logger();

    my $usage = "\nUsage: Arthropod build_tree alignment_mafft [options] <input_dir>
        
        Necessary input description:
        input_dir       <str>       This directory should contain one or more sra data.
                                    Files ended by .fa or .fasta or .fas or .faa

        Options:
        -h                          Print this usage page.
        -t              <int>       Number of threads(mafft). default=4.
        -c              <int>       Multi-threads number(multi_threads.pl). default=10.
        -o              <str>       Output directory. default=1_alignment_mafft_out
        ";

    die $usage if @ARGV!=1;
    die $usage if defined($opt_h);
    my ($input_dir)=@ARGV;

    # check the availability of mafft
    $logger->info("=========================== 1.Mafft alignment start ===========================");
    $logger->error("The mafft for alignment is not available!") unless check_command('mafft');
    my $mafft_path;
    if (check_command('mafft')){
        $mafft_path = check_command('mafft');
    } else {
        die "ERROR: The mafft for alignment is not available\n";
    }

    # create output directory
    my $out_dir = "1_alignment_mafft_out";
    $out_dir = $opt_o if defined($opt_o);
    system(qq(mkdir -p $out_dir)) unless -e $out_dir;
    # threads
    my $thread_num = 4;
    $thread_num = $opt_t if defined($opt_t);

    my $multi_thread_num = 10;
    $multi_thread_num = $opt_c if defined($opt_c);
    
    # 读取数据
    opendir(FASTA,$input_dir) || die("Error: can't open input data directory!\n");
    my @sample = readdir(FASTA);
    close FASTA;

    # 生成多线程shell脚本
    my $alignment_sh = "alignment_mafft.sh";
    $logger->info("Generate multi-thread mafft alignment script: $alignment_sh");
    open(SH,">$alignment_sh") || die("Error:can not write mafft alignment shell script:$alignment_sh\n");
    print STDOUT "The commands for alignment are:\n";
    foreach my $s(@sample){
        next if $s=~/^\./;
        my $input_file = $input_dir."/".$s;
        die("The file $input_file is not existed !") unless -e $input_file;
        my $OG_number = (split/\./,$s)[0];
        my $output_file = $out_dir."/".$OG_number.".alignment.fasta";
        print STDOUT "$mafft_path --thread $thread_num --quiet $input_file > $output_file \n";
        print SH "$mafft_path --thread $thread_num --quiet $input_file > $output_file \n";
    }
    close SH;

    # check multi_threads.pl
    my $multi_threads_path;
    if (check_command('multi_threads.pl')){
        $multi_threads_path=check_command('multi_threads.pl');
    } else {
        die "ERROR: multi_threads.pl was not found! \n";
    }
    $logger->info("The command of alignment multi-threads script:[ perl $multi_threads_path -c $multi_thread_num $alignment_sh ]");
    system(qq(perl $multi_threads_path -c $multi_thread_num $alignment_sh));
    $logger->info("=========================== Alignment finished ===========================");
}


sub trim{
    use strict;
    use warnings;
    use Getopt::Std;
    use Log::Log4perl;
    use vars qw($opt_h $opt_o $opt_t);
    getopts("ho:t:");

    Log::Log4perl->init("Build_tree.conf");     # 加载配置文件
    my $logger = Log::Log4perl->get_logger(); # 获取日志对象

    my $usage = "\nUsage: Arthropod build_tree trim [options] <input_dir>
        
        Necessary input description:
        input_dir       <str>       This directory of muscle or mafft alignment result.

        Options:
        -h                          Print this usage page.
        -t              <int>       Multi-threading number: default=10
        -o              <str>       Output directory. Default=2_trim_out
        ";

    die $usage if @ARGV!=1;
    die $usage if defined($opt_h);
    my ($input_dir)=@ARGV;

    # check the availability of trimal
    $logger->info("=========================== 2.Trimming start ===========================");
    $logger->error("ERROR: The trimal for trimming is not available") unless check_command('trimal');
    my $trimal_path;
    if (check_command('trimal')){
        $trimal_path = check_command('trimal');
    } else {
        die "ERROR: The trimal for trimming is not available\n";
    }

    # create output directory
    my $out_dir = "2_trim_out";
    $out_dir = $opt_o if defined($opt_o);
    if(-e $out_dir){
        die("Error: output directory \"$out_dir\" already exists. To avoid overwriting of existing files, we kindly request that the name of output directory should be changed.\n");
    } else{
        system(qq(mkdir -p $out_dir));
    }

    # threads
    my $thread_num = 10;
    $thread_num = $opt_t if defined($opt_t);
    
    # 读取数据
    opendir(FASTA,$input_dir) || die("Error: can't open input data directory!\n");
    my @sample = readdir(FASTA);
    close FASTA;
    # 生成多线程shell脚本
    my $trim_sh = "trimming.sh";
    $logger->info("Generate multi-thread smming cript: $trim_sh");
    open(SH,">$trim_sh") || die("Error:can not write trimming shell script:$trim_sh\n");
    print STDOUT "\nThe commands of $trim_sh are:\n";
    foreach my $s(@sample){
        next if $s=~/^\./;
        my $input_file = $input_dir."/".$s;
        die("The file $input_file is not existed !") unless -e $input_file;
        my $OG_number = (split/\./,$s)[0];
        my $output_file = $out_dir."/".$OG_number.".trim.fasta";
        print STDOUT "$trimal_path -in $input_file -out $output_file -automated1\n";
        print SH "$trimal_path -in $input_file -out $output_file -automated1\n";
    }
    close SH;

    # check multi_threads.pl
    my $multi_threads_path;
    if (check_command('multi_threads.pl')){
        $multi_threads_path=check_command('multi_threads.pl');
    } else {
        die "ERROR: multi_threads.pl was not found! \n";
    }
    $logger->info("The command of trimming multi-threads script:[ perl $multi_threads_path -c $thread_num $trim_sh ]\n");
    system(qq(perl $multi_threads_path -c $thread_num $trim_sh));
    $logger->info("=========================== Trimming finished ===========================");
} 


sub RAxMLtree{
    use strict;
    use warnings;
    use Getopt::Std;
    use Log::Log4perl;
    use vars qw($opt_h $opt_t $opt_o $opt_m $opt_N $opt_T);
    getopts("ho:t:m:N:T:");

    Log::Log4perl->init("Build_tree.conf");
    my $logger = Log::Log4perl->get_logger(); 

    my $usage = "\nUsage: Arthropod build_tree RAxMLtree [options] <input_dir>
        
        Necessary input description:
        input_dir       <str>       This directory of trimal trimming result.

        Options:
        -h                          Print this usage page.
        -t              <int>       Multi-threading number. default=4
        -o              <str>       Output directory. default=3_RAxMLtree_out
        -T              <str>       Specify the raxmlHPC-PTHREADS number of threads. default=4
        -m              <str>       Set the model of nucleotides substitution. default=GTRGAMMA
                                    if amino acid,you can use PROTGAMMAJTT
        -N             <int>        Set the number of bootstrap. default=100
        ";

    die $usage if @ARGV!=1;
    die $usage if defined($opt_h);
    my ($input_dir)=@ARGV;

    $logger->info("=========================== 3. Build tree start ===========================");
    $logger->info("=========================== raxmlHPC gene tree ===========================");
    # check the availability of raxml and astral
    $logger->error("The raxmlHPC-PTHREADS for build gene tree is not available") unless check_command('raxmlHPC-PTHREADS');
    my $raxml_path;
    my $astral_path;
    if (check_command('raxmlHPC-PTHREADS')){
        $raxml_path = check_command('raxmlHPC-PTHREADS');
    } else {
        die "ERROR: The raxmlHPC-PTHREADS for build gene tree is not available\n";
    }
    $logger->error("The astral for building the species tree is not available") unless check_command('astral');
    if (check_command('astral')){
        $astral_path = check_command('astral');
    } else {
        die "ERROR:The astral for building the species tree is not available\n";
    }
    # create output directory
    my $out_dir="3_RAxMLtree_out";
    $out_dir=$opt_o if defined($opt_o);
    system(qq(mkdir -p $out_dir)) unless -e $out_dir;
    # threads
    my $thread_num = 4;
    $thread_num=$opt_t if defined($opt_t);
    my $thread_raxml = 4;
    $thread_raxml=$opt_T if defined($opt_T);
    # set model and bootstrap
    my $model = "GTRGAMMA";
    $model=$opt_m if defined($opt_m);
    my $bootstrap = 100;
    $bootstrap=$opt_N if defined($opt_N);
    # 读取数据
    opendir(FASTA,$input_dir) || die("Error: can't open input data directory!\n");
    my @sample = readdir(FASTA);
    close FASTA;
    # 生成多线程shell脚本
    my $raxmlHPC_sh = "raxmlHPC.sh";
    $logger->info("Generate multi-thread raxmlHPC shell script: $raxmlHPC_sh");
    open(SH,">$raxmlHPC_sh") || die("Error:can not write raxmlHPC shell script:$raxmlHPC_sh\n");
    print STDOUT "\nThe commands of $raxmlHPC_sh are:\n";
    foreach my $s(@sample){
        next if $s=~/^\./;
        next if $s=~/reduced/;
        my $input_file = $input_dir."/".$s;
        die("The file $input_file is not existed !") unless -e $input_file;
        my $OG_number = (split/\./,$s)[0];
        my $output_dir = $out_dir."/raxmlHPC/".$OG_number;
        system(qq(mkdir -p $output_dir));
        my $output_suffix = $OG_number.".tre";
        print STDOUT "$raxml_path -s $input_file -m $model -p 12345 -n $output_suffix -f a -x 12345 -N $bootstrap -T $thread_raxml -w $output_dir\n";
        print SH "$raxml_path -s $input_file -m $model -p 12345 -n $output_suffix -f a -x 12345 -N $bootstrap -T $thread_raxml -w $output_dir\n";
    }
    close SH;

    # check multi_threads.pl
    my $multi_threads_path;
    if (check_command('multi_threads.pl')){
        $multi_threads_path=check_command('multi_threads.pl');
    } else {
        die "ERROR: multi_threads.pl was not found! \n";
    }
    $logger->info("The command of raxmlHPC-PTHREADS multi-thread script:[ perl $multi_threads_path -c $thread_num $raxmlHPC_sh ]");
    system(qq(perl $multi_threads_path -c $thread_num $raxmlHPC_sh));
    $logger->info("=========================== raxmlHPC finished ===========================");

    my $output_gene_tree_dir = $out_dir."/gene_tree";
    system(qq(mkdir -p $output_gene_tree_dir));
    my $all_OGs_bootstrap_pathfile = $out_dir."/gene_tree/all_bootstrap.txt";
    open(BOOTSTRAP,">$all_OGs_bootstrap_pathfile") || die("Error:can not write to $all_OGs_bootstrap_pathfile\n");
    foreach my $s(@sample){
        next if $s=~/^\./;
        next if $s=~/reduced/;
        my $OG_number = (split/\./,$s)[0];
        my $tree_file = $out_dir."/raxmlHPC/".$OG_number."/RAxML_bipartitions.$OG_number.tre";
        die("The tree file $tree_file is not existed !") unless -e $tree_file;
        my $bootstrap_file = $out_dir."/raxmlHPC/".$OG_number."/RAxML_bootstrap.$OG_number.tre";
        open(RAXMLTREE,"$bootstrap_file") || die("The bootstrap file $bootstrap_file is not existed !");
        my $bootstrap_file_tmp = $out_dir."/raxmlHPC/".$OG_number."/RAxML_bootstrap.$OG_number.tre.tmp";
        open(TMP,">$bootstrap_file_tmp") || die("The bootstrap file $bootstrap_file_tmp can not open !");
        while (<RAXMLTREE>){
            my $tmp = $_;
            $tmp=~s/_OG\d+//g;
            print TMP "$tmp";
        }
        close RAXMLTREE;
        close TMP;
        print BOOTSTRAP "$bootstrap_file_tmp\n";
        system(qq(cp -r $tree_file $output_gene_tree_dir));
    }
    close BOOTSTRAP;
    $logger->info("All OGs raxmlHPC-PTHREADS result bootstrap file path writes to $all_OGs_bootstrap_pathfile ");
    $logger->info("All gene trees were copied to $output_gene_tree_dir");

    $logger->info("========================== build species tree ==========================");
    # 读取数据
    opendir(TREE,$output_gene_tree_dir) || die("Error: can't open input data directory:$output_gene_tree_dir\n");
    my @sample1 = readdir(TREE);
    close TREE;
    my $output_astral_dir = $out_dir."/astral_out";
    system(qq(mkdir -p $output_astral_dir));
    # 生成树合并文件
    my $combine_tree = $output_astral_dir."/combine.tre";
    open(COMBINE,">$combine_tree") || die("Error:can not write to $combine_tree\n");
    foreach my $s(@sample1){
        next if $s=~/^\./;
        my $input_file = $output_gene_tree_dir."/".$s;
        if ($input_file=~/RAxML_bipartitions/){
            open(TREEFILE,"$input_file") || die("The input file $input_file cannot open !");
            while (<TREEFILE>){
                my $tree = $_;
                $tree=~s/_OG\d+:/:/g;
                print COMBINE "$tree";
            }
            close TREEFILE;
        }
    }
    close COMBINE;
    $logger->info("Generate combined gene tree:$combine_tree");
    
    $logger->info("========================== astral started ==========================");
    my $all_species_tree = $output_astral_dir."/astral_coalescent.result";
    unless (-e $combine_tree){
        die("Error: cannot find the combined gene tree:$combine_tree\n");
    } else {
        $logger->info("The command of building species is:[ $astral_path -i $combine_tree -o $all_species_tree -b $all_OGs_bootstrap_pathfile -r $bootstrap ]");
        system(qq($astral_path -i $combine_tree -o $all_species_tree -b $all_OGs_bootstrap_pathfile -r $bootstrap));
    }

    my $last_species_tree = $output_astral_dir."/Astral_coalescence.species.tre";
    unless (-e $all_species_tree){
        die("Error: cannot find the file: $all_species_tree");
    } else {
        system(qq(tail -n 1 $all_species_tree > $last_species_tree));
        $logger->info("Last species astral tree is $last_species_tree");
    }
    $logger->info("==================== Congratulations! All tasks finished!===================");
}


sub iqtree{
    use strict;
    use warnings;
    use Getopt::Std;
    use Log::Log4perl;
    use vars qw($opt_h $opt_t $opt_o $opt_m $opt_B $opt_T);
    getopts("ho:t:m:B:T:");

    Log::Log4perl->init("Build_tree.conf");
    my $logger = Log::Log4perl->get_logger(); 

    my $usage = "\nUsage: Arthropod build_tree iqtree [options] <input_dir>
        
        Necessary input description:
        input_dir       <str>       This directory of trimal trimming result.

        Options:
        -h                          Print this usage page.
        -t              <int>       Multi-threading number. default=4
        -o              <str>       Output directory. default=3_iqtree_out
        -T              <str>       Specify the iqtree number of threads. default=4
        -m              <str>       Set the model of DNA or Protein substitution. default=MFP
        -B              <int>       Replicates for ultrafast bootstrap (>=1000). default=1000
        ";

    die $usage if @ARGV!=1;
    die $usage if defined($opt_h);
    my ($input_dir)=@ARGV;

    $logger->info("=========================== 3. Build tree start ===========================");
    $logger->info("======================== Build gene tree with iqtree =========================");
    # check the availability of raxml and astral
    $logger->error("The iqtree for build gene tree is not available") unless check_command('iqtree');
    my $iqtree_path;
    my $astral_path;
    if (check_command('iqtree')){
        $iqtree_path = check_command('iqtree');
    } else {
        die "ERROR: The iqtree for build gene tree is not available\n";
    }
    $logger->error("The astral for building the species tree is not available") unless check_command('astral');
    if (check_command('astral')){
        $astral_path = check_command('astral');
    } else {
        die "ERROR:The astral for building the species tree is not available\n";
    }
    # create output directory
    my $out_dir="3_iqtree_out";
    $out_dir=$opt_o if defined($opt_o);
    system(qq(mkdir -p $out_dir)) unless -e $out_dir;
    # threads
    my $thread_num = 4;
    $thread_num = $opt_t if defined($opt_t);
    my $thread_iqtree = 4;
    $thread_iqtree = $opt_T if defined($opt_T);
    # set model and bootstrap
    my $model = "MFP";
    $model = $opt_m if defined($opt_m);
    my $bootstrap = 1000;
    $bootstrap = $opt_B if defined($opt_B);
    # 读取数据
    opendir(FASTA,$input_dir) || die("Error: can't open input data directory!\n");
    my @sample = readdir(FASTA);
    close FASTA;
    # 生成多线程shell脚本
    my $iqtree_sh = "iqtree.sh";
    $logger->info("Generate multi-thread iqtree shell script: $iqtree_sh");
    open(SH,">$iqtree_sh") || die("Error:can not write raxmlHPC shell script:$iqtree_sh\n");
    print STDOUT "\nThe commands of $iqtree_sh are:\n";
    foreach my $s(@sample){
        next if $s=~/^\./;
        next if $s=~/reduced/;
        my $input_file = $input_dir."/".$s;
        die("The file $input_file is not existed !") unless -e $input_file;
        my $OG_number = (split/\./,$s)[0];
        my $output_dir = $out_dir."/iqtree/".$OG_number;
        system(qq(mkdir -p $output_dir));
        my $output_suffix = $output_dir."/".$OG_number;
        print STDOUT "$iqtree_path -s $input_file --prefix $output_suffix -B $bootstrap -T $thread_iqtree -m $model --bnni\n";
        print SH "$iqtree_path -s $input_file --prefix $output_suffix -B $bootstrap -T $thread_iqtree -m $model --bnni\n";
    }
    close SH;

    # check multi_threads.pl
    my $multi_threads_path;
    if (check_command('multi_threads.pl')){
        $multi_threads_path=check_command('multi_threads.pl');
    } else {
        die "ERROR: multi_threads.pl was not found! \n";
    }
    $logger->info("The command of iqtree multi-thread shell script:[ perl $multi_threads_path -c $thread_num $iqtree_sh ]");
    system(qq(perl $multi_threads_path -c $thread_num $iqtree_sh));
    $logger->info("======================= Build gene tree with iqtree finished ======================");

    my $output_gene_tree_dir = $out_dir."/gene_tree";
    system(qq(mkdir -p $output_gene_tree_dir));
    foreach my $s(@sample){
        next if $s=~/^\./;
        next if $s=~/reduced/;
        my $OG_number = (split/\./,$s)[0];
        my $tree_file = $out_dir."/iqtree/$OG_number/$OG_number.treefile";
		if (-e $tree_file){
			system(qq(cp -r $tree_file $output_gene_tree_dir));
		} else {
			print "The tree file $tree_file is not existed !";
		}
    }
    $logger->info("All iqtree gene trees were copied to $output_gene_tree_dir");

    $logger->info("========================== build species tree ==========================");
    # 读取数据
    opendir(TREE,$output_gene_tree_dir) || die("Error: can't open input data directory:$output_gene_tree_dir\n");
    my @sample1 = readdir(TREE);
    close TREE;
    my $output_astral_dir = $out_dir."/astral_out";
    system(qq(mkdir -p $output_astral_dir));
    # 生成树合并文件
    my $combine_tree = $output_astral_dir."/combine.tre";
    open(COMBINE,">$combine_tree") || die("Error:can not write to $combine_tree\n");
    foreach my $s(@sample1){
        next if $s=~/^\./;
        my $input_file = $output_gene_tree_dir."/".$s;
        if ($input_file=~/treefile/){
            open(TREEFILE,"$input_file") || die("The input file $input_file cannot open !!!");
            while (<TREEFILE>){
                my $text = $_;
                $text=~s/_OG\d+:/:/g;
                print COMBINE "$text";
            }
            close TREEFILE;
        }
    }
    close COMBINE;
    $logger->info("Generate combined gene tree:$combine_tree");
    
    $logger->info("========================== astral started ==========================");
    my $all_species_tree = $output_astral_dir."/astral_coalescent.result";
    unless (-e $combine_tree){
        die("Error: cannot find the combined gene tree:$combine_tree\n");
    } else {
        $logger->info("The command of building species is:[ $astral_path -i $combine_tree -o $all_species_tree ]");
        system(qq($astral_path -i $combine_tree -o $all_species_tree));
    }

    $logger->info("Last species astral tree is $all_species_tree");
    $logger->info("================== Congratulations! All tasks of building tree finished !=================");
}



sub tree_plot{
    use strict;
    use warnings;
    use Getopt::Std;
    use Log::Log4perl;
    use vars qw($opt_h);
    getopts("h");

    Log::Log4perl->init("Build_tree.conf");
    my $logger = Log::Log4perl->get_logger(); 

    my $usage = "\nUsage: Arthropod build_tree tree_plot [options] <tree_file> <group_file> <output_prefix>
        
        Necessary input description:
        tree_file       <str>       The result tree file of raxml or astral.
        group_file      <str>       The group information of species in tree file (two columns).
                                    Example: 0=>black     1=>red
                                    Homo_sapiens    0
                                    Apis_cerana     1
        output_prefix   <str>       The output file prefix.

        Options:
        -h                          Print this usage page.
        ";

    die $usage if @ARGV!=3;
    die $usage if defined($opt_h);
    my ($tree_file, $group_file, $output_prefix)=@ARGV;

    $logger->error("Error: $tree_file was not exists !\n") unless (-e $tree_file);
    die("Error: $tree_file was not exists !\n") unless (-e $tree_file);
    $logger->error("Error: $group_file was not exists !\n") unless (-e $group_file);
    die("Error: $group_file was not exists !\n") unless (-e $group_file);

    $logger->info("Rscript tree_plot.r $tree_file $group_file $output_prefix");
    system(qq(Rscript tree_plot.r $tree_file $group_file $output_prefix));
}
1;

