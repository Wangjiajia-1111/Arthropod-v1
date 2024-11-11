#!/usr/bin/perl
package assembly;

# check software is available
sub check_command{
    my ($cmd) = @_;
    #return `which $cmd`;
    my $path = `which $cmd`;
    $path=~s/\n//g;
    return $path;
}

use Bio::Seq;
use Bio::SeqIO;
sub transCDS2PEP{
    my ($coding_seq) = @_;
    my $seq_obj = Bio::Seq->new(-seq => $coding_seq, -alphabet => 'dna');
    my $seq_trans = $seq_obj->translate();
    my $aa_sequence = $seq_trans->seq();
    $aa_sequence=~s/\*//g;
    return $aa_sequence;
}


sub assemble{
    my $usage="
        Usage: Arthropod assemble [commands] ...

        Commands:
        fastp               Quality control and data-filtering of FASTQ files.
        easy353             Based on reference sequences filter and assemble the transcriptome, the whole 
                            genome or genome skimming sequencing data to recover target genes.
        captus              Build phylogenomic datasets from multiple types of sequencing data based on 
                            reference sequences.
        transSeq            Convert coding sequence into amino acid sequence.
        ";
    
    die $usage if @ARGV<1;
    my $com=shift @ARGV;
    if ($com eq "fastp"){
        fastp(@ARGV);
    }
    elsif ($com eq "easy353"){
        easy353(@ARGV);
    }
    elsif ($com eq "captus"){
        captus(@ARGV);
    }
    elsif ($com eq "transSeq"){
        transSeq(@ARGV);
    }
    else {
        print STDERR "Unknown command: $com\n";
        die($usage);
        exit(-1);
    }
}


sub fastp{
    use strict;
    use warnings;
    use Log::Log4perl;
    use Getopt::Std;
    use vars qw($opt_h $opt_t);
    getopts("ht:");

    Log::Log4perl->init("Aseemble.conf");
    my $logger = Log::Log4perl->get_logger();

    my $usage = "\nUsage: Arthropod assemble fastp [options] <fq1> <fq2> <output_dir> <prefix>;

    Arthropod assemble fastp is used to quality control and data-filtering of FASTQ files.

    Necessary :
    fq1                                     The input files with paired-end reads, given in FASTQ format
    fq2                                     The input files with paired-end reads, given in FASTQ format
    output_dir                              The output files directory
    prefix                                  Prefix of output file, you can use the latin name of output data 
                                            eg: Homo_sapiens
    Options:
    -h                                      Print this usage page
    -t                      <int>           worker thread number, default=4
    ";

    die $usage if @ARGV!=4;
    die $usage if defined($opt_h);
    my ($fq1, $fq2, $output_dir, $prefix) = @ARGV;

    # check the availability of fastp
    $logger->info("=========================== fastp start ===========================");
    $logger->error("The fastp is not available!") unless check_command('fastp');
    my $fastp_path;
    if (check_command('fastp')){
        $fastp_path = check_command('fastp');
    } else {
        die("Error: cannot find fastp!");
    }

    ### read samples and process each sample
    die("Error: can not find the input file $fq1!\n") unless -e $fq1;
    die("Error: can not find the input file $fq2!\n") unless -e $fq2;
    unless (-e $output_dir){
        system(qq(mkdir -p $output_dir));
    }
    
    # work threads
    my $thread_num=4;
    $thread_num=$opt_t if defined($opt_t);
    # output file name
    my $output_fq1 = $output_dir."/".$prefix."_R1.fq.gz";
    my $output_fq2 = $output_dir."/".$prefix."_R2.fq.gz";

    $logger->info("The command of fastp script:[$fastp_path -i $fq1 -I $fq1 -o $output_fq1 -O $output_fq2 -w $thread_num]");
    system(qq($fastp_path -i $fq1 -I $fq1 -o $output_fq1 -O $output_fq2 -w $thread_num));
    $logger->info("=========================== fastp finished ===========================");
}


sub easy353{
    use strict;
    use warnings;
    use Log::Log4perl;
    use Getopt::Std;
    use vars qw($opt_h $opt_k $opt_K $opt_t $opt_T $opt_l);
    getopts("hk:K:t:T:l:");
    
    Log::Log4perl->init("Aseemble.conf");
    my $logger = Log::Log4perl->get_logger();

    my $usage = "\nUsage: Arthropod assemble easy353 [options] <fq1> <fq2> <reference_dir> <latin_name>;

    Arthropod assemble easy353 is used to filter and assemble reads to recover target genes.
    Necessary:
    fq1 fq2                             The input files with paired-end reads, given in FASTQ format.
    reference_dir                       The reference fasta directory.
    latin_name                          The latin name of output data. 
                                        output_dir: latin_name_easy353out, latin_name_easy353merge
    Options:
    -h                                  Print this usage page
    -k              <int>               K-mer length setting for filtering. Default:21
    -K              <int>               K-mer length setting for assembly. Default:31
    -t              <int>               the threads setting for reads filtering. Default:4
    -T              <int>               the threads setting for reads assembly. Default:8
    -l              <kmer_limit>        Set a limit of kmer count. Default:4
    ";

    die $usage if @ARGV!=4;
    die $usage if defined($opt_h);
    my ($fq1, $fq2, $reference_dir, $latin_name) = @ARGV;
    my $output_easy353_dir = $latin_name."_easy353out";

    # check the availability of easy353
    $logger->info("=========================== easy353 start ===========================");
    $logger->error("The easy353.py is not available!") unless check_command('easy353.py');
    my $easy353_path;
    if (check_command('easy353.py')){
        $easy353_path = check_command('easy353.py');
        $easy353_path=~s/\n//g;
    } else {
        die("Error: cannot find easy353.py!");
    }

    # read samples and process each sample
    die("Error: can not find the input file $fq1!\n") unless -e $fq1;
    die("Error: can not find the input file $fq2!\n") unless -e $fq2;
    die("Error: can not find the reference directory $reference_dir!\n") unless -e $reference_dir;
    unless (-e $output_easy353_dir){
        system(qq(mkdir -p $output_easy353_dir));
    }

    # kmer
    my $filter_kmer = 21;
    $filter_kmer = $opt_k if defined($opt_k);
    my $assemble_kmer = 31;
    $assemble_kmer = $opt_K if defined($opt_K);
    my $kmer_limit = 4;
    $kmer_limit = $opt_l if defined($opt_l);

    # threads
    my $filter_thread = 4;
    $filter_thread = $opt_t if defined($opt_t);
    my $assemble_thread = 8;
    $assemble_thread = $opt_T if defined($opt_T);

    $logger->info("$easy353_path -1 $fq1 -2 $fq2 -r $reference_dir -o $output_easy353_dir -fk $filter_kmer -ak $assemble_kmer -ft $filter_thread -at $assemble_thread -kmer_limit $kmer_limit");
    system(qq($easy353_path -1 $fq1 -2 $fq2 -r $reference_dir -o $output_easy353_dir -fk $filter_kmer -ak $assemble_kmer -ft $filter_thread -at $assemble_thread -kmer_limit $kmer_limit));
    $logger->info("========================== easy353 finished ===========================");

    $logger->info("========== Merge easy353 assemble fasta with reference OG fasta =========");
    use List::Util qw(min);
    opendir(REF,$reference_dir) || die("Error:can not open $reference_dir !");
    my @reference = readdir(REF);
    close REF;
    my $OutputDir = $latin_name."_easy353merge";
    unless (-e $OutputDir){
        system(qq(mkdir -p $OutputDir));
    }
    foreach my $r(@reference){
        next if $r=~/^\./;
        next if $r=~/captus\.faa/;
        my $OG_fasta = $reference_dir."/$r";
        my $OG_id = (split/\./,$r)[0];
        my $output_file = $OutputDir."/$OG_id.merge.fasta";
        open(OUTPUT,">$output_file");
        open(FASTA,"$OG_fasta");
        my @lens=();
        while (<FASTA>){
            chomp;
            print OUTPUT "$_\n";
            next if $_=~/^>/;
            my $len = length($_);
            push @lens,$len;
        }
        close FASTA;
        my $min_len = min(@lens);

        my $easy353_fasta = $output_easy353_dir."/assemble_out/$r";
		if (-e $easy353_fasta) {
			open(FASTA1,"$easy353_fasta");
			my $fasta1 = do { local $/; <FASTA1> };
			my @aa = split/>/,$fasta1;
			shift @aa;
			foreach my $a(@aa){
				my @bb=split/\n/,$a,2;
				my $seqence=$bb[1];
				$seqence=~s/\n//g;
				my $len = length($seqence);
				if ($len>=$min_len){
					print OUTPUT ">$latin_name\_$OG_id\n$seqence\n";
				}
			}
			close FASTA1;
			close OUTPUT;
		} else {
			system(qq(cp $OG_fasta $output_file));
		}
    }
    $logger->info("Merge fasta is under $OutputDir");
    $logger->info("All easy353 jobs are finished !!!");
}


sub captus{
    use strict;
    use warnings;
    use Log::Log4perl;
    use Getopt::Std;
    use vars qw($opt_h $opt_t $opt_T $opt_N $opt_d $opt_c $opt_R $opt_p);
    getopts("ht:T:N:d:c:R:p:");
    
    Log::Log4perl->init("Aseemble.conf");
    my $logger = Log::Log4perl->get_logger();

    my $usage = "\nUsage: Arthropod assemble captus [options] <FASTQ_files> <latin_name> <nuc_refs_dir>;
    
    Arthropod assemble captus is used to assemble of phylogenomic datasets from RNA-Seq, genome skimming and high-depth whole genome sequencing data.

    Necessary:
    FASTQ_files                             FASTQ files. Valid file name extensions are:.fq,.fastq,.fq.gz
                                            and .fastq.gz.The names must include the string '_R1' ('_R2'). 
                                            Everything before the string '_R1' will be used as sample name.
                                            There are a few ways to provide the FASTQ files:
                                            A directory = path to directory containing FASTQ files
                                            A list = file names separated by space
                                            A pattern = UNIX matching expression
    latin_name                              The latin name of output data. eg:Homo_sapiens
                                            Output directory: latin_name_CAPTUSout, latin_name_CAPTUSmerge
    nuc_refs_dir                            Provide a directory path cotain include FASTA files containing 
                                            your reference target protein sequences in either nucleotide or
                                            aminoacid. When the FASTA file is in nucleotides, '-N' will be 
                                            used to translate it to aminoacids.

    Options:
    -h                                      Print this usage page.
    -t                  <int>               The assemble threads. default=60
    -T                  <int>               The extract threads. default=30
    -N                  <int>               When the FASTA file is in nucleotides, '-N' will be used to 
                                            translate it to aminoacids. Genetic code table to translate 
                                            your nuclear proteins. Complete list of tables at: 
                                            https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
                                            default=1 :Standard
    -d                  <str>               Create the temporary directory 'captus_assembly_tmp' for 
                                            MEGAHIT assembly. Sometimes, when working on external hard 
                                            drives MEGAHIT will refuse to run unless this directory is 
                                            created in an internal hard drive. (default: current path)
    -c                  <int>               Multi-threading number of multi_threads.pl. default=4
    -R                  <int>               Maximum RAM in GB (e.g.: 4.5) dedicated to Captus,'auto' uses 
                                            99% of available RAM (default: auto)
    -p                  <PRESET>            The default preset is 'CAPSKIM', these settings work well with 
                                            either hybridization capture or genome skimming data 
                                            (or a combination of both). You can assemble RNA-Seq reads with
                                            the 'RNA' preset or high-coverage Whole Genome Sequencing reads
                                            with the 'WGS' preset,however, both presets require a minimum 
                                            of 8GB of RAM to work well.
                                            CAPSKIM = --k-list 31,39,47,63,79,95,111,127,143,159,175 
                                                    --min-count 2 --prune-level 2
                                            RNA = --k-list 27,47,67,87,107,127,147,167 
                                                    --min-count 2 --prune-level 2
                                            WGS = --k-list 31,39,49,69,89,109,129,149,169 
                                                    --min-count 3 --prune-level 2
    ";

    die $usage if @ARGV!=3;
    die $usage if defined($opt_h);
    my ($fastq_flie, $latin_name, $reference_dir) = @ARGV;
    
    my  $output_dir = $latin_name."_CAPTUSout";
    
    # check the availability of captus
    $logger->info("======================== Captus assemble start ========================");
    $logger->error("The captus_assembly is not available!") unless check_command('captus_assembly');
    my $captus_path;
    if (check_command('captus_assembly')){
        $captus_path = check_command('captus_assembly');
    } else {
        die("Error: cannot find captus_assembly !");
    }
    
    die("Error: can not find the fastq file $fastq_flie!\n") unless -e $fastq_flie;
    unless (-e $output_dir){
        system(qq(mkdir -p $output_dir));
    }

    # threads
    my $assemble_threads = 60;
    $assemble_threads = $opt_t if defined($opt_t);
    my $extract_threads = 30;
    $extract_threads = $opt_T if defined($opt_T);
    # RAM
    my $RAM = "auto";
    $RAM = $opt_R if defined($opt_R);

    # PRESET
    my $preset = "CAPSKIM";
    $preset = $opt_p if defined($opt_p);

    # temporary directory 
    use Cwd;
    my $temporary_directory = getcwd();
    $temporary_directory = $opt_d if defined($opt_d);

    my $assemble_out_dir = $output_dir."/assemble";
    if (-e $assemble_out_dir){
        system(qq(rm -rf $assemble_out_dir));
    }

    $logger->info("$captus_path assemble -r $fastq_flie -o $assemble_out_dir --threads $assemble_threads --tmp_dir $temporary_directory --preset $preset --ram $RAM --concurrent 1");
    system(qq($captus_path assemble -r $fastq_flie -o $assemble_out_dir --threads $assemble_threads --tmp_dir $temporary_directory --preset $preset --ram $RAM --concurrent 1));

    my $assemble_fasta = $assemble_out_dir."/$latin_name\__captus-asm/01_assembly/assembly.fasta";
    die("Captus assembly fasta $assemble_fasta is not found !!! Please view the captus_assembly log file for details.") unless -e $assemble_fasta;
    $logger->info("============================= Captus assemble finished ============================");

    
    $logger->info("============================= Captus extract start =============================");
    my $extract_out_dir = $output_dir."/extract";
    # # check the availability of multi_threads.pl
    my $multi_threads = 4;
    $multi_threads = $opt_c if defined($opt_c);
    my $multi_threads_path;
    if (check_command("multi_threads.pl")){
        $multi_threads_path = check_command("multi_threads.pl");
    } else {
        die("Error: the multi_threads script multi_threads.pl is not find !");
    }

    opendir(REF_DIR,"$reference_dir");
    my @reference = readdir(REF_DIR);
    closedir REF_DIR;
    my $nuc_transtable = 1;
    # 生成多线程shell脚本
    my $captus_extract_sh = "captus_extract.sh";
    $logger->info("Generate multi-thread captus_assembly extract script: $captus_extract_sh");
    open(SH,">$captus_extract_sh") || die("Error:can not write captus_assembly extract shell script:$captus_extract_sh\n");
    foreach my $r(@reference){
        next if $r=~/^\./;
        my $OG_id = (split/\./,$r)[0];
        my $reference_fasta = $reference_dir."/$r";
        if (defined($opt_N)){
            $nuc_transtable = $opt_N;
            print SH "$captus_path extract -a $assemble_out_dir -n $reference_fasta -o $extract_out_dir/$OG_id --threads $extract_threads --nuc_transtable $nuc_transtable\n";
        } else {
            print SH "$captus_path extract -a $assemble_out_dir -n $reference_fasta -o $extract_out_dir/$OG_id --threads $extract_threads\n";
        }
    }
    close SH;
    
    $logger->info("The command of captus_assembly extract multi-threads script:[ perl $multi_threads_path -c $multi_threads $captus_extract_sh ]");
    system(qq(perl $multi_threads_path -c $multi_threads $captus_extract_sh));
    $logger->info("============================== Captus extract finished ==============================");
    my $nuc_transtable_dir = $output_dir."/nuc_transtable";
    system(qq(mkdir -p $nuc_transtable_dir));
    foreach my $r(@reference){
        next if $r=~/^\./;
        my $OG_id = (split/\./,$r)[0];
        my $nuc_transtable_fasta = $reference_dir."/$OG_id.captus.faa";
        if (-e $nuc_transtable_fasta) {
            system(qq(mv $nuc_transtable_fasta $nuc_transtable_dir));
        }
    } 
    
    $logger->info("================== Merge captus assemble fasta with reference OG fasta =================");
    my $output_mergeFasta_dir = $latin_name."_CAPTUSmerge";
    system(qq(mkdir -p $output_mergeFasta_dir)) unless -e $output_mergeFasta_dir;
    foreach my $R(@reference){
        next if $R=~/^\./;
        my $reference_fasta = $reference_dir."/$R";
        my $OG_id = (split/\./,$R)[0];

        my $output_mergeFasta_file = $output_mergeFasta_dir."/$OG_id.merge.fasta";
        open(OUTPUT,">$output_mergeFasta_file");
        
        my @lens=();
	if (-e $reference_fasta){
	    open(FASTA,"$reference_fasta");
	} else {
		print "Warning: can not open $reference_fasta !!!\n";
	}
	#open(FASTA,"$reference_fasta") || die("Error:can not open $reference_fasta !");
        while (<FASTA>){
            chomp;
            print OUTPUT "$_\n";
            next if $_=~/^>/;
            my $len = length($_);
            push @lens,$len;
        }
        close FASTA;

        use List::Util qw(min);
        use List::Util qw(max);
        my $min_len = min(@lens);

        my $captus_fasta = "$output_dir/extract/$OG_id/$latin_name\__captus-ext/01_coding_NUC/NUC_coding_NT.fna";
        if (-e $captus_fasta){
            open(FASTA1,"$captus_fasta") || die("Error: $captus_fasta can not open !");
            my $fasta1 = do {local $/; <FASTA1>};
            my @aa = split/>/,$fasta1;
            shift @aa;
			my $max_seqence_len = 0;
			my $max_seqence = "";
            foreach my $a(@aa){
                my @bb=split/\n/,$a,2;
                my $seqence=$bb[1];
                $seqence=~s/\n//g;
				my $length = length($seqence);
                if ($length > $max_seqence_len) {
					$max_seqence_len = $length;
					$max_seqence = $seqence;
				}
            }
            if ($max_seqence_len >= $min_len){
                print OUTPUT ">$latin_name\_$OG_id\n$max_seqence\n";
            }
            close FASTA1;
            close OUTPUT;
        } else {
            $logger->info("$captus_fasta is not exist!!!");
        }
    }
    $logger->info("Merge captus assemble fastas are under $output_mergeFasta_dir");
    $logger->info("============================== All captus task finished ! =============================");
}


sub transSeq{
    use strict;
    use warnings;
    use Log::Log4perl;
    use Getopt::Std;
    use vars qw($opt_h $opt_o);
    getopts("ho:");

    Log::Log4perl->init("Aseemble.conf");
    my $logger = Log::Log4perl->get_logger();

    my $usage = "\nUsage: Arthropod assemble transSeq [options] <input_dir>;

    Arthropod assemble transSeq is used to convert coding sequence into amino acid sequence.
    
    Necessary:
    input_dir       <str>                   This directory should contain one or more coding sequence file.

    Options:
    -h                                      Print this usage page.
    -o              <str>                   Output directory. Defect:input_dir_pep
    ";

    die $usage if @ARGV!=1;
    die $usage if defined($opt_h);
    my ($input_dir) = @ARGV;

    # Output directory
    my $output_dir = $input_dir."_pep";
    $output_dir = $opt_o if defined($opt_o);
    system(qq(mkdir -p $output_dir));

    opendir(CDS_DIR,"$input_dir") || die("Error:can not open the directory $input_dir !");
    my @fastas = readdir(CDS_DIR);
    close CDS_DIR;
    foreach my $fasta(@fastas){
        next if $fasta=~/^\./;
        my $input_fasta = $input_dir."/$fasta";
        my $output_fasta = $output_dir."/$fasta";
        
        open(FASTA,$input_fasta);
        open(OUTPUT,">$output_fasta");
        my $contence = do { local $/; <FASTA> };
        my @sections = split/>/,$contence;
        shift @sections;
        foreach my $section(@sections){
            my @aa = split/\n/,$section,2;
            my $id = $aa[0];
            my $seq = $aa[1];
            $seq=~s/\n//g;
            my $trans_aa_sequence = transCDS2PEP($seq);
            print OUTPUT ">$id\n$trans_aa_sequence\n";
        }
        print "$input_fasta convert $output_fasta successfully.\n";
        close FASTA;
        close OUTPUT;
    }
}

1;

