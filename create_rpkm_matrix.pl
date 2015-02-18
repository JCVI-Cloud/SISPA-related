#!/usr/local/bin/perl -w
use warnings;
use strict;
$|++;

=head1 NAME

create_rpkm_matrix.pl

=head1 SYNOPSIS

    USAGE: create_rpkm_matrix.pl  -f </path/to/list_of_fasta_reads>
                                  -b </path/to/blast_db>
                                  -i <blast identity cutoff>
                                  -l <blast alignment length cutoff> 
                                  -c config.file
                                  [--no_blast <skips blast>
                                   --output <output directory>
                                   --help ]

=head1 OPTIONS

B<--fasta_list,f>          : File containing list of fasta files of reads and their associated genome names

B<--blast_db,b>            : Blast database to use

B<--identity,i>            : Blast identity cutoff to use for filtering results [ie. 99]

B<--alignment,l>           : Blast alignment cutoff to use for filtering results [ie. 14]

B<--no_blast>              : Will skip blast step. Will assume blast files named <genome>.blast are located in the output directory. 
                             Useful for when you want to create new rpkm matrices with different cutoffs, or when you are providing
                             already generated blast results (must be in -m8 format).

B<--config,c>              : Config file pointing to location of installed 3rd party tools

B<--help,-h>               : Display this help message.

B<--output,o>              : Output directory for blast files and log files
     
=head1  INPUT

--fasta_list : This script requires a list of fasta file locations and the associated genome name you'd like to give:
ie.UH0107.fasta<tab>UH0107

--config :  A generic config file used in JCVI's SISPA post processing scripts. These scripts require some 3rd party tools. 
The config allows the user to specify where these 3rd party tools are installed on their machines. Either find the config.txt
that was packaged with this script or copy the output below in a new text file.

raxmlHPC:/usr/local/bin/raxmlHPC
vcftools:/usr/local/bin/vcftools
formatdb: /usr/local/bin/formatdb
blastall:/usr/local/bin/blastall

=head1  DESCRIPTION

This script takes a list of Fasta file reads and run Blastn against a provided antibiotic resistance (abR) gene database. 
It then generates a matrix of RPKM values for each gene (row) and genome (column).                  

=head1  CONTACT

    Erin Beck
    ebeck@jcvi.org

=cut

use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use Pod::Usage;
use File::Slurp;
use Data::Dumper;
use File::Basename;
use Bio::AlignIO; 
use Cwd 'abs_path';
use Cwd;

my %opts;

#Arguments/Options
GetOptions( \%opts, 
	    'fasta_list|f=s',
	    'blast_db|b=s',
	    'identity|i=s',
	    'alignment|l=s',
	    'no_blast|n',
	    'output|o=s',
	    'config|c=s',
	    'help|h') || die "Error getting options! $!";
pod2usage( { -exitval => 1, -verbose => 2 } ) if $opts{help};

my ($output,$scripts) = check_params();

#Start logging
my $log_file = "$output/rpkm_matrix.log";
open( my $lfh, '>', $log_file ) || die "Can't open log file $log_file: $!";

#Wrapper to run grid array jobs
#Read fasta read input
my $files = parse_fasta_list($opts{fasta_list});

my $rpkm_hsh;
my @dbs;

#Format blast db
if(exists $scripts->{formatdb_exe}){
    my $formatdb_exe = $scripts->{formatdb_exe};
    my $formatcmd = "$formatdb_exe -p F -i $opts{blast_db}";
    &_log("Formatting blast db: $formatcmd");
    system($formatcmd) == 0 || die("ERROR: Failed running formatdb: $formatcmd") unless($opts{no_blast});
}else{
    die("ERROR: No location found in config file for formatdb");
}

#parsed db to calculate gene length
my $gene_lengths = parse_blast_db_for_length($opts{blast_db});

#Foreach fasta file preform actions
foreach my $file (keys %$files){

    $file =~ s/\s+$//;
    my $blast_file = "$output/$files->{$file}.blast";

    #Find the number of reads in fasta
    my $number_of_reads = `grep -c '>' $file`;
    $number_of_reads =~ s/\s+$//;

    #Blasts files locally
    blast_reads($file,$opts{blast_db},$files->{$file}); 
  
    #If blast file was found calculate rpkm
    if(-s $blast_file){

	&_log("Calculating rpkm values for $file");
	my $parsed_blast = parse_blast($blast_file,$opts{alignment},$opts{identity},$gene_lengths);
	$rpkm_hsh = calculate_rpkm($parsed_blast,$number_of_reads,$rpkm_hsh,$files->{$file},$gene_lengths);

    }else{

	die("ERROR: $blast_file not found");
    
    }
    
    #Keep track of genome names
    push(@dbs,$files->{$file});
}

#Print combined rpkm matix
&_log("Calculating rpkm matrix and printing");
print_rpkm($rpkm_hsh,\@dbs);

exit(0);

## Subs ##
sub parse_blast_db_for_length{
    
    my $file = shift;

    open(my $fh, "<", $file);

    my $current_gene;
    my $current_seq;
    my $hsh;

    while(<$fh>){
	my $line = $_;
	$line =~ s/\s+$//;

	if($line =~ /^>(.*)/){

	    if($current_gene){
		$current_gene =~ s/\s+$//;
		$current_seq =~ s/\s+//g;

		$hsh->{$current_gene} = length $current_seq;
		$current_seq = "";
	    }

	    $current_gene = $1;

	}else{
	    
	    $current_seq .= $line;

	}

    }

    #Push final gene to hash
    $current_gene =~ s/\s+$//;
    $current_seq =~ s/\s+//g;
    $hsh->{$current_gene} = length $current_seq;

    return $hsh;
}

sub parse_fasta_list{

    my $file = shift;

    my @lines = read_file($file);
    my $files;

    foreach my $line(@lines){

	$line =~ s/\s+$//;

	my ($path,$genome) = split(/\t/,$line);

	$files->{$path} = $genome;
    }

    return $files;
}
sub print_rpkm{

    my ($hsh,$dbs) = @_;
    my $out_name = "rpkm_matrix" . "_" . $opts{identity} . "_" . $opts{alignment} . ".txt";
    open(my $ofh, ">", "$output/$out_name");
    select $ofh;

    #print headers
    print "gene\t";
    print join("\t",sort @$dbs);
    print "\n";

    foreach my $gene (keys %$hsh){

	print $gene;

	foreach my $genome (sort @$dbs){
	    
	    if(exists $hsh->{$gene}->{$genome}){

		print "\t$hsh->{$gene}->{$genome}";

	    }else{

		print "\t";

	    }

	}

	print "\n";

    }

}
sub calculate_rpkm{
    
    my ($blast,$total_reads,$rpkm_hsh,$genome,$lengths) = @_;

    foreach my $gene (keys  %$blast){

	$gene =~ s/\s+/_/g;

	#Pull gene length
	my $gene_length = $lengths->{$gene};

	#Find number of reads that hit this gene
	my $num_reads = scalar keys %{$blast->{$gene}};
	
	#Calculate rpkm
	my $rpkm = ($num_reads * (10 ** 9)) / ($total_reads * $gene_length);
	$rpkm_hsh->{$gene}->{$genome} = sprintf("%.3f", $rpkm);
    }

    return $rpkm_hsh;

}
sub parse_blast{

    my ($blast_results,$alignment_cutoff,$identity_cutoff,$gene_lengths) = @_;

    open(my $bfh, "<", $blast_results);
    my $hsh;

    my $problem_genes;

    #go through each blast hit to filter based
    #on identity and alignment length
    while(<$bfh>){

	my $line = $_;
	$line =~ s/\s+$//;

	my @values = split(/\t/,$line);

	#Only keep hits that are above or equal to the set
	#identity cutoff and alignment length cutoff
	if(($values[2] >= $identity_cutoff) && ($values[3] >= $alignment_cutoff)){

	    $hsh->{$values[1]}->{$values[0]}->{identity} = $values[2];
	    $hsh->{$values[1]}->{$values[0]}->{alignment} = $values[3];

	    unless(exists $gene_lengths->{$values[1]}){

		$problem_genes->{$values[1]} = 1;

	    }

	}
    
    }

    if($problem_genes){
	
	my @genes = keys %$problem_genes;

	print "\nERROR: Some gene names were found in the blast result files that do not match the names found in the initial blast db. Please make sure there are no spaces in the gene names found in the blast database.\n";

	print join("\n",@genes);
	print "\n";

	exit;
    }

    return $hsh;
}
sub blast_reads{

    my ($file,$blastdb,$name) = @_;

    #Find the number of reads in fasta
    my $number_of_reads = `grep -c '>' $file`;
    $number_of_reads =~ s/\s+$//;

    #Run blast
    if(exists $scripts->{blastall_exe}){
	my $blast_exe = $scripts->{blastall_exe};
	my $blastcmd = "$blast_exe -p blastn";
	$blastcmd .= " -d $blastdb -i $file";
	$blastcmd .= " -m8 -F f -e .00001";
	$blastcmd .= " -o $name.blast";

	&_log("Running blast jobs locally: $blastcmd");
	system($blastcmd) == 0 || die("ERROR: Failed running blast: $blastcmd") unless($opts{no_blast});
    }else{
	die("ERROR: No location found in config file for blastall");
    }
    return("$name.blast",$number_of_reads,$name);
}
sub check_params{
    
    my $errors = "";
    my $scripts;

    unless($opts{identity} && $opts{fasta_list} && $opts{alignment} && $opts{blast_db}){
	$errors .= "\nUsage: ./create_rpkm_matrix.pl -f <FASTA reads> -b <AbR gene database> -i <identity cutoff> -l <alignmnt length cutoff> [--help]\n\n";
    }

    if($opts{config}){
	
	if(-s $opts{config}){

	    my @lines = read_file($opts{config});
	    
	    foreach my $line (@lines){

		$line =~ s/\s+$//;
		
		unless($line =~ /^#/){
		    my ($exe, $location) = split(/:/,$line);
		    $scripts->{$exe . "_exe"} = $location;
		}

	    }

	}else{
	    
	    $errors .= "File does not exist or is size zero: $opts{config}\n";
	}
       
    }else{

	$errors .= "Must provide config file\n";

    }

    my $output = $opts{output} // cwd;

    die($errors) if $errors;
    
    return($output,$scripts);
}
sub _log {

    my ( $msg, $lvl ) = ( @_ );

    chomp( $msg );
    print $lfh localtime . ":$msg\n";

}
