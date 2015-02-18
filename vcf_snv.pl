#!/usr/local/bin/perl -w
use warnings;
use strict;
$|++;

=head1 NAME

vcf_snv.pl - Wrapper to parse kSNP VCF files and run raxML

=head1 SYNOPSIS

    USAGE: vcf_snv.pl -v kSNP VCF File
                      -p recombinant region coordinate file
                      -t VCF File filtering threshold
                      -c config.file
                     [--output <output directory>
                      --help ]
=head1 OPTIONS

B<--vcf_file,v>      : kSNP VCF file output

B<--positions,p>     : List of coordinates for recombinant regions of a reference genomes

B<--threshold,t>     : Percentage thershold to use for SNV parsing. [Default: 90]

B<--config,c>        : Config file pointing to location of installed 3rd party tools

B<--help,-h>         : Display this help message.

B<--output,o>        : Output directory for blast files and log files

=head1  DESCRIPTION

 This scripts takes a VCF file from kSNP and does the following:

 1) Filter VCF file to remove recombinant regions
 2) Creat a SNV matrix (snv_matrix.txt), only retaining rows that have at least the threshold(ie. 90%) of genomes present
 3) Create a SNV Fasta file and Phylip file for each genome (vcf_snv.fasta, vcf_snv.phy)
 4) Run raxML 

 File will be placed in the current working directory. 

=head1  INPUT

--config :  A generic config file used in JCVI's SISPA post processing scripts. These scripts require some 3rd party tools. 
The config allows the user to specify where these 3rd party tools are installed on their machines. Either find the config.txt
that was packaged with this script or copy the output below in a new text file.

raxmlHPC:/usr/local/bin/raxmlHPC
vcftools:/usr/local/bin/vcftools
formatdb: /usr/local/bin/formatdb
blastall:/usr/local/bin/blastall

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

#GLOBAL VARS
my $THRESHOLD;

#Arguments/Options
GetOptions( \%opts, 'vcf_file|v=s',
	    'positions|p=s',
	    'threshold|t=s',
	    'config|c=s',
	    'output|o=s',
	    'help|h') || die "Error getting options! $!";
pod2usage( { -exitval => 1, -verbose => 2 } ) if $opts{help};

my($output,$scripts) = &check_params;

my $parsed_vcf_file = filter_vcf($opts{vcf_file},$opts{positions});
my ($vcf_info_hsh,$genomes,$pos_order) = parse_vcf($parsed_vcf_file);
my $snv_fasta = print_snv_matrix($vcf_info_hsh,$genomes,$pos_order);
my $fasta_file = print_snv_fasta($snv_fasta);
my $phylip = convert_fasta_phylip($fasta_file);
run_raxmlHPC($phylip);

### Subs ###

sub run_raxmlHPC{
    
    my $phy_file = shift;

    die("ERROR: No phylip file found, $phy_file") unless(-s $phy_file);
 
    my $current_dir = cwd;

    if(exists $scripts->{raxmlHPC_exe}){
	my $raxml_exe = $scripts->{raxmlHPC_exe};
	my $raxml_cmd = $raxml_exe . " -f a -s " . abs_path($phy_file) . " -n boot -m GTRCAT -x 200 -p 200 -# 100";
    
	my $std_err = "$current_dir/vcf_snv.stderr";
	my $std_out = "$current_dir/vcf_snv.stdout";

	system($raxml_cmd) == 0 || die ( "ERROR: Problem running raxml:$raxml_cmd");

    }else{
	die("ERROR: No location found in config file for raxmlHPC");
    }
}

sub convert_fasta_phylip{
    my $file = shift;

    die("ERROR: No fasta file found, $file") unless(-s $file);

    my $outfile = "vcf_snv.phy";

    my $in  = Bio::AlignIO->new(-file   => $file ,
				-format => 'fasta'
	);
    my $out = Bio::AlignIO->new(-file   => ">$outfile" ,
				-format => 'phylip');
    
    while ( my $aln = $in->next_aln() ) {
	$out->write_aln($aln);
    }

    return $outfile;
}

sub print_snv_fasta{

    my $fasta_hsh = shift;
    my $filename = "vcf_snv.fasta";

    open(my $fh, ">",$filename);

    #Print a multi fasta file containing the combined
    #sequence for all genomes
    foreach my $genome(keys %$fasta_hsh){

        my $print_length = 60;
        my $header = ">$genome\n";
	my $formatted_seq;

        while (my $chunk = substr($fasta_hsh->{$genome}, 0, $print_length, "")) {
	    $formatted_seq .= "$chunk\n";
        }

	print $fh $header;
	print $fh $formatted_seq;
	
    }
    
    return $filename;
}

sub print_snv_matrix{

    my ($hsh,$genome_header,$order) = @_;
    open(my $fh, ">", "snv_matrix.txt");

    my $genome_fasta;

    #print headers
    print $fh "POS\tREF\tNS\t"; 
    print $fh join("\t",@$genome_header);
    print $fh "\n";

    #print in same order as original vcf file
    foreach my $pos (@$pos_order){

	my($chrom,$position) = split(/:/,$pos);

	if(exists $hsh->{$chrom}->{$position}){
	
	    print $fh "$position\t$hsh->{$chrom}->{$position}->{ref}\t$hsh->{$chrom}->{$position}->{ns}";
	    
	    my @genome_values = @{$hsh->{$chrom}->{$position}->{genomes}};
	    my @alt_values =  split(",",$hsh->{$chrom}->{$position}->{alt});
	    
	    for (my $i = 0; $i < scalar @genome_values; $i++){
		my $base;

		if($genome_values[$i] eq '0'){
		    $base = $hsh->{$chrom}->{$position}->{ref};
		}elsif($genome_values[$i] eq '1'){
		    $base = $alt_values[0];
		}elsif($genome_values[$i] eq '2'){
		    $base = $alt_values[1];
		}elsif($genome_values[$i] eq '3'){
		    $base = $alt_values[2];
		}else{
		    $base =  "N";
		}

		print $fh "\t$base";
		$genome_fasta->{@$genome_header[$i]} .= $base;

	    }
	    
	    print $fh "\n";
	    
	}
	
    }
    
    return $genome_fasta;
}

sub parse_vcf{

    my $file = shift;
    my $hsh;
    my $header_index = 9; #Where the non set header columns(genomes) start
    my @genomes; #stores the genomes in analysis
    my @pos_order;

    open(my $fh, "<", $file);

    while(<$fh>){

	my $line = $_;
	$line =~ s/\s+$//; #remove trailing whitespace

	unless($line =~ /^##/){ #ignore first lines

	    if($line =~ /^#{1,1}/){ #pull header/genome information

		my @headers = split(/\t/,$line);

		#Store an array of the genomes seperately
		for(my $i=$header_index; $i < scalar @headers; $i++){
		    push(@genomes,$headers[$i]);
		}

	    }else{

		#Splice out only the genome data columns
		my @lvalues = split(/\t/,$line);
		my @genome_values = splice(@lvalues,9);

		#Determine the number of genomes with no data
		my $no_data_count = (join(",",@genome_values) =~ tr/.//);
		my $ns_number = scalar @genomes - $no_data_count;
		my $perc = ($ns_number / scalar @genomes) * 100;

		#Only store rows that have genome data for at least the
		#given threshold limit
		if($perc >= $THRESHOLD){
		    $hsh->{$lvalues[0]}->{$lvalues[1]}->{genomes} = \@genome_values;
		    $hsh->{$lvalues[0]}->{$lvalues[1]}->{ref} = $lvalues[3];
		    $hsh->{$lvalues[0]}->{$lvalues[1]}->{alt} = $lvalues[4];
		    $hsh->{$lvalues[0]}->{$lvalues[1]}->{ns} = $ns_number;
		}

		#Store order the rows were looped through to perserve order for output
		push(@pos_order, "$lvalues[0]:$lvalues[1]");

	    }	 
   
	}

    }

   return($hsh,\@genomes,\@pos_order);
}
    
sub filter_vcf{
    
    my ($file,$positions) = @_;

    my ($name,$path,$suffix) = fileparse($file);
    my $outfile = "parsed." . $name;

    if(exists $scripts->{vcftools_exe}){
	my $vcf_exe = $scripts->{vcftools_exe};
	my $cmd = $vcf_exe . " --vcf $file";
	$cmd .= " --positions $positions";
	$cmd .= " --out $outfile";
	$cmd .= " --recode";
	
	system($cmd) == 0 || die("ERROR running $cmd");
    }else{
	die("ERROR: No location found in config file for vcftools");
    }
    return "$outfile.recode.vcf";
}

sub check_params{

    my $errors = "";

    if(!($opts{vcf_file} || $opts{position} || $opts{threshold} || $opts{config})){

	$errors .= "USAGE: vcf_snv.pl -v kSNP VCF File -p recombinant region coordinate file  -t VCF File filtering threshold -c config.file [--help]\n";

    }else{
	
	$errors .= "File does not exist or is size zero: $opts{vcf_file}\n" unless(-s $opts{vcf_file});
	
	$THRESHOLD = $opts{threshold} // "90";
	
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
	
    }
    
    my $output = $opts{output} // cwd;

    die($errors) if $errors;

    return($output,$scripts);
}
