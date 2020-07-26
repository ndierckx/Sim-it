#!/usr/bin/env perl
######################################################
#         SOFTWARE COPYRIGHT NOTICE AGREEMENT        #
#  Copyright (C) {2019-2020}  {Nicolas Dierckxsens}  #
#              All Rights Reserved                   #
#         See file LICENSE for details.              #
######################################################
#           Sim-it - The Sequence Simulator
#           nicolasdierckxsens@hotmail.com
use strict;
use Getopt::Long;
use Time::HiRes qw(time);

print "\n\n-----------------------------------------------";
print "\nSim-it\n";
print "Version 1.1\n";
print "Author: Nicolas Dierckxsens, (c) 2020\n";
print "-----------------------------------------------\n\n";

my $reference = "";
my $ref_database = "";
my $project = "";
 #NC_000001.11
my $kmer = '30';
my $chromosome = "";
my $start_end = '250';
my $detection_limit = '30';
my $fraction = '2000';
my $left = '40';
my $right = '30';
my $reference_size = '0';
my $count_contigs = '1';
my $length_fasta_line = "";
my $SVs_whitin_1_line = "";
my $SVs_whitin_1_line_length = "";
my %SV;
undef %SV;

#error_profile-----
my $AAGA = '0';
my $AAGC = '0';
my $AAGT = '0';
my $AAGG = '0';
#------------------

#my $time_np = '0';
#my $time_np1 = '0';
#my $time_np2 = '0';
#my $time_np3 = '0';
#my $time_sv = '0';

my $VCF_input = "";
my %VCF_input;
undef %VCF_input;
my %VCF_input2;
undef %VCF_input2;
my %VCF_output;
undef %VCF_output;
my $next_pos = "";
my $next_type = "";
my $next_length = "";
my $next_hap = "";
my $next_seq = "";
my $VCF_input_now = "";
my $next_length_minus = '0';
my $next_length_minus_save = '0';
my $hap = "";
my $SEQ = "";
my %sequences;
undef %sequences;
my %sequences_final;
undef %sequences_final;
my %sequences_foreign;
undef %sequences_foreign;
my %foreign_contigs;
undef %foreign_contigs;
my $sequences_foreign = "";

my $DEL_input = "";
my $INS_input = "";
my $DUP_input = "";
my $INV_input = "";
my $TRA_input = "";
my $DEL_range = "";
my $INS_range = "";
my $INV_range = "";
my $DUP_range = "";
my $DUP_copies = "";
my $TRA_range = "";
my $CSUB_input = "";
my $CSUB_range = "";
my $IDUP_input = "";
my $IDUP_range = "";
my $heterozygosity = "";

my $haplotype1 = "";
my $haplotype2 = "";
my $haplotype1_print = "";
my $haplotype2_print = "";
my $finish_var = "";
my $finish_var_np = "";
my $variation_haplo = "";
my $haplo_merged = "";
my $haplo_merged_print = "";
my $ref_haplo = "";
my $pos_hap1 = '0';
my $pos_hap2 = '0';
my $pos_hap_merged = '0';
my $new_contig = "";
my $shortest_seq = '0';
my $shortest_seq1 = '0';
my $shortest_seq2 = '0';
my $seq_done = '0';
my $seq_done1 = '0';
my $seq_done2 = '0';

my $NP_coverage = '0';
my $NP_coverage1 = '0';
my $NP_coverage2 = '0';
my $NP_range = "";
my $NP_average = "";
my $NP_accuracy = "";
my $NP_error_profile = "";

my %graph_DEL;
undef %graph_DEL;
my %graph_INS;
undef %graph_INS;
my %graph_INV2;
undef %graph_INV2;

my $config = "";
my $output = "";
my $batch_file;
my $last_var = "";
my $test = "";
my $ambigious_replace = "";
my @nucs = ("A","C","T","G");

GetOptions (
            "c=s" => \$config,
            "o=s" => \$output,
            ) or die "Incorrect usage!\n";

open(CONFIG, $config) or die "Error: Can't open the configuration file, please check the manual!\n\nUsage: perl Sim-it.pl -c config.txt -o output_path\n\n";

my $last_output = substr $output, -1, 1;
if ($last_output ne "/" && $output ne "")
{
    $output .= "/";
}
while (my $line = <CONFIG>)
{
    chomp($line);
    if ($line =~ m/.*Project name\s+\=\s+(.*?)(Reference sequence.*)*$/)
    {
        $project = $1;
        chomp $project;
        my $project_tmp = $project;
        my $ggg;
        if ($project =~ m/batch\:(.*)/)
        {
            my $batch_file_tmp = $1;
            if ($batch_file eq "")
            {
                $batch_file = $batch_file_tmp;
                print "Batch file detected...\n\n";
                open(BATCH, $batch_file) or die "Error: $!\nCan't open the batch file, please check the manual!\n";
                $ggg = "yes";
            }
            while (my $line = <BATCH>)
            {
                $project = $line;
                chomp $project;
                last;
            }
            if ($project eq $project_tmp || $project eq "")
            {
                goto EXIT;
            }
            elsif ($ggg ne "yes")
            {
                    print "\n\n------------------------------\n------------------------------\n";
                    print "        NEXT SAMPLE\n";
                    print "------------------------------\n------------------------------\n\n\n";
            }
        }
    }
    if ($line =~ m/.*Reference sequence\s+\=\s+(.*?)(Replace ambiguous nts.*)*$/)
    {       
        $reference = $1;
        chomp $reference;
        if ($reference eq "batch" && $batch_file ne "")
        {
            while (my $line = <BATCH>)
            {
                $reference = $line;
                chomp $reference;
                last;
            }
        }   
    }
    if ($line =~ m/.*Replace ambiguous nts\(N\)\s+\=\s+(.*?)(Structural variation:.*)*$/)
    {
        $ambigious_replace = $1;
        chomp $ambigious_replace;
        if ($ambigious_replace eq "batch" && $batch_file ne "")
        {
            while (my $line = <BATCH>)
            {
                $ambigious_replace = $line;
                chomp $ambigious_replace;
                last;
            }
        }
    }
    if ($line =~ m/.*VCF input\s+\=\s+(.*?)(Foreign sequences.*)*$/)
    {
        $VCF_input = $1;
        chomp $VCF_input;
        if ($VCF_input eq "batch" && $batch_file ne "")
        {
            while (my $line = <BATCH>)
            {
                $VCF_input = $line;
                chomp $VCF_input;
                last;
            }
        }   
    }
    if ($line =~ m/.*Foreign sequences\s+\=\s+(.*?)(Deletions.*)*$/)
    {
        $sequences_foreign = $1;
        chomp $sequences_foreign;
        if ($sequences_foreign eq "batch" && $batch_file ne "")
        {
            while (my $line = <BATCH>)
            {
                $sequences_foreign = $line;
                chomp $sequences_foreign;
                last;
            }
        }   
    }
    if ($line =~ m/.*Deletions\s+\=\s+(.*?)(Length \(bp\).*)*$/)
    {
        $DEL_input = $1;
        chomp $DEL_input;
        $last_var = "DEL";
        if ($DEL_input eq "batch" && $batch_file ne "")
        {
            while (my $line = <BATCH>)
            {
                $DEL_input = $line;
                chomp $DEL_input;
                last;
            }
        }   
    }
    if ($line =~ m/.*Length \(bp\)\s+\=\s+(.*?)(Insertions.*)*$/ && $last_var eq "DEL")
    {
        $DEL_range = $1;
        chomp $DEL_range;
        if ($DEL_range eq "batch" && $batch_file ne "")
        {
            while (my $line = <BATCH>)
            {
                $DEL_range = $line;
                chomp $DEL_range;
                last;
            }
        }   
    }
    if ($line =~ m/.*Insertions\s+\=\s+(.*?)(Length \(bp\).*)*$/)
    {
        $INS_input = $1;
        chomp $INS_input;
        $last_var = "INS";
        if ($INS_input eq "batch" && $batch_file ne "")
        {
            while (my $line = <BATCH>)
            {
                $INS_input = $line;
                chomp $INS_input;
                last;
            }
        }   
    }
    if ($line =~ m/.*Length \(bp\)\s+\=\s+(.*?)(Tandem duplications.*)*$/ && $last_var eq "INS")
    {
        $INS_range = $1;
        chomp $INS_range;
        if ($INS_range eq "batch" && $batch_file ne "")
        {
            while (my $line = <BATCH>)
            {
                $INS_range = $line;
                chomp $INS_range;
                last;
            }
        }   
    }
    if ($line =~ m/.*Tandem duplications\s+\=\s+(.*?)(Length \(bp\).*)*$/)
    {
        $DUP_input = $1;
        chomp $DUP_input;
        $last_var = "DUP";
        if ($DUP_input eq "batch" && $batch_file ne "")
        {
            while (my $line = <BATCH>)
            {
                $DUP_input = $line;
                chomp $DUP_input;
                last;
            }
        }   
    }
    if ($line =~ m/.*Length \(bp\)\s+\=\s+(.*?)(Copies.*)*$/ && $last_var eq "DUP")
    {
        $DUP_range = $1;
        chomp $DUP_range;
        if ($DUP_range eq "batch" && $batch_file ne "")
        {
            while (my $line = <BATCH>)
            {
                $DUP_range = $line;
                chomp $DUP_range;
                last;
            }
        }   
    }
    if ($line =~ m/.*Copies\s+\=\s+(.*?)(Inversions.*)*$/)
    {
        $DUP_copies = $1;
        chomp $DUP_copies;
        if ($DUP_copies eq "batch" && $batch_file ne "")
        {
            while (my $line = <BATCH>)
            {
                $DUP_copies = $line;
                chomp $DUP_copies;
                last;
            }
        }   
    }
    if ($line =~ m/.*Inversions\s+\=\s+(.*?)(Length \(bp\).*)*$/)
    {
        $INV_input = $1;
        chomp $INV_input;
        $last_var = "INV";
        if ($INV_input eq "batch" && $batch_file ne "")
        {
            while (my $line = <BATCH>)
            {
                $INV_input = $line;
                chomp $INV_input;
                last;
            }
        }   
    }  
    if ($line =~ m/.*Length \(bp\)\s+\=\s+(.*?)(Complex substitutions.*)*$/ && $last_var eq "INV")
    {
        $INV_range = $1;
        chomp $INV_range;
        if ($INV_range eq "batch" && $batch_file ne "")
        {
            while (my $line = <BATCH>)
            {
                $INV_range = $line;
                chomp $INV_range;
                last;
            }
        }   
    }
    if ($line =~ m/.*Complex substitutions\s+\=\s+(.*?)(Length \(bp\).*)*$/)
    {
        $CSUB_input = $1;
        chomp $CSUB_input;
        $last_var = "CSUB";
        if ($CSUB_input eq "batch" && $batch_file ne "")
        {
            while (my $line = <BATCH>)
            {
                $CSUB_input = $line;
                chomp $CSUB_input;
                last;
            }
        }   
    }
    if ($line =~ m/.*Length \(bp\)\s+\=\s+(.*?)(Complex SVs.*)*$/ && $last_var eq "CSUB")
    {
        $CSUB_range = $1;
        chomp $CSUB_range;
        if ($CSUB_range eq "batch" && $batch_file ne "")
        {
            while (my $line = <BATCH>)
            {
                $CSUB_range = $line;
                chomp $CSUB_range;
                last;
            }
        }   
    }
    if ($line =~ m/.*Inverted duplications\s+\=\s+(.*?)(Length \(bp\).*)*$/)
    {
        $IDUP_input = $1;
        chomp $IDUP_input;
        $last_var = "IDUP";
        if ($IDUP_input eq "batch" && $batch_file ne "")
        {
            while (my $line = <BATCH>)
            {
                $IDUP_input = $line;
                chomp $IDUP_input;
                last;
            }
        }   
    }
    if ($line =~ m/.*Length \(bp\)\s+\=\s+(.*?)(Heterozygosity.*)*$/ && $last_var eq "IDUP")
    {
        $IDUP_range = $1;
        chomp $IDUP_range;
        if ($IDUP_range eq "batch" && $batch_file ne "")
        {
            while (my $line = <BATCH>)
            {
                $IDUP_range = $line;
                chomp $IDUP_range;
                last;
            }
        }   
    }
    if ($line =~ m/.*Heterozygosity\s+\=\s+(.*?)(Long Read simulation:.*)*$/)
    {
        $heterozygosity = $1;
        chomp $heterozygosity;
        if ($heterozygosity eq "batch" && $batch_file ne "")
        {
            while (my $line = <BATCH>)
            {
                $heterozygosity = $line;
                chomp $heterozygosity;
                last;
            }
        }   
    }
    
#Read simulation input----------------------------------------------------------------------------------------

    if ($line =~ m/.*Coverage\s+\=\s+(.*?)(Average length.*)*$/)
    {
        $NP_coverage = $1;
        chomp $NP_coverage;
        if ($NP_coverage eq "batch" && $batch_file ne "")
        {
            while (my $line = <BATCH>)
            {
                $NP_coverage = $line;
                chomp $NP_coverage;
                last;
            }
        }   
    }
    if ($line =~ m/.*Median length\s+\=\s+(.*?)(Length range.*)*$/)
    {
        $NP_average = $1;
        chomp $NP_average;
        if ($NP_average eq "batch" && $batch_file ne "")
        {
            while (my $line = <BATCH>)
            {
                $NP_average = $line;
                chomp $NP_average;
                last;
            }
        }   
    }
    if ($line =~ m/.*Length range\s+\=\s+(.*?)(Accuracy.*)*$/)
    {
        $NP_range = $1;
        chomp $NP_range;
        if ($NP_range eq "batch" && $batch_file ne "")
        {
            while (my $line = <BATCH>)
            {
                $NP_range = $line;
                chomp $NP_range;
                last;
            }
        }   
    }
    if ($line =~ m/.*Accuracy\s+\=\s+(.*?)(Error profile.*)*$/)
    {
        $NP_accuracy = $1;
        chomp $NP_accuracy;
        if ($NP_accuracy eq "batch" && $batch_file ne "")
        {
            while (my $line = <BATCH>)
            {
                $NP_accuracy = $line;
                chomp $NP_accuracy;
                last;
            }
        }   
    }
    if ($line =~ m/.*Error profile\s+\=\s+(.*?)$/)
    {
        $NP_error_profile = $1;
        chomp $NP_error_profile;
        if ($NP_error_profile eq "batch" && $batch_file ne "")
        {
            while (my $line = <BATCH>)
            {
                $NP_error_profile = $line;
                chomp $NP_error_profile;
                last;
            }
        }
        last;
    }
}
close CONFIG;

if ($reference eq "")
{
    die "\nYou need to give a reference sequence!\n";
}

#Check input parameters----------------------------------------------------------------------------

my $SV_input = "";
if ($VCF_input ne "" || $DEL_input > '0' || $INS_input > '0' || $INV_input > '0' || $TRA_input > '0' || $DUP_input > '0' || $CSUB_input > '0' || $IDUP_input > '0')
{
    $SV_input = "yes";
}
if ($SV_input eq "")
{
    $heterozygosity = "";
}
if ($SV_input ne "" && $heterozygosity eq "")
{
    $heterozygosity = '60%';
}

if ($DEL_input > 0 && $DEL_range ne "" && $DEL_range =~ m/^\d+\-\d+$/)
{}
elsif ($DEL_input > 0 && $DEL_range ne "")
{
    die "\n\nThe length parameter of the Deletions input is not in a correct format should be for example (30-150000), $!\n";
}
elsif ($DEL_input > 0)
{
    $DEL_range = "30-150000";
}
if ($INS_input > 0 && $INS_range ne "" && $INS_range =~ m/^\d+\-\d+$/)
{}
elsif ($INS_input > 0 && $INS_range ne "")
{
    die "\n\nThe length parameter of the Inserions input is not in a correct format should be for example (30-100000), $!\n";
}
elsif ($INS_input > 0)
{
    $INS_range = "30-100000";
}
if ($NP_error_profile eq "" && $NP_coverage > 0)
{
    die "\n\nAn error profile should be given$!\n";
}
#----------------------------------------------------------------------------------------------------

if ($SV_input ne "")
{
    my $output_DEL  = $output."graph_DEL_".$project.".txt";
    open(OUTPUT_DEL, ">" .$output_DEL) or die "\nCan't open file $output_DEL, $!\n";
    my $output_INS  = $output."graph_INS_".$project.".txt";
    open(OUTPUT_INS, ">" .$output_INS) or die "\nCan't open file $output_INS, $!\n";
    my $output_INV2  = $output."graph_INV_".$project.".txt";
    open(OUTPUT_INV2, ">" .$output_INV2) or die "\nCan't open file $output_INV2, $!\n";
}

my $output_LOG  = $output."log_".$project.".txt";
open(OUTPUT_LOG, ">" .$output_LOG) or die "\nCan't open file $output_LOG, $!\n";


my $output_vcf = $output.$project.".vcf";
my $output_ref = $output.$project.".fasta";
my $output_hap1 = $output.$project."_haplotype1.fasta";
my $output_hap2 = $output.$project."_haplotype2.fasta";


print "\n\n";
print "---------------------------------CONFIG INPUT-----------------------------------\n\n";
print "\nInput parameters from the configuration file:   *** Verify if everything is correct ***\n\n\n";
print "Project:\n";
print "-----------------------\n";
print "Project name             = ".$project."\n";
print "Reference sequence       = ".$reference."\n";
print "Replace ambiguous nts(N) = ".$ambigious_replace."\n\n\n";

print "Structural variation:\n";
print "-----------------------\n";
print "VCF input                = ".$VCF_input."\n";
print "Foreign sequences        = ".$sequences_foreign."\n\n";

print "Deletions                = ".$DEL_input."\n";
print "Length (bp)              = ".$DEL_range."\n\n";

print "Insertions               = ".$INS_input."\n";
print "Length (bp)              = ".$INS_range."\n\n";

print "Tandem duplications      = ".$DUP_input."\n";
print "Length (bp)              = ".$DUP_range."\n";
print "Copies                   = ".$DUP_copies."\n\n";

print "Inversions               = ".$INV_input."\n";
print "Length (bp)              = ".$INV_range."\n\n";

print "Complex substitutions    = ".$CSUB_input."\n";
print "Length (bp)              = ".$CSUB_range."\n\n";

print "Inverted duplications    = ".$IDUP_input."\n";
print "Length (bp)              = ".$IDUP_range."\n\n";

print "Heterozygosity           = ".$heterozygosity."\n\n\n";

print "Long Read simulation:\n";
print "-----------------------\n";
print "Coverage                 = ".$NP_coverage."\n";
print "Median length            = ".$NP_average."\n";
print "Length range             = ".$NP_range."\n";
print "Accuracy                 = ".$NP_accuracy."\n";
print "Error profile            = ".$NP_error_profile."\n\n\n";

print OUTPUT_LOG "\n\n";
print OUTPUT_LOG "---------------------------------CONFIG INPUT-----------------------------------\n\n";
print OUTPUT_LOG "\nInput parameters from the configuration file:   *** Verify if everything is correct ***\n\n\n";
print OUTPUT_LOG "Project:\n";
print OUTPUT_LOG "-----------------------\n";
print OUTPUT_LOG "Project name             = ".$project."\n";
print OUTPUT_LOG "Reference                = ".$reference."\n";
print OUTPUT_LOG "Replace ambiguous nts(N) = ".$ambigious_replace."\n\n\n";

print OUTPUT_LOG "Structural variation:\n";
print OUTPUT_LOG "-----------------------\n";
print OUTPUT_LOG "VCF input                = ".$VCF_input."\n";
print OUTPUT_LOG "Foreign sequences        = ".$sequences_foreign."\n\n";

print OUTPUT_LOG "Deletions                = ".$DEL_input."\n";
print OUTPUT_LOG "Length (bp)              = ".$DEL_range."\n\n";

print OUTPUT_LOG "Insertions               = ".$INS_input."\n";
print OUTPUT_LOG "Length (bp)              = ".$INS_range."\n\n";

print OUTPUT_LOG "Tandem duplications      = ".$DUP_input."\n";
print OUTPUT_LOG "Length (bp)              = ".$DUP_range."\n";
print OUTPUT_LOG "Copies                   = ".$DUP_copies."\n\n";

print OUTPUT_LOG "Inversions               = ".$INV_input."\n";
print OUTPUT_LOG "Length (bp)              = ".$INV_range."\n\n";

print OUTPUT_LOG "Complex substitutions    = ".$CSUB_input."\n";
print OUTPUT_LOG "Length (bp)              = ".$CSUB_range."\n\n";

print OUTPUT_LOG "Inverted duplications    = ".$IDUP_input."\n";
print OUTPUT_LOG "Length (bp)              = ".$IDUP_range."\n\n";

print OUTPUT_LOG "Heterozygosity           = ".$heterozygosity."\n\n\n";

print OUTPUT_LOG "Long Read simulation:\n";
print OUTPUT_LOG "-----------------------\n";
print OUTPUT_LOG "Coverage                 = ".$NP_coverage."\n";
print OUTPUT_LOG "Median length            = ".$NP_average."\n";
print OUTPUT_LOG "Length range             = ".$NP_range."\n";
print OUTPUT_LOG "Accuracy                 = ".$NP_accuracy."\n";
print OUTPUT_LOG "Error profile            = ".$NP_error_profile."\n\n";


if ($heterozygosity =~ m/^(\d+)%.*$/)
{
    $heterozygosity = $1;
    $heterozygosity /= 100;
}
else
{
    $heterozygosity = "no";
}
if ($heterozygosity eq '0')
{
    $heterozygosity = "no";
}

my $VCF_input_DEL = '0';
my $VCF_input_INS = '0';
my $VCF_input_INV = '0';
my $VCF_input_IDUP = '0';
my $VCF_input_DUP = '0';
my $VCF_input_TRA = '0';
my $VCF_input_CSUB = '0';
my $VCF_input_DEL_min = '100000000000000000';
my $VCF_input_INS_min = '100000000000000000';
my $VCF_input_INV_min = '100000000000000000';
my $VCF_input_IDUP_min = '100000000000000000';
my $VCF_input_DUP_min = '100000000000000000';
my $VCF_input_TRA_min = '100000000000000000';
my $VCF_input_CSUB_min = '100000000000000000';
my $VCF_input_DEL_max = '0';
my $VCF_input_INS_max = '0';
my $VCF_input_INV_max = '0';
my $VCF_input_IDUP_max = '0';
my $VCF_input_DUP_max = '0';
my $VCF_input_TRA_max = '0';
my $VCF_input_CSUB_max = '0';

my %TRA_hap1;
undef %TRA_hap1;
my %TRA_hap2;
undef %TRA_hap2;
my %TRA_hap_second;
undef %TRA_hap_second;
my %TRA_sequences;
undef %TRA_sequences;

if ($VCF_input ne "")
{
    if ($sequences_foreign ne "")
    {
        open (INPUT_SEQ, $sequences_foreign) or die "Can't open file $sequences_foreign, $!\n";
        my $id = "";
        while (my $line_seq = <INPUT_SEQ>)
        {
            chomp $line_seq;
            my $first_tmp = substr $line_seq, 0, 1, "";
            if ($first_tmp eq ">")
            {
                $id = $line_seq;
            }
            elsif ($line_seq ne "" && $line_seq ne " ")
            {
                $foreign_contigs{$id} = $line_seq;
            }  
        }
        close INPUT_SEQ;
    }
    open (INPUT_VCF, $VCF_input) or die "Can't open file $VCF_input, $!\n";
    while (my $line = <INPUT_VCF>)
    {
        if ($line =~ m/(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)*\s+(\S+)*.*/)
        {
            my $chr_tmp = $1;
            my $start_pos_tmp = $2;
            my $length_tmp = $3;
            my $type_tmp = $4;
            my $hap_tmp = $5;
            my $SEQ_tmp = $6;

            if ($SEQ_tmp ne "")
            {
                if (exists($foreign_contigs{$SEQ_tmp}))
                {}
                else
                {
                    $sequences{$SEQ_tmp} = undef;
                }
            }
            if ($type_tmp eq "DEL")
            {
                $VCF_input_DEL++;
                if ($length_tmp < $VCF_input_DEL_min)
                {
                    $VCF_input_DEL_min = $length_tmp;
                }
                if ($length_tmp > $VCF_input_DEL_max)
                {
                    $VCF_input_DEL_max = $length_tmp;
                }
            }
            if ($type_tmp eq "INS")
            {
                $VCF_input_INS++;
                if ($length_tmp < $VCF_input_INS_min)
                {
                    $VCF_input_INS_min = $length_tmp;
                }
                if ($length_tmp > $VCF_input_INS_max)
                {
                    $VCF_input_INS_max = $length_tmp;
                }
            }
            if ($type_tmp eq "DUP")
            {
                $VCF_input_DUP++;
                if ($length_tmp < $VCF_input_DUP_min)
                {
                    $VCF_input_DUP_min = $length_tmp;
                }
                if ($length_tmp > $VCF_input_DUP_max)
                {
                    $VCF_input_DUP_max = $length_tmp;
                }
            }
            if ($type_tmp eq "INV")
            {
                $VCF_input_INV++;
                if ($length_tmp < $VCF_input_INV_min)
                {
                    $VCF_input_INV_min = $length_tmp;
                }
                if ($length_tmp > $VCF_input_INV_max)
                {
                    $VCF_input_INV_max = $length_tmp;
                }
            }
            if ($type_tmp eq "CSUB")
            {
                $VCF_input_CSUB++;
                if ($length_tmp < $VCF_input_CSUB_min)
                {
                    $VCF_input_CSUB_min = $length_tmp;
                }
                if ($length_tmp > $VCF_input_CSUB_max)
                {
                    $VCF_input_CSUB_max = $length_tmp;
                }
            }
            if ($type_tmp eq "TRA")
            {
                $VCF_input_TRA++;
                if ($length_tmp < $VCF_input_TRA_min)
                {
                    $VCF_input_TRA_min = $length_tmp;
                }
                if ($length_tmp > $VCF_input_TRA_max)
                {
                    $VCF_input_TRA_max = $length_tmp;
                }
                my @chr_tmp = split /\:/, $chr_tmp;
                if (exists($TRA_hap1{$chr_tmp[1]}))
                {
                    die "\n\nERROR: ".$line." You can only have one translocation per chromosome\n";
                }
                if (exists($TRA_hap2{$chr_tmp[1]}))
                {
                    die "\n\nERROR: ".$line." You can only have one translocation per chromosome\n";
                }               
                if ($hap_tmp ne '1')
                {               
                    $TRA_hap1{$chr_tmp[1]}{$chr_tmp[0]} = $length_tmp;
                    $TRA_hap1{$chr_tmp[0]}{$chr_tmp[1]} = $length_tmp;
                    $TRA_hap_second{$chr_tmp[1]} = undef;
                }
                if ($hap_tmp ne '0')
                {
                    $TRA_hap2{$chr_tmp[1]}{$chr_tmp[0]} = $length_tmp;
                    $TRA_hap2{$chr_tmp[0]}{$chr_tmp[1]} = $length_tmp;
                    $TRA_hap_second{$chr_tmp[1]} = undef;
                }
                
            }
            if ($type_tmp eq "IDUP" || $type_tmp eq "INV-DUP")
            {
                $VCF_input_IDUP++;
                if ($length_tmp < $VCF_input_IDUP_min)
                {
                    $VCF_input_IDUP_min = $length_tmp;
                }
                if ($length_tmp > $VCF_input_IDUP_max)
                {
                    $VCF_input_IDUP_max = $length_tmp;
                }
            }
        }
    }
    print "\n----------------------------------VCF INPUT------------------------------------\n\n";
    print OUTPUT_LOG "\n----------------------------------VCF INPUT------------------------------------\n\n";
    if ($VCF_input_DEL > 0)
    {
        print "Deletions             = ".$VCF_input_DEL."\n";
        print "Length (bp)           = ".$VCF_input_DEL_min."-".$VCF_input_DEL_max."\n\n";
        print OUTPUT_LOG "Deletions             = ".$VCF_input_DEL."\n";
        print OUTPUT_LOG "Length (bp)           = ".$VCF_input_DEL_min."-".$VCF_input_DEL_max."\n\n";
    }
    if ($VCF_input_INS > 0)
    {
        print "Insertions            = ".$VCF_input_INS."\n";
        print "Length (bp)           = ".$VCF_input_INS_min."-".$VCF_input_INS_max."\n\n";
        print OUTPUT_LOG "Insertions            = ".$VCF_input_INS."\n";
        print OUTPUT_LOG "Length (bp)           = ".$VCF_input_INS_min."-".$VCF_input_INS_max."\n\n";
    }
    if ($VCF_input_DUP > 0)
    {
        print "Tandem duplications   = ".$VCF_input_DUP."\n";
        print "Length (bp)           = ".$VCF_input_DUP_min."-".$VCF_input_DUP_max."\n\n";
        print OUTPUT_LOG "Tandem duplications   = ".$VCF_input_DUP."\n";
        print OUTPUT_LOG "Length (bp)           = ".$VCF_input_DUP_min."-".$VCF_input_DUP_max."\n\n";
    }
    if ($VCF_input_INV > 0)
    {
        print "Inversions            = ".$VCF_input_INV."\n";
        print "Length (bp)           = ".$VCF_input_INV_min."-".$VCF_input_INV_max."\n\n";
        print OUTPUT_LOG "Inversions            = ".$VCF_input_INV."\n";
        print OUTPUT_LOG "Length (bp)           = ".$VCF_input_INV_min."-".$VCF_input_INV_max."\n\n";
    }if ($VCF_input_TRA > 0)
    {
        print "Translocations        = ".$VCF_input_TRA."\n";
        print "Length (bp)           = ".$VCF_input_TRA_min."-".$VCF_input_TRA_max."\n\n";
        print OUTPUT_LOG "Translocations        = ".$VCF_input_TRA."\n";
        print OUTPUT_LOG "Length (bp)           = ".$VCF_input_TRA_min."-".$VCF_input_TRA_max."\n\n";    
    }
    if ($VCF_input_CSUB > 0)
    {
        print "Complex substitutions = ".$VCF_input_CSUB."\n";
        print "Length (bp)           = ".$VCF_input_CSUB_min."-".$VCF_input_CSUB_max."\n\n";
        print OUTPUT_LOG "Complex substitutions = ".$VCF_input_CSUB."\n";
        print OUTPUT_LOG "Length (bp)           = ".$VCF_input_CSUB_min."-".$VCF_input_CSUB_max."\n\n";
    }
    if ($VCF_input_IDUP > 0)
    {
        print "Inverted duplications = ".$VCF_input_IDUP."\n";
        print "Length (bp)           = ".$VCF_input_IDUP_min."-".$VCF_input_IDUP_max."\n\n";
        print OUTPUT_LOG "Inverted duplications = ".$VCF_input_IDUP."\n";
        print OUTPUT_LOG "Length (bp)           = ".$VCF_input_IDUP_min."-".$VCF_input_IDUP_max."\n\n";
    }
    close INPUT_VCF;
}

print "--------------------------------------------------------------------------------\n\n";
print OUTPUT_LOG "--------------------------------------------------------------------------------\n\n";

#Generate vcf file---------------------------------------------------------------------------------------------------------------------------------------
if ($SV_input eq "yes")
{
    open(OUTPUT_VCF, ">" .$output_vcf) or die "Can't open variance file $output_vcf, $!\n";
    
    my ($wday, $mon, $mday, $hour, $min, $sec, $year) = localtime;
    my @localtime = split / /, localtime;
    my %mon2num = qw(
    Jan 01  Feb 02  Mar 03  Apr 04  May 05  Jun 06
    Jul 07  Aug 08  Sep 09  Oct 10 Nov 11 Dec 12
    );
    my $month = $localtime[1];
    if (exists($mon2num{$localtime[1]}))
    {
       $month = $mon2num{$localtime[1]};
    }
    print OUTPUT_VCF "##fileformat=VCFv4.0\n";
    print OUTPUT_VCF "##fileDate=".$localtime[4].$month.$localtime[2]."\n";
    print OUTPUT_VCF "##reference=".$reference."\n";
    #print OUTPUT_VCF "##INFO=<ID=tPOS,Number=1,Type=Float,Description=\"Allele Frequency\">\n";
    print OUTPUT_VCF "#CHR\tPOS\tSVLENGTH\tTYPE\tVARHAP\tVARSEQ\n";

#--------------------------------------------------------------------------------------------------------------------------------------------------------
    open(OUTPUT_REF, ">" .$output_ref) or die "Can't open fasta output file $output_ref, $!\n";
    if ($heterozygosity ne "no" || $VCF_input ne "")
    {
        open(OUTPUT_HAP1, ">" .$output_hap1) or die "Can't open haplotype 1 fasta output file $output_hap1, $!\n";
        open(OUTPUT_HAP2, ">" .$output_hap2) or die "Can't open haplotype 2 fasta output file $output_hap2, $!\n";
    }
}

select(STDERR);
$| = 1;
select(STDOUT); # default
$| = 1;
print "\nScan reference sequence...";

my %contigs;
undef %contigs;
my %contigs_chr;
undef %contigs_chr;

my @TRA_range = split /-/, $TRA_range;
my $TRA_range_low = $TRA_range[0];
my $TRA_range_high = $TRA_range[1];
#my $TRA_interval = int($reference_size/($TRA_input+$VCF_input_TRA+1));
#my $TRA_interval_tmp = $TRA_interval; 
#my $TRA_input_tmp = $TRA_input;
my $TRA_range2 = $TRA_range_high-$TRA_range_low;
my $random_length_TRA = "";
my $translocating = "";
my $TRA_contig_select = "";
my $TRA_contig_select2 = "";
my $TRA_set = '0';
my $TRA_set2 = "";

my $check_zip = substr $reference, -2;
my $check_zip2 = substr $reference, -3;
my $firstLine = "";
my $secondLine = "";
my $thirdLine = "";
my $FILE_REF;
if ($check_zip eq "gz")
{
    open ($FILE_REF, '-|', 'gzip', '-dc', $reference) or die "Can't open file $reference, $!\n";
    $firstLine = <$FILE_REF>;
    chomp $firstLine;
    $secondLine = <$FILE_REF>;
    chomp $secondLine;
}
elsif ($check_zip2 eq "bz2")
{
    open ($FILE_REF, '-|', 'bzip2', '-dc', $reference) or die "Can't open file $reference, $!\n";
    $firstLine = <$FILE_REF>;
    chomp $firstLine;
    $secondLine = <$FILE_REF>;
    chomp $secondLine;
}
else
{
    open($FILE_REF, $reference) or die "\n\nCan't open reference file $reference, $!\n";
    $firstLine = <$FILE_REF>;
    chomp $firstLine;
    $secondLine = <$FILE_REF>;
    chomp $secondLine;
}
close $FILE_REF;
if ($check_zip eq "gz")
{
    open ($FILE_REF, '-|', 'gzip', '-dc', $reference) or die "Can't open file $reference, $!\n";
}
elsif ($check_zip2 eq "bz2")
{
    open ($FILE_REF, '-|', 'bzip2', '-dc', $reference) or die "Can't open file $reference, $!\n";
}
else
{
    open($FILE_REF, $reference) or die "\n\nCan't open reference file $reference, $!\n";
}

my $value_ref2 = "";
my $s1000 = "";
my $pos_last_contig = '0';
my $size_last_contig = '0';
my $last_contig_id = "";
my $chromosome_tmp = "";
my $chr_tmp_prev = "";
my $current_contig_tmp = '0';
my %sequences_chr;
undef %sequences_chr;
my $start_pos_next = "";
my $end_pos_next = "";
my $insert_seq = "";
my $insert_seq_check = "";
my $sequence_final = "";


while (my $line = <$FILE_REF>)
{     
    chomp $line;
    if ($reference_size < 1)
    {
        $last_contig_id = $line;
        $reference_size++;
        $current_contig_tmp++;
        if ($line =~ m/>(\d+|X|Y|MT)\s+dna:chromosome.*/)
        {
            $chromosome_tmp = $1;
        }
        elsif ($line =~ m/>chr(\d+|X|Y|MT)$/)
        {
            $chromosome_tmp = $1;
        }
        elsif ($line =~ m/>.*chromosome\s(\d+|X|Y|MT).*/)
        {
            $chromosome_tmp = $1;
        }
        else
        {
            $chromosome_tmp = substr $line, 1;
        }
        foreach my $seqs (keys %sequences)
        {
            my @sequence = split /,/, $seqs;
            my $pos_seq = $sequence[1];
            my $ch_seq = $sequence[0];
            my $line_chr = substr $line, 1;
            if ($ch_seq eq $chromosome_tmp || $ch_seq eq $line_chr)
            {
                my @pos_seq = split /-/, $pos_seq;
                $sequences_chr{$pos_seq[0]} = $pos_seq[1]."*".$seqs;
            }
        }
        next;
    }
    elsif ($length_fasta_line eq "")
    {
        $length_fasta_line = length($line);
    }
    if ($sequence_final eq "")
    {
        foreach my $start_pos_tmpie (sort {$a <=> $b} keys %sequences_chr)
        {
            $start_pos_next = $start_pos_tmpie;
            my @sequences_chr = split /\*/, $sequences_chr{$start_pos_tmpie};
            $end_pos_next = $sequences_chr[0];
            $sequence_final = $sequences_chr[1];
            last;
        }
    }
    
    my $first = substr $line, 0, 1;
    my $line3 = "";  
    my $end_contig = "";
    if ($first eq '>' || $first eq '@')
    {
        $contigs_chr{$chromosome_tmp} = $size_last_contig;
        if ($line =~ m/>(\d+|X|Y|MT)\s+dna:chromosome.*/)
        {
            $chromosome_tmp = $1;
        }
        elsif ($line =~ m/>chr(\d+|X|Y|MT)$/)
        {
            $chromosome_tmp = $1;
        }
        elsif ($line =~ m/>.*chromosome\s(\d+|X|Y|MT).*/)
        {
            $chromosome_tmp = $1;
        }
        else
        {
            $chromosome_tmp = substr $line, 1;
        }
        $end_contig = "yes";
        $count_contigs++;
        $size_last_contig = $reference_size-$pos_last_contig;
        $pos_last_contig = $reference_size;

        $contigs{$last_contig_id} = $size_last_contig;
        $last_contig_id = $line;
        $current_contig_tmp = '1';
        
        foreach my $seqs (keys %sequences)
        {
            my @sequence = split /,/, $seqs;
            my $pos_seq = $sequence[1];
            my $ch_seq = $sequence[0];
            my $line_chr = substr $line, 1;
            if ($ch_seq eq $chromosome_tmp || $ch_seq eq $line_chr)
            {
                my @pos_seq = split /-/, $pos_seq;
                $sequences_chr{$pos_seq[0]} = $pos_seq[1]."*".$seqs;
            }
        }
        #$line3 = $value_ref2."NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN";
    }
    else
    {
#Make hash from sequences file------------------------------------------------------------------------------------------------------------------        
        my $line_tmp = $line;
        if ($start_pos_next <= $current_contig_tmp+length($line) && $start_pos_next >= $current_contig_tmp)
        { 
            if ($start_pos_next-$current_contig_tmp > 0)
            {
                substr $line_tmp, 0, $start_pos_next-$current_contig_tmp, "";
            }
            $insert_seq_check = "yes";
        }
        if ($current_contig_tmp+length($line) >= $end_pos_next+1 && $insert_seq_check eq "yes")
        {
            if ((($current_contig_tmp+length($line))-$end_pos_next-1) > 0)
            {
                substr $line_tmp, -(($current_contig_tmp+length($line))-$end_pos_next-1), (($current_contig_tmp+length($line))-$end_pos_next-1), "";
            }
            $insert_seq_check = "yes2";
        }
        if ($insert_seq_check ne "")
        {
            $insert_seq .= $line_tmp;
            if ($insert_seq_check eq "yes2")
            {
                $insert_seq_check = "";
                $sequences_final{$sequence_final} = $insert_seq;
                $insert_seq = "";
                delete $sequences_chr{$start_pos_next};
                $start_pos_next = "";
                $end_pos_next = "";
                $sequence_final = "";
            }
        }
#------------------------------------------------------------------------------------------------------------------------------------------------
        $reference_size += length($line);
        $current_contig_tmp += length($line);
    }
}
$contigs{$last_contig_id} = $reference_size-$pos_last_contig;
$contigs_chr{$chromosome_tmp} = $reference_size-$pos_last_contig;
close $FILE_REF;


#preprep TRA seq-------------------------------------------------------------------------------------------------------
my $reference_size_tmp = '0';
if ($VCF_input_TRA > 0)
{
    my $TRA_length_tmp = "";
    my $TRA_seq_tmp = "";
    
    if ($check_zip eq "gz")
    {
        open ($FILE_REF, '-|', 'gzip', '-dc', $reference) or die "Can't open file $reference, $!\n";
    }
    elsif ($check_zip2 eq "bz2")
    {
        open ($FILE_REF, '-|', 'bzip2', '-dc', $reference) or die "Can't open file $reference, $!\n";
    }
    else
    {
        open($FILE_REF, $reference) or die "\n\nCan't open reference file $reference, $!\n";
    }
    
    while (my $line = <$FILE_REF>)
    {     
        chomp $line;
        if ($reference_size_tmp < 1)
        {
            $last_contig_id = $line;
            $reference_size_tmp++;
            $current_contig_tmp++;
            if ($line =~ m/>(\d+|X|Y|MT)\s+dna:chromosome.*/)
            {
                $chromosome_tmp = $1;
            }
            elsif ($line =~ m/>chr(\d+|X|Y|MT)$/)
            {
                $chromosome_tmp = $1;
            }
            elsif ($line =~ m/>.*chromosome\s(\d+|X|Y|MT).*/)
            {
                $chromosome_tmp = $1;
            }
            else
            {
                $chromosome_tmp = substr $line, 1;
            }

            $TRA_length_tmp = "";
            $TRA_seq_tmp = "";
            if (exists($TRA_hap1{$chromosome_tmp}))
            {          
                foreach my $line_tmp2 (sort {$a <=> $b} keys %{$TRA_hap1{$chromosome_tmp}})
                {
                    $TRA_length_tmp = $TRA_hap1{$chromosome_tmp}{$line_tmp2};
                }
            }
            if (exists($TRA_hap2{$chromosome_tmp}))
            {             
                foreach my $line_tmp2 (sort {$a <=> $b} keys %{$TRA_hap2{$chromosome_tmp}})
                {
                    $TRA_length_tmp = $TRA_hap2{$chromosome_tmp}{$line_tmp2};
                }
            }
            next;
        }
        elsif ($length_fasta_line eq "")
        {
            $length_fasta_line = length($line);
        }
        
        my $first = substr $line, 0, 1;
        my $line3 = "";  
        my $end_contig = "";
        if ($first eq '>' || $first eq '@')
        {
            $chr_tmp_prev = $chromosome_tmp;
            if ($line =~ m/>(\d+|X|Y|MT)\s+dna:chromosome.*/)
            {
                $chromosome_tmp = $1;
            }
            elsif ($line =~ m/>chr(\d+|X|Y|MT)$/)
            {
                $chromosome_tmp = $1;
            }
            elsif ($line =~ m/>.*chromosome\s(\d+|X|Y|MT).*/)
            {
                $chromosome_tmp = $1;
            }
            else
            {
                $chromosome_tmp = substr $line, 1;
            }
            $last_contig_id = $line;
            $current_contig_tmp = '1';

            if ($TRA_seq_tmp ne "")
            {
                $TRA_sequences{$chr_tmp_prev} = $TRA_seq_tmp;
            }
            $TRA_length_tmp = "";
            $TRA_seq_tmp = "";
            if (exists($TRA_hap1{$chromosome_tmp}))
            {
                foreach my $line_tmp2 (sort {$a <=> $b} keys %{$TRA_hap1{$chromosome_tmp}})
                {
                    $TRA_length_tmp = $TRA_hap1{$chromosome_tmp}{$line_tmp2};
                    if ($TRA_length_tmp > $contigs{$last_contig_id})
                    {
                        die "\n\nERROR: ".$TRA_length_tmp." bp. The translocation length can not be longer than the length of the chromosomes\n";
                    }
                }
            }
            if (exists($TRA_hap2{$chromosome_tmp}))
            {
                foreach my $line_tmp2 (sort {$a <=> $b} keys %{$TRA_hap2{$chromosome_tmp}})
                {
                    $TRA_length_tmp = $TRA_hap2{$chromosome_tmp}{$line_tmp2};
                    if ($TRA_length_tmp > $contigs{$last_contig_id})
                    {
                        die "\n\nERROR: ".$TRA_length_tmp." bp. The translocation length can not be longer than the length of the chromosomes\n";
                    }
                }
            }
        }
        else
        {
            my $line_tmp = $line;

            if ($TRA_length_tmp ne "")
            {
                if ($contigs{$last_contig_id}-$current_contig_tmp <= $TRA_length_tmp)
                {
                    my $first = substr $line, 0, 1;
                    my $last_nuc = "";
                    if ($ambigious_replace eq "yes")
                    {
                        $last_nuc = substr $line, -1, 1;
                    }
                    if (($first eq "N" && $last_nuc eq "N") && $ambigious_replace eq "yes")
                    {
                        my $line_tmp = "";
                        while (length($line_tmp) < length($line))
                        {
                            my $chance_tmp = int(rand(4));
                            my $nuc_tmp = $nucs[$chance_tmp];
                            $line_tmp .= $nuc_tmp; 
                        }
                        $line = $line_tmp;
                    }
                    elsif (($first eq "N" || $last_nuc eq "N") && $ambigious_replace eq "yes")
                    {
                        my $chance_tmp = int(rand(4));
                        my $nuc_tmp = $nucs[$chance_tmp];
                        my $next_tmp = '1';
                        while ($next_tmp eq '1')
                        {
                            $next_tmp = $line =~ s/N/$nuc_tmp/;
                            $chance_tmp = int(rand(4));
                            $nuc_tmp = $nucs[$chance_tmp];
                        }
                    }
                    
                    $TRA_seq_tmp .= $line;
                }
            }         
            $current_contig_tmp += length($line);
        }
    }        
    if ($TRA_seq_tmp ne "")
    {
        $TRA_sequences{$chromosome_tmp} = $TRA_seq_tmp;        
        $TRA_seq_tmp = "";  
    }
    close $FILE_REF;
    
    foreach my $TRA_hap1_tmp (keys %TRA_hap1)
    {
        if (exists($contigs{$TRA_hap1_tmp}))
        {       
        }
        elsif (exists($contigs_chr{$TRA_hap1_tmp}))
        {       
        }
        else
        {
            die "\n\nERROR: This chromosome (".$TRA_hap1_tmp.") of the translocation input is not recognizied\n";
        }
    }
    foreach my $TRA_hap2_tmp (keys %TRA_hap2)
    {
        if (exists($contigs{$TRA_hap2_tmp}))
        {       
        }
        elsif (exists($contigs_chr{$TRA_hap2_tmp}))
        {       
        }
        else
        {
            die "\n\nERROR: This chromosome (".$TRA_hap2_tmp.") of the translocation input is not recognizied\n";
        }
    }
}

#------------------------------------------------------------------------------------------------------------------------------------------------------

print "...OK\n\n\n";
print "Reference sequence\n";
print "-------------------\n";
print "Contigs          : ".$count_contigs."\n";
print "Reference Length : ".$reference_size." bp\n\n\n";
print OUTPUT_LOG "Reference sequence\n";
print OUTPUT_LOG "-------------------\n";
print OUTPUT_LOG "Contigs          : ".$count_contigs."\n";
print OUTPUT_LOG "Reference Length : ".$reference_size." bp\n\n\n";


#Read Error profiles-----------------------------------------------------------------------------------------------------------------------------
my %NP_error_profile_mismatch;
undef %NP_error_profile_mismatch;
my %NP_error_profile_mismatch2;
undef %NP_error_profile_mismatch2;
my %NP_error_profile_ins;
undef %NP_error_profile_ins;
my %NP_error_profile_ins2;
undef %NP_error_profile_ins2;
my %NP_error_profile_del;
undef %NP_error_profile_del;
my %NP_del_length;
undef %NP_del_length;
my %NP_ins_length;
undef %NP_ins_length;
my $NP_mismatch_multi2 = "";
my $NP_ins_multi = "";
my $NP_ins_multi2 = "";
my $NP_del_multi = "";
my $NP_mismatch_ratio = "";
my $NP_ins_ratio = "";
my $NP_del_ratio = "";
my $mismatch_percentage = "";
my $ins_percentage = "";
my $del_percentage = "";
    
if ($NP_coverage > 0)
{
    my $NP_ERROR;
    open($NP_ERROR, $NP_error_profile) or die "\n\nCan't open the error profile file $NP_error_profile, $!\n";
    
    my $error_type = "";
    my $lowest_del = '100';
    my $lowest_ins = '100';
    while (my $line = <$NP_ERROR>)
    {
        chomp($line);
        if ($line eq "DELETIONS")
        {
            $error_type = "del";
            next;
        }
        elsif ($line eq "INSERTIONS")
        {
            $error_type = "ins";
            next;
        }
        if ($error_type eq "del")
        {
            my @del = split /\:/, $line;
            if ($del[1] < $lowest_del && $del[1] > 0)
            {
                $lowest_del = $del[1];  
            }
        }     
        elsif ($error_type eq "ins")
        {
            my @ins = split /\:/, $line;
            if ($ins[1] < $lowest_ins && $ins[1] > 0)
            {
                $lowest_ins = $ins[1];  
            }
        }
    }
    my $ins_above_one = '1';
    my $del_above_one = '1';
    if ($lowest_ins < 1 && $lowest_ins > 0)
    {
        $ins_above_one = (1/$lowest_ins)+0.1;
    }
    if ($lowest_del < 1 && $lowest_del > 0)
    {
        $del_above_one = (1/$lowest_del)+0.1;
    }
    close $NP_ERROR;
    
    open($NP_ERROR, $NP_error_profile) or die "\n\nCan't open the error profile file $NP_error_profile, $!\n";
    
    $error_type = "";
    my $highest_mismatch2 = '0';
    my $highest_ins = '0';
    my $highest_ins2 = '0';
    my $highest_del = '0';
    my $l = '0';

    while (my $line = <$NP_ERROR>)
    {
        chomp($line);
        if ($l eq '0')
        {
            $mismatch_percentage = $line;
            $l++;
            next;
        }
        elsif ($l eq '1')
        {
            $del_percentage = $line;
            $l++;
            next;
        }
        elsif ($l eq '2')
        {
            $ins_percentage = $line;
            $l++;
            next;
        }
        elsif ($line eq "DEL_LENGTH")
        {
            $error_type = "del_length";
            next;
        }
        elsif ($line eq "INS_LENGTH")
        {
            $error_type = "ins_length";
            next;
        }
        elsif ($line eq "MISMATCHES")
        {
            $error_type = "mismatch";
            next;
        }
        elsif ($line eq "MISMATCHES2")
        {
            $error_type = "mismatch2";
            next;
        }
        elsif ($line eq "DELETIONS")
        {
            $error_type = "del";
            next;
        }
        elsif ($line eq "INSERTIONS")
        {
            $error_type = "ins";
            next;
        }
        elsif ($line eq "INSERTIONS2")
        {
            $error_type = "ins2";
            next;
        }
        if ($error_type eq "del_length")
        {
            my @del_length = split /\:/, $line;
            $NP_del_length{$del_length[0]} = $del_length[1];
        }
        elsif ($error_type eq "ins_length")
        {
            my @ins_length = split /\:/, $line;
            $NP_ins_length{$ins_length[0]} = $ins_length[1];
        }
        elsif ($error_type eq "mismatch")
        {
            my @mismatch = split /\:/, $line;
            $NP_error_profile_mismatch{$mismatch[0]} = $mismatch[1];
        }
        elsif ($error_type eq "mismatch2")
        {
            my @mismatch2 = split /\:/, $line;
            my $tmp_rate = $mismatch2[1];
            $NP_error_profile_mismatch2{$mismatch2[0]} = $tmp_rate;
            if ($tmp_rate > $highest_mismatch2)
            {
                $highest_mismatch2 = $tmp_rate;  
            }
        }      
        elsif ($error_type eq "del")
        {
            my @del = split /\:/, $line;
            my $tmpie = ($del[1]*$del_above_one)**2;
            if ($del[0] eq "TTT" || $del[0] eq "CCC" || $del[0] eq "AAA" ||$del[0] eq "GGG")
            {
                $tmpie *= 1.35
            }
            $NP_error_profile_del{$del[0]} = $tmpie;
            if ($tmpie > $highest_del)
            {
                $highest_del = $tmpie;  
            }
        }     
        elsif ($error_type eq "ins")
        {
            my @ins = split /\:/, $line;
            my $tmpie = ($ins[1]*$ins_above_one)**4;
            #print $tmpie."\n";
            $NP_error_profile_ins{$ins[0]} = $tmpie;
            if ($tmpie > $highest_ins)
            {
                $highest_ins = $tmpie;  
            }
        }
        elsif ($error_type eq "ins2")
        {
            my @ins2 = split /\:/, $line;
            $NP_error_profile_ins2{$ins2[0]} = $ins2[1];
            if ($ins2[1] > $highest_ins2)
            {
                $highest_ins2 = $ins2[1];  
            }
        }
    }
    $NP_mismatch_multi2 = 100/$highest_mismatch2;
    $NP_ins_multi = 100/$highest_ins;
    $NP_ins_multi2 = 100/$highest_ins2;
    $NP_del_multi = 100/$highest_del;
    
    $NP_mismatch_ratio = $mismatch_percentage/($mismatch_percentage+$ins_percentage+$del_percentage);
    $NP_ins_ratio = $ins_percentage/($mismatch_percentage+$ins_percentage+$del_percentage);
    $NP_del_ratio = $del_percentage/($mismatch_percentage+$ins_percentage+$del_percentage);
    
    foreach my $tmpie2 (keys %NP_error_profile_mismatch2)
    {
        my $total_chance = '0';
        foreach my $tmpie1 (keys %NP_error_profile_mismatch)
        {
            my $tmpie1_min = $tmpie1;
            substr $tmpie1_min, 2, 1, "";
            if ($tmpie1_min eq $tmpie2)
            {
                $total_chance += $NP_error_profile_mismatch{$tmpie1};
            }
        }
        my $multie_tmp = 100/$total_chance;
        foreach my $tmpie1 (keys %NP_error_profile_mismatch)
        {
            my $tmpie1_min = $tmpie1;
            substr $tmpie1_min, 2, 1, "";
            if ($tmpie1_min eq $tmpie2)
            {
                my $rate_tmp = $NP_error_profile_mismatch{$tmpie1}*$multie_tmp;
                $NP_error_profile_mismatch{$tmpie1} = $rate_tmp;
            }
        }
    }
    
    close $NP_ERROR;
}
#-----------------------------------------------------------------------------------------------------------------------------------------------------------


my $TOTAL_input = $DEL_input+$INS_input+$INV_input+$DUP_input+$TRA_input+$CSUB_input+$IDUP_input+$VCF_input_DEL+$VCF_input_INS+$VCF_input_DUP+$VCF_input_INV+$VCF_input_TRA+$VCF_input_CSUB+$VCF_input_IDUP;
my $TOTAL_interval = int($reference_size/(($TOTAL_input*1.05)+1));
my $TOTAL_interval_tmp = $TOTAL_interval;
my $random_length_interval = $TOTAL_interval_tmp;

my $DEL_interval = int($reference_size/($DEL_input+$VCF_input_DEL+1));
my $DEL_interval_tmp = $DEL_interval;
my @DEL_range = split /-/, $DEL_range;
my $DEL_range_low = $DEL_range[0];
my $DEL_range_high = $DEL_range[1];
my $DEL_input_tmp = $DEL_input;
my $DEL_range2 = $DEL_range_high-$DEL_range_low;
my $DEL_range2_short = 350-$DEL_range_low;
my $deleting = "";
my $random_length_DEL = "";

my $IDUP_interval = int($reference_size/($IDUP_input+$VCF_input_IDUP+1));
my $IDUP_interval_tmp = $IDUP_interval;
my @IDUP_range = split /-/, $IDUP_range;
my $IDUP_range_low = $IDUP_range[0];
my $IDUP_range_high = $IDUP_range[1];
my $IDUP_range2 = $IDUP_range_high-$IDUP_range_low;
my $IDUP_input_tmp = $IDUP_input;
my $random_length_IDUP = "";
my $IDUP_sequence = "";

my $INS_interval = int($reference_size/($INS_input+1+$VCF_input_INS));
my $INS_interval_tmp = $INS_interval;
my @INS_range = split /-/, $INS_range;
my $INS_range_low = $INS_range[0];
my $INS_range_high = $INS_range[1];
my $INS_input_tmp = $INS_input;
my $INS_range2 = $INS_range_high-$INS_range_low;
my $INS_range2_short = 350-$INS_range_low;
my $random_length_INS = "";

my $INV_interval = int($reference_size/($INV_input+$VCF_input_INV+1));
my $INV_interval_tmp = $INV_interval;
my @INV_range = split /-/, $INV_range;
my $INV_range_low = $INV_range[0];
my $INV_range_high = $INV_range[1];
my $INV_input_tmp = $INV_input;
my $random_length_INV = "";
my $inverting = "";
my $inversion_seq = "";
my $INV_graph_interval = int(($INV_range_high-$INV_range_low)/($INV_input+1/50));

my $CSUB_interval = int($reference_size/($CSUB_input+$VCF_input_CSUB+1));
my $CSUB_interval_tmp = $CSUB_interval;
my @CSUB_range = split /-/, $CSUB_range;
my $CSUB_range_low = $CSUB_range[0];
my $CSUB_range_high = $CSUB_range[1];
my $CSUB_input_tmp = $CSUB_input;
my $csubbing = "";

my $DUP_interval = int($reference_size/($DUP_input+$VCF_input_DUP+1));
my $DUP_interval_tmp = $DUP_interval;
my @DUP_range = split /-/, $DUP_range;
my $DUP_range_low = $DUP_range[0];
my $DUP_range_high = $DUP_range[1];
my @DUP_copies = split /-/, $DUP_copies;
my $DUP_copies_low = $DUP_copies[0];
my $DUP_copies_high = $DUP_copies[1];
my $DUP_input_tmp = $DUP_input;
my $DUP_range2 = $DUP_range_high-$DUP_range_low;
my $DUP_range2_short = 350-$DUP_range_low;
my $DUP_copies2 = $DUP_copies_high-$DUP_copies_low;
my $DUP_dis = '0';
my $random_length_DUP = "";
my $random_copies_DUP = "";
my $duplicating = "";
my $duplication_seq = "";

my $NEXT_SV = "";
my $reference_size2 = '0';
my $current_contig = "";
my $size_current_contig = '0';

#Subroutines------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------------

#Define the haplotype-------------------------------------------------------
sub hap
{
    my $c = rand(1);
    my $rand_h = rand(1);
    if ($rand_h < $heterozygosity)
    {    
        if ($c > 0.5)
        {
            $hap = "1/0";
        }
        else
        {
            $hap = "0/1";
        }
    }
    else
    {
        $hap = "1/1";
    }
}
#Resets values for next SV & next SV is determined (only for VCF input)-------------------------------------------------------------
sub VCF_input_sub
{
    my ($next_length_start_tmp) = (@_);
    $hap = $next_hap;
    $SEQ = $next_seq;
    $next_pos = "";
    $next_length = "";
    $next_type = "";
    $next_seq = "";
    $VCF_input_now = "";
    
    if ($next_length_start_tmp ne "")
    {
        $variation_haplo .= $next_length_start_tmp;
    }

SVCF:foreach my $start_pos_tmp (sort {$a <=> $b} keys %VCF_input)
    {
        my @vcf_input = split /\*/, $VCF_input{$start_pos_tmp};
        $next_pos = $start_pos_tmp;
        $next_type = $vcf_input[0];
        $next_length = $vcf_input[1];
        $next_hap = $vcf_input[2];
        $next_seq = $vcf_input[3];
        delete $VCF_input{$start_pos_tmp};  
        last SVCF;
    }
}
#------------------------------------------------------------------------------
sub insert_seq
{
    my @value = @_;
    my $random_length_tmp = $value[0];
    my $insert_tmp = $value[1];
    my $nuc_range = '4';
    
    if ($insert_tmp eq "")
    {
        while (length($insert_tmp) < $random_length_tmp)
        {
            my $random_nuc = int(rand($nuc_range));
            $insert_tmp .= $nucs[$random_nuc-1];
        }
    }

    while (length($insert_tmp) > 0)
    {
        my $insert_part = substr $insert_tmp, 0, $length_fasta_line, "";
        $variation_haplo .= $insert_part;
    }
    return $insert_tmp;
}
#------------------------------------------------------------------------------------------------------------------------------------
#Read simulation parameters---------------------------------------------------------------------------------------------

my @NP_range = split /-/, $NP_range;
my $NP_range_low = $NP_range[0];
my $NP_range_high = $NP_range[1];
my $output_NP = $output."Nanopore_".$project.".fasta";
my $NP_read_count = '0';
my %NP_seq_length;
undef %NP_seq_length;
my %NP_seq_length1;
undef %NP_seq_length1;
my %NP_seq_length2;
undef %NP_seq_length2;

if ($NP_accuracy =~ m/^(\d+)%.*$/)
{
    $NP_accuracy = $1;
}
else
{
    $NP_accuracy = 100-($mismatch_percentage+$ins_percentage+$del_percentage);
}
my $NP_error_rate = 100-$NP_accuracy;
my $NP_min_read_length = '10000000000000000';
my $NP_max_read_length = '0';
my $NP_total_length = '0';

my $start_read_lq = '15';
#if ($NP_chemistry eq "R7")
#{
    #$start_read_lq = '60';
#}

my $NP_extra_long_reads = '0';
my %NP_extra_long_reads;
undef %NP_extra_long_reads;
if ($NP_range_high > $NP_average*2+5000)
{
    my $read_count_estimation = ($reference_size*15)/$NP_average;      
    $NP_extra_long_reads = int($read_count_estimation*0.005);
    
    my $n = '0';
    while ($n < $NP_extra_long_reads)
    {
        my $pos_tmp = int(rand($read_count_estimation));
        $NP_extra_long_reads{$pos_tmp} = undef;
        $n++;
    }
}
if ($NP_coverage > 0)
{
    open(OUTPUT_NP, ">" .$output_NP) or die "Can't open Nanopore file $output_NP, $!\n";
}

select(STDERR);
$| = 1;
select(STDOUT); # default
$| = 1;
print "\nStart simulation...     ";

$chromosome = "";
my $chr_prev = "";
if ($check_zip eq "gz")
{
    open ($FILE_REF, '-|', 'gzip', '-dc', $reference) or die "Can't open file $reference, $!\n";
}
elsif ($check_zip2 eq "bz2")
{
    open ($FILE_REF, '-|', 'bzip2', '-dc', $reference) or die "Can't open file $reference, $!\n";
}
else
{
    open($FILE_REF, $reference) or die "\n\nCan't open reference file $reference, $!\n";
}

if ($VCF_input ne "")
{
    open (INPUT_VCF, $VCF_input) or die "Can't open file $VCF_input, $!\n";
}


my $progress_print_prev = '0';
my $progress_print = "";

if ($heterozygosity > 0 && $NP_coverage > 0)
{
    $NP_coverage1 = ($NP_coverage)/2;
    if ($NP_coverage1 % 2 == 0)
    {     
    }
    else
    {
        $NP_coverage1 = ($NP_coverage-1)/2;
        $NP_coverage1 += int(rand(1));
    }
    $NP_coverage2 = $NP_coverage-$NP_coverage1;
}

#----------------------------------------------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------START SIMULATION-------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------------------------------

REF:while (my $line = <$FILE_REF>)
{
    chomp $line;
#my $time_start_tmp = time;
    my $progress = ($reference_size2/$reference_size)*100;
    $progress += 1;
    $progress_print = int($progress)."%";
    if (int($progress) > 100)
    {
        $progress_print = "100%";
    }
    if ($progress_print > $progress_print_prev)
    {
        $|=1;
        print "\b" x length($progress_print);
        print ' ' x length($progress_print);
        print "\b" x length($progress_print);
        print $progress_print;
    }
    $progress_print_prev = $progress_print;

    if ($reference_size2 < 1)
    {
        $reference_size2++;
          
        if ($line =~ m/>(\d+|X|Y|MT)\s+dna:chromosome.*/)
        {
            $chromosome = $1;
        }
        elsif ($line =~ m/>chr(\d+|X|Y|MT)$/)
        {
            $chromosome = $1;
        }
        elsif ($line =~ m/>.*chromosome\s(\d+|X|Y|MT).*/)
        {
            $chromosome = $1;
        }
        else
        {
            $chromosome = substr $line, 1;
        }

        $current_contig = $line;
        print OUTPUT_REF $line."\n";
        print OUTPUT_HAP1 $line."\n";
        print OUTPUT_HAP2 $line."\n";

        if ($VCF_input ne "")
        {
            while (my $line = <INPUT_VCF>)
            {
                if ($line =~ m/(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)*\s+(\S+)*.*/)
                {
                    my $chr_tmp = $1;
                    my $start_pos_tmp = $2;
                    my $length_tmp = $3;
                    my $type_tmp = $4;
                    my $hap_tmp = $5;
                    my $SEQ_tmp = $6;
                  
                    if ($SEQ_tmp ne "" && $type_tmp ne "TRA")
                    {
                        if (exists($foreign_contigs{$SEQ_tmp}))
                        {
                            $sequences_foreign{$start_pos_tmp} = $foreign_contigs{$SEQ_tmp};
                        }
                    }
                    if ($chr_tmp eq $chromosome && $type_tmp ne "TRA")
                    {
                        $VCF_input{$start_pos_tmp} = $type_tmp."\*".$length_tmp."\*".$hap_tmp."\*".$SEQ_tmp;
                    }
                }
            }
            foreach my $start_pos_tmp (sort {$a <=> $b} keys %VCF_input)
            {
               my @vcf_input = split /\*/, $VCF_input{$start_pos_tmp};
               $next_pos = $start_pos_tmp;
               $next_type = $vcf_input[0];
               $next_length = $vcf_input[1];
               $next_hap = $vcf_input[2];
               $next_seq = $vcf_input[3];
               delete $VCF_input{$start_pos_tmp};
               last;
            }
        }
        if (exists($TRA_hap1{$chromosome}))
        {
            foreach my $line_tmp2 (sort {$a <=> $b} keys %{$TRA_hap1{$chromosome}})
            {
                $random_length_TRA = $TRA_hap1{$chromosome}{$line_tmp2};
                $TRA_contig_select = $chromosome;
                $TRA_contig_select2 = $line_tmp2;
            }
        }
        if (exists($TRA_hap2{$chromosome}))
        {
            foreach my $line_tmp2 (sort {$a <=> $b} keys %{$TRA_hap2{$chromosome}})
            {
                $random_length_TRA = $TRA_hap2{$chromosome}{$line_tmp2};
                $TRA_contig_select = $chromosome;
                $TRA_contig_select2 = $line_tmp2;
            }          
        }        
        next REF;
    }     
    #$line =~ tr/actgn/ACTGN/;
    my $first = substr $line, 0, 1;
    my $line3 = "";
NEW_CONTIG0:
    if ($first eq '>' || $first eq '@')
    {
#print chromosome to haplo files-------------------------------------------------------------------------------------------------    
        
        if ($translocating eq "yes")
        {
            $translocating = "";
            $TRA_contig_select = "";
            $TRA_contig_select2 = "";
            $random_length_TRA = "";
        }
        if ($chromosome ne "")
        {
            $chr_prev = $chromosome;
        }
        
        if ($line =~ m/>(\d+|X|Y|MT)\s+dna:chromosome.*/)
        {
            $chromosome = $1;
        }
        elsif ($line =~ m/>chr(\d+|X|Y|MT)$/)
        {
            $chromosome = $1;
        }
        elsif ($line =~ m/>.*chromosome\s(\d+|X|Y|MT).*/)
        {
            $chromosome = $1;
        }
        else
        {
            $chromosome = substr $line, 1;
        }
        
        while ($pos_hap1 < length($haplotype1_print))
        {
            my $print_line = substr $haplotype1_print, $pos_hap1, $length_fasta_line;
            print OUTPUT_HAP1 $print_line."\n";
            $pos_hap1 += length($print_line);
        }
        while ($pos_hap2 < length($haplotype2_print))
        {
            my $print_line = substr $haplotype2_print, $pos_hap2, $length_fasta_line;
            print OUTPUT_HAP2 $print_line."\n";
            $pos_hap2 += length($print_line);
        }
        while ($pos_hap_merged < length($haplo_merged_print))
        {
            my $print_line = substr $haplo_merged_print, $pos_hap_merged, $length_fasta_line;
            print OUTPUT_REF $print_line."\n";
            $pos_hap_merged += length($print_line);
        }
        $haplotype1_print = "";
        $haplotype2_print = "";
        $haplo_merged_print = "";
        $pos_hap1 = '0';
        $pos_hap2 = '0';
        $pos_hap_merged = '0';

        if ($NP_coverage > 0 && $new_contig eq "" && $translocating eq "")
        {
            $new_contig = "yes";
            $chromosome = "";
            goto NEW_CONTIG;
        }
        $new_contig = "";
#-------------------------------------------------------------------------------------------------------------------------------   
       
        if (exists($TRA_hap1{$chromosome}))
        {
            foreach my $line_tmp2 (sort {$a <=> $b} keys %{$TRA_hap1{$chromosome}})
            {
                $random_length_TRA = $TRA_hap1{$chromosome}{$line_tmp2};
                $TRA_contig_select = $chromosome;
                $TRA_contig_select2 = $line_tmp2;
            }
        }
        if (exists($TRA_hap2{$chromosome}))
        {
            foreach my $line_tmp2 (sort {$a <=> $b} keys %{$TRA_hap2{$chromosome}})
            {
                $random_length_TRA = $TRA_hap2{$chromosome}{$line_tmp2};
                $TRA_contig_select = $chromosome;
                $TRA_contig_select2 = $line_tmp2;
            }          
        }        
        
        print OUTPUT_REF $line."\n";
        print OUTPUT_HAP1 $line."\n";
        print OUTPUT_HAP2 $line."\n";
        
        $current_contig = $line;
        $size_current_contig = '0';       
        
        if ($VCF_input ne "")
        {
            #print OUTPUT_LOG $chromosome." Chr\n";

            foreach my $VCF_input (sort {$a <=> $b} keys %VCF_input2)
            {
                if (exists($VCF_output{$VCF_input}))
                {
                }
                else
                {
                    print OUTPUT_LOG $VCF_input." NOT\n";
                    print OUTPUT_LOG "Skipped SV: ".$chr_prev."\t".$VCF_input."\n";
                }
            }
            undef %VCF_input;
            undef %VCF_input2;
            undef %VCF_output;
            close INPUT_VCF;
            open (INPUT_VCF, $VCF_input) or die "Can't open file $VCF_input, $!\n";
            while (my $line = <INPUT_VCF>)
            {
                if ($line =~ m/(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)*\s+(\S+)*.*/)
                {
                    my $chr_tmp = $1;
                    my $start_pos_tmp = $2;
                    my $length_tmp = $3;
                    my $type_tmp = $4;
                    my $hap_tmp = $5;
                    my $SEQ_tmp = $6;
                    if ($chr_tmp eq $chromosome && $type_tmp ne "TRA")
                    {
                        $VCF_input{$start_pos_tmp} = $type_tmp."\*".$length_tmp."\*".$hap_tmp."\*".$SEQ_tmp;
                        $VCF_input2{$start_pos_tmp} = undef;
                    }
                }
            }
            foreach my $start_pos_tmp (sort {$a <=> $b} keys %VCF_input)
            {
               my @vcf_input = split /\*/, $VCF_input{$start_pos_tmp};
               $next_pos = $start_pos_tmp;
               $next_type = $vcf_input[0];
               $next_length = $vcf_input[1];
               $next_hap = $vcf_input[2];
               $next_seq = $vcf_input[3];
               delete $VCF_input{$start_pos_tmp};
               last;
            }
        }
        next REF;
    }
#--------   
    if ($translocating eq "yes")
    {
        if (eof)
        {
            $new_contig = "yes";
            goto NEW_CONTIG;
        }
        next REF;
    }

    my $last_nuc = "";
    if ($ambigious_replace eq "yes")
    {
        $last_nuc = substr $line, -1, 1;
    }
    if (($first eq "N" && $last_nuc eq "N") && $ambigious_replace eq "yes")
    {
        my $line_tmp = "";
        while (length($line_tmp) < length($line))
        {
            my $chance_tmp = int(rand(4));
            my $nuc_tmp = $nucs[$chance_tmp];
            $line_tmp .= $nuc_tmp; 
        }
        $line = $line_tmp;
    }
    elsif (($first eq "N" || $last_nuc eq "N") && $ambigious_replace eq "yes")
    {
        my $chance_tmp = int(rand(4));
        my $nuc_tmp = $nucs[$chance_tmp];
        my $next_tmp = '1';
        while ($next_tmp eq '1')
        {
            $next_tmp = $line =~ s/N/$nuc_tmp/;
            $chance_tmp = int(rand(4));
            $nuc_tmp = $nucs[$chance_tmp];
        }
    }
#Translocating shit----------------------------------------------------------------------------------------------     

    if ($first eq "N" && $next_pos >= $size_current_contig+length($line) && $deleting eq "" && $inverting eq "" && $duplicating eq "")
    {

        $reference_size2 += length($line);
        $size_current_contig += length($line);
        $haplotype1 .= $line;
        $haplotype2 .= $line;
        $haplo_merged .= $line;
        $haplotype1_print .= $line;
        $haplotype2_print .= $line;
        $haplo_merged_print .= $line;
        next REF;
    }
   
    if ($TRA_contig_select eq $chromosome && $contigs{$current_contig}-$size_current_contig < $random_length_TRA && $translocating eq "")
    {
        my $tra_seq = $TRA_sequences{$TRA_contig_select2};
        my $check1 = "";
        my $check2 = "";
        my $hap_tmp = "1/1";
        if (exists($TRA_hap1{$chromosome}))
        {
            $check1 = "yes";
        }
        if (exists($TRA_hap2{$chromosome}))
        {
            $check2 = "yes";
        }
        
        while (length($tra_seq) > 0)
        {
            my $print_line = substr $tra_seq, 0, $length_fasta_line, "";
            $haplo_merged_print .= $print_line;
            $reference_size2 += length($print_line);
            if ($check1 eq "yes")
            {
                if ($NP_coverage > 0)
                {
                    $haplotype1 .= $print_line;
                }
                $haplotype1_print .= $print_line;
            }
            if ($check2 eq "yes")
            {
                if ($NP_coverage > 0)
                {
                    $haplotype2 .= $print_line;
                }
                $haplotype2_print .= $print_line;
            }
            if ($NP_coverage > 0)
            {
                $haplo_merged .= $print_line;
            }           
        }
        if ($check1 ne "yes" || $check2 ne "yes")
        {
            my $tra_seq2 = $TRA_sequences{$TRA_contig_select};
            while (length($tra_seq2) > 0)
            {
                my $print_line = substr $tra_seq2, 0, $length_fasta_line, "";
                if ($check1 ne "yes")
                {
                    if ($NP_coverage > 0)
                    {
                        $haplotype1 .= $print_line;
                    }
                    $haplotype1_print .= $print_line;
                    $hap_tmp = "1";
                }
                if ($check2 ne "yes")
                {
                    if ($NP_coverage > 0)
                    {
                        $haplotype2 .= $print_line;
                    }
                    $haplotype2_print .= $print_line;
                    $hap_tmp = "0";
                }
            }
        }

        $translocating = "yes";        
        my $chromosome2 = $TRA_contig_select2;
        
        if (exists($TRA_hap_second{$chromosome}))
        {}
        else
        {
            print OUTPUT_VCF $chromosome.":".$chromosome2."\t".$size_current_contig."\t".$random_length_TRA."\tTRA\t".$hap_tmp."\n";
        }
        next REF;
    }
SELECT_SV:

#Select type SV----------------------------------------------------------------------------------------------------------------        
    if ($SV_input ne "")
    {
        if ($inverting eq "")
        {
            $inversion_seq = "";
        }

        if ($DEL_input_tmp eq '0')
        {
            $DEL_interval_tmp = $reference_size;
        }
        if ($INS_input_tmp eq '0')
        {
            $INS_interval_tmp = $reference_size;
        }
        if ($DUP_input_tmp eq '0')
        {
            $DUP_interval_tmp = $reference_size;
        }
        if ($INV_input_tmp eq '0')
        {
            $INV_interval_tmp = $reference_size;
        }
        if ($IDUP_input_tmp eq '0')
        {
            $IDUP_interval_tmp = $reference_size;
        }
        if ($CSUB_input_tmp eq '0')
        {
            $CSUB_interval_tmp = $reference_size;
        }
        if ($next_pos >= $size_current_contig && $next_pos < $size_current_contig+length($line) && $VCF_input ne "" && $next_type ne "")
        {
            $NEXT_SV = $next_type;

            $VCF_input_now = "yes";
            $next_length_minus = $next_pos-$size_current_contig;
            $next_length_minus_save = $next_length_minus;
            
            if ($inverting ne "" || $duplicating ne "" || $deleting ne "")
            {
                #print $NEXT_SV." ALARM\n";
                #print $inverting." INV\n";
                #print $duplicating." DUP\n";
                #print $deleting." DEL\n";           
                $SVs_whitin_1_line = "yes";
            }
            if ($size_current_contig eq "160654020")
                {
                    #print OUTPUT_LOG $next_length_minus." NLM\n";
                }
        }
        elsif ($next_pos < $size_current_contig && $VCF_input ne "" && $inverting eq "" && $duplicating eq "" && $deleting eq "" && $next_pos ne "")
        {
            print OUTPUT_LOG $NEXT_SV." ALARM\n";
            print OUTPUT_LOG $next_type." next_type\n";
            print OUTPUT_LOG $next_pos." next_pos\n";
            print OUTPUT_LOG $size_current_contig." CURRENT_C\n";
            VCF_input_sub;
        }
        elsif ($NEXT_SV eq "" && $reference_size2 > $TOTAL_interval_tmp && $deleting eq "" && $inverting eq "" && $duplicating eq "")
        {
            if ($DEL_input_tmp > 0 && $DEL_interval_tmp <= $INS_interval_tmp && $DEL_interval_tmp <= $INV_interval_tmp && $DEL_interval_tmp <= $DUP_interval_tmp && $DEL_interval_tmp <= $CSUB_interval_tmp && $DEL_interval_tmp <= $IDUP_interval_tmp)
            {
                $NEXT_SV = "DEL";
            }
            elsif ($INS_input_tmp > 0 && $INS_interval_tmp <= $DEL_interval_tmp && $INS_interval_tmp <= $INV_interval_tmp && $INS_interval_tmp <= $DUP_interval_tmp && $INS_interval_tmp <= $CSUB_interval_tmp && $INS_interval_tmp <= $IDUP_interval_tmp)
            {
                $NEXT_SV = "INS";
            }
            elsif ($DUP_input_tmp > 0 && $DUP_interval_tmp <= $DEL_interval_tmp && $DUP_interval_tmp <= $INV_interval_tmp && $DUP_interval_tmp <= $INS_interval_tmp && $DUP_interval_tmp <= $CSUB_interval_tmp && $DUP_interval_tmp <= $IDUP_interval_tmp)
            {
                $NEXT_SV = "DUP";
            }
            elsif ($INV_input_tmp > 0 && $INV_interval_tmp <= $DEL_interval_tmp && $INV_interval_tmp <= $INS_interval_tmp && $INV_interval_tmp <= $DUP_interval_tmp && $INV_interval_tmp <= $CSUB_interval_tmp && $INV_interval_tmp <= $IDUP_interval_tmp)
            {
                $NEXT_SV = "INV";
            }
            elsif ($IDUP_input_tmp > 0 && $IDUP_interval_tmp <= $DEL_interval_tmp && $IDUP_interval_tmp <= $INS_interval_tmp && $IDUP_interval_tmp <= $DUP_interval_tmp && $IDUP_interval_tmp <= $INV_interval_tmp && $IDUP_interval_tmp <= $CSUB_interval_tmp)
            {
                $NEXT_SV = "IDUP";
            }
            elsif ($CSUB_input_tmp > 0 && $CSUB_interval_tmp <= $DEL_interval_tmp && $CSUB_interval_tmp <= $INS_interval_tmp && $CSUB_interval_tmp <= $DUP_interval_tmp && $CSUB_interval_tmp <= $INV_interval_tmp && $CSUB_interval_tmp <= $IDUP_interval_tmp)
            {
                $NEXT_SV = "CSUB";
            }
        }
#Input SV------------------------------------------------------------------------------------------------------------------------                
        if ($NEXT_SV ne "DEL" && $NEXT_SV ne "INV" && $NEXT_SV ne "INS" && $NEXT_SV ne "DUP" && $NEXT_SV ne "TRA" && $NEXT_SV ne "CSUB" && $NEXT_SV ne "IDUP" && $NEXT_SV ne "")
        {
            die "\n".$NEXT_SV.": This type of SV is not supported\n";
        }
#Deletions & Complex substitutions--------------------------------------------------------------------------------------------------------------------------------
        if ((($NEXT_SV eq "DEL" && (($reference_size2 > $random_length_interval && $DEL_input_tmp > 0 && $deleting eq "") || $VCF_input_now eq "yes")) ||
            ($NEXT_SV eq "CSUB" && (($reference_size2 > $random_length_interval && $CSUB_input_tmp > 0 && $deleting eq "") || $VCF_input_now eq "yes"))) && $SVs_whitin_1_line ne "yes")
        {
            if ($VCF_input_now eq "yes")
            {
                goto VCF_INPUT_DEL;
            }
            my $w=0;
            $w =rand(15);
            my $distribution = int($w); 
            my $random_length = "50";
            
            if ($distribution > -1 && $distribution < 9)
            {
DEL_RANGE2:
                my $t=0;
                $t+=rand() for(2..11);
                #my $random_length = int($t*60+1);
                $random_length = int((2**($t*1.3))-10);
                if ($random_length > 1000 || $random_length < 30)
                {
                    goto DEL_RANGE2;
                }
                elsif ($random_length < $DEL_range_low)
                {
                    goto DEL_RANGE2;
                }
            }
            elsif ($distribution > 8 && $distribution < 10)
            {
DEL_RANGE2b:
                my $s=0;
                $s+=rand() for(-1.5..1.5);
                $random_length = int($s*30+300);
                if ($random_length < $DEL_range_low)
                {
                    goto DEL_RANGE2b;
                }
            }
            else
            {
DEL_RANGE3:                
                my $s=0;
                $s+=rand() for(0..4);
                $random_length = int($s*($DEL_range_high/2.5))-$DEL_range_high+15;
                if ($random_length < $DEL_range_low)
                {
                    goto DEL_RANGE3;
                }           
            }
            $random_length_DEL = $random_length;
VCF_INPUT_DEL:

            my $pos_tmp = $size_current_contig;
            
            my $line_tmp = $line;
            if ($VCF_input_now eq "yes")
            {
                $random_length_DEL = $next_length;
                $pos_tmp = $next_pos;
                
                my $next_length_start = "";
                
                if ($next_length_minus > 0)
                {
                    $next_length_start = substr $line_tmp, 0, $next_length_minus, "";
                }
                VCF_input_sub $next_length_start; 
            }
            if ($hap eq "" && $heterozygosity ne "no")
            {
                hap
            }
     
            if ($random_length_DEL <= length($line)-$next_length_minus)
            {
                substr $line_tmp, 0, $random_length_DEL, "";

                if ($NEXT_SV eq "CSUB")
                {
                    if (exists($sequences_foreign{$SEQ}))
                    {
                        insert_seq ($random_length_DEL,$sequences_foreign{$SEQ});
                    }
                    elsif (exists($sequences_final{$SEQ}))
                    {
                        insert_seq ($random_length_DEL,$sequences_final{$SEQ});
                    }
                    else
                    {
                        insert_seq ($random_length_DEL);
                    }                  
                }

                $variation_haplo .= $line_tmp;
                $finish_var = "yes";
            }
            else
            {
                $deleting = length($line)-$next_length_minus;
                if ($NEXT_SV eq "CSUB")
                {
                    $csubbing = "yes";
                }
            }
            $next_length_minus = '0';
            
            if ($NEXT_SV eq "DEL")
            {
                $DEL_input_tmp--;
                $DEL_interval_tmp = $DEL_interval+$reference_size2;
                $random_length_interval = int(rand($TOTAL_interval-2000-$random_length_DEL)) + 2000 + $random_length_DEL + $reference_size2;
                print OUTPUT_VCF $chromosome."\t".$pos_tmp."\t".$random_length_DEL."\tDEL\t".$hap."\t".$SEQ."\n";
                $VCF_output{$pos_tmp} = undef;
                
                my $range_graph = int(int($random_length_DEL)/15);
                if (exists($graph_DEL{$range_graph}))
                { 
                    my $f = $graph_DEL{$range_graph}+1;
                    $graph_DEL{$range_graph} = $f;
                }
                else
                {
                    $graph_DEL{$range_graph} = 1;
                }
            }
            else
            {
                $CSUB_input_tmp--;
                $CSUB_interval_tmp = $CSUB_interval+$reference_size2;
                $random_length_interval = int(rand($TOTAL_interval-2000)) + 2000 + $reference_size2;
                print OUTPUT_VCF $chromosome."\t".$pos_tmp."\t".$random_length_DEL."\tCSUB\t".$hap."\t".$SEQ."\n";
                $VCF_output{$pos_tmp} = undef;
            }
            $TOTAL_interval_tmp = $TOTAL_interval+$reference_size2;
            $NEXT_SV = "";
            $ref_haplo .= $line;
        }
        elsif ($deleting < $random_length_DEL && $deleting ne "")
        { 
            $ref_haplo .= $line;
            if ($random_length_DEL-$deleting <= length($line))
            {
                my $line_tmp = $line;
                substr $line_tmp, 0, $random_length_DEL-$deleting, "";
                
                if ($csubbing eq "yes")
                {
                    $csubbing = "";
                    if (exists($sequences_foreign{$SEQ}))
                    {
                        insert_seq ($random_length_DEL,$sequences_foreign{$SEQ});
                    }
                    elsif (exists($sequences_final{$SEQ}))
                    {
                        insert_seq ($random_length_DEL,$sequences_final{$SEQ});
                    }
                    else
                    {
                        insert_seq ($random_length_DEL);
                    }
                }
                $deleting = "";
                $finish_var = "yes";
                if ($SVs_whitin_1_line eq "yes")
                {
                    $SVs_whitin_1_line_length = length($line);
                    $line = $line_tmp;
                    goto TRA_SKIP;
                }

                $variation_haplo .= $line_tmp;  
            }
            else
            {
                $deleting += length($line);
            } 
        }
#Insertions and inverted duplications----------------------------------------------------------------------------------------------
        elsif ((($NEXT_SV eq "INS" && (($reference_size2 > $random_length_interval && $INS_input_tmp > 0) || $VCF_input_now eq "yes")) ||
               ($NEXT_SV eq "IDUP" && $VCF_input_now eq "yes")) && $SVs_whitin_1_line ne "yes")
        { 
            my $insert = "";          
            $random_length_INS = "50";

            my $next_length_end = "";
            if ($VCF_input_now eq "yes")
            {
                $random_length_INS = $next_length;
    
                $next_length_end = $line;

                my $next_length_start = "";
                if ($next_length_minus > 0)
                {
                    $next_length_start = substr $next_length_end, 0, $next_length_minus, "";
                }

                $variation_haplo .= $next_length_start;

                print OUTPUT_VCF $chromosome."\t".$next_pos."\t".$random_length_INS."\t".$next_type."\t".$next_hap."\t".$next_seq."\n";
                $VCF_output{$next_pos} = undef;
                
                if (exists($sequences_foreign{$next_seq}))
                {
                    $insert = $sequences_foreign{$next_seq};
                }
                elsif (exists($sequences_final{$next_seq}))
                {
                   $insert = $sequences_final{$next_seq};
                }
                my $insert_tmp = $insert;
                if ($NEXT_SV eq "IDUP")
                {
                    $insert = reverse($insert_tmp);
                }
            }
            if ($VCF_input_now eq "")
            {
                $INS_input_tmp--;
                my $w=0;
                $w =rand(15);
                my $distribution = int($w);
                            
                if ($distribution > -1 && $distribution < 9)
                {
INS_RANGE2:
                    my $t=0;
                    $t+=rand() for(2..11);
                    #my $random_length = int($t*60+1);
                    $random_length_INS = int((2**($t*1.3))-10);
                    if ($random_length_INS > 1500 || $random_length_INS < 30)
                    {
                        goto INS_RANGE2;
                    }
                    elsif ($random_length_INS < $INS_range_low)
                    {
                        goto INS_RANGE2;
                    }
                }
                elsif ($distribution > 8 && $distribution < 10)
                {
INS_RANGE2b:
                    my $s=0;
                    $s+=rand() for(-1.5..1.5);
                    $random_length_INS = int($s*30+300);
                    if ($random_length_INS < $INS_range_low)
                    {
                        goto INS_RANGE2b;
                    }
                }
                else
                {
INS_RANGE3:                
                    my $s=0;
                    $s+=rand() for(0..4);
                    $random_length_INS = int($s*($INS_range_high/2.5))-$INS_range_high+15;
                    if ($random_length_INS < $INS_range_low)
                    {
                        goto INS_RANGE3;
                    }           
                }
                if ($hap eq "" && $heterozygosity ne "no")
                {
                    hap
                }
                print OUTPUT_VCF $chromosome."\t".$size_current_contig."\t".$random_length_INS."\tINS\t".$hap."\t".$next_seq."\n";
                $VCF_output{$size_current_contig} = undef;
            }
            
            if ($NEXT_SV eq "INS")
            {
                my $range_graph = int($random_length_INS/15);
                if (exists($graph_INS{$range_graph}))
                { 
                    my $f = $graph_INS{$range_graph}+1;
                    $graph_INS{$range_graph} = $f;
                }
                else
                {
                    $graph_INS{$range_graph} = 1;
                }
            }            

            insert_seq ($random_length_INS, $insert);
    
            my $line_tmp = $line;
            if ($VCF_input_now eq "yes")
            {
                $line_tmp = $next_length_end;
                VCF_input_sub;
            }

            $variation_haplo .= $line_tmp;
            $ref_haplo .= $line;
            $finish_var = "yes";

            $next_length_minus = '0';
            $TOTAL_interval_tmp = $TOTAL_interval+$reference_size2;
            if ($NEXT_SV eq "INS")
            {
                $INS_interval_tmp = $INS_interval+$reference_size2;
            }
            else
            {
                $IDUP_interval_tmp = $IDUP_interval+$reference_size2;
            }
            $random_length_interval = int(rand($TOTAL_interval-2000)) + 2000 + $reference_size2;
            $NEXT_SV = "";
        }        
#Duplications-and IDUP from random input---------------------------------------------------------------------------------------------
        elsif (($NEXT_SV eq "DUP" && (($reference_size2 > $random_length_interval && $DUP_input_tmp > 0 && $duplicating eq "") || $VCF_input_now eq "yes")) ||
               ($NEXT_SV eq "IDUP" && $reference_size2 > $random_length_interval && $IDUP_input_tmp > 0) && $SVs_whitin_1_line ne "yes")
        {
            if ($VCF_input_now eq "yes")
            {
                goto VCF_INPUT_DUP;
            }
            my $random_length = int(rand($DUP_range2)) + $DUP_range_low;
            my $random_length_short = int(rand($DUP_range2_short)) + $DUP_range_low;
            $random_length_DUP = $random_length_short;
            $random_copies_DUP = int(rand($DUP_copies2)) + $DUP_copies_low;
            if ($NEXT_SV eq "IDUP")
            {
                $random_copies_DUP = '0';
                $random_length_DUP = int(rand($IDUP_range2)) + $IDUP_range_low;
            }
            else
            {
                $DUP_dis++;
            }

            if ($DUP_dis > 3 && $NEXT_SV eq "DUP")
            {
                $DUP_dis = '0';
                $random_length_DUP = $random_length;
            }
VCF_INPUT_DUP:

            my $pos_tmp = $size_current_contig;
            my $line_tmp = $line;

            if ($VCF_input_now eq "yes")
            {
                my @next_length = split /x/, $next_length;
                $random_length_DUP = $next_length[0];
                $random_copies_DUP = $next_length[1];
                $pos_tmp = $next_pos;
         
                my $next_length_start = "";
                if ($next_length_minus > '0')
                {
                    $next_length_start = substr $line_tmp, 0, $next_length_minus, "";
                }
                VCF_input_sub $next_length_start;        
            }
            if ($hap eq "" && $heterozygosity ne "no")
            {
                hap
            }
            if ($NEXT_SV eq "DUP")
            {
                print OUTPUT_VCF $chromosome."\t".$pos_tmp."\t".$random_length_DUP."x".$random_copies_DUP."\t".$NEXT_SV."\t".$hap."\t".$SEQ."\n";
            }
            else
            {
                print OUTPUT_VCF $chromosome."\t".$pos_tmp."\t".$random_length_DUP."\t".$NEXT_SV."\t".$hap."\t".$SEQ."\n";
            }
            $VCF_output{$pos_tmp} = undef;
            
            if (exists($sequences_final{$SEQ}) || exists($sequences_foreign{$SEQ}))
            {
                my $duplication_tmp = $sequences_foreign{$SEQ};
                if (exists($sequences_final{$SEQ}))
                {
                    $duplication_tmp = $sequences_final{$SEQ};
                }
                my $b = '1';
                my $duplication = $duplication_tmp;
                while ($b < $random_copies_DUP)
                {
                    $duplication .= $duplication_tmp;
                    $b++;
                }
                my $duplication_tmp2 = $duplication;
                if ($NEXT_SV eq "IDUP")
                {
                    $duplication = reverse($duplication_tmp2);
                }

                $variation_haplo .= $duplication;
                $variation_haplo .= $line_tmp;
                $finish_var = "yes";
            }
            elsif ($random_length_DUP <= length($line)-$next_length_minus)
            {
                my $duplication_tmp = substr $line_tmp, 0, $random_length_DUP, "";
                
                my $b = '0';
                my $duplication = $duplication_tmp;
                while ($b < $random_copies_DUP)
                {
                    $duplication .= $duplication_tmp;
                    $b++;
                }
                my $duplication_tmp2 = $duplication;
                if ($NEXT_SV eq "IDUP")
                {
                    $duplication = reverse($duplication_tmp2);
                }
                
                $variation_haplo .= $duplication;
                $variation_haplo .= $line_tmp;
                $finish_var = "yes";
            }
            else
            {
                $duplicating = length($line_tmp);
                $duplication_seq .= $line_tmp;
            }
            $ref_haplo .= $line;
                 
            $next_length_minus = '0';
            if ($NEXT_SV eq "IDUP")
            {
                $IDUP_input_tmp--;
                $IDUP_interval_tmp = $IDUP_interval+$reference_size2;
                $random_length_interval = int(rand($TOTAL_interval-2000-$random_length_DUP)) + $random_length_DUP + 2000 + $reference_size2;
            }
            else
            {
                $DUP_input_tmp--;
                $DUP_interval_tmp = $DUP_interval+$reference_size2;
                $random_length_interval = int(rand($TOTAL_interval-2000-($random_length_DUP*$random_copies_DUP))) + ($random_length_DUP*$random_copies_DUP) + 2000 + $reference_size2;
            }             
            
            $TOTAL_interval_tmp = $TOTAL_interval+$reference_size2;               
            $NEXT_SV = "";
        }
        elsif ($duplicating < $random_length_DUP && $duplicating ne "")
        {             
            $ref_haplo .= $line;
            if ($random_length_DUP-$duplicating <= length($line))
            {
                my $line_tmp = $line;
                my $duplication_tmp = substr $line_tmp, 0, $random_length_DUP-$duplicating, "";
                $duplication_seq .= $duplication_tmp;

                my $b = '1';
                my $duplication = $duplication_seq;
                while ($b < $random_copies_DUP)
                {
                    $duplication .= $duplication_seq;
                    $b++;
                }
                
                $duplicating = "";
                $duplication_seq = "";
                $finish_var = "yes";
                
                my $duplication_tmp2 = $duplication;
                if ($NEXT_SV eq "IDUP")
                {
                    $duplication = reverse($duplication_tmp2);
                }
                
                if ($SVs_whitin_1_line eq "yes")
                {
                    $SVs_whitin_1_line_length = length($line);
                    $line = $line_tmp;
      
                    $variation_haplo .= $duplication;
                    goto TRA_SKIP;
                }
                
                $variation_haplo .= $duplication;
                $variation_haplo .= $line_tmp;
            }
            else
            {
               $duplicating += length($line);
               $duplication_seq .= $line;
            }
        }
#Inversions-------------------------------------------------------------------------------------------------------------------
        elsif ($NEXT_SV eq "INV" && (($reference_size2 > $random_length_interval && $INV_input_tmp > 0 && $inverting eq "") || $VCF_input_now eq "yes") && $SVs_whitin_1_line ne "yes")
        {
            if ($VCF_input_now eq "yes")
            {
                goto VCF_INPUT_INV;
            }
INV_RANGE:
            my $t=0;
            $t+=rand() for(-1..1);
            #my $random_length = int($t*60+1);
            $random_length_INV = 23*int(2**($t*6));
            #print $contigs{$current_contig}." CC\n";
            #print $size_current_contig." CN\n";
 
            if ($random_length_INV > $INV_range_high || $random_length_INV < $INV_range_low || $random_length_INV > $contigs{$current_contig}-$size_current_contig)
            {
                goto INV_RANGE;
            }

           # my $s=0;
            #$s+=rand() for(0..9);
           # $random_length_INV = int($s*(($INV_range_high*1.3)/5))-($INV_range_high*1.3);
           # if ($random_length_INV < $INV_range_low)
           # {
           #     goto INV_RANGE;
           # } 
            
            if (exists($contigs{$current_contig}) && $contigs{$current_contig}-$size_current_contig < $random_length_INV)
            {
                if ($contigs{$current_contig}-$size_current_contig < 10000)
                {
                    $random_length_INV = $contigs{$current_contig}-$size_current_contig-500 
                }
                else
                {
                    goto INV_RANGE;
                }
            }
VCF_INPUT_INV:
            my $pos_tmp = $size_current_contig;
            my $line_tmp = $line;
            if ($VCF_input_now eq "yes")
            {
                $random_length_INV = $next_length;
                $pos_tmp = $next_pos;
                
                my $next_length_start = "";
                if ($next_length_minus > '0')
                {
                    $next_length_start = substr $line_tmp, 0, $next_length_minus, "";
                }
                VCF_input_sub $next_length_start;
            }
            if ($hap eq "" && $heterozygosity ne "no")
            {
                hap
            }
            print OUTPUT_VCF $chromosome."\t".$pos_tmp."\t".$random_length_INV."\tINV\t".$hap."\t".$SEQ."\n";
            $VCF_output{$pos_tmp} = undef;
                      
            if (exists($graph_INV2{$random_length_INV}))
            { 
                my $f = $graph_INV2{$random_length_INV}+1;
                $graph_INV2{$random_length_INV} = $f;
            }
            else
            {
                $graph_INV2{$random_length_INV} = 1;
            }
            
            if ($random_length_INV <= length($line))
            {
                my $inversion_tmp = substr $line_tmp, 0, $random_length_INV, "";
                my $inversion = reverse($inversion_tmp);

                $variation_haplo .= $inversion;
                $variation_haplo .= $line_tmp;
                $finish_var = "yes";
            }
            else
            {
                $inverting = length($line_tmp);
                $inversion_seq .= $line_tmp;
            }

            $next_length_minus = '0';
            $INV_input_tmp--;
            $TOTAL_interval_tmp = $TOTAL_interval+$reference_size2;
            $INV_interval_tmp = $INV_interval+$reference_size2;
            $random_length_interval = int(rand($TOTAL_interval-2000-$random_length_INV)) + $random_length_INV + 2000 + $reference_size2;
            $NEXT_SV = "";
            $ref_haplo .= $line;
        }
        elsif ($inverting < $random_length_INV && $inverting ne "")
        {             
            $ref_haplo .= $line;
            if ($random_length_INV-$inverting <= length($line))
            {
                my $line_tmp = $line;
                my $inversion_tmp = substr $line_tmp, 0, $random_length_INV-$inverting, "";
                $inversion_seq .= $inversion_tmp;
                my $inversion = reverse($inversion_seq);
                $inverting = "";
                $finish_var = "yes";

                if ($SVs_whitin_1_line eq "yes")
                {
                    $SVs_whitin_1_line_length = length($line);
                    $line = $line_tmp;
                    $variation_haplo .= $inversion;
                    goto TRA_SKIP;
                }
                $variation_haplo .= $inversion;
                $variation_haplo .= $line_tmp;      
            }
            else
            {
               $inverting += length($line);
               $inversion_seq .= $line;
            }
        }
        else
        {        
            $haplotype1 .= $line;
            $haplotype2 .= $line;
            $haplo_merged .= $line;
            $haplotype1_print .= $line;
            $haplotype2_print .= $line;
            $haplo_merged_print .= $line;
        }
    }
    else
    {
        $haplotype1 .= $line;
        $haplotype2 .= $line;
        $haplo_merged .= $line;
        $haplotype1_print .= $line;
        $haplotype2_print .= $line;
        $haplo_merged_print .= $line;
    }
    if ($next_pos < $size_current_contig+length($line) && $next_pos ne "")
    {
        if ($next_length_minus_save > 0)
        {
            substr $line, 0, $next_length_minus_save, "";
        }
        $reference_size2 += $next_length_minus_save;
        $size_current_contig += $next_length_minus_save;
        goto SELECT_SV;
    }
TRA_SKIP:
    if ($SVs_whitin_1_line eq "yes")
    {
        $reference_size2 += $SVs_whitin_1_line_length-length($line);
        $size_current_contig += $SVs_whitin_1_line_length-length($line);
    }
    else
    {
        $SVs_whitin_1_line = "";
        $reference_size2 += length($line);
        $size_current_contig += length($line);
    }
#Select haplotype----------------------------------------------------------------------------------------------------------------        
    if ($finish_var eq "yes")
    {       
        $haplo_merged .= $variation_haplo;
        $haplo_merged_print .= $variation_haplo;
        if ($heterozygosity eq "no")
        {
            $hap = "1/1";
        }
        if ($hap eq "1/1")
        {
            $haplotype1 .= $variation_haplo;
            $haplotype2 .= $variation_haplo;
            $haplotype1_print .= $variation_haplo;
            $haplotype2_print .= $variation_haplo;
        }
        elsif ($hap eq "1/0")
        {
            $haplotype1 .= $variation_haplo;
            $haplotype2 .= $ref_haplo;
            $haplotype1_print .= $variation_haplo;
            $haplotype2_print .= $ref_haplo;
        }
        elsif ($hap eq "0/1")
        {
            $haplotype2 .= $variation_haplo;
            $haplotype1 .= $ref_haplo;
            $haplotype2_print .= $variation_haplo;
            $haplotype1_print .= $ref_haplo;
        }
        else
        {
            print $hap." HAP\n";
            print  $size_current_contig." CURRENT_CONTIG\n";
            print  $next_pos." NEXT_POS\n";
            print  length($line)." LL\n";
            die "\nHeterozygosity parameter is not filled in correctly: give a percentage (for example: 60%). If no heterozygosity is desired, give 0%\n";
        }

        $hap = "";
        $variation_haplo = "";
        $ref_haplo = "";
        $finish_var = "";
        $finish_var_np = "yes";
        
        if ($SVs_whitin_1_line eq "yes")
        {
            $SVs_whitin_1_line = "yes2";
            #print  $size_current_contig." CURRENT_CONTIG\n";
           # print  $next_pos." NEXT_POS\n";
            #print  length($line)." LL\n";
            goto SELECT_SV;
        }
    }
NEW_CONTIG:

#print chromosome to haplo files2------------------------------------------------------------------------------------------------- 
    if (eof || $reference_size2 eq $reference_size)
    {    
        while ($pos_hap1 < length($haplotype1_print))
        {
            my $print_line = substr $haplotype1_print, $pos_hap1, $length_fasta_line;
            print OUTPUT_HAP1 $print_line."\n";
            $pos_hap1 += length($print_line);
        }
        while ($pos_hap2 < length($haplotype2_print))
        {
            my $print_line = substr $haplotype2_print, $pos_hap2, $length_fasta_line;
            print OUTPUT_HAP2 $print_line."\n";
            $pos_hap2 += length($print_line);
        }      
        while ($pos_hap_merged < length($haplo_merged_print))
        {
            my $print_line = substr $haplo_merged_print, $pos_hap_merged, $length_fasta_line;
            print OUTPUT_REF $print_line."\n";
            $pos_hap_merged += length($print_line);
        }
        if (eof && $NP_coverage > 0 && $new_contig eq "")
        {
            $new_contig = "yes";
        }
    }
    
#----------------------------------------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------Simulate Nanopore reads----------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------------------------------

    #my $time_np_start = time;
    #$time_sv += $time_np_start-$time_start_tmp;   
    

    if ($NP_coverage > 0)
    {
        my $check_haplo1_start = "";       

        if ($heterozygosity > 0 && $SV_input ne "" && $NP_coverage1 > 0 && length($haplotype1) > 0)
        {                     
            while (((length($haplotype1) >= $NP_range_high && ($finish_var_np eq "yes" || $reference_size2 < $random_length_interval || $SV_input eq ""))
                    || ($new_contig eq "yes" && length($haplotype1) > 0)))
            {
                my $o = '0';
                my $d = '0';    
                my $count_seqs = keys %NP_seq_length1;
                $check_haplo1_start = "yes";

                while ($o <= $NP_coverage1-$count_seqs)
                {
                    $NP_read_count++;
NP_LENGTH1:      
                    my $random_NP_length = '0';
                    if (exists($NP_extra_long_reads{$NP_read_count}))
                    {
                         my $range_tmp = $NP_range_high-($NP_average*2);
                         $random_NP_length = int(rand($range_tmp)) + ($NP_average*2);
                    }
                    else
                    {
                        my $j = (($NP_average*2)+400)/3;
                        my $t=0;
                        $t+=rand() for(1..3);
                        $random_NP_length = int($t*$j+1);
                    }
                    if ($random_NP_length < $NP_min_read_length)
                    {
                       $NP_min_read_length = $random_NP_length; 
                    }
                    if ($random_NP_length > $NP_max_read_length)
                    {
                        $NP_max_read_length = $random_NP_length;
                    }
                    $NP_total_length += $random_NP_length;
        
                    if (exists($NP_seq_length1{$random_NP_length+$shortest_seq1}))
                    {
                        goto NP_LENGTH1;
                    }
                    my $chance_sense = int(rand(2));
                    my $seq_tmp = "";
                    if ($chance_sense eq '0')
                    {
                        $seq_tmp = reverse(substr $haplotype1, 0, $random_NP_length-$start_read_lq);
                        $seq_tmp =~ tr/ACTG/TGAC/;
                    }
                    else
                    {
                        $seq_tmp = substr $haplotype1, 0, $random_NP_length-$start_read_lq;
                    }
                    my $seq_tmp2 = $seq_tmp;
                    my $length_read = length($seq_tmp2);

#input sequencing errors---------------------------------------------------------------------------------------------------------------------              

                    my $length_tmp = length($seq_tmp);                      
                    my $j = 1.2;
                    my $t= rand(1.5);
                    $t+= 1;
                    my $random_error_rate = int(($NP_error_rate/($t*$j))+$NP_error_rate-($NP_error_rate/2.2));
                    my $extra_error = rand(($NP_error_rate/12)/100);

                    my $NP_deletion_rate = $extra_error+($random_error_rate*($del_percentage/($del_percentage+$ins_percentage+$mismatch_percentage)))/100;
                    my $NP_insertion_rate = $extra_error+($random_error_rate*($ins_percentage/($del_percentage+$ins_percentage+$mismatch_percentage)))/100;
                    my $NP_mismatch_rate = $extra_error+($random_error_rate*($mismatch_percentage/($del_percentage+$ins_percentage+$mismatch_percentage)))/100;
                    
                    if ($NP_error_rate eq '0')
                    {
                        $NP_deletion_rate = '0';
                        $NP_insertion_rate = '0';
                        $NP_mismatch_rate = '0';
                    }
                                          
#mismatch-----------------------------------------------------------------------------------------------          
                    my $mismatch_nuc_count = int($length_tmp*$NP_mismatch_rate);
                    my $f = '0';
                    my %mismatches;
                    undef %mismatches;

                    while ($f < $mismatch_nuc_count)
                    {                           
                        my $pos = int(rand($length_read-3));
                        if (exists($mismatches{$pos}))
                        {
                            next;
                        }
                        if (exists($mismatches{$pos+1}))
                        {
                            next;
                        }
                        if (exists($mismatches{$pos-1}))
                        {
                            next;
                        }
                 
                        my $nuc = substr $seq_tmp2, $pos, 1;
                        my $nuc_prev = substr $seq_tmp2, $pos-1, 1;
                        my $nuc_next = substr $seq_tmp2, $pos+1, 1;
                        my $mismatch_nuc = "";
                        my $no_extra = "";
                                
                        my $three_nuc = $nuc_prev.$nuc.$nuc_next;
                        my $chance = '1';
                        my $random_chance = rand(1);
                   
                        if (exists($NP_error_profile_mismatch2{$three_nuc}))
                        {
                            $chance = ($NP_error_profile_mismatch2{$three_nuc}*$NP_mismatch_multi2)/100;     
                        }
                        if ($random_chance > $chance)
                        {
                            next;
                        }
                     
                        my $four_nuc = $nuc_prev.$nuc."A".$nuc_next;
                        if (exists($NP_error_profile_mismatch{$four_nuc}))
                        {
                            my $random_chance2 = rand(1);                              
                            my $chance2 = $NP_error_profile_mismatch{$four_nuc}/100;
                            if ($random_chance2 <= $chance2 && $nuc ne "A")
                            {
                               $mismatch_nuc = "A"; 
                            }
                            else
                            {
                                $four_nuc = $nuc_prev.$nuc."C".$nuc_next;
                                my $chance3 = $NP_error_profile_mismatch{$four_nuc}/100;
                                if ($random_chance2 <= $chance2+$chance3 && $nuc ne "C")
                                {
                                   $mismatch_nuc = "C"; 
                                }
                                else
                                {
                                    $four_nuc = $nuc_prev.$nuc."T".$nuc_next;
                                    my $chance4 = $NP_error_profile_mismatch{$four_nuc}/100;
                                    if ($random_chance2 <= $chance2+$chance3+$chance4 && $nuc ne "T")
                                    {
                                       $mismatch_nuc = "T";
                                    }
                                    elsif ($nuc ne "G")
                                    {
                                        $mismatch_nuc = "G";
                                    }
                                }
                                $mismatches{$pos} = undef;
                            }
                        }
                        else
                        {
                            my $range = int(rand(4));
                            $mismatch_nuc = $nucs[$range];
                            $mismatches{$pos} = undef;
                            $no_extra = "yes";
                        }
                    
                        my $q = '1';
                        
                        while ($pos+$q < $length_tmp && $no_extra eq "")
                        {
                            $nuc = $nuc_next;
                            $nuc_next = substr $seq_tmp2, $pos+1+$q, 1;
                            $nuc_prev = $nuc;
                            
                            $three_nuc = $nuc_prev.$nuc.$nuc_next;
                            $random_chance = rand(1);
                            
                            if (exists($NP_error_profile_mismatch2{$three_nuc}))
                            {
                                $chance = ($NP_error_profile_mismatch2{$three_nuc}*$NP_mismatch_multi2)/100;
                            }
                            if ($random_chance > $chance)
                            {
                                last;
                            }
                            
                            my $random_chance2 = rand(1);
                            my $four_nuc = $nuc_prev.$nuc."A".$nuc_next;
                            my $chance2 = $NP_error_profile_mismatch{$four_nuc}/100;
                            if ($random_chance2 <= $chance2 && $nuc ne "A")
                            {
                               $mismatch_nuc .= "A"; 
                            }
                            else
                            {
                                $four_nuc = $nuc_prev.$nuc."C".$nuc_next;
                                my $chance3 = $NP_error_profile_mismatch{$four_nuc}/100;
                                if ($random_chance2 <= $chance2+$chance3 && $nuc ne "C")
                                {
                                   $mismatch_nuc .= "C"; 
                                }
                                else
                                {
                                    $four_nuc = $nuc_prev.$nuc."T".$nuc_next;
                                    my $chance4 = $NP_error_profile_mismatch{$four_nuc}/100;
                                    if ($random_chance2 <= $chance2+$chance3+$chance4 && $nuc ne "T")
                                    {
                                       $mismatch_nuc .= "T";
                                    }
                                    elsif ($nuc ne "G")
                                    {
                                        $mismatch_nuc .= "G";
                                    }
                                }
                                $mismatches{$pos+$q} = undef;
                            }
                            $q++;
                        }                     
                                                   
                        if ($mismatch_nuc ne "")
                        {
                            $f += length($mismatch_nuc);                         
                            substr $seq_tmp, $pos, length($mismatch_nuc), $mismatch_nuc;
                        }                         
                    }
#indel------------------------------------------------------------------------------------------------------           
                    my %indels;
                    undef %indels;
                    my $delete_nuc_count = int($length_tmp*$NP_deletion_rate);
                    my $s = '0';

                    while ($s < $delete_nuc_count)
                    {
                        my $pos = int(rand($length_read-$delete_nuc_count));
                        if (exists($mismatches{$pos}))
                        {
                            next;
                        }
                        my $three_nuc = substr $seq_tmp2, $pos-2, 3;
                        
                        if (exists($NP_error_profile_del{$three_nuc}))
                        {                   
                            my $random_chance = rand(1);
                            if ($random_chance <= (($NP_error_profile_del{$three_nuc}*$NP_del_multi)/100))
                            {       
                            }
                            else
                            {
                                next;
                            }
                        }

                        my $extra_del = rand(1);
                        my $length_tmp2 = '1';
                        my $threshold = 0.7;
                        if (exists ($NP_del_length{$length_tmp2}))
                        {
                            $threshold = $NP_del_length{$length_tmp2};
                        }
                        while ($extra_del > $threshold && $pos < length($seq_tmp2)+$length_tmp2 && $length_tmp2 < 30)
                        {
                            $length_tmp2++;
                            $extra_del = rand(1);
                            if (exists ($NP_del_length{$length_tmp2}))
                            {
                                $threshold = $NP_del_length{$length_tmp2};
                            }
                            $s++;
                        }
                        $indels{$pos-1} = $length_tmp2;
                        $s++;
                    }

                    my $r = '0';

                    my $insert_nuc_count = int($length_tmp*$NP_insertion_rate);
                    while ($r < $insert_nuc_count)
                    {                                          
                        my $pos = rand($length_read);
                        if (exists($mismatches{$pos}))
                        {
                            next;
                        }
                        my $nuc_prev = substr $seq_tmp2, $pos, 1;
                        my $nuc_next = substr $seq_tmp2, $pos+1, 1;
                        my $insert_seq = "";
                        
                        
                        my $two_nuc = $nuc_prev.$nuc_next;
                        my $chance = '1';
                        my $random_chance = rand(1);
                    
                        if (exists($NP_error_profile_ins2{$two_nuc}))
                        {
                            $chance = ($NP_error_profile_ins2{$two_nuc}*$NP_ins_multi2)/100;
                        }
                        if ($random_chance > $chance)
                        {
                            next;
                        }
                        elsif (exists($NP_error_profile_ins2{$two_nuc}))
                        {
                            my $random_chance2 = rand(1);
                            my $three_nuc = $nuc_prev."A".$nuc_next;
                            my $three_nuc2 = $nuc_prev."C".$nuc_next;
                            my $three_nuc3 = $nuc_prev."T".$nuc_next;
                            my $three_nuc4 = $nuc_prev."G".$nuc_next;
                            my $chance2 = ($NP_error_profile_ins{$three_nuc})/
                            ($NP_error_profile_ins{$three_nuc}+$NP_error_profile_ins{$three_nuc2}+$NP_error_profile_ins{$three_nuc3}+$NP_error_profile_ins{$three_nuc4});

                            if ($random_chance2 <= $chance2)
                            {
                               $insert_seq = "A"; 
                            }
                            else
                            {
                                my $chance3 = ($NP_error_profile_ins{$three_nuc2})/
                            ($NP_error_profile_ins{$three_nuc}+$NP_error_profile_ins{$three_nuc2}+$NP_error_profile_ins{$three_nuc3}+$NP_error_profile_ins{$three_nuc4});
                                if ($random_chance2 <= $chance2+$chance3)
                                {
                                   $insert_seq = "C"; 
                                }
                                else
                                {
                                    my $chance4 = ($NP_error_profile_ins{$three_nuc3})/
                            ($NP_error_profile_ins{$three_nuc}+$NP_error_profile_ins{$three_nuc2}+$NP_error_profile_ins{$three_nuc3}+$NP_error_profile_ins{$three_nuc4});
                                    if ($random_chance2 <= $chance2+$chance3+$chance4)
                                    {
                                       $insert_seq = "T";
                                    }
                                    else
                                    {
                                        $insert_seq = "G";
                                    }
                                }
                            }
                        }
                        else
                        {
                            my $range = int(rand(4));
                            $insert_seq = $nucs[$range];
                            $mismatches{$pos} = undef;
                        }
                        
                        my $extra_ins = rand(1);
                        my $threshold = 0.8;
                        if (exists ($NP_ins_length{length($insert_seq)}))
                        {
                            $threshold = $NP_ins_length{length($insert_seq)};
                        }
                        while ($extra_ins > $threshold && $pos < length($seq_tmp2))
                        {
                            my $random_nuc2 = int(rand(4));
                            $insert_seq .= $nucs[$random_nuc2-1];
                            $extra_ins = rand(1);
                            if (exists ($NP_ins_length{length($insert_seq)}))
                            {
                                $threshold = $NP_ins_length{length($insert_seq)};
                            }
                            $r++;
                        }
                        $indels{$pos} = $insert_seq;
                        $r++;    
                    }
#-----------------------------------------------------------------------------------------------------------                                   
                    
                    foreach my $indel_tmp (sort {$b <=> $a} keys %indels)
                    {                        
                        if ($indels{$indel_tmp} > 0)
                        {
                            substr $seq_tmp, $indel_tmp, $indels{$indel_tmp}, "";
                        }
                        else
                        {
                            my $one = substr $seq_tmp, 0, $indel_tmp;
                            my $two = substr $seq_tmp, $indel_tmp;
                            $seq_tmp = $one.$indels{$indel_tmp}.$two;
                        }
                    }

                    print OUTPUT_NP ">".$project."_".$NP_read_count."_length=".length($seq_tmp)."_HAP0\n";
                    print OUTPUT_NP $seq_tmp."\n";

                    $NP_seq_length1{$shortest_seq1+$random_NP_length} = undef;
                    $o++;
                }
                
                if ($d eq '0')
                {
                    foreach my $seq_length (sort {$a <=> $b} keys %NP_seq_length1)
                    {           
                        $shortest_seq1 = $seq_length;
                        delete $NP_seq_length1{$seq_length};           
                        $d++;
                        last;
                    }
                }
                    
                my $shortest_seq_tmp = '0';

                $shortest_seq_tmp = $shortest_seq1-$seq_done1;
                $seq_done1 += $shortest_seq_tmp;
                
                substr $haplotype1, 0, $shortest_seq_tmp, "";
            }
            if ($new_contig ne "")
            {
                $new_contig = "yes2";
            }
        }  
        if ($check_haplo1_start eq "yes")
        {         
            if ($NP_coverage2 > 0 && length($haplotype2) > 0)
            {
                while (((length($haplotype2) >= $NP_range_high && ($finish_var_np eq "yes" || $reference_size2 < $random_length_interval || $SV_input eq ""))
                        || ($new_contig eq "yes2" && length($haplotype2) > 0)))
                {
                    my $o = '0';
                    my $d = '0';
                    
                    my $count_seqs = keys %NP_seq_length2;
           
                    while ($o <= $NP_coverage2-$count_seqs)
                    {
                        $NP_read_count++;
NP_LENGTH2:      
                        my $random_NP_length = '0';
                        
                        if (exists($NP_extra_long_reads{$NP_read_count}))
                        {
                             my $range_tmp = $NP_range_high-($NP_average*2);
                             $random_NP_length = int(rand($range_tmp)) + ($NP_average*2);
                        }
                        else
                        {
                            my $j = (($NP_average*2)+400)/3;
                            my $t=0;
                            $t+=rand() for(1..3);
                            $random_NP_length = int($t*$j+1);
                        }
                        
                        if ($random_NP_length < $NP_min_read_length)
                        {
                           $NP_min_read_length = $random_NP_length; 
                        }
                        if ($random_NP_length > $NP_max_read_length)
                        {
                            $NP_max_read_length = $random_NP_length;
                        }
                        $NP_total_length += $random_NP_length;
            
                        if (exists($NP_seq_length2{$random_NP_length+$shortest_seq2}))
                        {
                            goto NP_LENGTH2;
                        }

                        my $chance_sense = int(rand(2));
                        my $seq_tmp = "";
                        if ($chance_sense eq '0')
                        {
                            $seq_tmp = reverse(substr $haplotype2, 0, $random_NP_length-$start_read_lq);
                            $seq_tmp =~ tr/ACTG/TGAC/;
                        }
                        else
                        {
                            $seq_tmp = substr $haplotype2, 0, $random_NP_length-$start_read_lq;
                        }
                        my $seq_tmp2 = $seq_tmp;
                        my $length_read = length($seq_tmp2);

#input sequencing errors---------------------------------------------------------------------------------------------------------------------              

                        my $length_tmp = length($seq_tmp);                      
                        my $j = 1.2;
                        my $t= rand(1.5);
                        $t+= 1;
                        my $random_error_rate = int(($NP_error_rate/($t*$j))+$NP_error_rate-($NP_error_rate/2.2));
                        my $extra_error = rand(($NP_error_rate/12)/100);
                        
                        my $NP_deletion_rate = $extra_error+($random_error_rate*($del_percentage/($del_percentage+$ins_percentage+$mismatch_percentage)))/100;
                        my $NP_insertion_rate = $extra_error+($random_error_rate*($ins_percentage/($del_percentage+$ins_percentage+$mismatch_percentage)))/100;
                        my $NP_mismatch_rate = $extra_error+($random_error_rate*($mismatch_percentage/($del_percentage+$ins_percentage+$mismatch_percentage)))/100;             
                        
                        if ($NP_error_rate eq '0')
                        {
                            $NP_deletion_rate = '0';
                            $NP_insertion_rate = '0';
                            $NP_mismatch_rate = '0';
                        }
                        
#mismatch-----------------------------------------------------------------------------------------------          
                        my $mismatch_nuc_count = int($length_tmp*$NP_mismatch_rate);
                        my $f = '0';
                        my %mismatches;
                        undef %mismatches;

                        while ($f < $mismatch_nuc_count)
                        {                           
                            my $pos = int(rand($length_read-3));
                            if (exists($mismatches{$pos}))
                            {
                                next;
                            }
                            if (exists($mismatches{$pos+1}))
                            {
                                next;
                            }
                            if (exists($mismatches{$pos-1}))
                            {
                                next;
                            }
                     
                            my $nuc = substr $seq_tmp2, $pos, 1;
                            my $nuc_prev = substr $seq_tmp2, $pos-1, 1;
                            my $nuc_next = substr $seq_tmp2, $pos+1, 1;
                            my $mismatch_nuc = "";
                            my $no_extra = "";
                                    
                            my $three_nuc = $nuc_prev.$nuc.$nuc_next;
                            my $chance = '1';
                            my $random_chance = rand(1);
                       
                            if (exists($NP_error_profile_mismatch2{$three_nuc}))
                            {
                                $chance = ($NP_error_profile_mismatch2{$three_nuc}*$NP_mismatch_multi2)/100;     
                            }
                            if ($random_chance > $chance)
                            {
                                next;
                            }
                         
                            my $four_nuc = $nuc_prev.$nuc."A".$nuc_next;
                            if (exists($NP_error_profile_mismatch{$four_nuc}))
                            {
                                my $random_chance2 = rand(1);                              
                                my $chance2 = $NP_error_profile_mismatch{$four_nuc}/100;
                                if ($random_chance2 <= $chance2 && $nuc ne "A")
                                {
                                   $mismatch_nuc = "A"; 
                                }
                                else
                                {
                                    $four_nuc = $nuc_prev.$nuc."C".$nuc_next;
                                    my $chance3 = $NP_error_profile_mismatch{$four_nuc}/100;
                                    if ($random_chance2 <= $chance2+$chance3 && $nuc ne "C")
                                    {
                                       $mismatch_nuc = "C"; 
                                    }
                                    else
                                    {
                                        $four_nuc = $nuc_prev.$nuc."T".$nuc_next;
                                        my $chance4 = $NP_error_profile_mismatch{$four_nuc}/100;
                                        if ($random_chance2 <= $chance2+$chance3+$chance4 && $nuc ne "T")
                                        {
                                           $mismatch_nuc = "T";
                                        }
                                        elsif ($nuc ne "G")
                                        {
                                            $mismatch_nuc = "G";
                                        }
                                    }
                                    $mismatches{$pos} = undef;
                                }
                            }
                            else
                            {
                                my $range = int(rand(4));
                                $mismatch_nuc = $nucs[$range];
                                $mismatches{$pos} = undef;
                                $no_extra = "yes";
                            }
                        
                            my $q = '1';
                            
                            while ($pos+$q < $length_tmp && $no_extra eq "")
                            {
                                $nuc = $nuc_next;
                                $nuc_next = substr $seq_tmp2, $pos+1+$q, 1;
                                $nuc_prev = $nuc;
                                
                                $three_nuc = $nuc_prev.$nuc.$nuc_next;
                                $random_chance = rand(1);
                                
                                if (exists($NP_error_profile_mismatch2{$three_nuc}))
                                {
                                    $chance = ($NP_error_profile_mismatch2{$three_nuc}*$NP_mismatch_multi2)/100;
                                }
                                if ($random_chance > $chance)
                                {
                                    last;
                                }
                                
                                my $random_chance2 = rand(1);
                                my $four_nuc = $nuc_prev.$nuc."A".$nuc_next;
                                my $chance2 = $NP_error_profile_mismatch{$four_nuc}/100;
                                if ($random_chance2 <= $chance2 && $nuc ne "A")
                                {
                                   $mismatch_nuc .= "A"; 
                                }
                                else
                                {
                                    $four_nuc = $nuc_prev.$nuc."C".$nuc_next;
                                    my $chance3 = $NP_error_profile_mismatch{$four_nuc}/100;
                                    if ($random_chance2 <= $chance2+$chance3 && $nuc ne "C")
                                    {
                                       $mismatch_nuc .= "C"; 
                                    }
                                    else
                                    {
                                        $four_nuc = $nuc_prev.$nuc."T".$nuc_next;
                                        my $chance4 = $NP_error_profile_mismatch{$four_nuc}/100;
                                        if ($random_chance2 <= $chance2+$chance3+$chance4 && $nuc ne "T")
                                        {
                                           $mismatch_nuc .= "T";
                                        }
                                        elsif ($nuc ne "G")
                                        {
                                            $mismatch_nuc .= "G";
                                        }
                                    }
                                    $mismatches{$pos+$q} = undef;
                                }
                                $q++;
                            }                     
                                                       
                            if ($mismatch_nuc ne "")
                            {
                                $f += length($mismatch_nuc);                         
                                substr $seq_tmp, $pos, length($mismatch_nuc), $mismatch_nuc;
                            }                         
                        }
#indel------------------------------------------------------------------------------------------------------           
                        my %indels;
                        undef %indels;
                        my @indels;
                        undef @indels;
                        my $delete_nuc_count = int($length_tmp*$NP_deletion_rate);
                        my $s = '0';

                        while ($s < $delete_nuc_count)
                        {
                            my $pos = int(rand($length_read-$delete_nuc_count));
                            if (exists($mismatches{$pos}))
                            {
                                next;
                            }
                            my $three_nuc = substr $seq_tmp2, $pos-2, 3;
                            
                            if (exists($NP_error_profile_del{$three_nuc}))
                            {                   
                                my $random_chance = rand(1);
                                if ($random_chance <= (($NP_error_profile_del{$three_nuc}*$NP_del_multi)/100))
                                {       
                                }
                                else
                                {
                                    next;
                                }
                            }

                            my $extra_del = rand(1);
                            my $length_tmp2 = '1';
                            my $threshold = 0.7;
                            if (exists ($NP_del_length{$length_tmp2}))
                            {
                                $threshold = $NP_del_length{$length_tmp2};
                            }
                            while ($extra_del > $threshold && $pos < length($seq_tmp2)+$length_tmp2 && $length_tmp2 < 30)
                            {
                                $length_tmp2++;
                                $extra_del = rand(1);
                                if (exists ($NP_del_length{$length_tmp2}))
                                {
                                    $threshold = $NP_del_length{$length_tmp2};
                                }
                                $s++;
                            }
                            $indels{$pos-1} = $length_tmp2;
                           
                            $s++;
                        }

                        my $r = '0';

                        my $insert_nuc_count = int($length_tmp*$NP_insertion_rate);
                        while ($r < $insert_nuc_count)
                        {                                          
                            my $pos = rand($length_read);
                            if (exists($mismatches{$pos}))
                            {
                                next;
                            }
                            my $nuc_prev = substr $seq_tmp2, $pos, 1;
                            my $nuc_next = substr $seq_tmp2, $pos+1, 1;
                            my $insert_seq = "";
                                                     
                            my $two_nuc = $nuc_prev.$nuc_next;
                            my $chance = '1';
                            my $random_chance = rand(1);
                        
                            if (exists($NP_error_profile_ins2{$two_nuc}))
                            {
                                $chance = ($NP_error_profile_ins2{$two_nuc}*$NP_ins_multi2)/100;
                            }
                            if ($random_chance > $chance)
                            {
                                next;
                            }
                            elsif (exists($NP_error_profile_ins2{$two_nuc}))
                            {
                                my $random_chance2 = rand(1);
                                my $three_nuc = $nuc_prev."A".$nuc_next;
                                my $three_nuc2 = $nuc_prev."C".$nuc_next;
                                my $three_nuc3 = $nuc_prev."T".$nuc_next;
                                my $three_nuc4 = $nuc_prev."G".$nuc_next;
                                my $chance2 = ($NP_error_profile_ins{$three_nuc})/
                                ($NP_error_profile_ins{$three_nuc}+$NP_error_profile_ins{$three_nuc2}+$NP_error_profile_ins{$three_nuc3}+$NP_error_profile_ins{$three_nuc4});

                                if ($random_chance2 <= $chance2)
                                {
                                   $insert_seq = "A"; 
                                }
                                else
                                {
                                    my $chance3 = ($NP_error_profile_ins{$three_nuc2})/
                                ($NP_error_profile_ins{$three_nuc}+$NP_error_profile_ins{$three_nuc2}+$NP_error_profile_ins{$three_nuc3}+$NP_error_profile_ins{$three_nuc4});
                                    if ($random_chance2 <= $chance2+$chance3)
                                    {
                                       $insert_seq = "C"; 
                                    }
                                    else
                                    {
                                        my $chance4 = ($NP_error_profile_ins{$three_nuc3})/
                                ($NP_error_profile_ins{$three_nuc}+$NP_error_profile_ins{$three_nuc2}+$NP_error_profile_ins{$three_nuc3}+$NP_error_profile_ins{$three_nuc4});
                                        if ($random_chance2 <= $chance2+$chance3+$chance4)
                                        {
                                           $insert_seq = "T";
                                        }
                                        else
                                        {
                                            $insert_seq = "G";
                                        }
                                    }
                                }
                            }
                            else
                            {
                                my $range = int(rand(4));
                                $insert_seq = $nucs[$range];
                                $mismatches{$pos} = undef;
                            }
                            
                            my $extra_ins = rand(1);
                            my $threshold = 0.8;
                            if (exists ($NP_ins_length{length($insert_seq)}))
                            {
                                $threshold = $NP_ins_length{length($insert_seq)};
                            }
                            while ($extra_ins > $threshold && $pos < length($seq_tmp2))
                            {
                                my $random_nuc2 = int(rand(4));
                                $insert_seq .= $nucs[$random_nuc2-1];
                                $extra_ins = rand(1);
                                if (exists ($NP_ins_length{length($insert_seq)}))
                                {
                                    $threshold = $NP_ins_length{length($insert_seq)};
                                }
                                $r++;
                            }
                            $indels{$pos} = $insert_seq;
                            $r++;    
                        }
#-----------------------------------------------------------------------------------------------------------                                   
                        
                        foreach my $indel_tmp (sort {$b <=> $a} keys %indels)
                        {                        
                            if ($indels{$indel_tmp} > 0)
                            {
                                substr $seq_tmp, $indel_tmp, $indels{$indel_tmp}, "";
                            }
                            else
                            {
                                my $one = substr $seq_tmp, 0, $indel_tmp;
                                my $two = substr $seq_tmp, $indel_tmp;
                                $seq_tmp = $one.$indels{$indel_tmp}.$two;
                            }
                        }
 
                        print OUTPUT_NP ">".$project."_".$NP_read_count."_length=".length($seq_tmp)."_HAP1\n";
                        print OUTPUT_NP $seq_tmp."\n";

                        $NP_seq_length2{$shortest_seq2+$random_NP_length} = undef;
                        $o++;
                    }

                    if ($d eq '0')
                    {
                        foreach my $seq_length (sort {$a <=> $b} keys %NP_seq_length2)
                        {                         
                            $shortest_seq2 = $seq_length;
                            delete $NP_seq_length2{$seq_length};
                            last;
                        }
                        $d++;
                    }
            
                    my $shortest_seq_tmp = '0';

                    $shortest_seq_tmp = $shortest_seq2-$seq_done2;
                    $seq_done2 += $shortest_seq_tmp;
                    
                    substr $haplotype2, 0, $shortest_seq_tmp, "";
                }
                if ($new_contig ne "")
                {
                    $new_contig = "yes2";
                }
            }
        }
#No heterozygosity-------------------------------------------------------
#------------------------------------------------------------------------
        if ($heterozygosity eq "no")
        {                   
            if ($NP_coverage > 0 && length($haplo_merged) > 0)
            {        
                while (((length($haplo_merged) >= $NP_range_high && ($finish_var_np eq "yes" || $reference_size2 < $random_length_interval || $SV_input eq ""))
                        || ($new_contig eq "yes" && length($haplo_merged) > 0)))
                {
                    my $o = '0';
                    my $d = '0';
                    
                    my $count_seqs = keys %NP_seq_length;
           
                    while ($o <= $NP_coverage-$count_seqs)
                    {
                        $NP_read_count++;
NP_LENGTH:      
                        my $random_NP_length = '0';
                        
                        if (exists($NP_extra_long_reads{$NP_read_count}))
                        {
                             my $range_tmp = $NP_range_high-($NP_average*2);
                             $random_NP_length = int(rand($range_tmp)) + ($NP_average*2);
                        }
                        else
                        {
                            my $j = (($NP_average*2)+400)/3;
                            my $t=0;
                            $t+=rand() for(1..3);
                            $random_NP_length = int($t*$j+1);
                        }
                        
                        if ($random_NP_length < $NP_min_read_length)
                        {
                           $NP_min_read_length = $random_NP_length; 
                        }
                        if ($random_NP_length > $NP_max_read_length)
                        {
                            $NP_max_read_length = $random_NP_length;
                        }
                        $NP_total_length += $random_NP_length;
            
                        if (exists($NP_seq_length{$random_NP_length+$shortest_seq}))
                        {
                            goto NP_LENGTH;
                        }

                        #print OUTPUT_NP $random_NP_length." RANDOM_LENGTH\n";
                        
                        #my $NP_start = substr $haplo_merged, 0, $start_read_lq;
                        #my $NP_read_tmp = substr $haplo_merged, $start_read_lq, $random_NP_length-$start_read_lq;
                        #my $NP_start_check = "";
                        #my $seq_tmp = $NP_start;
                        #my $seq_tmp2 = $NP_start;
                        my $chance_sense = int(rand(2));
                        my $seq_tmp = "";
                        if ($chance_sense eq '0')
                        {
                            $seq_tmp = reverse(substr $haplo_merged, 0, $random_NP_length-$start_read_lq);
                            $seq_tmp =~ tr/ACTG/TGAC/;
                        }
                        else
                        {
                            $seq_tmp = substr $haplo_merged, 0, $random_NP_length-$start_read_lq;
                        }
                        my $seq_tmp2 = $seq_tmp;
                        my $length_read = length($seq_tmp2);
                        #my $extra_error = '0';
                        #my $extra_error = '30';
                        #if ($start_read_lq < 50)
                        #{
                            #$extra_error = '10';
                        #}
                        
                        #my @NP_read_tmp = split //, $seq_tmp;
                               
#input sequencing errors---------------------------------------------------------------------------------------------------------------------
NP_FIRST_50:                
                        #if ($NP_start_check eq "yes")
                        #{
                            #$seq_tmp = $NP_read_tmp;
                            #$seq_tmp2 = $NP_read_tmp;
                           # $extra_error = '0';
                        #}
                        my $length_tmp = length($seq_tmp);                        
                        my $j = 1.2;
                        my $t= rand(1.5);
                        $t+= 1;
                        my $random_error_rate = int(($NP_error_rate/($t*$j))+$NP_error_rate-($NP_error_rate/2.2));
                        my $extra_error = rand(($NP_error_rate/12)/100);

                        my $NP_deletion_rate = $extra_error+($random_error_rate*($del_percentage/($del_percentage+$ins_percentage+$mismatch_percentage)))/100;
                        my $NP_insertion_rate = $extra_error+($random_error_rate*($ins_percentage/($del_percentage+$ins_percentage+$mismatch_percentage)))/100;
                        my $NP_mismatch_rate = $extra_error+($random_error_rate*($mismatch_percentage/($del_percentage+$ins_percentage+$mismatch_percentage)))/100;          
                        
                        if ($NP_error_rate eq '0')
                        {
                            $NP_deletion_rate = '0';
                            $NP_insertion_rate = '0';
                            $NP_mismatch_rate = '0';
                        }
                     
#mismatch-----------------------------------------------------------------------------------------------          
                        my $mismatch_nuc_count = int($length_tmp*$NP_mismatch_rate);
                        my $f = '0';
                        my %mismatches;
                        undef %mismatches;
    
#my $time_np_start_loop = time;
                        while ($f < $mismatch_nuc_count)
                        {                           
                            my $pos = int(rand($length_read-3));
                            if (exists($mismatches{$pos}))
                            {
                                next;
                            }
                            if (exists($mismatches{$pos+1}))
                            {
                                next;
                            }
                            if (exists($mismatches{$pos-1}))
                            {
                                next;
                            }
                     
                            my $nuc = substr $seq_tmp2, $pos, 1;
                            my $nuc_prev = substr $seq_tmp2, $pos-1, 1;
                            my $nuc_next = substr $seq_tmp2, $pos+1, 1;
                            my $mismatch_nuc = "";
                            my $no_extra = "";
                                    
                            my $three_nuc = $nuc_prev.$nuc.$nuc_next;
                            my $chance = '1';
                            my $random_chance = rand(1);
                       
                            if (exists($NP_error_profile_mismatch2{$three_nuc}))
                            {
                                $chance = ($NP_error_profile_mismatch2{$three_nuc}*$NP_mismatch_multi2)/100;     
                            }
                            if ($random_chance > $chance)
                            {
                                next;
                            }
                         
                            my $four_nuc = $nuc_prev.$nuc."A".$nuc_next;
                            if (exists($NP_error_profile_mismatch{$four_nuc}))
                            {
                                my $random_chance2 = rand(1);                              
                                my $chance2 = $NP_error_profile_mismatch{$four_nuc}/100;
                                if ($random_chance2 <= $chance2 && $nuc ne "A")
                                {
                                   $mismatch_nuc = "A"; 
                                }
                                else
                                {
                                    $four_nuc = $nuc_prev.$nuc."C".$nuc_next;
                                    my $chance3 = $NP_error_profile_mismatch{$four_nuc}/100;
                                    if ($random_chance2 <= $chance2+$chance3 && $nuc ne "C")
                                    {
                                       $mismatch_nuc = "C"; 
                                    }
                                    else
                                    {
                                        $four_nuc = $nuc_prev.$nuc."T".$nuc_next;
                                        my $chance4 = $NP_error_profile_mismatch{$four_nuc}/100;
                                        if ($random_chance2 <= $chance2+$chance3+$chance4 && $nuc ne "T")
                                        {
                                           $mismatch_nuc = "T";
                                        }
                                        elsif ($nuc ne "G")
                                        {
                                            $mismatch_nuc = "G";
                                        }
                                    }
                                    $mismatches{$pos} = undef;
                                }
                            }
                            else
                            {
                                my $range = int(rand(4));
                                $mismatch_nuc = $nucs[$range];
                                $mismatches{$pos} = undef;
                                $no_extra = "yes";
                            }
                        
                            my $q = '1';
                            
                            while ($pos+$q < $length_tmp && $no_extra eq "")
                            {
                                $nuc = $nuc_next;
                                $nuc_next = substr $seq_tmp2, $pos+1+$q, 1;
                                $nuc_prev = $nuc;
                                
                                $three_nuc = $nuc_prev.$nuc.$nuc_next;
                                $random_chance = rand(1);
                                
                                if (exists($NP_error_profile_mismatch2{$three_nuc}))
                                {
                                    $chance = ($NP_error_profile_mismatch2{$three_nuc}*$NP_mismatch_multi2)/100;
                                }
                                if ($random_chance > $chance)
                                {
                                    last;
                                }
                                
                                my $random_chance2 = rand(1);
                                my $four_nuc = $nuc_prev.$nuc."A".$nuc_next;
                                my $chance2 = $NP_error_profile_mismatch{$four_nuc}/100;
                                if ($random_chance2 <= $chance2 && $nuc ne "A")
                                {
                                   $mismatch_nuc .= "A"; 
                                }
                                else
                                {
                                    $four_nuc = $nuc_prev.$nuc."C".$nuc_next;
                                    my $chance3 = $NP_error_profile_mismatch{$four_nuc}/100;
                                    if ($random_chance2 <= $chance2+$chance3 && $nuc ne "C")
                                    {
                                       $mismatch_nuc .= "C"; 
                                    }
                                    else
                                    {
                                        $four_nuc = $nuc_prev.$nuc."T".$nuc_next;
                                        my $chance4 = $NP_error_profile_mismatch{$four_nuc}/100;
                                        if ($random_chance2 <= $chance2+$chance3+$chance4 && $nuc ne "T")
                                        {
                                           $mismatch_nuc .= "T";
                                        }
                                        elsif ($nuc ne "G")
                                        {
                                            $mismatch_nuc .= "G";
                                        }
                                    }
                                    $mismatches{$pos+$q} = undef;
                                }
                                #$extra_mismatch = rand(1);
                                $q++;
                            }                     
                                                       
                            if ($mismatch_nuc ne "")
                            {
                                $f += length($mismatch_nuc);                         
                                substr $seq_tmp, $pos, length($mismatch_nuc), $mismatch_nuc;
                            }                         
                        }
          #$time_np1 += time-$time_np_start_loop;               
                        #foreach my $mismatch (sort {$a <=> $b} keys %mismatches)
                        #{
                            #substr $seq_tmp, $mismatch, 1, $mismatches{$mismatch};
                        #}                  
                 
                #my $time_np_indel = time;
#indel------------------------------------------------------------------------------------------------------           
                        my %indels;
                        undef %indels;
                        my %dels;
                        undef %dels;
                        my %inserts;
                        undef %inserts;
                        my @indels;
                        undef @indels;
                        my $delete_nuc_count = int($length_tmp*$NP_deletion_rate);
                        my $s = '0';

                        while ($s < $delete_nuc_count)
                        {
                            my $pos = int(rand($length_read-$delete_nuc_count));
                            if (exists($mismatches{$pos}))
                            {
                                next;
                            }
                            my $three_nuc = substr $seq_tmp2, $pos-2, 3;
                            
                            if (exists($NP_error_profile_del{$three_nuc}))
                            {                   
                                my $random_chance = rand(1);
                                if ($random_chance <= (($NP_error_profile_del{$three_nuc}*$NP_del_multi)/100))
                                {       
                                }
                                else
                                {
                                    next;
                                }
                            }

                            my $extra_del = rand(1);
                            my $length_tmp2 = '1';
                            my $threshold = 0.7;
                            if (exists ($NP_del_length{$length_tmp2}))
                            {
                                $threshold = $NP_del_length{$length_tmp2};
                            }
                            while ($extra_del > $threshold && $pos < length($seq_tmp2)+$length_tmp2 && $length_tmp2 < 30)
                            {
                                $length_tmp2++;
                                $extra_del = rand(1);
                                if (exists ($NP_del_length{$length_tmp2}))
                                {
                                    $threshold = $NP_del_length{$length_tmp2};
                                }
                                $s++;
                            }
                            #substr $seq_tmp, $pos, $length_tmp2, "";
                            #$indels{$pos-1} = $length_tmp2;
                            #splice @NP_read_tmp, $pos, 1+$extra;
                            push @indels, $pos-1;
                            $dels{$pos-1} = $length_tmp2;
                            $s++;
                        }

                        my $r = '0';
                #my $time_np_last = time;

                        my $insert_nuc_count = int($length_tmp*$NP_insertion_rate);
                        while ($r < $insert_nuc_count)
                        {                                          
                            my $pos = rand($length_read);
                            if (exists($mismatches{$pos}))
                            {
                                next;
                            }
                            my $nuc_prev = substr $seq_tmp2, $pos, 1;
                            my $nuc_next = substr $seq_tmp2, $pos+1, 1;
                            my $insert_seq_tmp = "";                          
                            
                            my $two_nuc = $nuc_prev.$nuc_next;
                            my $chance = '1';
                            my $random_chance = rand(1);
                        
                            if (exists($NP_error_profile_ins2{$two_nuc}))
                            {
                                $chance = ($NP_error_profile_ins2{$two_nuc}*$NP_ins_multi2)/100;
                            }
                            if ($random_chance > $chance)
                            {
                                next;
                            }
                            elsif (exists($NP_error_profile_ins2{$two_nuc}))
                            {
                                my $random_chance2 = rand(1);
                                my $three_nuc = $nuc_prev."A".$nuc_next;
                                my $three_nuc2 = $nuc_prev."C".$nuc_next;
                                my $three_nuc3 = $nuc_prev."T".$nuc_next;
                                my $three_nuc4 = $nuc_prev."G".$nuc_next;
                                my $chance2 = ($NP_error_profile_ins{$three_nuc})/
                                ($NP_error_profile_ins{$three_nuc}+$NP_error_profile_ins{$three_nuc2}+$NP_error_profile_ins{$three_nuc3}+$NP_error_profile_ins{$three_nuc4});

                                if ($random_chance2 <= $chance2)
                                {
                                   $insert_seq_tmp = "A"; 
                                }
                                else
                                {
                                    my $chance3 = ($NP_error_profile_ins{$three_nuc2})/
                                ($NP_error_profile_ins{$three_nuc}+$NP_error_profile_ins{$three_nuc2}+$NP_error_profile_ins{$three_nuc3}+$NP_error_profile_ins{$three_nuc4});
                                    if ($random_chance2 <= $chance2+$chance3)
                                    {
                                       $insert_seq_tmp = "C"; 
                                    }
                                    else
                                    {
                                        my $chance4 = ($NP_error_profile_ins{$three_nuc3})/
                                ($NP_error_profile_ins{$three_nuc}+$NP_error_profile_ins{$three_nuc2}+$NP_error_profile_ins{$three_nuc3}+$NP_error_profile_ins{$three_nuc4});
                                        if ($random_chance2 <= $chance2+$chance3+$chance4)
                                        {
                                           $insert_seq_tmp = "T";
                                        }
                                        else
                                        {
                                            $insert_seq_tmp = "G";
                                        }
                                    }
                                }
                            }
                            else
                            {
                                my $range = int(rand(4));
                                $insert_seq_tmp = $nucs[$range];
                                $mismatches{$pos} = undef;
                            }
                               
                            my $extra_ins = rand(1);
                            my $threshold = 0.8;
                            if (exists ($NP_ins_length{length($insert_seq_tmp)}))
                            {
                                $threshold = $NP_ins_length{length($insert_seq_tmp)};
                            }
                            while ($extra_ins > $threshold && $pos < length($seq_tmp2))
                            {
                                my $random_nuc2 = int(rand(4));
                                $insert_seq_tmp .= $nucs[$random_nuc2-1];
                                $extra_ins = rand(1);
                                if (exists ($NP_ins_length{length($insert_seq_tmp)}))
                                {
                                    $threshold = $NP_ins_length{length($insert_seq_tmp)};
                                }
                                $r++;
                            }
                            #my $one = substr $seq_tmp, 0, $pos;
                            #my $two = substr $seq_tmp, $pos;
                            #$seq_tmp = $one.$insert_seq.$two;
                            #substr $seq_tmp, $pos, 0, $insert_seq;
                            #$indels{$pos} = $insert_seq_tmp;
                            push @indels, $pos;
                            $inserts{$pos} = $insert_seq_tmp;
                            $r++;    
                        }
                #$time_np3 += time-$time_np_last;
                #$time_np2 += time-$time_np_indel;
#-----------------------------------------------------------------------------------------------------------                                   
                        #if ($NP_start_check eq "")
                        #{
                            #$NP_start = $seq_tmp;
                            #$NP_start_check = "yes";
                           # goto NP_FIRST_50;
                        #}
                        foreach my $indel_tmp (sort {$b <=> $a} @indels)
                        {                        
                            if (exists($dels{$indel_tmp}))
                            {
                                substr $seq_tmp, $indel_tmp, $dels{$indel_tmp}, "";
                            }
                            elsif (exists($inserts{$indel_tmp}))
                            {
                                my $one = substr $seq_tmp, 0, $indel_tmp;
                                my $two = substr $seq_tmp, $indel_tmp;
                                $seq_tmp = $one.$inserts{$indel_tmp}.$two;
                            }
                        }
                        
                        #my $NP_read = $NP_start.$seq_tmp;
                        
                        #print OUTPUT_NP ">".$project."_".$NP_read_count."_length=".length($NP_read)." ".$two_haplos." ".$NP_coverage_tmp." ".length($haplo_tmp)." ".$o." ".$NP_range_high."\n";
                        print OUTPUT_NP ">".$project."_".$NP_read_count."_length=".length($seq_tmp)."\n";
                        print OUTPUT_NP $seq_tmp."\n";

                        $NP_seq_length{$shortest_seq+$random_NP_length} = undef;
                        $o++;  
                    }
                    if ($d eq '0')
                    {
                        foreach my $seq_length (sort {$a <=> $b} keys %NP_seq_length)
                        {
                            $shortest_seq = $seq_length;
                            delete $NP_seq_length{$seq_length};
                            $d++;
                            last;
                        }
                    }
                    my $shortest_seq_tmp = '0';

                    $shortest_seq_tmp = $shortest_seq-$seq_done;
                    $seq_done += $shortest_seq_tmp;
                    
                    #print OUTPUT_NP $shortest_seq_tmp." SHORTEST_SEQ\n";
                    substr $haplo_merged, 0, $shortest_seq_tmp, "";
                }
                if ($new_contig ne "")
                {
                    $new_contig = "yes2";
                }
            }
        }  
        $finish_var_np = "";
    }
    if (eof)
    {}
    elsif ($new_contig ne "")
    {
        goto NEW_CONTIG0;
    }
    #if ($reference_size2 > 185000)
    #{
        #last;
    #}
    #$time_np += time-$time_np_start;
}

close $FILE_REF;

print "\n\n";
#print STATS------------------------------------------------------------------------------------------------------------------------
#print "TIME_SV = ".$time_sv."\n";
#print "TIME_NP = ".$time_np."\n";
#print "TIME_MISMATCH = ".$time_np1."\n";
#print "TIME_INDEL = ".$time_np2."\n";
#print "TIME_DEL = ".$time_np3."\n";

if ($NP_coverage > 0)
{
    my $NP_average_read_length = $NP_total_length/$NP_read_count;
    print "\nLong read simulation\n";
    print "---------------------\n";
    print "Read count          : ".$NP_read_count."\n";
    print "Average read length : ".int($NP_average_read_length)." bp\n";
    print "Shortest read       : ".$NP_min_read_length." bp\n";
    print "Longest read        : ".$NP_max_read_length." bp\n\n";
    
    print OUTPUT_LOG "\nLong read simulation\n";
    print OUTPUT_LOG "---------------------\n";
    print OUTPUT_LOG "Read count          : ".$NP_read_count."\n";
    print OUTPUT_LOG "Average read length : ".int($NP_average_read_length)." bp\n";
    print OUTPUT_LOG "Shortest read       : ".$NP_min_read_length." bp\n";
    print OUTPUT_LOG "Longest read        : ".$NP_max_read_length." bp\n\n";
}

#print Graphs--------------------------------------------------------------------------------------------------
if ($SV_input ne "")
{
    foreach my $graph_DEL (sort {$a <=> $b} keys %graph_DEL)
    {
        my $graph_DEL2 = $graph_DEL."0";
        my $graph_DEL3 = int($graph_DEL2*1.5);
        print OUTPUT_DEL $graph_DEL3."\t".$graph_DEL{$graph_DEL}."\n";
    }
    foreach my $graph_INS (sort {$a <=> $b} keys %graph_INS)
    {
        my $graph_INS2 = $graph_INS."0";
        my $graph_INS3 = int($graph_INS2*1.5);
        print OUTPUT_INS $graph_INS3."\t".$graph_INS{$graph_INS}."\n";
    }
    my $kb1 = '0';
    my $kb10 = '0';
    my $kb100 = '0';
    my $Mb1 = '0';
    my $Mb1plus = '0';
    
    foreach my $graph_INV2 (sort {$a <=> $b} keys %graph_INV2)
    {
        if ($graph_INV2 > 0 && $graph_INV2 < 1000)
        {
            $kb1 += $graph_INV2{$graph_INV2};
        }
        if ($graph_INV2 > 1000 && $graph_INV2 < 10000)
        {
            $kb10 += $graph_INV2{$graph_INV2};
        }
        if ($graph_INV2 > 10000 && $graph_INV2 < 100000)
        {
            $kb100 += $graph_INV2{$graph_INV2};
        }
        if ($graph_INV2 > 100000 && $graph_INV2 < 1000000)
        {
            $Mb1 += $graph_INV2{$graph_INV2};
        }
        if ($graph_INV2 > 1000000)
        {
            $Mb1plus += $graph_INV2{$graph_INV2};
        }
    }
    print OUTPUT_INV2 "1\t".$kb1."\n";
    print OUTPUT_INV2 "2\t".$kb10."\n";
    print OUTPUT_INV2 "3\t".$kb100."\n";
    print OUTPUT_INV2 "4\t".$Mb1."\n";
    print OUTPUT_INV2 "5\t".$Mb1plus."\n";
}

close OUTPUT_NP;
close OUTPUT_VCF;
close OUTPUT_REF;
close OUTPUT_HAP1;
close OUTPUT_HAP2;
close OUTPUT_DEL;
close OUTPUT_INS;
close OUTPUT_INV2;
close INPUT_VCF;
close OUTPUT_LOG;
