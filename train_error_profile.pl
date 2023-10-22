#!/usr/bin/env perl
use strict;
use Getopt::Long;
use MCE::Child;
use MCE::Channel;

print "\n\n-----------------------------------------------";
print "\nSim-it - error profile training\n";
print "Version 1.0\n";
print "Author: Nicolas Dierckxsens, (c) 2020\n";
print "-----------------------------------------------\n\n";

my $ref_file = "";
my $reads_file = "";

GetOptions (
            "ref=s" => \$ref_file,
            "reads=s" => \$reads_file,
            ) or die "Incorrect usage!\n";

my $USAGE = "\nUsage: perl train_error_profile.pl -ref reference.fasta -reads ONT_or_PacBio.fasta\n\n";

print $USAGE and exit if !$ref_file or !$reads_file;

my $output_file = "error_profile_figure.txt";
open(OUTPUT, ">" .$output_file) or die "Can't open file $output_file, $!\n";

my $FILE;
open($FILE, $reads_file) or die "\n\nCan't open the reads file $reads_file, $!\n";

my $error_profile = "error_profile_1.txt";
open(ERROR_PROFILE, ">" .$error_profile) or die "\n\nCan't open error profile file $error_profile, $!\n";

my $chnl;  
$chnl = MCE::Channel->new( impl => 'Simple' );

#spin up worker early before creating big hash---------
mce_child
{
    local $SIG{__WARN__} = sub {};
    while ( my ($cmd, @args) = $chnl->recv ) {
        local ($?, $!);
        system($cmd, @args);
        $chnl->send2($?, $!);
    }
};

sub syscmd {
    my $cmd = shift;
    return unless $cmd;

    $chnl->send($cmd, @_);
    my ($status, $errmsg) = $chnl->recv2;
    
    if ($status == -1) {
        print "SYSTEM: failed to execute ($cmd): $errmsg\n";
    }
    elsif ($status & 127) {
        printf "SYSTEM: $cmd died with signal %s, %s coredump\n",
            ($status & 127), ($status & 128) ? 'with' : 'without';
    }
    else {
        #printf "SYSTEM: $cmd exited with status %d\n", $status >> 8;
    }
}   
  
  
my $total_SV = "";
my $count_deletions_final = "";
my $count_insertions_final = "";
my $count_inversions_final = "";
my $total_matches_total = '0';
my $total_mismatches_total = '0';
my $total_del_read = '0';
my %del_length;
my $total_ins_read = '0';
my %ins_length;
my $length_read_total = '0';
my $length_ref = '0';
my $A_to_C = '0';
my $A_to_T = '0';
my $A_to_G = '0';
my $C_to_A = '0';
my $C_to_T = '0';
my $C_to_G = '0';
my $T_to_A = '0';
my $T_to_C = '0';
my $T_to_G = '0';
my $G_to_A = '0';
my $G_to_C = '0';
my $G_to_T = '0';
my $A_to_A = '0';
my $C_to_C = '0';
my $T_to_T = '0';
my $G_to_G = '0';
my $reference_length = '7000';

my %mismatches;
my %mismatches2;
my %mismatches0;
my %indel0;
my $mismatch = "";
my $mismatch2 = "";
my %ins;
my %ins2;
my $ins = "";
my $ins2 = "";
my %del;
my $del = "";
my $total_del_pos = '0';
my $total_ins_pos = '0';

my $count_reads = '0';
my $line_content = "";
my $read_tmp = "";
my $blast_now = "";
BLAST: while (my $line = <$FILE>)
{
    chomp($line);
    my $first_nuc = substr $line, 0, 1;
    if (($first_nuc eq ">" || $first_nuc eq "@") && $line_content eq "")
    {
        $line_content = "id";
    }
    elsif ($first_nuc eq ">" || $first_nuc eq "@")
    {
        $blast_now = "yes";
    }
    elsif ($first_nuc eq "+")
    {
        next;
    }
    elsif ($first_nuc eq "A" || $first_nuc eq "C" || $first_nuc eq "T" || $first_nuc eq "G" || $first_nuc eq "N" || $first_nuc eq "a" || $first_nuc eq "c" || $first_nuc eq "t" || $first_nuc eq "g" || $first_nuc eq "n")
    {
        $read_tmp .= $line;
    }
    else
    {
        next;
    }
    my $total_matches_read = '0';
    my $total_mismatches_read = '0';
    my $length_ref_read = '0';
    my $total_ins_read_single = '0';
    
    if ($blast_now eq "yes" || eof)
    {
        $count_reads++;
        #$length_read = length($line);
        my $output_file1  = "sequence_tmp.fasta";
        open(OUTPUT1, ">" .$output_file1) or die "\nCan't open file $output_file1, $!\n";
        print OUTPUT1 ">".$count_reads."\n";
        print OUTPUT1 $read_tmp."\n";
        my $command = "blastn -query ".$output_file1." -subject ".$ref_file." -out blast_tmp.txt -outfmt 4 -strand plus -evalue 1e-50 -word_size 60 -culling_limit 1";                
        syscmd($command);
        close OUTPUT1;
        
        my $input_fle2  = "blast_tmp.txt";
        open(INPUT2, $input_fle2) or die "\n\nCan't open blast file $input_fle2, $!\n";
        
        my $reference0 = "";
        my $multiple_hits0 = "";
        my $seq_ref = "";
        my $seq_ref_indel = "";
        $read_tmp = "";
        $blast_now = "";

        while (my $line2 = <INPUT2>)
        {
            chomp($line2);
            my $part2 = substr $line2, 0, 7;
                 
            if ($part2 eq "Query_1")
            {
                $multiple_hits0 = "";
            }
            if ($part2 eq "Subject" && $multiple_hits0 eq "")
            {
                my @split = split /\s+/, $line2;
                $reference0 = $split[2];
                $reference0 =~ tr/actgn/ACTGN/;         
                my @reference0 =  split //, $reference0;
                foreach my $ref0 (@reference0)
                {
                    if ($ref0 eq "A" || $ref0 eq "C" || $ref0 eq "T" || $ref0 eq "G")
                    {
                        $seq_ref .= $ref0;
                        $seq_ref_indel .= $ref0;
                        if (length($seq_ref) eq '3')
                        {
                            my $count_tmp = '0';
                            if (exists($mismatches0{$seq_ref}))
                            {
                                $count_tmp = $mismatches0{$seq_ref};
                            }
                            $count_tmp++;
                            $mismatches0{$seq_ref} = $count_tmp;
                            substr $seq_ref, 0, 1, "";
                        }
                        if (length($seq_ref_indel) eq '2')
                        {
                            my $count_tmp = '0';
                            if (exists($indel0{$seq_ref_indel}))
                            {
                                $count_tmp = $indel0{$seq_ref_indel};
                            }
                            $count_tmp++;
                            $indel0{$seq_ref_indel} = $count_tmp;
                            substr $seq_ref_indel, 0, 1, "";
                        }
                    }
                }
                $multiple_hits0 = "yes";
            }
        }
        close INPUT2;
 #---------------------------------       
        
        open(INPUT2, $input_fle2) or die "\n\nCan't open blast file $input_fle2, $!\n";
        my $multiple_hits = "";
        my $read = "";
        my $reference = "";
        my @read;
        my @reference;
        undef @read;
        undef @reference;
        
        my $prev_nuc = "";
            
        while (my $line2 = <INPUT2>)
        {
            chomp($line2);
            my $part2 = substr $line2, 0, 7;
            
            if ($part2 eq "Query_1")
            {
                my @split = split /\s+/, $line2;
                $read = $split[2];
                $read =~ tr/actgn/ACTGN/;
                @read =  split //, $read;
                $multiple_hits = "";
            }
            if ($part2 eq "Subject" && $multiple_hits eq "")
            {
                my @split = split /\s+/, $line2;            
                $reference = $split[2];
                $reference =~ tr/actgn/ACTGN/;
                @reference =  split //, $reference;
                $multiple_hits = "yes";
                
                my $u = '0';
                my $del_length = '0';
                my $no_del = "yes";
                my $ins_length = '0';
                my $no_ins = "yes";
                
                while ($u < length($read))
                {
                    if ($mismatch ne "" && $reference[$u] =~ m/A|C|T|G/)
                    { 
                        $mismatch2 = $mismatch;
                        chop($mismatch2);
                        $mismatch .= $reference[$u];
                        $mismatch2 .= $reference[$u];
                        
                        my $count_tmp = '0';
                        if (exists($mismatches{$mismatch}))
                        {
                            $count_tmp = $mismatches{$mismatch};
                        }
                        $count_tmp++;
                        $mismatches{$mismatch} = $count_tmp;
                        $mismatch = "";
                        
                        my $count_tmp2 = '0';
                        if (exists($mismatches2{$mismatch2}))
                        {
                            $count_tmp2 = $mismatches2{$mismatch2};
                        }
                        $count_tmp2++;
                        $mismatches2{$mismatch2} = $count_tmp2;
                        $mismatch2 = "";     
                    }
                    if ($del ne "" && $reference[$u] =~ m/A|C|T|G/)
                    {
                        $del .= $reference[$u];
                        my $count_tmp = '0';
                        if (exists($del{$del}))
                        {
                            $count_tmp = $del{$del};
                        }
                        $count_tmp++;
                        $del{$del} = $count_tmp;
                        $del = "";
                    }
                    if ($ins ne "" && $reference[$u] =~ m/A|C|T|G/)
                    {
                        $ins2 = $ins;
                        chop($ins2);
                        $ins .= $reference[$u];
                        $ins2 .= $reference[$u];
                        
                        my $count_tmp = '0';
                        if (exists($ins{$ins}))
                        {
                            $count_tmp = $ins{$ins};
                        }
                        $count_tmp++;
                        $ins{$ins} = $count_tmp;
                        $ins = "";
                        
                        my $count_tmp2 = '0';
                        if (exists($ins2{$ins2}))
                        {
                            $count_tmp2 = $ins2{$ins2};
                        }
                        $count_tmp2++;
                        $ins2{$ins2} = $count_tmp2;
                        $ins2 = "";
                    }
                    
                    if ($read[$u] eq $reference[$u] && $read[$u] ne "-")
                    {
                        $total_matches_read++;
                        $total_matches_total++;
                        $length_ref_read++;
                        $length_read_total++;
                        $length_ref++;
                        $no_del = "yes";
                        $no_ins = "yes";
                        
                        if ($reference[$u] eq "A" && $read[$u] eq "A")
                        {
                            $A_to_A++;
                        }
                        elsif ($reference[$u] eq "C" && $read[$u] eq "C")
                        {
                            $C_to_C++;
                        }
                        elsif ($reference[$u] eq "T" && $read[$u] eq "T")
                        {
                            $T_to_T++;
                        }
                        elsif ($reference[$u] eq "G" && $read[$u] eq "G")
                        {
                            $G_to_G++;
                        }
                        $prev_nuc = $reference[$u];
                    }
                    elsif (($read[$u] eq "A" || $read[$u] eq "C" || $read[$u] eq "T" || $read[$u] eq "G" || $read[$u] eq "N") && ($reference[$u] eq "A" || $reference[$u] eq "C" || $reference[$u] eq "T" || $reference[$u] eq "G" || $reference[$u] eq "N"))
                    {
                        $total_mismatches_read++;
                        $total_mismatches_total++;
                        
                        $mismatch = $prev_nuc.$reference[$u].$read[$u];
                        
                        if ($reference[$u] eq "A" && $read[$u] eq "C")
                        {
                            $A_to_C++;
                        }
                        elsif ($reference[$u] eq "A" && $read[$u] eq "T")
                        {
                            $A_to_T++;
                        }
                        elsif ($reference[$u] eq "A" && $read[$u] eq "G")
                        {
                            $A_to_G++;
                        }
                        elsif ($reference[$u] eq "C" && $read[$u] eq "A")
                        {
                            $C_to_A++;
                        }
                        elsif ($reference[$u] eq "C" && $read[$u] eq "T")
                        {
                            $C_to_T++;
                        }
                        elsif ($reference[$u] eq "C" && $read[$u] eq "G")
                        {
                            $C_to_G++;
                        }
                        elsif ($reference[$u] eq "T" && $read[$u] eq "A")
                        {
                            $T_to_A++;
                        }
                        elsif ($reference[$u] eq "T" && $read[$u] eq "C")
                        {
                            $T_to_C++;
                        }
                        elsif ($reference[$u] eq "T" && $read[$u] eq "G")
                        {
                            $T_to_G++;
                        }
                        elsif ($reference[$u] eq "G" && $read[$u] eq "A")
                        {
                            $G_to_A++;
                        }
                        elsif ($reference[$u] eq "G" && $read[$u] eq "C")
                        {
                            $G_to_C++;
                        }
                        elsif ($reference[$u] eq "G" && $read[$u] eq "T")
                        {
                            $G_to_T++;
                        }
                        $length_ref_read++;
                        $length_read_total++;
                        $length_ref++;
                        $no_del = "yes";
                        $no_ins = "yes";
                        $prev_nuc = $reference[$u];
                    }
                    elsif ($read[$u] eq "-" && $reference[$u] ne "-")
                    {
                        $total_del_read++;
                        $del_length++;
                        $length_ref++;
                        $length_ref_read++;
                        $no_del = "";
                        $no_ins = "yes";
                        
                        $del = $prev_nuc.$reference[$u];
                        $prev_nuc = $reference[$u];
                    }
                    elsif ($read[$u] ne "-" && $reference[$u] eq "-")
                    {
                        $total_ins_read++;
                        $total_ins_read_single++;
                        $length_read_total++;
                        $ins_length++;
                        $no_del = "yes";
                        $no_ins = "";
                        
                        $ins = $prev_nuc.$read[$u];
                    }
                    if ($no_del eq "yes" && $del_length > 0)
                    {
                         if (exists($del_length{$del_length}))
                         {
                            my $count = $del_length{$del_length}+1;
                            $del_length{$del_length} = $count;
                         }
                         else
                         {
                            $del_length{$del_length} = '1';
                         }
                         $total_del_pos++;
                         $del_length = '0';
                    }
                    if ($no_ins eq "yes" && $ins_length > 0)
                    {
                         if (exists($ins_length{$ins_length}))
                         {
                            my $count = $ins_length{$ins_length}+1;
                            $ins_length{$ins_length} = $count;
                         }
                         else
                         {
                            $ins_length{$ins_length} = '1';
                         }
                         $total_ins_pos++;
                         $ins_length = '0';
                    }
                    $u++;                  
                }
                $read = "";
                $reference = "";
                undef @read;
                undef @reference;
            }
        }
        close INPUT2;
    }
    else
    {
        next BLAST;
    }
    my $accuracy_read = '0';
    if ($length_ref_read > 0)
    {
        my $accuracy_read = sprintf("%.3g",($total_matches_read-$total_ins_read_single)/$length_ref_read)*100;
        if ($accuracy_read >= 180 && $accuracy_read < 82)
        {
            print $count_reads."\n";
            print $total_matches_read." MAT\n";
            print $total_ins_read_single." INS\n";
            print $length_ref_read." REF\n";
            print "Accuracy                        : ".$accuracy_read."%\n";
        }
    }
    my $total_matches = '0';
    my $t = '0';
    my $t2 = '0';
}
my $coverage = $length_read_total/$reference_length;
my $accuracy_total = sprintf("%.3g",($total_matches_total-$total_ins_read)/$length_ref)*100;
my $del_read = sprintf("%.3g",$total_del_read/$length_read_total)*100;
my $ins_read = sprintf("%.3g",$total_ins_read/$length_read_total)*100;
my $mismatches_percentage = sprintf("%.3g", $total_mismatches_total/$length_ref)*100;
my $A_mm =  sprintf("%.3g",($A_to_C+$A_to_T+$A_to_G)/($A_to_A+$A_to_C+$A_to_T+$A_to_G))*100;
my $C_mm =  sprintf("%.3g",($C_to_A+$C_to_T+$C_to_G)/($C_to_A+$C_to_C+$C_to_T+$C_to_G))*100;
my $T_mm =  sprintf("%.3g",($T_to_C+$T_to_A+$T_to_G)/($T_to_A+$T_to_C+$T_to_T+$T_to_G))*100;
my $G_mm =  sprintf("%.3g",($G_to_C+$G_to_T+$G_to_A)/($G_to_A+$G_to_C+$G_to_T+$G_to_G))*100;

print "Mismatches                      : ".$mismatches_percentage." %\n";
print "A                               : ".$A_mm."%\n";
print "A to A                          : ".$A_to_A."\n";
print "A to C                          : ".$A_to_C."\n";
print "A to T                          : ".$A_to_T."\n";
print "A to G                          : ".$A_to_G."\n";
print "C                               : ".$C_mm."%\n";
print "C to C                          : ".$C_to_C."\n";
print "C to A                          : ".$C_to_A."\n";
print "C to T                          : ".$C_to_T."\n";
print "C to G                          : ".$C_to_G."\n";
print "T                               : ".$T_mm."%\n";
print "T to T                          : ".$T_to_T."\n";
print "T to A                          : ".$T_to_A."\n";
print "T to C                          : ".$T_to_C."\n";
print "T to G                          : ".$T_to_G."\n";
print "G                               : ".$G_mm."%\n";
print "G to G                          : ".$G_to_G."\n";
print "G to A                          : ".$G_to_A."\n";
print "G to C                          : ".$G_to_C."\n";
print "G to T                          : ".$G_to_T."\n\n";
print "Deletions                       : ".$del_read."%\n";

print ERROR_PROFILE $mismatches_percentage."\n";
print ERROR_PROFILE $del_read."\n";
print ERROR_PROFILE $ins_read."\n";
print ERROR_PROFILE "DEL_LENGTH\n";
my $fraction_del_length = '0';
foreach my $del_length_tmp (sort {$a <=> $b} keys %del_length)
{
    print $del_length_tmp." bp                            : ".$del_length{$del_length_tmp}."\n";
    my $fraction_tmp = $del_length{$del_length_tmp}/($total_del_pos-$fraction_del_length);
    print ERROR_PROFILE $del_length_tmp.":".$fraction_tmp."\n";
    $fraction_del_length += $del_length{$del_length_tmp};
}

print ERROR_PROFILE "INS_LENGTH\n";
my $fraction_ins_length = '0';
print "\nInsertions                      : ".$ins_read."%\n";
foreach my $ins_length_tmp (sort {$a <=> $b} keys %ins_length)
{
    print $ins_length_tmp." bp                            : ".$ins_length{$ins_length_tmp}."\n";
    my $fraction_tmp = $ins_length{$ins_length_tmp}/($total_ins_pos-$fraction_ins_length);
    print ERROR_PROFILE $ins_length_tmp.":".$fraction_tmp."\n";
    $fraction_ins_length += $ins_length{$ins_length_tmp};
}

print "\nCoverage                        : ".$coverage."\n";
print "\nAccuracy                        : ".$accuracy_total."%\n\n";
print "\n\nMISMATCH:  BEFORE/REFERENCE/READ/AFTER\n";

print ERROR_PROFILE "MISMATCHES\n";
print OUTPUT "MISMATCHES\n\n";
my $count_tmp = '0';
my $count_tmp2 = '0';

foreach my $mismatch_tmp (sort keys %mismatches)
{
    if (length($mismatch_tmp) eq '4')
    {
        print $mismatch_tmp." : ".$mismatches{$mismatch_tmp}."\n";
        my $mismatch_ref = $mismatch_tmp;
        substr $mismatch_ref, 2, 1, "";
        if (exists($mismatches0{$mismatch_ref}))
        {
            my $per = sprintf("%.3g",$mismatches{$mismatch_tmp}/$mismatches0{$mismatch_ref})*100;
            print $mismatch_tmp." : ".$per." %\n";
            print ERROR_PROFILE $mismatch_tmp.":".$per."\n";
            
            print OUTPUT $per."\t";
            $count_tmp++;
            $count_tmp2++;
            if ($count_tmp2 eq '12')
            {
                print OUTPUT "\n";
                $count_tmp = '0';
                $count_tmp2 = '0';
            }
            elsif ($count_tmp eq '4')
            {
                print OUTPUT "\t";
                $count_tmp = '0';
            }
        }
    }
}

print ERROR_PROFILE "MISMATCHES2\n";
foreach my $mismatch_tmp2 (sort keys %mismatches2)
{
    if (length($mismatch_tmp2) eq '3')
    {
        if (exists($mismatches0{$mismatch_tmp2}))
        {
            my $per = sprintf("%.3g",$mismatches2{$mismatch_tmp2}/$mismatches0{$mismatch_tmp2})*100;
            print ERROR_PROFILE $mismatch_tmp2.":".$per."\n";
        }
    }
}

$count_tmp = '0';
$count_tmp2 = '0';
print ERROR_PROFILE "DELETIONS\n";
print OUTPUT "\n\nDELETIONS\n\n";
print "\n\nDEL:  BEFORE/DELETION/AFTER\n";
foreach my $del_tmp (sort keys %del)
{
    if (length($del_tmp) eq '3')
    {
        print $del_tmp." : ".$del{$del_tmp}."\n";
        if (exists($mismatches0{$del_tmp}))
        {
            my $per = sprintf("%.3g",$del{$del_tmp}/$mismatches0{$del_tmp})*100;
            print $del_tmp." : ".$per." %\n";
            print ERROR_PROFILE $del_tmp.":".$per."\n";
            
            print OUTPUT $per."\t";
            $count_tmp++;
            $count_tmp2++;
            if ($count_tmp2 eq '16')
            {
                print OUTPUT "\n";
                $count_tmp = '0';
                $count_tmp2 = '0';
            }
            elsif ($count_tmp eq '4')
            {
                print OUTPUT "\t";
                $count_tmp = '0';
            }
        }
    }
}

$count_tmp = '0';
$count_tmp2 = '0';
print ERROR_PROFILE "INSERTIONS\n";
print OUTPUT "\n\nINSERTIONS\n\n";
print "\n\nINS:  BEFORE/INSERTION/AFTER\n";
foreach my $ins_tmp (sort keys %ins)
{ 
    if (length($ins_tmp) eq '3')
    {
        print $ins_tmp." : ".$ins{$ins_tmp}."\n";
        my $ins_tmp_ref = $ins_tmp;
        substr $ins_tmp_ref, 1, 1, "";
        if (exists($indel0{$ins_tmp_ref}))
        {
            my $per = sprintf("%.3g",$ins{$ins_tmp}/$indel0{$ins_tmp_ref})*100;
            print $ins_tmp." : ".$per." %\n";
            print ERROR_PROFILE $ins_tmp.":".$per."\n";
            
            print OUTPUT $per."\t";
            $count_tmp++;
            $count_tmp2++;
            if ($count_tmp2 eq '16')
            {
                print OUTPUT "\n";
                $count_tmp = '0';
                $count_tmp2 = '0';
            }
            elsif ($count_tmp eq '4')
            {
                print OUTPUT "\t";
                $count_tmp = '0';
            }
        }
    }
}
print ERROR_PROFILE "INSERTIONS2\n";
foreach my $ins_tmp2 (sort keys %ins2)
{
    if (length($ins_tmp2) eq '2')
    {
        if (exists($indel0{$ins_tmp2}))
        {
            my $per = sprintf("%.3g",$ins2{$ins_tmp2}/$indel0{$ins_tmp2})*100;
            print ERROR_PROFILE $ins_tmp2.":".$per."\n";
        }
    }
}
close $FILE;
close ERROR_PROFILE;
close OUTPUT1;
close OUTPUT;
