# Sim-it: A structural variance and Nanopore/PacBio sequencing reads simulator

## Getting help

Any issues/requests/problems/comments that are not yet addressed on this page can be posted on [Github issues](https://github.com/ndierckx/Sim-it/issues) and I will try to reply the same day.

Or you can contact me directly through the following email address:

nicolasdierckxsens at hotmail dot com 

## Prerequisites

Perl

From version 1.3.2 on, additional Perl modules need to be installed:

<code>cpan Parallel::ForkManager</code>

## Instructions

**!A more complete manual can found under the wiki section!**
https://github.com/ndierckx/Sim-it/wiki


Usage:

<code>perl Sim-it1.3.9.pl -c config_Sim-it.txt -o output/directory/path</code>


----------------------------------------------------------------------------------------------------------


## Configuration file

This is an example of a configuration file.
To make the simulation work, your configuration file has to have the exact same structure.
(Make sure there is always a space after the equals sign and every parameter is captured in one single line)

**1. Example of configuration file:**
<pre>

Project:
-----------------------------
Project name             = Test
Reference sequence       = /path/to/reference.fasta
Replace ambiguous nts(N) = 
Max threads              = 8


Structural variation:
-----------------------------
VCF input                =
Foreign sequences        =

Deletions                = 0
Length (bp)              = 30-150000

Insertions               = 0
Length (bp)              = 30-100000

Tandem duplications      = 0
Length (bp)              = 50-10000
Copies                   = 1-20

Inversions               = 0
Length (bp)              = 150-1300000

Complex substitutions    = 0
Length (bp)              = 30-100000

Inverted duplications    = 0
Length (bp)              = 150-350000

Heterozygosity           = 60%


Long Read simulation:
-----------------------------
Sequencing depth         = 0
Median length            = 25000
Length range             = 500-150000
Accuracy                 = 94%
Error profile            = error_profile_ONT.txt
</pre>

**2. Explanation parameters:**
<pre>

Project:
-----------------------------
Project name             = Choose a name for your project, it will be used for the output files.
Reference sequence       = /path/to/reference.fasta  #it can be a gzipped file
Replace ambiguous nts(N) = If the reference contains regions with ambiguous nucleotides, these can be replaced by 
                           random nucleotides to avoid reads consisting out of Ns (yes/no) (Default: no)
Max threads              = You can choose the mamximum parallel threads that can be used (Default: no)


Structural variation:
-----------------------------
VCF input                = A list of SV positions can be given as input, look at the wiki section for the format of 
                           the VCF input. This option can be combined with random SV inputs
Foreign sequences        = If foreign sequences have to be inserted in the reference, they need to listed in a fasta file 

Deletions                = Deletions count (Default: 0)
Length (bp)              = The range for the length of the Deletons (Default: 30-150000)

Insertions               = Insertions count (Default: 0)
Length (bp)              = The range for the length of the Insertions (Default: 30-100000)

Tandem duplications      = Tandem duplications count (Default: 0)
Length (bp)              = The range for the length of a single duplication (Default: 50-10000)
Copies                   = The range of how many times the duplicated sequences is repeated (Default: 1-20)

Inversions               = Inversions count (Default: 0)
Length (bp)              = The range for the length of the Inversions (Default: 150-1300000)

Complex substitutions    = Complex substitutions count (Default: 0)
Length (bp)              = The range for the length of the Complex substitutions (Default: 30-100000)

Inverted duplications    = Inverted duplications count (Default: 0)
Length (bp)              = The range for the length of the Inverted duplications (Default: 150-350000)

Heterozygosity           = The percentage of heterozygous SVs (Default: 60%)


Read simulation:
-----------------------------
Sequencing depth         = The genome coverage of the long reads (Default: 0). 
                           Besides a fixed value, the path to a sequencing depth profile obtained with samtools can be given.
Median length            = The average length of the reads (Default: 25000)
Length range             = The range for the lengths of the reads, it will give a normal distribution around the 
                           median length (Default: 500-150000)
Accuracy                 = The average accuracy of the reads (Default: 88%)
Error profile            = This is a file that contains the error profile of the sequencing technology, the github page 
                           provides ONT, PacBio RS2, Sequel II & Sequel CCS/HiFi. These error profiles will provide the 
                           ratios, which will kept when adjusting the average accuracy to your liking
</pre>
