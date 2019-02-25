
#Main author of this script: Dr Marc-Olivier Duceppe: marc-olivier.duceppe@canada.ca
#https://github.com/duceppemo
#Co-author Dr Emilie Tremblay: emilie.tremblay@canada.ca


#!/bin/perl

use strict;
use warnings;
use diagnostics;
use Data::Dumper;  # debug
use List::MoreUtils qw(indexes);  # to get indexes of a match


###########
#         #
#   I/O   #
#         #
###########


# Must respect the order for the arguments
my $query_in = $ARGV[0];  # text file with word to query. One per line
my $otu_in = $ARGV[1];  # OTU table from QIIME
my $meta_in = $ARGV[2];  # Metadata for/from QIIME
my $fasta_in = $ARGV[3];  # new_refseqs.fna from QIIME
my $table_taxon_out = $ARGV[4];  # Output table by taxon
my $table_sample_out = $ARGV[5];  # Output table by sample
my $fasta_out_folder = $ARGV[6];

# test if all required arguments are present
if (scalar(@ARGV) != 7)  # "scalar" give the length of an array
{
    print "Usage: perl metaResultExtractor.pl <query.txt> <otu_in.tsv> <meta_in.tsv> <new_refseqs.fna> <table_taxon_out.tsv> <table_sample_out.tsv> <fasta_out_folder>\n";
    exit;  # Stop script execution
}


########################
#                      #
#   Parse query file   #
#                      #
########################


#Open query file
open (my $query_in_fh, "<", $query_in) or die "Cannot open $query_in: $!";  # Create file handle

# Create array (list) to store words to search (query)
my @queries;

# Parse query file
while (my $line = <$query_in_fh>)  # Read line by line
{
    chomp($line); #remove carriage return at end of line
    next if $line eq ""; #skip empty lines

    push( @queries, $line);  # add "$line" value to the end of "queries" array
}

# Close query file handle
close ($query_in_fh);

# Debug
# print ( Dumper (\@queries));  # print the content of "@queries" on screen in a nice readable way


###########################
#                         #
#   Parse metadata file   #
#                         #
###########################


open (my $meta_in_fh, "<", $meta_in) or die "Cannot open $meta_in: $!";

#Block comment (multiline) to show the structure of the file we're going to parse
=pod
row.names   BarcodeSequence LinkerPrimerSequence    InputFileName   TrapType    Lure    CollectionDate  Province    City    Longitude   Latitude    Description
EM13S01-ITS2    CTAAGGTAAC  GATGCTGCGTTCTTCATCGATGC EM13S01-ITS2_BC01-Run1.fasta    Insect  alpha_pinene    2016-07-10  British-Colombia    Vernon  50.281538   -119.275796 EM13S01
=cut

my %meta;  # Hash to parse the metadata
my $dummy = <$meta_in_fh>;  # Burn header line
while (my $line = <$meta_in_fh>)  # Read line by line
{
    chomp($line); #remove carriage return at end of line
    next if $line eq ""; #skip empty lines

    # Split the line at the TAB characters ("\t") and assign only the ones that we want to variable
    my ($sample, $traptype, $lure, $date, $province, $city) = (split(/\t/, $line))[0, 4..8];

    # Store the desired information from the metadata file into a hash. Will be used to produce the outputs
    $meta{$sample}{'trap'} = $traptype;
    $meta{$sample}{'lure'} = $lure;
    $meta{$sample}{'date'} = $date;
    $meta{$sample}{'prov'} = $province;
    $meta{$sample}{'city'} = $city;
}

# Close metadata file handle
close ($meta_in_fh);


########################
#                      #
#   Parse fasta file   #
#                      #
########################


open (my $fasta_in_fh, "<", $fasta_in) or die "Cannot open $fasta_in: $!";

my (%fastas, $id);  # data structure to store data

while (my $line = <$fasta_in_fh>)  # Read line by line
{
    chomp($line); #remove carriage return at end of line
    next if $line eq "";  # skip empty lines

    if ($line =~ /^>/)  # If the header line
    {
        $line =~ s/>//;  # remove the ">" to make sure the ID matches the ones from the other input files
        $id = (split(/ /, $line))[0];  # only keep the part of the header before the first space
    }
    else  # it's a sequence line
    {
        push ( @{ $fastas{$id} }, $line);  # we use a push because sequence can be extended over more than one line
    }
}

close ($fasta_in_fh);


#######################
#                     #
#   Parse OTU table   #
#                     #
#######################


open (my $otu_in_fh, "<", $otu_in) or die "Cannot open $otu_in: $!";

# block comment
=pod
row.names   EM13S01-ITS2    EM13S02-ITS2 taxonomy
SH210217.07FU_AM231334_reps 0   0   k__Fungi; p__unidentified; c__unidentified; o__unidentified; f__unidentified; g__unidentified; s__Fungi sp
=cut

# Create data structures (hash) to store OTU table (and other files) data
my (%taxons, %samples);

my $first_line = <$otu_in_fh>;  # Read first line (header)
chomp($first_line);  # Remove carriage return

#Process header
my @header = split(/\t/, $first_line);  # Split header at TAB character and store in array

while (my $line = <$otu_in_fh>)  # Read line by line
{
    chomp($line); #remove carriage return at end of line
    next if $line eq ""; #skip empty lines

    my @fields = split(/\t/, $line);  # put each field of the line in an array
    my $seq = $fields[0];  # sequence header information (same as in rep_set)
    my $taxo = $fields[- 1];  # Last field

    # Search for words in query fine in the taxonomy column
    foreach my $word (@queries)
    {
        if (index($taxo, $word) != - 1)  # if word is found in taxonomy
        {
            #find the index of the samples with positive counts
            my @index = indexes { $_ > 0 } (@fields[1..($#fields - 1)]);  # the line iwthout first and last fields

            foreach my $i (@index)  #for all the indexes found
            {
                # Get the sample name
                my $sample = $header[$i + 1];  # +1 because the index omits the first column (non numerical)

                # Get the count value
                my $cnt = $fields[$i + 1];

                # Store desired information in hash

                # By taxon
                $taxons{$word}{$sample}{'count'} += $cnt;  # Count total reads
                push ( @{ $taxons{$word}{$sample}{'seq'} }, $fastas{$seq});  # push because can have more than one seq

                # Add metadata info
                $taxons{$word}{$sample}{'trap'} = $meta{$sample}{'trap'};
                $taxons{$word}{$sample}{'lure'} = $meta{$sample}{'lure'};
                $taxons{$word}{$sample}{'date'} = $meta{$sample}{'date'};
                $taxons{$word}{$sample}{'prov'} = $meta{$sample}{'prov'};
                $taxons{$word}{$sample}{'city'} = $meta{$sample}{'city'};

                # By sample
                $samples{$sample}{$word}{'count'} += $cnt;  # Count total reads
                push ( @{ $samples{$sample}{$word}{'seq'} }, $fastas{$seq});  # push because can have more than one seq

                # Add metadata info
                $samples{$sample}{$word}{'trap'} = $meta{$sample}{'trap'};
                $samples{$sample}{$word}{'lure'} = $meta{$sample}{'lure'};
                $samples{$sample}{$word}{'date'} = $meta{$sample}{'date'};
                $samples{$sample}{$word}{'prov'} = $meta{$sample}{'prov'};
                $samples{$sample}{$word}{'city'} = $meta{$sample}{'city'};
            }

            last;  # Don't look for other words if found because can't happen
        }
    }
}

# Close otu table handle
close ($otu_in_fh);


#############################
#                           #
#   Output table by taxon   #
#                           #
#############################


open(my $table_taxon_out_fh, '>', $table_taxon_out) or die "Cannot write to $table_taxon_out: $!";

#print header in tab-separated format
print { $table_taxon_out_fh } ("Taxon\tSample\tCity\tProvince\tCollection_Date\tTrap_type\tLure\tCount\tSequence\n");


#Loop through the hash to be able to print it's content to a file
foreach my $word (sort(keys(%taxons)))
{
    foreach my $sample (sort keys %{ $taxons{$word} } )
    {
        #populate output file in tab-separated format
        print { $table_taxon_out_fh } ($word
            . "\t" . $sample
            . "\t" . $taxons{$word}{$sample}{'city'}
            . "\t" . $taxons{$word}{$sample}{'prov'}
            . "\t" . $taxons{$word}{$sample}{'date'}
            . "\t" . $taxons{$word}{$sample}{'trap'}
            . "\t" . $taxons{$word}{$sample}{'lure'}
            . "\t" . $taxons{$word}{$sample}{'count'});

        # Print all the sequences that were associated to a specific "word" (or taxon)
        foreach my $s (@{ $taxons{$word}{$sample}{'seq'} })
        {
            foreach my $x (@{ $s })
            {
                print { $table_taxon_out_fh } ("\t" . $x);
            }
        }

        print { $table_taxon_out_fh } ("\n");  # Add a carriage return
    }
}


##############################
#                            #
#   Output table by sample   #
#                            #
##############################


#output table by sample
open(my $table_sample_out_fh, '>', $table_sample_out) or die "Cannot write to $table_sample_out: $!";

#print header
print { $table_sample_out_fh } ("Sample\tTaxon\tCity\tProvince\tCollection_Date\tTrap_type\tLure\tCount\tSequence\n");

foreach my $sample (sort(keys(%samples)))
{
    # Create one output file for each sample
    my $fasta_taxon_out = "$fasta_out_folder" . "/" . $sample . ".fasta";
    open(my $fasta_taxon_out_fh, '>', $fasta_taxon_out) or die "Cannot write to $fasta_taxon_out: $!";

    foreach my $taxon (sort keys %{ $samples{$sample} } )
    {
        my $city = $samples{$sample}{$taxon}{'city'};
        my $prov = $samples{$sample}{$taxon}{'prov'};
        my $date = $samples{$sample}{$taxon}{'date'};
        my $trap = $samples{$sample}{$taxon}{'trap'};
        my $lure = $samples{$sample}{$taxon}{'lure'};
        my $count = $samples{$sample}{$taxon}{'count'};

        print { $table_sample_out_fh } ($sample
            . "\t" . $taxon
            . "\t" . $city
            . "\t" . $prov
            . "\t" . $date
            . "\t" . $trap
            . "\t" . $lure
            . "\t" . $count);

        foreach my $s (@{ $samples{$sample}{$taxon}{'seq'} })
        {
            my $seqCounter = 0;  # to track how many sequences liked to a taxon
            foreach my $x (@{ $s })
            {
                print { $table_sample_out_fh } ("\t" . $x);

                # print fasta out
                $seqCounter += 1;
                my $fastaHeader = ">" . join('_', $sample, $taxon, $city, $prov,
                    $date, $trap, $lure, $count, $seqCounter) . "\n";
                print { $fasta_taxon_out_fh } ($fastaHeader . $x . "\n");
            }
        }

        print { $table_sample_out_fh } ("\n");  # End of line
    }
}
