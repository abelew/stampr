#!/usr/bin/env perl
use warnings;

use File::Basename;
use FileHandle;
use IO::Handle;
use Modern::Perl;
use autodie qw":all";
use diagnostics;
use warnings qw"all";
use Bio::Seq;
use Bio::SeqIO;
use Cwd;
use FileHandle;
use File::Basename qw"dirname basename";
use File::Path qw"make_path remove_tree";
use File::Which qw"which";

use strict;
use 5.008_005;
our $VERSION = '0.01';

=head1 NAME

count_idx.pl - A quick and dirty tag counter.

=head1 SYNOPSIS

count_idx.pl <inputfile.fastq.whatever>

=cut

## Handle arbitrary compression of the fastq input.
## Normally I would use Getopt::Long and make this pretty, but no.
my $inputted = FileHandle->new("less $ARGV[0] |");
my $hist = FileHandle->new(">$ARGV[1]");
my $in = Bio::SeqIO->new(-fh => $inputted, -format => 'Fastq');

=item C<BODY>

This body of this script does the following:

1.  Creates the hash 'indices'.

2.  iterates over every fastq tag, pulls the first 9 nucleotides (this also
should be a parameter), and checks indices if it already exists.  If so,
increment that number, if not, set it to 1.

3.  Write the indices as a 2 column tsv file.

=cut
my $indices = {};
READS: while (my $in_seq = $in->next_dataset()) {
    my $id = $in_seq->{'-descriptor'};
    my $sequence = $in_seq->{'-seq'};
    my $qual = $in_seq->{'-raw_quality'};
    ## Reverse it so that when we pull the first 9 nucleotides, we get the ones
    ## closest to the adapter.
    $sequence = substr(reverse($sequence), 0, 9);
    ## $sequence = substr($sequence, 0, 9);
    if (defined($indices->{$sequence})) {
        $indices->{$sequence}++;
    } else {
        $indices->{$sequence} = 1;
    }
}

foreach my $k (keys %{$indices}) {
    print $hist qq"${k}\t$indices->{$k}\n";
}
close($hist);
close($inputted);
