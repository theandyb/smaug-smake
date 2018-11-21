#!/usr/local/bin/perl

use strict;
use warnings;
use POSIX;
use File::Basename;
use File::Path qw(make_path);
use FindBin;
use feature 'say';

my @macs = ("common", "singletons");
use lib "$FindBin::Bin/../lib";
use SmaugFunctions qw(forkExecWait getRef);

my $relpath = $FindBin::Bin; # Location of this script
my $analysisdir = dirname(dirname(dirname($relpath))); # root of the smaug-smake project
my $vcftoolsdir = "$analysisdir/vcftools";

my $vcffile = $ARGV[0];
if (not defined $vcffile) {
    die "Need to provide a vcffile input";
}

for my $mac (@macs){
    my $filename = fileparse($vcffile);
    my @parts = split(/\.vcf.gz/, $filename);
    my $basename = $parts[0];
    my @nameparts = split(/\./, $basename);
    my $i = $nameparts[0];
    $i =~ s/chr//g;

    my $newvcf = "$analysisdir/vcfs/$basename.ma.aa.$mac.vcf.gz";
    my $ancestral = "$analysisdir/reference_data/human_ancestor_GRCh37_e59/human_ancestor_$i.fa.gz";
    my $fasta = "$analysisdir/reference_data/human_g1k_v37/chr$i.fasta.gz";

    # Extract singletons with filter PASS
    my $maparse = "perl $relpath/ma_parse.pl --i $vcffile";

    # Fill ancestral allele to AA field
    my $aaparse = "perl $vcftoolsdir/perl/fill-aa -a $ancestral";

    # Add Motif and Category info fields
    my $infoparse = "perl $relpath/fill_motif.pl -a $fasta";

    # pipe commands and execute
    my $pipe;
    if ($mac eq "singletons"){
        $pipe = "$maparse | $aaparse | $infoparse | bgzip -c > $newvcf";
    } elsif($mac eq "common"){
        my $filter = "bcftools view -i 'AC>=10' -f PASS";
        $pipe = "$filter $vcffile | $infoparse | bgzip -c > $newvcf";
    }
    print "$pipe\n";
    # forkExecWait($pipe);
}
