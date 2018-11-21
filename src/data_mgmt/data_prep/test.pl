#!/usr/local/bin/perl

use strict;
use warnings;
use POSIX;
use File::Basename;
use FindBin;

my $relpath = $FindBin::Bin;
print "This is FindBin: $relpath \n";
print "This is 0: $0 \n";
print "Can we get lower?\n";
my $lowerdir = dirname(dirname(dirname($relpath)));
print "Yes! $lowerdir\n";
