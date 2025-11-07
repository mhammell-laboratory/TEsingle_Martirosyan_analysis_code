#!/bin/env perl

use strict;
use warnings;

die "Usage: annotate_barcode.pl [library annotation] [zipped barcode] > [csv output]\n" if scalar @ARGV < 2;

my $annotFile = shift @ARGV;
open my $ifh, "<", $annotFile or die $!;

my %annot;
while(my $line = <$ifh>){
  chomp $line;
  my @tmp = split "\t", $line;
  $annot{$tmp[0]} = join(",",@tmp[1..$#tmp]);
}
close $ifh or die $!;

my $barcodeFile = shift @ARGV;
open $ifh, "-|", "gunzip", "-cf", $barcodeFile or die $!;

while(my $line = <$ifh>){
  chomp $line;
  my ($bc,$lib) = split "-", $line;
  print $line, ",", $annot{$lib}, "\n";
}
close $ifh or die $!;

