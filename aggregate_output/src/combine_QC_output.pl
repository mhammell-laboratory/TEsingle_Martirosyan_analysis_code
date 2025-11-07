#!/bin/env perl

use strict;
use warnings;

die "Usage: perl combine_ED_DF_output.pl [lib order] [preprocessing output ...] > [combined output CSV]\n" if scalar @ARGV < 2;

my %lib;

my $infile = shift @ARGV;
open my $ifh, "<", $infile or die $!;
while(my $line = <$ifh>){
  chomp $line;
  my @tmp = split "\t", $line;
  $lib{$tmp[0]} = $tmp[1];
}
close $ifh or die $!;

print "orig.ident,nCount_RNA,nFeature_RNA,EmptyDropsFDR,pANN,DoubletFinder\n";
foreach $infile (@ARGV){
  open $ifh, "<", $infile or die $!;
  while(my $line = <$ifh>){
    next if $line =~ /^orig/;
    chomp $line;
    my @tmp = split "\t", $line;
    if(! exists $lib{$tmp[1]}){
      warn "Library $tmp[1] not found. Entry skipped\n";
      next;
    }
    if(scalar @tmp < 7){
      push @tmp, "0";
      push @tmp, "Singlet";
    }
    $tmp[0] .= "-" . $lib{$tmp[1]};
    print join(",",@tmp), "\n";
  }
  close $ifh or die $!;
}


