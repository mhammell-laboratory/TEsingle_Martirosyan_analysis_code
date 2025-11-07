#!/bin/env perl

use strict;
use warnings;
use File::Basename;

usage() if scalar @ARGV < 2;

my $outprefix = shift @ARGV;
if (! -d $outprefix){
  mkdir $outprefix, 0755;
}
else{
  print STDERR "Folder $outprefix exists. Files within will be overridden\n";
}

my $outbc = "${outprefix}/barcodes.tsv.gz";
my $outfeat = "${outprefix}/features.tsv.gz";
my $outtmp = "${outprefix}/matrix.tmp";
my $outmtx = "${outprefix}/matrix.mtx.gz";
my $outmdata = "${outprefix}/libOrder.txt";

open my $bfh, "|-", "gzip -c > $outbc" or die $!;
open my $tfh, ">", $outtmp or die $!;
open my $ofh, ">", $outmdata or die $!;

my $file_counter = 0;
my $total_bc = 0;
my $total_count = 0;
my $total_entries = 0;
my $total_feature = 0;
my $infile;
my $ifh;
foreach my $aFile (@ARGV){
  $file_counter ++;
  my $libID = basename($aFile);
  $libID =~ s/\.gz//;
  $libID =~ s/\.annots//;
  print {$ofh} "$libID\t$file_counter\n";
  if($file_counter == 1){
    $infile = $aFile;
    open $ifh, "-|", "gunzip", "-cf", $infile or die $!;
    open my $ffh, "|-", "gzip -c > $outfeat" or die $!;
    while(my $line = <$ifh>){
      $total_feature ++;
      print {$ffh} $line;
    }
    close $ifh or die;
    close $ffh or die;
  }
  $infile = $aFile;
  $infile =~ s/annots/cbcs/;
  open $ifh, "-|", "gunzip", "-cf", $infile or die $!;
  my $bc_count = 0;
  while(my $line = <$ifh>){
    $bc_count ++;
    chomp $line;
    print {$bfh} "${line}-${file_counter}\n";
  }
  close $ifh or die $!;
  $infile = $aFile;
  $infile =~ s/annots/mtx/;
  open $ifh, "-|", "gunzip", "-cf", $infile or die $!;
  while(my $line = <$ifh>){
    next if $. < 4;
    chomp $line;
    my @tmp = split " ", $line;
    $total_entries ++;
    $total_count += $tmp[2];
    my $current_bc = $tmp[1] + $total_bc;
    print {$tfh} "$tmp[0] $current_bc $tmp[2]\n";
  }
  close $ifh or die $!;
  $total_bc += $bc_count;
}
close $bfh or die $!;

open my $mfh, "|-", "gzip -c > $outmtx" or die $!;
print {$mfh} '%%', "MatrixMarket matrix coordinate integer general\n";
print {$mfh} '%', "\n";
print {$mfh} "$total_feature $total_bc $total_entries\n";
open $ifh, "-|", "sort", "-k2,2n", "-k1,1n", $outtmp or die $!;
while(my $line = <$ifh>){
  print {$mfh} $line;
}
close $ifh or die $!;
close $mfh or die $!;
unlink $outtmp or die "Can't delete $outtmp: $!\n";

sub usage{
  print STDERR <<EOF
  Usage: perl TEsingle_aggregate.pl [outprefix] [TEsingle annots files ...]
    Assumes cbcs and matrix files are located in the same folder.

EOF
    ;
  exit 1;
}
