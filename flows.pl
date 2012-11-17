#!/usr/bin/perl
use strict;
use warnings;

use lib::stat_mod qw(:DEFAULT); 

my $vol  = $ARGV[0];
my $beta = $ARGV[1];
my $mass = $ARGV[2];

my @w0 = ();
my @t0 = ();
my $base_name = "4flav_${vol}/BlockedWflow_low_${vol}_${beta}_-0.25_${mass}.";

my @files = <$base_name*>;

foreach my $f (@files) {
  open(IN, "<", "$f") or die "cannot open < $f: $!";
  chomp(my @in = <IN>);
  close IN;
  my @lines = grep{$_=~/^WFLOW/} @in;
  my $tfirst = 0;
  my $wfirst = 0;
  foreach my $l (@lines) {
    my @split = split(/\s+/, $l);
    if (($split[4] >= 0.3) && ($tfirst == 0)) {
      push(@t0, $split[1]);
      $tfirst=1;
    }
    if (($split[5] >= 0.3) && ($wfirst == 0)) {
      push(@w0, sqrt($split[1]));
      $wfirst = 1;
    }
    if (($wfirst == 1) && ($tfirst == 1)) {
      last;
    }
  }
  $tfirst = 0;
  $wfirst = 0;
  my @backlines = reverse(@lines);
  foreach my $l (@backlines) {
    my @split = split(/\s+/, $l);
    if (($split[4] <= 0.3) && ($tfirst == 0)) {
      push(@t0, $split[1]);
      $tfirst=1;
    }
    if (($split[5] <= 0.3) && ($wfirst == 0)) {
      push(@w0, sqrt($split[1]));
      $wfirst = 1;
    }
    if (($wfirst == 1) && ($tfirst == 1)) {
      last;
    }
  }
}
print "$#t0\t$#w0\n";
my $t0     = stat_mod::avg(@t0);
my $t0_err = stat_mod::stdev(@t0);
my $w0     = stat_mod::avg(@w0);
my $w0_err = stat_mod::stdev(@w0);
my $a_from_t0 = sqrt(0.096/$t0);
my $a_from_w0 = sqrt(.1755/$w0);

print "*\t*\t*\nBeta:  ${beta}\nMass:  ${mass}\n*\t*\t*\nWILSON:  $t0 $t0_err $a_from_t0\nWINTER:  $w0 $w0_err $a_from_w0\n";
