#!/usr/bin/perl
use strict;
use warnings;

use lib '/Users/gregpetrop/perl5/lib/perl5/darwin-thread-multi-2level';
use lib '/Users/gregpetrop/perl5/lib/perl5';
use Chart::Gnuplot;
use lib::stat_mod qw(:DEFAULT); 

my $vol = $ARGV[0];
my $beta = $ARGV[1];
my $mass = $ARGV[2];

my (%t2E, %tdt2E) = ();
my (@avg_t2E, @avg_tdt2E, @std_t2E, @std_tdt2E, @t) = ();
my $base_name = "4flav_${vol}/BlockedWflow_low_${vol}_${beta}_-0.25_${mass}.";
my @files = <$base_name*>;
foreach my $f (@files) {
  print"on file $f\n";
  open(IN, "<", "$f") or die "cannot open < $f: $!";
  chomp(my @in = <IN>);
  close IN;
  my @lines = grep{$_=~/^WFLOW/} @in;
  foreach my $l (@lines) {
    my @split = split(/\s+/, $l);
    push(@{$t2E{$split[1]}}, $split[4]);
    push(@{$tdt2E{$split[1]}}, $split[5]);
    if (!grep{$_ eq $split[1]} @t) {
      push(@t, $split[1]);
    }
  }
}
#plotting
foreach my $tt (@t) {
  push(@avg_t2E, stat_mod::avg(@{$t2E{$tt}}));
  push(@std_t2E, stat_mod::stdev(@{$t2E{$tt}}));
  push(@avg_tdt2E, stat_mod::avg(@{$tdt2E{$tt}}));
  push(@std_tdt2E, stat_mod::stdev(@{$tdt2E{$tt}}));
}
my $chart = Chart::Gnuplot->new(                    #Create chart object 
  output => "Plots/${vol}_${beta}_${mass}.png",
  title  => "Flows",
  xlabel => "t",
);
my $dataSet0 = Chart::Gnuplot::DataSet->new(        #Create dataset object for small volumes
  xdata => \@t,
  ydata => [\@avg_t2E, \@std_t2E],
  title => "Wilson Flow",
  style => "yerrorbars",
);
my $dataSet1 = Chart::Gnuplot::DataSet->new(        #Create dataset object for large volume
  xdata => \@t,
  ydata => [\@avg_tdt2E, \@std_tdt2E],
  title => "Winter Flow",
  style => "yerrorbars",
);
$chart->plot2d($dataSet0, $dataSet1);
