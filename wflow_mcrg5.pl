#!/usr/bin/perl
use strict;
use warnings;

use lib '/Users/gregpetrop/perl5/lib/perl5/darwin-thread-multi-2level';
use lib '/Users/gregpetrop/perl5/lib/perl5';
use lib::stat_mod qw(:DEFAULT);
use PDL;
use PDL::Fit::Polynomial;
use Math::Polynomial::Solve qw(poly_roots);
use Chart::Gnuplot;

my $NF     = $ARGV[0];
my $SmallV = $ARGV[1];                                                            # reads in command line arguments
my $SmallM = $ARGV[2];
my $LargeV = $ARGV[3];
my $LargeM = $ARGV[4];
my %Mass = (
  "$SmallV" => "$SmallM",
  "$LargeV" => "$LargeM",
);

my (%val, %avg, %err, %Beta, %block, %Full_delta_beta) = ();
my @smearingt = ();

print"Reading File:\n";                                                           # gets data
foreach my $vol ($SmallV, $LargeV) {
  my $base_name = "${NF}flav_${vol}/BlockedWflow4_";
  my @files = grep { /BlockedWflow4_(high|low)_${vol}_[0-9].[0-9]_-0.25_$Mass{$vol}/ } glob( "$base_name*" );           # globs for file names
  foreach my $f (@files) {                                                        # loops through matching files
    print"... $f\n";
    my @split = split(/_/, $f);
    my $beta = $split[4];
    if (!grep{$_ eq $beta} @{$Beta{$vol}}) {                                      # constructs Beta hash
      push(@{$Beta{$vol}}, $split[4]);
    }
    open(IN, "<", "$f") or die "cannot open < $f: $!";                            # reads in the file
    chomp(my @in = <IN>);
    close IN;
    my @lines = grep{$_=~/^LOOPS/} @in;                                           # greps for lines with the header LOOPS
    foreach my $l (@lines) {                                                      # loops through matched lines
      my @split = split(/\s+/, $l);                                               # splits matched lines
      if (($split[1] < 0.01) || ($split[1] > 0.4)) {next;}                                                # some of the files have legacy mistake t's that are too large
      push(@{$val{$vol}{$beta}{$split[1]}{$split[2]}{$split[4]}}, $split[6]);     # reads data into hash
      if (!grep{$_ eq $split[1]} @smearingt) {                                    # fills the smearing time array
        push(@smearingt, $split[1]);
      }
      if (!grep{$_ eq $split[4]} @{$block{$vol}}) {                               # fills the number of blockings
        push(@{$block{$vol}}, $split[4]);
      }
    }
  }
  foreach my $beta (@{$Beta{$vol}}) {
    foreach my $loop (0,1,2,3,4) {                                                # loops over observables
      foreach my $b (@{$block{$vol}}) {                                           # beta values
        foreach my $t (@smearingt) {                                              # and smearing times
          $avg{$vol}{$beta}{$t}{$loop}{$b} = stat_mod::avg(@{$val{$vol}{$beta}{$t}{$loop}{$b}});     # to find statistics
          $err{$vol}{$beta}{$t}{$loop}{$b} = stat_mod::stdev(@{$val{$vol}{$beta}{$t}{$loop}{$b}});
        }
      }
    }
  }
}
print"File Read in Complete!\n";

print"Finding Delta Beta:\n";
foreach my $block (@{$block{$LargeV}}) {                                    
  if ($block == 0) {next;}                                                        # skips no matching on the large volume
  foreach my $t (@smearingt) {
    print"... Large blocking:  $block\tSmearing Time: $t\n";
    foreach my $loop (0,1,2,3,4) {
      foreach my $largeb (@{$Beta{$LargeV}}) {                                    # loops over large volume beta
        my $lv_value = $avg{$LargeV}{$largeb}{$t}{$loop}{$block};                 # large volume value at large volume beta
        my @x1 = @{$Beta{$SmallV}};                                               # arrays for plots
        my @y1 = ();
        my @e1 = ();
        my @x2 = @{$Beta{$LargeV}};
        my @y2 = ();
        my @e2 = ();
        my $index=0;
        my $count=0;
        my $n_diff;
        my $o_diff=100;
        foreach my $smallb (@{$Beta{$SmallV}}) {                                  # finds area near intersection
          $n_diff = $avg{$SmallV}{$smallb}{$t}{$loop}{$block-1} - $lv_value;
          if (($n_diff <= $o_diff) && ($n_diff >= 0)) {
            $o_diff = $n_diff;
            $index  = $count;
          }
          $count++;
          push(@y1, $avg{$SmallV}{$smallb}{$t}{$loop}{$block-1});
          push(@e1, $err{$SmallV}{$smallb}{$t}{$loop}{$block-1});
        }
        foreach my $lv_beta (@{$Beta{$LargeV}}) {
          push(@y2, $avg{$LargeV}{$largeb}{$t}{$loop}{$block});
          push(@e2, $err{$LargeV}{$largeb}{$t}{$loop}{$block});
        }
        my (@x,@y,@e) =();                                                        # declares fitting arrays
        if ($index < 2) {
          @x = @{$Beta{$SmallV}}[0..3];
          @y = @y1[0..3];
          @e = @e1[0..3];
        }
        elsif ($index > (@{$Beta{$SmallV}}-3)) {
          @x = @{$Beta{$SmallV}}[(@{$Beta{$SmallV}}-4)..(@{$Beta{$SmallV}}-1)];
          @y = @y1[(@y1-4)..(@y1-1)];
          @e = @e1[(@e1-4)..(@e1-1)];
        }
        else{
          @x = @{$Beta{$SmallV}}[($index-1)..($index+2)];
          @y = @y1[($index-1)..($index+2)];
          @e = @e1[($index-1)..($index+2)];
        }
        my $x=pdl(@x);                                                            # puts small volume data into piddle for fitting
        my $y=pdl(@y);
        my $e=pdl(@e);
        (my $fit ,my $coeffs)=fitpoly1d $x, $y, 2;                                # fits the small volumes
        my $a=$coeffs->at(1);                                                     # extracts out the coefficients
        my $b=$coeffs->at(0);
        #my $a=$coeffs->at(3);                                                     # extracts out the coefficients
        #my $b=$coeffs->at(2);
        #my $c=$coeffs->at(1);
        #my $d=$coeffs->at(0);
        my @roots=poly_roots($a,($b-$lv_value));                                  # solves for the difference between the fit and the large mass value
        my $beta_diff;
        #my $hasroot = 0;
        foreach my $r (@roots){
          #if (($r < 4.0)||($r =~ /i$/)){next;}          #skips imaginary roots
          #if (($r < $Beta{$SmallV}[0])||($r =~ /i$/)){next;}          #skips imaginary roots
          $Full_delta_beta{$block}{$largeb}{$loop}{$t} = $largeb - $r;            # gets the right root
          $beta_diff = $r;                                                        # right now only written for linear fit
          #$hasroot++;
        }
        #if ($hasroot == 0) {
          #$Full_delta_beta{$block}{$largeb}{$loop}{$t}="NaN";
        #}
        my $chart = Chart::Gnuplot->new(                                          # Create chart object 
          output => "Plots_${NF}flav5/deltabeta/${largeb}_${block}_${t}_${loop}_full.png",
          title  => "Deltabeta for beta ${largeb}, matching with ${LargeV} blocked ${block} after ${t} smearing",
          xlabel => "Beta",
          ylabel => "Expectation Value",
        );
        $chart->command("set obj 1 rectangle behind from screen 0,0 to screen 1,1");
        $chart->command("set obj 1 fillstyle solid 1.0 fillcolor rgbcolor \"white\"");
        $chart->command("set label 1 \"Delta Beta:  $Full_delta_beta{$block}{$largeb}{$loop}{$t}\"");
        $chart->command("set label 1 at graph 0.02, 0.85 tc lt 3");
        #if ($hasroot > 0) {
          $chart->command("set arrow from $largeb,$lv_value to $beta_diff,$lv_value");
        #}
        my $dataSet0 = Chart::Gnuplot::DataSet->new(                              # Create dataset object for small volumes
          xdata => \@x1,
          ydata => [\@y1, \@e1],
          title => "Small Volume Observable: ${loop}",
          style => "yerrorbars",
        );
        my $dataSet1 = Chart::Gnuplot::DataSet->new(                              # Create dataset object for large volume
          xdata => \@x2,
          ydata => [\@y2, \@e2],
          title => "Large Volume Observable: ${loop}",
          style => "yerrorbars",
        );
        my $dataSet2 = Chart::Gnuplot::DataSet->new(                              # Create dataset object for the fit
          #func => "$a*x**3+$b*x**2+$c*x+$d",
          func => "$a*x+$b",
          title => "Fit to Small Volume",
        );
        $chart->plot2d($dataSet0, $dataSet1, $dataSet2);                          # plots the chart
      }
    }
  }
}
print"Finding Delta Beta Complete!\n";

my (%T_optimal, %Delta_beta_optimal) = ();

print"Finding Optimal Smearing Time:\n";
foreach my $largeb (@{$Beta{$LargeV}}) {
  print"... Beta:  $largeb\n";
  my (@t_opt, @b_opt) = ();
  foreach my $loop (0,1,2,3,4) {
    my @x = @smearingt;
    my @y1 = ();
    my @y2 = ();
    my @y3 = ();
    my $count = 0;
    my $index = 0;
    my $die = 0;
    foreach my $t (@smearingt) {
      my $diff;
      if (($Full_delta_beta{2}{$largeb}{$loop}{$t} =~ "NaN") || ($Full_delta_beta{3}{$largeb}{$loop}{$t} =~ "NaN")){$diff="NaN";}
      $diff = $Full_delta_beta{2}{$largeb}{$loop}{$t} - $Full_delta_beta{3}{$largeb}{$loop}{$t};
      if (($diff > 0) && ($die == 0)) {$index = $count; $die = 1}
      $count ++;
      push(@y1,$Full_delta_beta{1}{$largeb}{$loop}{$t});
      push(@y2,$Full_delta_beta{2}{$largeb}{$loop}{$t});
      push(@y3,$Full_delta_beta{3}{$largeb}{$loop}{$t});
    }

    my @xf=@smearingt[($index-1)..$index];
    my @yf1=@y2[($index-1)..$index];
    my @yf2=@y3[($index-1)..$index];
    my $x1=pdl(@xf);                                                              # puts data into a piddle for fitting
    my $y1=pdl(@yf1);
    my $x2=pdl(@xf);                                                              # puts data into a piddle for fitting
    my $y2=pdl(@yf2);
    (my $fit1,my $coeffs1)=fitpoly1d $x1, $y1, 2;                                 # fits the data
    (my $fit2,my $coeffs2)=fitpoly1d $x2, $y2, 2;
    my $b1=$coeffs1->at(0);                                                       # extracts out the coefficients
    my $a1=$coeffs1->at(1);
    my $b2=$coeffs2->at(0);      
    my $a2=$coeffs2->at(1);

    my @roots=poly_roots($a1-$a2,$b1-$b2);                                        # finds intersection i.e. optimal smearing time
    foreach my $r (@roots) {
      if ($r =~ /i$/){next;}                                                      # skips imaginary roots
      $T_optimal{$loop}{$largeb}=$r;
      $Delta_beta_optimal{$loop}{$largeb}=$a1*$r+$b1;
      print"$loop\t$r\t$Delta_beta_optimal{$loop}{$largeb}\n";
      push(@t_opt,$r);
      push(@b_opt, $a1*$r+$b1);
    }

    my $chart = Chart::Gnuplot->new(                                              # Create chart object 
      output => "Plots_${NF}flav5/smearing_time/${largeb}_${loop}.png",
      title  => "Delta Beta as a function of smearing time for Beta: ${largeb}",
      xlabel => "Smearing Time",
      ylabel => "Delta Beta",
      #yrange => [0,1]
    );
    $chart->command("set obj 1 rectangle behind from screen 0,0 to screen 1,1");
    $chart->command("set obj 1 fillstyle solid 1.0 fillcolor rgbcolor \"white\"");
    $chart->command("set label 1 \"Optimal Smearing Time:  $T_optimal{$loop}{$largeb}\"");
    $chart->command("set label 1 at graph 0.02, 0.85 tc lt 3");
    $chart->command("set label 2 \"Optimal Delta Beta:     $Delta_beta_optimal{$loop}{$largeb}\"");
    $chart->command("set label 2 at graph 0.02, 0.75 tc lt 3");
    my $dataSet1 = Chart::Gnuplot::DataSet->new(                                  # Create dataset object for small volumes
      xdata => \@x,
      ydata => \@y1,
      title => "Large volume blocked once: ${loop}",
    );
    my $dataSet2 = Chart::Gnuplot::DataSet->new(                                  # Create dataset object for large volume
      xdata => \@x,
      ydata => \@y2,
      title => "Large volume blocked twice: ${loop}",
    );
    my $dataSet3 = Chart::Gnuplot::DataSet->new(                                  # Create dataset object for the fit
      xdata => \@x,
      ydata => \@y3,
      title => "Large volume blocked thrice: ${loop}",
    );
    my $dataSet4 = Chart::Gnuplot::DataSet->new(                                  # Create dataset object for the fit
      func => "$a1*x+$b1",
      title => "Fit: Blocked Twice",
    );
    my $dataSet5 = Chart::Gnuplot::DataSet->new(                                  # Create dataset object for the fit
      func => "$a2*x+$b2",
      title => "Fit: Blocked Thrice",
    );
    $chart->plot2d($dataSet1, $dataSet2, $dataSet3, $dataSet4, $dataSet5);
  }

  my $b_avg = stat_mod::avg(@b_opt);
  my $b_std = stat_mod::stdev(@b_opt);
  print"AVG $b_avg\tSTDEV $b_std\n";
}
print"Finding Optimal Smearing Time Complete!\n";
