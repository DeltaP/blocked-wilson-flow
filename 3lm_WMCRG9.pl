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
my $MaxBlock = $ARGV[1];
my $SmallV = $ARGV[2];                                                            # reads in command line arguments
my $SmallM = $ARGV[3];
my $MediumV = $ARGV[4];                                                            # reads in command line arguments
my $MediumM = $ARGV[5];
my $LargeV = $ARGV[6];
my $LargeM = $ARGV[7];
my %Mass = (
  "$SmallV" => "$SmallM",
  "$MediumV"=> "$MediumM",
  "$LargeV" => "$LargeM",
);

sub zup {
  join "\n", map { my $i = $_; join ' ', map $_->[ $i ], @_ } 0 .. $#{
+ $_[0] }
}

my (%val, %avg, %err, %Beta, %block, %Full_delta_beta) = ();
my @smearingt = ();

print"Reading File:\n";                                                           # gets data
foreach my $vol ($SmallV, $MediumV, $LargeV) {
  my $base_name = "${NF}flav_${vol}/WMCRG9_";
  my @files = grep { /WMCRG9_(high|low|mix)_${vol}_[0-9].[0-9]_-0.25_$Mass{$vol}/ } glob( "$base_name*" );           # globs for file names
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
      push(@{$val{$vol}{$beta}{$split[1]}{$split[2]}{$split[4]}}, $split[6]);     # reads data into hash
      if (!grep{$_ eq $split[1]} @smearingt) {                                    # fills the smearing time array
        push(@smearingt, $split[1]);
      }
      if (!grep{$_ eq $split[4]} @{$block{$vol}}) {                               # fills the number of blockings
        push(@{$block{$vol}}, $split[4]);
      }
    }
  }
  my @sorted_beta = sort { $a <=> $b } @{$Beta{$vol}};
  $Beta{$vol} = [@sorted_beta];
  print"$vol\t@{$Beta{$vol}}\n";
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
my $block = $MaxBlock;                                  
foreach my $t (@smearingt) {
  print"... Large blocking:  $block\tSmearing Time: $t\n";
  foreach my $loop (0,1,2,3,4) {
    foreach my $largeb (@{$Beta{$LargeV}}) {                                    # loops over large volume beta
      my $lv_value = $avg{$LargeV}{$largeb}{$t}{$loop}{$block};                 # large volume value at large volume beta
      my @x1 = @{$Beta{$MediumV}};                                              # arrays for plots
      my @y1 = ();
      my @e1 = ();
      my @x2 = @{$Beta{$LargeV}};
      my @y2 = ();
      my @e2 = ();

      foreach my $smallb (@{$Beta{$MediumV}}) {                                 # finds area near intersection
        push(@y1, $avg{$MediumV}{$smallb}{$t}{$loop}{$block-1});
        push(@e1, $err{$MediumV}{$smallb}{$t}{$loop}{$block-1});
      }
      foreach my $lv_beta (@{$Beta{$LargeV}}) {
        push(@y2, $avg{$LargeV}{$lv_beta}{$t}{$loop}{$block});
        push(@e2, $err{$LargeV}{$lv_beta}{$t}{$loop}{$block});
      }

      my $x1=pdl(@x1);                                                          # puts small volume data into piddle for fitting
      my $y1=pdl(@y1);
      my $e1=pdl(@e1);
      (my $fit1 , my $coeffs1)=fitpoly1d $x1, $y1, 4;                           # fits the small volumes
      my $a1=$coeffs1->at(3);                                                   # extracts out the coefficients
      my $b1=$coeffs1->at(2);
      my $c1=$coeffs1->at(1);
      my $d1=$coeffs1->at(0);
      my $temp1=pdl(@e1);
      my $tt = $fit1-$y1;
      $temp1=(($fit1-$y1)**2/$e1**2);
      my $chi1=sum $temp1;
      $chi1/=(@x1-4-1);
      #print"CHI^2 small volume:  $chi1\n";

      my $x2=pdl(@x2);                                                          # puts small volume data into piddle for fitting
      my $y2=pdl(@y2);
      my $e2=pdl(@e2);
      (my $fit2 , my $coeffs2)=fitpoly1d $x2, $y2, 4;                           # fits the small volumes
      my $a2=$coeffs2->at(3);                                                   # extracts out the coefficients
      my $b2=$coeffs2->at(2);
      my $c2=$coeffs2->at(1);
      my $d2=$coeffs2->at(0);
      my $temp2=pdl(($fit2-$y2)**2/$e2**2);
      my $chi2=sum $temp2;
      $chi2/=(@x1-4-1);
      #print"CHI^2 large volume:  $chi2\n";

      my $lv_fit_value = $a2*$largeb**3+$b2*$largeb**2+$c2*$largeb+$d2;
      my @roots=poly_roots(($a1),($b1),($c1),($d1-$lv_fit_value));              # solves for the difference between the fit and the large mass value
      my $beta_diff;
      my $hasroot = 0;
      foreach my $r (@roots){
        if ($r =~ /i$/){next;}
        if (($r < $x1[0])||($r>$x1[$#x1])){next;}                               #skips imaginary roots
        $Full_delta_beta{$block}{$largeb}{$loop}{$t} = $largeb - $r;            # gets the right root
        $beta_diff = $r;                                                        # right now only written for linear fit
        $hasroot++;
      }
      if ($hasroot != 1) {
        $Full_delta_beta{$block}{$largeb}{$loop}{$t}="NaN";
        #print"$hasroot number of roots found:  $t\t$loop\t$largeb\n";
        next;
      }
      my $chart = Chart::Gnuplot->new(                                          # Create chart object 
        output => "Plots_${NF}_WMCRG9/deltabeta/${largeb}_${block}_${t}_${loop}_full.png",
        title  => "Deltabeta for beta ${largeb}, matching with ${LargeV} blocked ${block} after ${t} smearing",
        xlabel => "Beta",
        ylabel => "Expectation Value",
      );
      $chart->command("set obj 1 rectangle behind from screen 0,0 to screen 1,1");
      $chart->command("set obj 1 fillstyle solid 1.0 fillcolor rgbcolor \"white\"");
      $chart->command("set label 1 \"Delta Beta:  $Full_delta_beta{$block}{$largeb}{$loop}{$t}\"");
      $chart->command("set label 1 at graph 0.02, 0.85 tc lt 3");
      $chart->command("set label 2 \"Small V chi^2/dof:  $chi1\"");
      $chart->command("set label 2 at graph 0.02, 0.75 tc lt 3");
      $chart->command("set label 3 \"Large V chi^2/dof:  $chi2\"");
      $chart->command("set label 3 at graph 0.02, 0.65 tc lt 3");
      #if ($hasroot > 0) {
        $chart->command("set arrow from $largeb,$lv_fit_value to $beta_diff,$lv_fit_value");
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
        func => "$a1*x**3+$b1*x**2+$c1*x+$d1",
        title => "Fit to Small Volume",
      );
      my $dataSet3 = Chart::Gnuplot::DataSet->new(                              # Create dataset object for the fit
        func => "$a2*x**3+$b2*x**2+$c2*x+$d2",
        title => "Fit to Large Volume",
      );
      $chart->plot2d($dataSet0, $dataSet1, $dataSet2, $dataSet3);               # plots the chart
    }
  }
}
$block--;                                  
foreach my $t (@smearingt) {
  print"... Large blocking:  $block\tSmearing Time: $t\n";
  foreach my $loop (0,1,2,3,4) {
    foreach my $largeb (@{$Beta{$LargeV}}) {                                    # loops over large volume beta
      my $lv_value = $avg{$MediumV}{$largeb}{$t}{$loop}{$block};                # large volume value at large volume beta
      my @x1 = @{$Beta{$SmallV}};                                               # arrays for plots
      my @y1 = ();
      my @e1 = ();
      my @x2 = @{$Beta{$MediumV}};
      my @y2 = ();
      my @e2 = ();

      foreach my $smallb (@{$Beta{$SmallV}}) {                                  # finds area near intersection
        push(@y1, $avg{$SmallV}{$smallb}{$t}{$loop}{$block-1});
        push(@e1, $err{$SmallV}{$smallb}{$t}{$loop}{$block-1});
      }
      foreach my $lv_beta (@{$Beta{$MediumV}}) {
        push(@y2, $avg{$MediumV}{$lv_beta}{$t}{$loop}{$block});
        push(@e2, $err{$MediumV}{$lv_beta}{$t}{$loop}{$block});
      }

      my $x1=pdl(@x1);                                                          # puts small volume data into piddle for fitting
      my $y1=pdl(@y1);
      my $e1=pdl(@e1);
      (my $fit1 , my $coeffs1)=fitpoly1d $x1, $y1, 4;                           # fits the small volumes
      my $a1=$coeffs1->at(3);                                                   # extracts out the coefficients
      my $b1=$coeffs1->at(2);
      my $c1=$coeffs1->at(1);
      my $d1=$coeffs1->at(0);
      my $temp1=pdl(@e1);
      my $tt = $fit1-$y1;
      $temp1=(($fit1-$y1)**2/$e1**2);
      my $chi1=sum $temp1;
      $chi1/=(@x1-4-1);
      #print"CHI^2 small volume:  $chi1\n";

      my $x2=pdl(@x2);                                                          # puts small volume data into piddle for fitting
      my $y2=pdl(@y2);
      my $e2=pdl(@e2);
      (my $fit2 , my $coeffs2)=fitpoly1d $x2, $y2, 4;                           # fits the small volumes
      my $a2=$coeffs2->at(3);                                                   # extracts out the coefficients
      my $b2=$coeffs2->at(2);
      my $c2=$coeffs2->at(1);
      my $d2=$coeffs2->at(0);
      my $temp2=pdl(($fit2-$y2)**2/$e2**2);
      my $chi2=sum $temp2;
      $chi2/=(@x1-4-1);
      #print"CHI^2 large volume:  $chi2\n";

      my $lv_fit_value = $a2*$largeb**3+$b2*$largeb**2+$c2*$largeb+$d2;
      my @roots=poly_roots(($a1),($b1),($c1),($d1-$lv_fit_value));              # solves for the difference between the fit and the large mass value
      my $beta_diff;
      my $hasroot = 0;
      foreach my $r (@roots){
        if ($r =~ /i$/){next;}
        if (($r < $x1[0])||($r>$x1[$#x1])){next;}                               # skips imaginary roots
        $Full_delta_beta{$block}{$largeb}{$loop}{$t} = $largeb - $r;            # gets the right root
        $beta_diff = $r;                                                        # right now only written for linear fit
        $hasroot++;
      }
      if ($hasroot != 1) {
        $Full_delta_beta{$block}{$largeb}{$loop}{$t}="NaN";
        #print"$hasroot number of roots found:  $t\t$loop\t$largeb\n";
        next;
      }
      my $chart = Chart::Gnuplot->new(                                          # Create chart object 
        output => "Plots_${NF}_WMCRG9/deltabeta/${largeb}_${block}_${t}_${loop}_full.png",
        title  => "Deltabeta for beta ${largeb}, matching with ${MediumV} blocked ${block} after ${t} smearing",
        xlabel => "Beta",
        ylabel => "Expectation Value",
      );
      $chart->command("set obj 1 rectangle behind from screen 0,0 to screen 1,1");
      $chart->command("set obj 1 fillstyle solid 1.0 fillcolor rgbcolor \"white\"");
      $chart->command("set label 1 \"Delta Beta:  $Full_delta_beta{$block}{$largeb}{$loop}{$t}\"");
      $chart->command("set label 1 at graph 0.02, 0.85 tc lt 3");
      $chart->command("set label 2 \"Small V chi^2/dof:  $chi1\"");
      $chart->command("set label 2 at graph 0.02, 0.75 tc lt 3");
      $chart->command("set label 3 \"Large V chi^2/dof:  $chi2\"");
      $chart->command("set label 3 at graph 0.02, 0.65 tc lt 3");
      #if ($hasroot > 0) {
        $chart->command("set arrow from $largeb,$lv_fit_value to $beta_diff,$lv_fit_value");
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
        func => "$a1*x**3+$b1*x**2+$c1*x+$d1",
        title => "Fit to Small Volume",
      );
      my $dataSet3 = Chart::Gnuplot::DataSet->new(                              # Create dataset object for the fit
        func => "$a2*x**3+$b2*x**2+$c2*x+$d2",
        title => "Fit to Large Volume",
      );
      $chart->plot2d($dataSet0, $dataSet1, $dataSet2, $dataSet3);                          # plots the chart
    }
  }
}
print"Finding Delta Beta Complete!\n";

my (%T_optimal, %Delta_beta_optimal, %T_optimal_avg) = ();

#averages over different loops
print"Interpolating Delta Beta Averaging over Loops First:\n";
my (%Avg_t_optimal,%Avg_delta_beta_optimal) = ();
my (%Avg_t_optimal_1,%Avg_delta_beta_optimal_1) = ();
my (%Avg_t_optimal_2,%Avg_delta_beta_optimal_2) = ();
my (%Avg_t_optimal_3,%Avg_delta_beta_optimal_3) = ();
my (%Avg_t_optimal_4,%Avg_delta_beta_optimal_4) = ();
foreach my $largeb (@{$Beta{$LargeV}}) {
  print"... Beta:  $largeb\n";
  my (@db2_avg, @db2_err, @db3_avg, @db3_err, @plott) = ();
  foreach my $t (@smearingt) {
    my (@db2, @db3) = ();
    foreach my $loop (0,1,2,3,4) {
      if ($Full_delta_beta{$block}{$largeb}{$loop}{$t} ne 'NaN') {
        push(@db2, $Full_delta_beta{$block}{$largeb}{$loop}{$t});
      }
      if ($Full_delta_beta{$block+1}{$largeb}{$loop}{$t} ne 'NaN') {
        push(@db3, $Full_delta_beta{$block+1}{$largeb}{$loop}{$t});
      }
    }
    if ((@db2 < 2)||(@db3 < 2)) {next;}
    push(@db2_avg,stat_mod::avg(@db2));
    push(@db2_err,stat_mod::stdev(@db2));
    push(@db3_avg,stat_mod::avg(@db3));
    push(@db3_err,stat_mod::stdev(@db3));
    push(@plott,$t);
  }
  if ((@db2_avg < 2)||(@db3_avg < 2)) {next;}

  my $count = 0;
  my $index = 0;
  my $die = 0;
  foreach my $t (@plott) {
    my $diff;
    $diff = $db2_avg[$count] - $db3_avg[$count];
    if (($diff > 0) && ($die == 0)) {$index = $count; $die = 1}
    $count ++;
  }

  my $x=pdl(@plott[$index-1..$index+1]);                                                              # puts data into a piddle for fitting
  my $y2=pdl(@db2_avg[$index-1..$index+1]);
  my $y2_err=pdl(@db2_err[$index-1..$index+1]);
  my $y3=pdl(@db3_avg[$index-1..$index+1]);
  my $y3_err=pdl(@db3_err[$index-1..$index+1]);
  my $y2_plus=$y2+$y2_err;
  my $y2_minus=$y2-$y2_err;
  my $y3_plus=$y3+$y3_err;
  my $y3_minus=$y3-$y3_err;

  (my $fit2,my $coeffs2)=fitpoly1d $x, $y2, 2;
  my $b2=$coeffs2->at(0);      
  my $a2=$coeffs2->at(1);
  (my $fit2_plus,my $coeffs2_plus)=fitpoly1d $x, $y2_plus, 2;
  my $b2_plus=$coeffs2_plus->at(0);      
  my $a2_plus=$coeffs2_plus->at(1);
  (my $fit2_minus,my $coeffs2_minus)=fitpoly1d $x, $y2_minus, 2;
  my $b2_minus=$coeffs2_minus->at(0);      
  my $a2_minus=$coeffs2_minus->at(1);
  (my $fit3,my $coeffs3)=fitpoly1d $x, $y3, 2;
  my $b3=$coeffs3->at(0);      
  my $a3=$coeffs3->at(1);
  (my $fit3_plus,my $coeffs3_plus)=fitpoly1d $x, $y3_plus, 2;
  my $b3_plus=$coeffs3_plus->at(0);      
  my $a3_plus=$coeffs3_plus->at(1);
  (my $fit3_minus,my $coeffs3_minus)=fitpoly1d $x, $y3_minus, 2;
  my $b3_minus=$coeffs3_minus->at(0);      
  my $a3_minus=$coeffs3_minus->at(1);

  my @roots=poly_roots($a2-$a3,$b2-$b3);                #finds intersection i.e. alpha optimal
  foreach my $r (@roots){
    $Avg_t_optimal{$largeb}=$r;
    $Avg_delta_beta_optimal{$largeb}=$a2*$r+$b2;
    print"$largeb\t$r\t$Avg_delta_beta_optimal{$largeb}\n";
  }

  @roots=poly_roots($a2_plus-$a3_plus,$b2_plus-$b3_plus);                #finds intersection i.e. alpha optimal
  foreach my $r (@roots){
    $Avg_t_optimal_1{$largeb}=$r;
    $Avg_delta_beta_optimal_1{$largeb}=$a2_plus*$r+$b2_plus;
    print"$largeb\t$r\t$Avg_delta_beta_optimal_1{$largeb}\n";
  }
  @roots=poly_roots($a2_minus-$a3_minus,$b2_minus-$b3_minus);                #finds intersection i.e. alpha optimal
  foreach my $r (@roots){
    $Avg_t_optimal_2{$largeb}=$r;
    $Avg_delta_beta_optimal_2{$largeb}=$a2_minus*$r+$b2_minus;
    print"$largeb\t$r\t$Avg_delta_beta_optimal_2{$largeb}\n";
  }
  @roots=poly_roots($a2_minus-$a3_plus,$b2_minus-$b3_plus);                #finds intersection i.e. alpha optimal
  foreach my $r (@roots){
    $Avg_t_optimal_3{$largeb}=$r;
    $Avg_delta_beta_optimal_3{$largeb}=$a2_minus*$r+$b2_minus;
    print"$largeb\t$r\t$Avg_delta_beta_optimal_3{$largeb}\n";
  }
  @roots=poly_roots($a2_plus-$a3_minus,$b2_plus-$b3_minus);                #finds intersection i.e. alpha optimal
  foreach my $r (@roots){
    $Avg_t_optimal_4{$largeb}=$r;
    $Avg_delta_beta_optimal_4{$largeb}=$a2_plus*$r+$b2_plus;
    print"$largeb\t$r\t$Avg_delta_beta_optimal_4{$largeb}\n";
  }
  print"PLOT ";
  my $high_err=($Avg_delta_beta_optimal_1{$largeb}+$Avg_delta_beta_optimal_3{$largeb})/2;
  my $low_err=($Avg_delta_beta_optimal_2{$largeb}+$Avg_delta_beta_optimal_4{$largeb})/2;
  print"$largeb $Avg_t_optimal_1{$largeb} $Avg_delta_beta_optimal{$largeb} $low_err $high_err\n";


  my $chart = Chart::Gnuplot->new(                                              # Create chart object 
    output => "Plots_${NF}_WMCRG9/avg_smearing_time/avg_${largeb}.png",
    title  => "Delta Beta as a function of smearing time for Beta: ${largeb}",
    xlabel => "Smearing",
    ylabel => "Delta Beta",
  );
    $chart->command("set obj 1 rectangle behind from screen 0,0 to screen 1,1");
    $chart->command("set obj 1 fillstyle solid 1.0 fillcolor rgbcolor \"white\"");
    $chart->command("set label 1 \"Optimal Smearing Time:  $Avg_t_optimal{$largeb}\"");
    $chart->command("set label 1 at graph 0.02, 0.85 tc lt 3");
    $chart->command("set label 2 \"Optimal Delta Beta:  $Avg_delta_beta_optimal{$largeb}\"");
    $chart->command("set label 2 at graph 0.02, 0.75 tc lt 3");
    my $dataSet1 = Chart::Gnuplot::DataSet->new(                                  # Create dataset object for large volume
      xdata => \@plott,
      ydata => [\@db2_avg,\@db2_err],
      title => "Large volume blocked twice",
      style => "yerrorbars",
    );
    my $dataSet2 = Chart::Gnuplot::DataSet->new(                                  # Create dataset object for large volume
      xdata => \@plott,
      ydata => [\@db3_avg,\@db3_err],
      title => "Large volume blocked thrice",
      style => "yerrorbars",
    );
    my $dataSet3 = Chart::Gnuplot::DataSet->new(                                  # Create dataset object for the fit
      func => "$a2_plus*x+$b2_plus",
      title => "Fit: Blocked 2+error",
    );
    my $dataSet4 = Chart::Gnuplot::DataSet->new(                                  # Create dataset object for the fit
      func => "$a2_minus*x+$b2_minus",
      title => "Fit: Blocked 2-error",
    );
    my $dataSet5 = Chart::Gnuplot::DataSet->new(                                  # Create dataset object for the fit
      func => "$a2*x+$b2",
      title => "Fit: Blocked 2",
    );
    my $dataSet6 = Chart::Gnuplot::DataSet->new(                                  # Create dataset object for the fit
      func => "$a3_plus*x+$b3_plus",
      title => "Fit: Blocked 3+error",
    );
    my $dataSet7 = Chart::Gnuplot::DataSet->new(                                  # Create dataset object for the fit
      func => "$a3_minus*x+$b3_minus",
      title => "Fit: Blocked 3-error",
    );
    my $dataSet8 = Chart::Gnuplot::DataSet->new(                                  # Create dataset object for the fit
      func => "$a3*x+$b3",
      title => "Fit: Blocked 3",
    );
    $chart->plot2d($dataSet1, $dataSet2, $dataSet3, $dataSet4, $dataSet5, $dataSet6, $dataSet7, $dataSet8);
    open FILE, ">", "Plots_${NF}_WMCRG9/out/avg_${largeb}" or die $!;
    print FILE zup \(@plott, @db2_avg, @db2_err, @db3_avg, @db3_err); 
    print FILE "\n"; 
    close FILE;
}
print"Extrapolating Delta Beta Averaging over Loops First Complete!\n";
