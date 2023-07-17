#!/usr/bin/env perl

use warnings;
use strict;
use Util::H2O::More qw/ddd Getopt2h2o/;

# build and load subroutines
use Alien::OpenMP;
use OpenMP::Environment;
use Inline ( 
  C    => 'DATA',
  with => qw/Alien::OpenMP/,
);

my $P           = 5;
my $V           = undef;
my $ABSOLUTEMIN = 1;

my $o = Getopt2h2o \@ARGV, {
    eps     => 0.1,
    f       => q{./fort.14},
    minPts  => 5,
    maxSize => undef,
    threads => 2,
  },
  qw/eps=s f=s minPts=i maxSize=i threads=i/;

my $env = OpenMP::Environment->new();

$env->omp_num_threads($o->threads);


sub show_help {
    print <<END;

Reference implementation of row based grid partitioning for ADCIRC grid files. 

Input

  ADCIRC grid file.

Output

  ADCIRC partmesh.txt, used by adcprep >= 47.x to define subdomains.

Usage

  adcpart-smart.pl -d <numdomains> [-f ./fort.14] [-p <precision>] [-m <min_nodes_per_bin>] > partmesh.txt 

 ... then run adcprep, skipping the option (1) that creates partmesh.txt.

Options        Description
   -n             Number of subdomains
   -f             Specifies file (default is ./fort.14)
END
    exit;
}

# number of significant figures to use when counting nodes as being in the same longitude
my $precision = 5;
my %bins      = ();

# this "higher order" function reads over each node in the fort.14 and applies
# $sub_ref (user provided subroutine reference) to the longitudinal position, $y
sub process_f14 {
    my ($o)      = @_;
    my $NODES    = {};
    my $ELEMENTS = {};

    # open file - put in initial set of bins
    open my $F14, q{<}, $o->f or die $!;

    my $mesh_name = <$F14>;
    chomp $mesh_name;

    my $line = <$F14>;
    $line =~ m/(\d+)\s+(\d+)/;
    my $NE = $1;
    my $NN = $2;

    $o->maxSize(4000) if not $o->maxSize();

    # collect node positions
    for ( my $i = 1; $i <= $NN; ++$i ) {
        my $line = <$F14>;
        $line =~ m/^\s*(\d+)\s+([+\-]*\d+\.\d+)\s+([+\-]*\d+\.\d+)/;
        my $node = $1;
        my $x    = $2;
        my $y    = $3;
        $NODES->{$node} = { x => $x, y => $y };
    }

    # collect element members
    my $sum = 0;
    for ( my $i = 1; $i <= $NE; ++$i ) {
        my $line = <$F14>;
        $line =~ m/^\s*(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)/;
        my $element = $1;
        $ELEMENTS->{$element} = { nodes => [ $3, $4, $5 ] };
       foreach my $n ( $3, $4, $5 ) {
            $NODES->{$n}->{e}->{$element} = undef;
            map { $NODES->{$n}->{k}->{$_} = undef } ( $3, $4, $5 );

# stab at estimating epsilon?
            my $x = $NODES->{$3}->{x} - $NODES->{$4}->{x};
            my $y = $NODES->{$3}->{y} - $NODES->{$4}->{y};
            my $d1= sqrt($x**2 + $y**2);
            $x = $NODES->{$3}->{x} - $NODES->{$5}->{x};
            $y = $NODES->{$3}->{y} - $NODES->{$5}->{y};
            my $d2= sqrt($x**2 + $y**2);
            $x = $NODES->{$4}->{x} - $NODES->{$5}->{x};
            $y = $NODES->{$4}->{y} - $NODES->{$5}->{y};
            my $d3   = sqrt($x**2 + $y**2);
            $sum += $d1 + $d2 + $d3;

        }
    }

    #$o->eps($sum / (3*$NE));


    close $F14;

    my $c = 1;
    my $S = DBSCAN( $NODES, \&distance, $o->eps, $o->minPts );
    foreach my $C (sort { $a <=> $b } keys %$S) {
      foreach my $n (@{$S->{$C}}) {
        my $x = $NODES->{$n}->{x};
        my $y = $NODES->{$n}->{y};
        printf qq{%d %s %s %d\n}, $c++, $x, $y, $C;
      }
    }
}

sub DBSCAN {
    my ( $NODES, $distFunc, $eps, $minPts ) = @_;
    my $C        = 0;                                             # initialize cluster counter
    my $labels   = {};
    my $Clusters = {};
  NODE:
    foreach my $P ( sort {$a <=> $b} keys %$NODES ) {
        #printf qq{Processing Node % 10d\n}, $P;
        next NODE if $labels->{$P};                             # continue to next $P if node has been labeled
        my $N = rangeQuery( $NODES, $distFunc, $P, $eps );      # get neighbors of $P
        #printf qq{\tP Node %d has %d neighbors ...\n}, $P, scalar @$N if @$N;
        if ( @$N < $minPts ) {                                  # density check 
            #printf qq{\tMarking P Node %d as "noise" ...\n}, $P;
            $labels->{$P} = q{noise};                           # mark $P as "noise"
            next NODE;                                          # continue to next $P, not a "core" point
        }
        #if ( $C == 0 or ( $Clusters->{$C} and @{$Clusters->{$C}} >= $o->maxSize )) {
          $C = $C + 1;                                          # increment cluster number
          printf STDERR qq{Cluster %d hash been initialized ... \n}, $C;
        #}
        $labels->{$P} = $C;                                     # add $P to cluster $C
        #printf qq{\tP Node %d added to Cluster %d ... Processing %d's Neighbors ...}, $P, $C, $P;
        my @S = @$N;                                            # use $P's neighbors ($N) as next candidate
      SEEDS:
        foreach my $Q (@S) {                                    # iterate over each neighbor of $P, in @N
          #printf qq{\tProcessing % 10d's Neighbor Node % 10d\n}, $P, $Q;
          if ($labels->{$Q} and $labels->{$Q} eq q{noise}) {    # if previously marked as "noise" update $Q to be in $C
            #printf qq{\t\tQ Node %d has been relabeled form 'noise' to be in Cluster %d ...\n}, $Q, $C;
            $labels->{$Q} = $C;
            push @{$Clusters->{$C}}, $Q;                        # capture nodes by $C
          }
          next SEEDS if defined $labels->{$Q};                  # continue to next $Q if label is now defined (previously processed)
          $labels->{$Q} = $C;                                   # set undefined label for $Q to be $C (current cluster)  
          push @{$Clusters->{$C}}, $Q;                          # capture nodes by $C
          #printf STDERR qq{\t\tQ Node %d has been added to Cluster %d ...\n}, $Q, $C;
          my $N = rangeQuery( $NODES, $distFunc, $Q, $eps );    # get neighbors of $Q
          
          if (@$N and @$N >= $minPts ) {
            #printf qq{\t\tQ Node %d has %d neighbors ...\n}, $P, scalar @$N;
            push @S, @$N;                                       # if @N >= min number of nodes, push it onto @S for continued processing
          }
        }
    }
    return $Clusters;
}

sub rangeQuery {
    my ( $NODES, $distFunc, $A, $epsilon ) = @_;
    my $_k = omp_rangeQuery($NODES, $NODES->{$A}, $epsilon);
    return $_k;
}

sub distance {
    my ( $x, $y ) = @_;
    return sqrt( $x**2 + $y**2 );
}

process_f14($o);

__DATA__
__C__

#define _POMP_ENV_UPDATE_NUM_THREADS_ char *num = getenv("OMP_NUM_THREADS"); omp_set_num_threads(atoi(num));

AV* omp_rangeQuery(SV* NODES_ref, SV* A, float eps) {

  /*
   * set up return array ref
  */
  _POMP_ENV_UPDATE_NUM_THREADS_

  AV* ret = newAV();
  sv_2mortal((SV*)ret);

  HV *NODES_hash;
  HE *NODES_hash_entry;
  HV *A_hash;
  int num_keys, i;
  int a = 0;

  SV **NODE_ref;
  HV *NODE_hash;
 
  NODES_hash = (HV*)SvRV(NODES_ref);
  num_keys = hv_iterinit(NODES_hash);

  A_hash = (HV*)SvRV(A);
  SV **node_x = hv_fetch(A_hash, "x", 1, FALSE);
  SV **node_y = hv_fetch(A_hash, "y", 1, FALSE);
  float A_X = SvNV(*node_x);
  float A_Y = SvNV(*node_y);

  int accumulator[num_keys];
  int sum = 0;
  #pragma omp parallel reduce(+:sum)
  {
    sum = omp_get_thread_num();

    #pragma omp for
    for (i = 0; i < num_keys; i++) {
      int length = floor(log10(abs(i+1))) + 1;
      char key[length]; 
      sprintf(key,"%d",i+1);
      NODE_ref = hv_fetch(NODES_hash, key, length, FALSE);
      NODE_hash = (HV*)SvRV(*NODE_ref);
      SV **node_x = hv_fetch(NODE_hash, "x", 1, FALSE); // hv_fetch returns a SV** -see deref below
      SV **node_y = hv_fetch(NODE_hash, "y", 1, FALSE); // hv_fetch returns a SV** -see deref below
      float X = SvNV(*node_x);
      float Y = SvNV(*node_y);
  
      float dist = sqrt(pow(A_X-X,2)+pow(A_Y-Y,2));
      //printf("distance between (%f,%f) and (%f,%f) is %f ...\n", A_X, A_Y, X, Y, dist);
  
      // add to arrayref that gets returned
      if (dist <= eps) {
        accumulator[i] =  1;
      }
    }
    #pragma omp single
    for (int i = 0; i < num_keys; i++) {
      int tid = omp_get_thread_num();
      if (accumulator[i] == 1) { 
        av_push(ret, newSViv(i+1));
      }
    }
   
  }
  return ret;
}
