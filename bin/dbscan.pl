#!/usr/bin/env perl

use warnings;
use strict;
use Util::H2O::More qw/ddd Getopt2h2o/;

my $P           = 5;
my $V           = undef;
my $ABSOLUTEMIN = 1;

my $o = Getopt2h2o \@ARGV, {
    eps     => 0.1,
    f       => q{./fort.14},
    minPts  => 5,
    minSize => undef,
  },
  qw/eps=s f=s minPts=i/;

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

    $o->minSize(4000) if not $o->minSize();

    # collect node positions
    for ( my $i = 1; $i <= $NN; ++$i ) {
        my $line = <$F14>;
        $line =~ m/^\s*(\d+)\s+([+\-]*\d+\.\d+)\s+([+\-]*\d+\.\d+)/;
        my $node = $1;
        my $x    = $2;
        my $y    = $3;
        $NODES->{$node} = { x => $x, y => $y };
    }

#    # collect element members
#    for ( my $i = 1; $i <= $NE; ++$i ) {
#        my $line = <$F14>;
#        $line =~ m/^\s*(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)/;
#        my $element = $1;
#        $ELEMENTS->{$element} = { nodes => [ $3, $4, $5 ] };
#        foreach my $n ( $3, $4, $5 ) {
#            $NODES->{$n}->{e}->{$element} = undef;
#            map { $NODES->{$n}->{k}->{$_} = undef } ( $3, $4, $5 );
#        }
#    }
    close $F14;
#
# DBSCAN!! ... HDBSCAN!! ??
#
    my $S = DBSCAN( $NODES, \&distance, $o->eps, $o->minPts );
    ddd $S;
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
        printf qq{\tP Node %d has %d neighbors ...\n}, $P, scalar @$N if @$N;
        if ( @$N < $minPts ) {                                  # density check 
            printf qq{\tMarking P Node %d as "noise" ...\n}, $P;
            $labels->{$P} = q{noise};                           # mark $P as "noise"
            next NODE;                                          # continue to next $P, not a "core" point
        }
        if ( $C == 0 or ( $Clusters->{$C} and @{$Clusters->{$C}} >= $o->minSize )) {
print $o->minSize;
ddd $Clusters if $C > 0;
die if $C > 0;
          $C = $C + 1;                                          # increment cluster number
        }
        $labels->{$P} = $C;                                     # add $P to cluster $C
        #printf qq{\tP Node %d added to Cluster %d ... Processing %d's Neighbors ...}, $P, $C, $P;
        my @S = @$N;                                            # use $P's neighbors ($N) as next candidate
      SEEDS:
        foreach my $Q (@S) {                                    # iterate over each neighbor of $P, in @N
          #printf qq{\tProcessing % 10d's Neighbor Node % 10d\n}, $P, $Q;
          if ($labels->{$Q} and $labels->{$Q} eq q{noise}) {    # if previously marked as "noise" update $Q to be in $C
            printf qq{\t\tQ Node %d has been relabeled form 'noise' to be in Cluster %d ...\n}, $Q, $C;
            $labels->{$Q} = $C;
            push @{$Clusters->{$C}}, $Q;                        # capture nodes by $C
          }
          next SEEDS if defined $labels->{$Q};                  # continue to next $Q if label is now defined (previously processed)
          $labels->{$Q} = $C;                                   # set undefined label for $Q to be $C (current cluster)  
          push @{$Clusters->{$C}}, $Q;                          # capture nodes by $C
          printf qq{\t\tQ Node %d has been added to Cluster %d ...\n}, $Q, $C;
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
    my $k = [];
  NODE:
    foreach my $B ( keys %$NODES ) {
        next NODE if $A == $B;
        my $x = $NODES->{$A}->{x} - $NODES->{$B}->{x};
        my $y = $NODES->{$A}->{y} - $NODES->{$B}->{y};
        my $dist = $distFunc->( $x, $y ); 
        if ( $distFunc->( $x, $y ) <= $epsilon ) {
            push @$k, $B;
        }
    }
    return $k;
}

sub distance {
    my ( $x, $y ) = @_;
    return sqrt( $x**2 + $y**2 );
}

process_f14($o);
