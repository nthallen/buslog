#! /usr/bin/perl -w
use strict;

my %stop;
my @stop;
my $stop_num = 1;
my $dir_num = 1;
my $path_num = 1;
my $route = "route";

while (<>) {
  if ( m/^<route tag="([^"]+)" .* latMin="(-?[\d.]+)" latMax="(-?[\d.]+)" lonMin="(-?[\d.]+)" lonMax="(-?[\d.]+)">/ ) {
    $route = "route_$1";
    print "$route.latMin = $2;\n$route.latMax = $3;\n$route.lonMin = $4;\n$route.lonMax = $5;\n";
  } elsif ( m|^<stop tag="([^"]+)" title="([^"]+)" lat="(-?[\d.]+)" lon="(-?[\d.]+)"(?: stopId="[^"]+")?/>| ) {
    die "Stop $1 redefined\n" if $stop{$1};
    $stop{$1} = { title => $2, lat => $3, lon => $4, idx => $stop_num++ };
    push @stop, $stop{$1};
  } elsif ( m|<direction tag="([^"]+)" title="([^"]+)" name="([^"]+)" useForUI="([^"]+)">| ) {
    print "$route.direction($dir_num) = struct('tag','$1','title','$2','name','$3','stops', [\n";
    while (<>) {
      if ( m|^\s*<stop tag="([^"]+)" />| ) {
        print "  $stop{$1}->{idx} $stop{$1}->{lat} $stop{$1}->{lon}\n";
      } elsif ( m|^</direction>| ) {
        print "  ] );\n";
        last;
      }
    }
    ++$dir_num;
  } elsif ( m|^<path>| ) {
    print "$route.path($path_num).points = [\n";
    while (<>) {
      if ( m|^<point lat="(-?[\d.]+)" lon="(-?[\d.]+)"/>| ) {
        print "  $1 $2\n";
      } elsif ( m|^</path>| ) {
        print "];\n";
        ++$path_num;
        last;
      }
    }
  }
}
print "$route.stops = [\n";
foreach my $i ( 0 .. $#stop ) {
  print "  struct('title', '$stop[$i]->{title}', 'lat', $stop[$i]->{lat}, 'lon', $stop[$i]->{lon} )\n";
}
print "];\n";

