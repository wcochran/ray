#!/usr/bin/perl

use strict;
use Math::Trig;

my $SLICES = 128;
my $DICES = 128;

my @verts;

my $dphi = pi/$SLICES;
my $dtheta = 2*pi/$DICES;

push @verts, (0, 0, 1);

for (my $j = 1; $j < $DICES; $j++) {
	my $phi = $j*$dphi;
	for (my $i = 0; $i < $SLICES; $i++) {
		my $theta = $i*$dtheta;
		my $x = cos($theta) * sin($phi);
		my $y = sin($theta) * sin($phi);
		my $z = cos($phi);
		push @verts, ($x, $y, $z);
	}
}

push @verts, (0, 0, -1);

my $N = @verts / 3;

open(my $fout, '>:raw', "sphere.ply") or die "opening for write: $!\n";

print $fout <<"EOH";
ply
format binary_little_endian 1.0
comment sphere.pl generated
element vertex $N
property float x
property float y
property float z
property uchar red
property uchar green
property uchar blue
property uchar alpha
property float nx
property float ny
property float nz
element face 0
property list uchar int vertex_indices
end_header
EOH

for (my $i = 0; $i < $N; $i++) {
   my $idx = 3*$i;
   my $x = $verts[$idx++];
   my $y = $verts[$idx++];
   my $z = $verts[$idx++];
   print $fout pack("f<", $x);
   print $fout pack("f<", $y);
   print $fout pack("f<", $z);
   print $fout pack("C", 255);
   print $fout pack("C", 255);
   print $fout pack("C", 255);
   print $fout pack("C", 255);
   print $fout pack("f<", $x);
   print $fout pack("f<", $y);
   print $fout pack("f<", $z);
}

close $fout;



