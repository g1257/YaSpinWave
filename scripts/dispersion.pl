#!/usr/bin/perl
use warnings;
use strict;

my ($mfile) = @ARGV;
defined($mfile) or die "USAGE: $0 mfile\n";
my $total = 100;
my $b = 0.5;
my $largestX = -100;

for (my $i = 0; $i < $total; ++$i) {
	my $kx = -1+$i/$total;
	my $exe = "./spaceToReciprocal -m $mfile -q $kx,0";
	my $w = `$exe`;
	printResult($kx,$w);
	$largestX = $kx if ($kx > $largestX);
}

my $offset = $largestX + 1.0/$total;
for (my $i = 0; $i < $total; ++$i) {
	my $ky = $i/$total;
	my $exe = "./spaceToReciprocal -m $mfile -q 0,$ky";
	my $w = `$exe`;
	$_ = $ky * $b + $offset;
	printResult($_,$w);
	$largestX = $_ if ($_ > $largestX);
}

$offset = $largestX + 1.0/$total;
for (my $i = 0; $i < $total; ++$i) {
	my $kx = $i/$total;
	my $ky = 0;
	my $exe = "./spaceToReciprocal -m $mfile -q $kx,$ky";
	my $w = `$exe`;
	$_ = $kx + $offset;
	printResult($_ ,$w);
	$largestX = $_ if ($_ > $largestX);
}

$offset = $largestX + 1.0/$total;
for (my $i = 0; $i < $total; ++$i) {
	my $kx = $i/$total;
	my $ky = $kx;
	my $exe = "./spaceToReciprocal -m $mfile -q $kx,$ky";
	my $w = `$exe`;
	$_ = $kx + $offset;
	printResult($_, $w);
	$largestX = $_ if ($_ > $largestX);
}

sub printResult
{
	my ($k,$val) = @_;
	chomp($val);
	my @val= split/ /,$val;
	return unless (scalar(@val) == 2);
	print "$k $val[0] $val[1]\n";
}



