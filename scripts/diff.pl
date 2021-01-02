#!/usr/bin/perl

use strict;
use warnings;
use utf8;

my ($file1, $file2, $upto, $eps) = @ARGV;
defined($file2) or die "USAGE: file1 file2 [upto eps=1e-5]\n";
defined($eps) or $eps = 1e-5;
my @a1 = loadFile($file1);
my @a2 = loadFile($file2);

diffOf(\@a1, \@a2);

sub loadFile
{
	my ($file) = @_;
	my @a;
	open(FILE, "<", "$file") or die "$0: Cannot open file $file : $!\n";
	my $counter = 0;
	while (<FILE>) {
		my @temp = split;
		$a[$counter++] = \@temp;
	}

	close(FILE);
	return @a;
}

sub diffOf
{
	my ($a1, $a2) = @_;
	my $n1 = scalar(@$a1);
	my $n2 = scalar(@$a2);
	($n1 == $n2) or die "$0: Number of rows different\n";

	my @a;
	for (my $i = 1; $i < $n1; ++$i) {
		my $ptr1 = $a1->[$i];
		my $ptr2 = $a2->[$i];
		my $m1 = scalar(@$ptr1);
		my $m2 = scalar(@$ptr2);
		my $m = ($m1 < $m2) ? $m1 : $m2;
		$m = $upto if (defined($upto));
		#($m1 == $m2) or die "$0: Number of cols different in row $i\n";
		for (my $j = 0; $j < $m; ++$j) {
			my $val = abs($ptr1->[$j]) - abs($ptr2->[$j]);
			$val = 0 if ($val < $eps);

			print "$val ";
		}

		print "\n";
	}
}

