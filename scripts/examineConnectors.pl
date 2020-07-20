#!/usr/bin/perl

use strict;
use warnings;
use utf8;

my ($file) = @ARGV;

defined($file) or die "USAGE: $0 filename\n";

examineFile($file);

sub examineFile
{
	my ($file) = @_;
	open(FILE, "<", "$file") or die "$0: Cannot open $file : $!\n";
	my $row = 0;
	while (<FILE>) {
		my @temp = split;
		my $n = scalar(@temp);
		if ($n < 5) {
			$row = 0;
			print;
			next;
		}

		my $ret = procLine(\@temp, $row++);
		print "\n" if ($ret);
	}

	close(FILE);
}

sub procLine
{
	my ($array, $row) = @_;
	my $n = scalar(@$array);
	my $ret = 0;
	for (my $i = 0; $i < $n; ++$i) {
		my $val = $array->[$i];
		if ($val ne "(0,0)") {
			print "$row,$i,$val ";
			$ret = 1;
		}
	}

	return $ret;
}
