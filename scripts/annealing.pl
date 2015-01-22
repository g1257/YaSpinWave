#!/usr/bin/perl

use strict;
use warnings;

my ($template,@betas) = @ARGV;
my $n = scalar(@betas);
(defined($template) and $n > 0) or die "USAGE: $0 template.txt beta1 beta2 ...\n";

my ($rootInput,$rootOut) = ("Input","Output");
my $prevOutput;

for (my $i = 0; $i < $n; ++$i) {
	my $input = "${rootInput}$i.inp";
	my $output = "${rootOut}$i.txt";
	open(FIN,$template) or die "$0: Cannot open $template : $!\n";
	createInput($input,$output,$betas[$i],\*FIN,$prevOutput);
	runInput($input);
	getLast("AcceptancePercentage spin",$output);
	$prevOutput = $output;
}

printFinalAngles($prevOutput);

sub createInput
{
	my ($input,$output,$beta,$fh,$prevOutput) = @_;
	my $random = int(10000*rand());
	open(FOUT,">$input") or die "$0: Open > failed for $input : $!\n";

	while (<$fh>) {
		next if (/^#/);
		if (/\$([a-zA-Z0-9\[\]]+)/) {
			my $name = $1;
			my $str = "\$".$name;
			my $val = eval "$str";
			defined($val) or die "$0: Undefined substitution for $name\n";
			s/\$\Q$name/$val/g;
		}

		if (/MonteCarloStartTypeOrFile/ && defined($prevOutput)) {
			print FOUT "MonteCarloStartTypeOrFile=$prevOutput\n";
		}

		print FOUT;
	}

	close(FOUT);
}

sub runInput
{
	my ($input) = @_;
	my $cmd = "./monteCarlo -f $input &> /dev/null";
	print STDERR "$0: Executing $cmd\n";
	system($cmd);
}

sub printFinalAngles2
{
	my ($theta,$phi) = @_;
	my $n = scalar(@$theta);
	die "$0: Theta and Phi not of same size\n" unless ($n == scalar(@$phi));
	print "Angles\n";
	print "$n 2\n";
	for (my $i = 0; $i < $n; ++$i) {
		print "$theta->[$i] $phi->[$i]\n";
	}
}

sub printFinalAngles
{
	my ($file) = @_;
	my @theta;
	my @phi;
	my ($total,$total2);
	my ($flag,$flag2) = (0,0);
	open(FILE,"$file") or die "$0: Cannot open $file : $!\n";
	my $e;
	while (<FILE>) {
		next if (/^#/);
		chomp;
		if (/TotalEnergy=(.*$)/) {
			$e = $1;
			next;
		}

		if (/^Theta/) {
			$flag++;
			next;
		}

		if (/^Phi/) {
			$flag2++;
			next;
		}

		if ($flag == 1) {
			$total = $_;
			$flag++;
			next;
		}

		if ($flag2 == 1) {
			$total2 = $_;
			$flag2++;
			next;
		}

		if ($flag >= 2 and defined($total) and $flag-2 < $total) {
			$theta[$flag-2] = $_;
			$flag++;
		}

		if ($flag2 >= 2 and defined($total2) and $flag2-2 < $total2) {
			$phi[$flag2-2] = $_;
			$flag2++;
		}
	}

	close(FILE);

	print STDERR "LastEnergy=$e\n" if defined($e);
	printFinalAngles2(\@theta,\@phi);
}

sub getLast
{
	my ($label,$file) = @_;
	my $e;
	open(FILE,"$file") or die "$0: Cannot open $file : $!\n";
	while (<FILE>) {
		chomp;
		if (/$label=(.*$)/) {
			$e = $1;
			next;
		}
	}

	close(FILE);

	print STDERR "Last$label=$e\n" if defined($e);
}
