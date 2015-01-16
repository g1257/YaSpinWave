#!/usr/bin/perl

use strict;
use warnings;

use lib '../../PsimagLite/scripts';
use Make;

my @drivers = ("findAngles","anglesToMatrix","spaceToReciprocal");

my $lapack = " ";
#Make::findLapack();
Make::backupMakefile();
writeMakefile();
make();

sub make
{
	system("make");
}

sub writeMakefile
{
	open(my $fh,">Makefile") or die "Cannot open Makefile for writing: $!\n";

	my $libs = "$lapack -L../../PsimagLite/lib   -lm  -lpthread -lgsl -lgslcblas -lpsimaglite";
	my $cxx = "g++ -O3 -DNDEBUG -DNO_LAPACK -DUSE_GSL";
	my $cppflags = " -IEngine  -I../../PsimagLite/src -I../../PsimagLite";
	Make::make($fh,\@drivers,"PsimagLite","Linux",0,$libs,$cxx,$cppflags,"true"," "," ");

	close($fh);
	print "$0: Done writing Makefile\n";
}

