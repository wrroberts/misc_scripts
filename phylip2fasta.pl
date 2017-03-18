#1 /usr/bin/perl -w
use strict;
use warnings;

my $usage = "Usage: perl phylip2fasta.pl <inputPhylipFile> <outputFastaFile>\n";
my $infile = shift or die $usage;
my $outfile = shift or die $usage;
my $unixfile = $infile."_unix";
my $count = 0;

system ("tr '\r' '\n' <$infile >$unixfile");
open(IN, $unixfile) || die "Can't open $unixfile: $!\n";
open(OUT, ">$outfile") || die "Can't open out $outfile: $!\n";

while(my $line = <IN>) {
	chomp $line;
	next if (!$line || $line =~ /\d+\s+\d+/);
	if($line =~ /^(.*)\s+(\S+)\s*$/) {
		my $name = $1;
		my $sequence = $2;
		$name =~ s/\s+$//;
		print OUT ">$name\n$sequence\n";
		$count++;
	}
}

close IN;
close OUT;
unlink $unixfile;
print "There are a total of $count sequences\n";
