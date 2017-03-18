#!/usr/bin/perl
#contig_stats.pl
#AUTHOR: Wade R. Roberts (modified 27 Feb 2016)
#
#This program calculates and outputs the following statistics for a de novo transcriptome assembly.
#	1. number of entries
#	2. max length
#	3. min length
#	4. mean length
#	5. total bases
#	6. length quantile1
#	7. median length
#	8. length quantile3
#	9. N50
#	10. N90
#	11. N50 > 500
#	12. N90 > 500
#	13. contigs > 500
#	14. bases in contigs > 500
#	15. dup-mer-21
#	16. dup-mer-count
#	17. dup-pct-21

$max_mers = 40000000;

if (@ARGV > 1) {
for (@ARGV) {
if ( (stat($_))[9] > (stat("$_.stats"))[9] ) {
if (!fork) {
system("grun '$0 $_ > $_.stats'");
exit(0);
}
}
}
while (wait != -1) {};
exit(0);
}

$in = $ARGV[0];

open(IN, $in) || die "$in:$!\n";

my $seq;
while (<IN>) {
if (/^>/) {
summ();
$ln=0;
$seq='';
next;
} else {
chomp;
$seq .= $_;
$ln += length($_);
}
}

summ() if $ln;

print "entries\t$gc\n";
print "len max\t$mx\n";
print "len min\t$mn\n";
printf "len mean\t%d\n", $tl/$gc;
print "total bases\t$tl\n";

@ln = sort {$a <=> $b} @ln;

printf "len q1\t%d\n", quantile(\@ln, .25);
printf "len median\t%d\n", quantile(\@ln, .50);
printf "len q3\t%d\n", quantile(\@ln, .75);

@gt500 = grep {$_>500} @ln;
$bgt500 = 0; for (@gt500) { $bgt500 += $_ };

for (@ln) {
$sum += $_;
if (!$n50 && $sum >= $tl*.5) {
$n50=$_;
}
if (!$n90 && $sum >= $tl*.9) {
$n90=$_;
}
}

for (@ln) {
if ($_ >= 500) {
$xsum += $_;
if (!$xn50 && $xsum >= $bgt500*.5) {
$xn50=$_;
}
if (!$xn90 && $xsum >= $bgt500*.9) {
$xn90=$_;
}
}
}

printf "N50\t%d\n", $n50;
printf "N90\t%d\n", $n90;
printf "N50 > 500\t%d\n", $xn50;
printf "N90 > 500\t%d\n", $xn90;

printf "contigs > 500\t%d\n", scalar @gt500;
printf "bases in contigs > 500\t%d\n", $bgt500;

my $mer21 = 0;
my $pct21 = 0;

for (values(%mer)) {$mer21+=1 if $_>1};
for (values(%mer)) {$pct21+=$_ if $_>1};
$mer21 = 100*$mer21/scalar(keys(%mer));
$pct21 = 100*$pct21/$tl;
printf "dup-mer-21\t%2.2f\n", $mer21;
printf "dup-mer-cnt\t%d\n", $mers;
printf "dup-pct-21\t%2.2f\n", $pct21;


sub summ {
if ($ln > 0) {
++$gc;
$tl+=$ln;
$mx = $ln if $ln>$mx;
$mn = $ln if $ln<$mn || !$mn;
for ($i=0;$i<length($seq)-21;++$i) {
$tmp = substr($seq,$i,21);
if ($mers < $max_mers) {
if (!exists $mer{$tmp}) {
++$mers;
++$mer{$tmp};
} else {
++$mer{$tmp};
}
}
}
}
push @ln, $ln;
}

sub quantile {
my ($a,$p) = @_;
my $l = scalar(@{$a});
my $t = ($l-1)*$p;
my $v=$a->[int($t)];
if ($t > int($t)) {
return $v + $p * ($a->[int($t)+1] - $v);
} else {
return $v;
}
}
