#!/usr/bin/perl

my %depth1;
my %depth5;
my %depth20;
my %totalcount;
#my $count=0;

while (<>) {
#$count++;
	my @line=split(/\s/,$_);
        my $depth=$line[2];
        my $chrom=$line[0];
        my $chromtype="UNK";
#       if ($chrom=~/WKLF/) {$chrom="nDi";}
#       if ($chrom=~/NJBR02|HQ/) {next;}
#       if ($chrom=~/NJBR02/) {$chrom="wWb";}
#       if ($chrom=~/OM/) {$chromtype="nuclear";}
        if ($chrom=~/OO/) {$chromtype="contig";}
        if ($chrom=~/OM1/) {$chromtype="autosomal";}
        if ($chrom=~/OM3/) {$chromtype="autosomal";}
        if ($chrom=~/OM4/) {$chromtype="autosomal";}
        if ($chrom=~/OM2/) {$chromtype="se";}
#       $depth1{$chrom}=$depth1{$chrom}+$depth;
#       $depth1{$chromtype}=$depth1{$chromtype}+$depth;
        if ($depth >=1) { $depth1{$chrom}++; $depth1{$chromtype}++; }
        if ($depth >= 5) { $depth5{$chrom}++; $depth5{$chromtype}++; }
        if ($depth >= 20) { $depth20{$chrom}++; $depth20{$chromtype}++; }
        $totalcount{$chrom}++;
        $totalcount{$chromtype}++;
#print "$count\t$totalcount{nWb}\n";
#print "\t$depth5{$chrom}\t$depth20{$chrom}\t$totalcount{$chrom}\t";

}

print "\n\nChromosome\tNum_bp_genome\tNum_bp_covered\tAverage\t%Cov_1\t%_Cov_5\t%_Cov_20";

while (my ($k,$c) = each(%totalcount)) {
        my $percent5=100*$depth5{$k}/$c;
        my $percent20=100*$depth20{$k}/$c;
        my $percent1=100*$depth1{$k}/$c;
        my $averagedepth=$depth1{$k}/$c;
        print "\n$k\t$c\t$depth1{$k}\t$averagedepth\t$percent1\t$percent5\t$percent20";
}

print "\n\n";
#print "om2 $depth1{se} \n";
my $se="se";
my $au="autosomal";
#print "om2 $depth1{$se}\n";
my $ratio=$depth1{$se}/$depth1{$au};
print "The ratio of OM2 to autosomes is $ratio\n\n";
end;