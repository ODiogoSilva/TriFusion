#!/usr/bin/perl

use strict;

my $blastFile = shift(@ARGV);
my $fastaFilesDir = shift(@ARGV);

usage() unless $blastFile;

usage() unless $fastaFilesDir;
 
opendir(DIR, $fastaFilesDir) || die "Can't open fasta directory '$fastaFilesDir'\n";
my @fastaFiles = readdir(DIR);
closedir(DIR);

my $genes = getGenesFromFasta($fastaFilesDir, @fastaFiles);

open(F,$blastFile) || die "can't open BLAST file '$blastFile'\n";

=pod
query_name, hitname, 
$pcid, len, 
mismatches, ngaps, 
start('query'), end('query'),
start('hit'), end('hit'),
evalue, bits
=cut

my $prevSubjectId = 'blah';
my $prevQueryId = 'blah';
my $subject;  # hash to hold subject info
my $queryShorter;

while(<F>) {
    chomp;
    my ($queryId, $subjectId, $percentIdentity, $length, $mismatches, $ngaps, $queryStart, $queryEnd, $subjectStart, $subjectEnd, $evalue, $bits) = split;

    if ($queryId ne $prevQueryId || $subjectId ne $prevSubjectId) {

	# print previous subject
	printPreviousSubject($subject) if $subject;

	# initialize new one from first HSP
	$prevSubjectId = $subjectId;
        $prevQueryId = $queryId;

	$subject = {}; 
	$subject->{queryId} = $queryId;
	$subject->{subjectId} = $subjectId;
	$subject->{queryShorter} = getTaxonAndLength($subject, $genes);
	
	($subject->{evalueMant}, $subject->{evalueExp})
	    = formatEvalue($evalue); # from first hsp
    }

    # get additional info from subsequent HSPs
    my $hspspan = [$subjectStart, $subjectEnd];
    $hspspan = [$queryStart, $queryEnd] if $subject->{queryShorter};
    push(@{$subject->{hspspans}}, $hspspan);
    $subject->{totalIdentities} += $percentIdentity * $length;
    $subject->{totalLength} += $length;
}
printPreviousSubject($subject);

########################################################################################

sub getGenesFromFasta {
    my $fastaFilesDir = shift(@_);
    my (@fastaFiles) = @_;

    my $genes;
    foreach my $fastaFile (@fastaFiles) {
	next if $fastaFile =~ /^\./;
	print STDERR "acquiring genes from $fastaFile\n";
	$fastaFile =~ /(\w+).fasta/ || die "'$fastaFile' is not in 'taxon.fasta' format\n";
	my $taxon = $1;
	open(FF,"$fastaFilesDir/$fastaFile") || die "can't open fasta file '$fastaFilesDir/$fastaFile'";
	my $gene;
	my $length;
	while (<FF>) {
	    chomp;
	    next if /^\s*$/;
	    if (/\>(\S+)/) {
		$genes->{$gene}->{length} = $length if $gene;
		$genes->{$gene}->{taxon} = $taxon if $gene;
		$gene = $1;
		$length = 0;
	    } else {
		$length += length($_);
	    }
	}
	$genes->{$gene}->{length} = $length if $gene;
	$genes->{$gene}->{taxon} = $taxon if $gene;
	close(FF);
    }
    return $genes;
}

sub getTaxonAndLength {
    my ($subject, $genes) = @_;
    $subject->{queryTaxon} = $genes->{$subject->{queryId}}->{taxon};
    $subject->{subjectTaxon} = $genes->{$subject->{subjectId}}->{taxon};
    $subject->{queryLength} = $genes->{$subject->{queryId}}->{length};
    $subject->{subjectLength} = $genes->{$subject->{subjectId}}->{length};
    die "couldn't find taxon for gene '$subject->{subjectId}'" unless $subject->{subjectTaxon};
    die "couldn't find taxon for gene '$subject->{queryId}'" unless $subject->{queryTaxon};
    return $subject->{queryLength} < $subject->{subjectLength};
}

sub printPreviousSubject {
    my ($subject) = @_;

    my $nonOverlapMatchLen = computeNonOverlappingMatchLength($subject);

    my $percentIdent =
	int($subject->{totalIdentities} / $subject->{totalLength} * 10 + .5)/10;
    my $shorterLength = $subject->{queryShorter}? $subject->{queryLength} : $subject->{subjectLength};
    my $percentMatch = int($nonOverlapMatchLen / $shorterLength * 1000 + .5) / 10;
    print "$subject->{queryId}\t$subject->{subjectId}\t$subject->{queryTaxon}\t$subject->{subjectTaxon}\t$subject->{evalueMant}\t$subject->{evalueExp}\t$percentIdent\t$percentMatch\n";
}

# this (corrected) version of formatEvalue provided by Robson de Souza
sub formatEvalue {
    my ($evalue) = @_;
    $evalue = '1' . $evalue if ($evalue =~ /^e/);
    $evalue = sprintf("%.3e",$evalue);
    my ($evalue_mant, $evalue_exp) = split(/e/, $evalue);
    $evalue_mant = sprintf("%.2f",$evalue_mant);
    $evalue_mant =~ s/\.0+$//;
    $evalue_exp =~ s/\+//;
    $evalue_exp = 0 if ($evalue_exp eq '00');
    return ($evalue_mant, $evalue_exp);
}

sub computeNonOverlappingMatchLength {
    my ($subject) = @_;

    my @hsps = sort {$a->[0] <=> $b->[0]} @{$subject->{hspspans}};
    my $first = shift @hsps;
    return 0 unless $first;
    my ($start, $end) = getStartEnd($first);
    my $len = 0;
    foreach my $h (@hsps){
	my ($hspStart,$hspEnd) = getStartEnd($h);
	
	next if $hspEnd <= $end; ##does not extend
	if ($hspStart <= $end) {  ##overlaps
	    $end = $hspEnd;  #extend end ... already dealt with if new end is less
	} else {  ##there is a gap in between ..
	    $len += $end - $start + 1;
	    $start = $hspStart;
	    $end = $hspEnd;
	}
    }
    $len += $end - $start + 1; # deal with the last one 
    return $len
}

#  flip orientation if nec.
sub getStartEnd {
    my ($h) = @_;
    my $hspStart = $h->[0];
    my $hspEnd = $h->[1];
    if ($hspStart > $hspEnd) {
	$hspEnd = $h->[0];
	$hspStart = $h->[1];
    }
    return($hspStart,$hspEnd);
}

sub usage {
  print STDERR "

orthomclBlastParser blast_file fasta_files_dir

where:
  blast_file:       BLAST output in m8 format.
  fasta_files_dir:  a directory of compliant fasta files as produced by
                    orthomclAdjustFasta 

  
m8 format has these columns:
  query_name, hitname, pcid, len, mismatches, ngaps, start('query'), 
  end('query'), start('hit'), end('hit'), evalue, bits

output:
  tab delimited text file, with one row per query-subject match. the columns are:
     query_id, subject_id, query_taxon, subject_taxon,
     evalue_mant, evalue_exp, percent_ident, percent_match

(percent_match is computed by counting the number of bases or amino acids in the shorter sequence that are matched in any hsp, and dividing by the length of that shorter sequence)

EXAMPLE: orthomclSoftware/bin/orthomclBlastParser my_blast_results my_orthomcl_dir/compliantFasta >> my_orthomcl_dir/similarSequences.txt


";


  exit(1);
}
