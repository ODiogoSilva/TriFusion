#!/usr/bin/perl

my ($prefix, $startId) = @ARGV;

&usage unless ($prefix && ($startId =~ /\d+/));


while (<STDIN>) {
  s/\t/ /g;
  print "$prefix$startId: $_";
  $startId++;
}

sub usage {
  print STDERR "
mclOutput2groupsFile prefix starting_id_num

create an orthomcl groups file from an mcl output file. just generate a group ID for each group, and prepend it to that group's line.

where:
 prefix           a prefix to use when generating group ids.  For example OG2_
 starting_id_num  a number to start the id generating with.  For example 1000

std input:  mcl output file (label mode)
std output: orthomcl groups file

an orthomcl group file has one line per group and looks like this:

OG2_1009: osa|ENS1222992 pfa|PF11_0844


";
  exit(1);
}
