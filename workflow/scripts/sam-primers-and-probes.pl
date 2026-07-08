#!/usr/bin/env perl

use warnings;
use strict;

use FindBin;                     # locate this script
use lib "$FindBin::RealBin/lib";  # use the lib directory
use biointsam;
use biointbasics;
use YAML 'LoadFile';
use Data::Dumper;

#===DESCRIPTION=================================================================

my $description = 
    "Description:\n\t" .
    "A tool to find PCR amplicons from a SAM file and YAML PCR test specification file.\n";
my $usage = 
    "Usage:\n\t$0 [OPTIONS]  <SAM file> <YAML file>\n";
my $options = 
    "Options:\n" .
    "\t-h | --help\n\t\tPrint the help message; ignore other arguments.\n" .
    "\n";
my $info = {
    description => $description,
    usage => $usage,
    options => $options,
};

#===MAIN========================================================================

my %ref;
my %query;

# my @keep;
# for (@ARGV) {
# }
# @ARGV = @keep;

biointbasics::print_help(\@ARGV, $info, "ERROR: This program requires exactly two arguments.\n\n") unless @ARGV == 2;

my ($sam, $confile) = @ARGV;

my $debug = 0;

# Print help if needed
biointbasics::print_help(\@ARGV, $info);


my $config = LoadFile($confile);
# Store test info and mapping of IDs to tests
my %info = %$config;
for my $test (keys %$config) {
    for my $type (keys %{ $config->{$test} }) {
		my $id = $config->{$test}->{$type};
		$info{$id} = {"test", $test, "type", $type};
		my $pair;
		if ($type eq 'F') {
			$pair = $config->{$test}->{'R'};
		} elsif ($type eq 'R') {
				$pair = $config->{$test}->{'F'};
		}
		if ($pair) {
			$info{$id}->{'pair'} = $pair;
		}
    }
}

print STDERR Dumper(\%info) if $debug;

# Remember read sequence until QNAME changes
my $qname = "";
my $mem = "";
my $memflag = 0;
my $in;
if ($sam eq '-') {
    $in = *STDIN;
} else {
    open($in, '<', $sam) || die $!;
}
my %map;
my @order;
print STDERR "Start looking for hits\n" if $debug;
while(<$in>) {
    my %hit;
    biointsam::parse_sam($_, \%ref, \%hit);
    if (%hit) {
		unless ($hit{'RNAME'} eq "*") {
			$hit{'aln'} = biointsam::parse_cigar($hit{'CIGAR'}, $hit{'FLAG'}, $hit{'SEQ'});
			push @order, $hit{'QNAME'} unless $map{ $hit{'QNAME'} };
			push @{ $map{ $hit{'QNAME'} } }, \%hit; 
		}
    }
}

print join("\t", "test", "loc", "inner amplicon loc", "amplicon size", "total missmatches", "individual missmatches"), "\n";

for my $q (@order) {
    my @hits = sort{ $a->{'aln'}->{'start'} <=> $b->{'aln'}->{'start'} } @{ $map{$q} };
    my @worked;
    while (@hits > 1) {
		my $a = shift @hits;
		# need positive forward or a positive reverse
		# 1. check orientation
		# 2. check reference and config
		unless ($a->{'FLAG'} & 16) {
			my $id = $a->{'RNAME'};
			next unless $info{$id}; # Skip if not in specified tests
			my $type = $info{$id}->{'type'};
			my $test = $info{$id}->{'test'};
			unless ($type eq 'P') {
				my $stop = $info{$id}->{'pair'};
				print STDERR "Found start on '$q' for $test $type and looking for $stop\n" if $debug;
				# need a probe if the set has one
				my $probe = $info{$test}->{'P'};
				my $signal;
				# Here we assume that the hits are sorted based on position
				for (@hits) {
					$signal = $_ if ($probe && $_->{'RNAME'} eq $probe);
					# need the matching end (negative reverse or negative forward)
					if ($_->{'RNAME'} eq $stop && $_->{'FLAG'} & 16) {
						my $b = $_;
						print STDERR "\tFound stop\n" if $debug;
						if ( ($probe && $signal) || ! $probe) {
							print STDERR "\tit is working\n" if $debug;
							# Save all three hits
							push @worked, [$a, $signal, $b];
						}
						last;
					}
				}
			}
		}
    }
    for (@worked) {
		# three hit objects
		my ($a, $signal, $b) = @$_;
		# Get Test name
		my $test = $info{ $a->{'RNAME'} }->{'test'};
		my $ori = "+";
		if ($info{ $a->{'RNAME'} }->{'type'} eq 'R') {
			$ori = "-";
		}
		my $start = $a->{'aln'}->{'start'};

		# Get region
		my $id = $a->{'QNAME'};
		my $aln = $b->{'aln'};
		my $end = $aln->{'start'} + $aln->{'length'} - $aln->{'deletion'} - 1;
		my $region = $start . ".." . $end;
		my $amplicon = "" . ($start + $a->{'aln'}->{'length'} - $a->{'aln'}->{'deletion'}) . ".." . ($aln->{'start'} - 1) . "";
		if ($ori eq "-") {
			$region = "complement($region)";
			$amplicon = "complement($amplicon)";
		}
		# Get mismatches
		my ($sum, @counts) = &get_mismatches($_, \%ref);
		# my $sum = 0;
		# for (@counts) {
		#     $sum += $_;
		# }
		my $len = $end - $start + 1;
		print join("\t", $test, $id . ":" . $region, $id . ":" . $amplicon, $len, $sum, join(";", @counts)), "\n";
    }
}

sub get_mismatches {
    my ($list, $ref) = @_;
    my @counts;
    my $sum = 0;
    for (@$list) {
		next unless $_;
		my $name = $_->{'RNAME'};
		my $aln = $_->{'aln'};
		my $len = $ref->{$name};
		my $covered = $aln->{'length'} - $aln->{'insertion'};
		my $miss = $_->{'NM:i'} + $len - $covered;
		push @counts, $name . ":" . $miss;
		$sum += $miss;
    }
    return $sum, @counts;
}
