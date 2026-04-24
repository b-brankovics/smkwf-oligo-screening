#!/usr/bin/perl -w
use strict;
use Pod::Usage;

#===DESCRIPTION=================================================================
# A tool to extract sequence regions from a multi fasta file


#===HELP========================================================================

if (grep{/^--?h(elp)?$/} @ARGV) {
    pod2usage({-verbose=>2, -output=>\*STDOUT, -noperldoc => 1, -exitval=>0});
}

#===MAIN========================================================================

# Hash to store the sequences: key (id and defintion) and value (sequence) 
my %fas_data;
# Array to store the order of the sequences
my @ids;

my @requests;
my $tabfile;
my $exact;
my @fastas;
for (@ARGV) {
    if (/\.fas$/ || /\.f[nas]?a$/ || /\.fasta$/ || /^-$/) {
		push @fastas, $_;
    } elsif (/^--?tab=(.+)$/) {
		$tabfile = $1;
		$tabfile =~ s/^~/$ENV{"HOME"}/;
    } else {
		push @requests, $_;
    }
}
@ARGV = @fastas;

unless (@requests || $tabfile) {
    print STDERR "ERROR: There ware no regions specified\n\n";
    pod2usage({-verbose => 2, -output => \*STDERR, -exitval => 1, -noperldoc => 1});
}


# If there were files specified then read the sequence from them
# OR read from STDIN
if (@ARGV) {
    for (@ARGV) {
		&read_fasta(\%fas_data, \@ids, $_);
    }
} else {
    &read_fasta(\%fas_data, \@ids);
}

my %tags;
if ($tabfile) {
    my $intab;
    if ($tabfile eq "-") {
		# print STDERR "Reading from STDIN\n";
		$intab = *STDIN;
    } else {
		open($intab, '<', $tabfile) || die $!;
    }
    for (<$intab>) {
		s/\R+//;
		if (/^$/) {
			print STDERR "WARNING: skipping empty lines in the tabfile (line #$.)\n";
			next;
		}
		# my ($id, $tag) = split/\t/;
		# #    print STDERR "$id -> $tag\n";
		# push @requests, $id;
		# $tags{$id} = $tag;
		push @requests, $_;
    }
    close $intab;
}

# Go through the requests and extract the target regions
# Print the sequences to STDOUT
my $output;
# Store the output in a string and print it only when no errors were found
for my $req (@requests) {
	my ($pos_string, $tag) = split/\t/, $req;
	my ($seq_id, $pos) = split/:/, $pos_string;
	my $rev;
	if ($pos =~ /^complement\(([\d,.]+)\)$/) {
	    $rev = 1;
	    $pos = $1;
	}
	if ($pos !~ /^[\d,.]+$/) {
	    die "ERROR: invalid position string '$pos' in request '$req'\n";
	}
	unless ($fas_data{$seq_id}) {
	    die "ERROR: sequence id '$seq_id' in request '$req' not found in the fasta data\n";
	}
	# Extract the region(s)
	my $chunk = "";
	for my $p (split/,/, $pos) {
	    my ($start, $end) = split/\.\./, $p;
		$chunk .= substr($fas_data{$seq_id}, $start-1, $end-$start+1);
	}
	$chunk = &reverse_seq($chunk) if $rev;
	my $seq_name = $pos_string;
	if ($tag) {
	    $seq_name .= " " . $tag;
	}
	$output .= &to_fasta($seq_name, $chunk);
}

print $output;

#===SUBROUTINES=================================================================

sub to_fasta {
    # Return a fasta formated string
    my ($seq_name, $seq, $len) = @_;
    # default to 60 characters of sequence per line
    $len = 60 unless $len;
    my $formatted_seq = ">$seq_name\n";
    while (my $chunk = substr($seq, 0, $len, "")) {
	$formatted_seq .= "$chunk\n";
    }
    return $formatted_seq;
}

sub read_fasta {
    # This loads the sequences into a hash and an array
    my ($hash, $list, $file) = @_;
    my $in;
    $file = undef if ($file && $file eq "-");
    if ($file) {
	open $in, '<', $file || die $!;
    } else {
	$in = *STDIN;
    }
    if ($file) {
	unless (-e $file) {
	    die "ERROR: no such file '$file'\n";
	}
	if (-z $file) {
	    print STDERR "WARNING: '$file' is empty\n";
	}
    }
    my $seq_id;
    for (<$in>) {
	# Skip empty lines
	next if /^\s*$/;
	s/\R+//g;
	# Check wheter it is an id line
	if (/>(\S+)/) {
	    # Save the id and the definition and store it in the array
	    $seq_id = $1;
	    print STDERR "WARNING: <" . $seq_id . "> is present in multiple copies\n" if $hash->{$seq_id};
	    push @$list, $seq_id;
	} else {
	    # If there was no id lines before this then throw an error
	    unless (defined $seq_id) {
		print "Format error! Check the file!\n";
		last;
	    }
	    # Remove lineendings and white space
	    s/\R//g;
	    s/\s+//g;
	    # Add to the sequence
	    $hash->{$seq_id} .= $_;
	}
    }
    close $in;
}

sub print_fasta {
    # Print all the sequences to STDOUT
    my ($hash, $list, $tag) = @_;
    for (@$list) {
	# /^(\S+)/;
	my $id = $_;
	for my $regex (keys %$tag) {
	    next unless $tag->{$regex};
	    if ($id =~ /$regex/) {
		# my $t = $tag->{$regex};
		# $id =~ s/($regex)/$1 $t/;
		$id .= " " . $tag->{$regex};
	    }
	}
	# if ($tag->{$id}) {
	# #    print STDERR "Tag: $_ -> " . $tag->{$id} . "\n";
	#     $id = $_ . " " . $tag->{$id};
	# } else {
	#     $id = $_;
	# }
	print &to_fasta($id, $hash->{$_});
    }
}

sub reverse_seq {
    # Reverse complements the sequences
    my ($seq) = @_;
    # Reverse the sequnce
    my $complement = reverse $seq;
    # Complement the sequence
    $complement =~ tr/ACGTacgtWwMmRrSsKkYyBbVvDdHh/TGCAtgcaWwKkYySsMmRrVvBbHhDd/;
    return $complement;
}

##############END OF SCRIPT#########

__END__
# POD: Plain Old Documentation block


=pod

=head3 fasta_grep

This is a script to extract sequence regions from a multi FASTA file.

=head4 Synopsis

    fasta_get_regions [-tab=<tabular file>] [<input.fas>]...

=head4 Options

=over 12

=item B<-h> | B<--help>

Print the help message; ignore other arguments.

=item B<-tab=FILE> | B<--tab=FILE>

Reads the region position values from the file. It splits the lines by tabs and
uses the first part as NCBI style possition (SeqID:start..end) and the second one
to add as tag when printing the output.

=back

=head4 Input

FASTA format file containing sequences that will be filtered.
There can be multiple sequence files specified and the script can also read from STDIN.
If no FASTA files are specified then it will automatically read from STDIN.
Also if B<-> is included as argument, then STDIN will be read.

Files with the following extensions are recognized as FASTA:
'*.fas' | '*.fa' | '*.fsa' | '*.fasta'

=head4 Output

FASTA format file is printed to STDOUT after filtering.

=cut
