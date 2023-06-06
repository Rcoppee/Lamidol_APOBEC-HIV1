#!/usr/bin/perl

use strict;
use warnings;

################################################################################

############
# PACKAGES #
############

use Cwd;
use Getopt::Long;

################################################################################

#############
# FUNCTIONS #
#############

sub usage {

	print STDERR <<EOF;

Usage:
	extract_info_from_pileup.pl -p PILEUP [-o OUTDIR] [-v]
	extract_info_from_pileup.pl [-h]

	extract_info_from_pileup.pl --pileup PILEUP [--outdir OUTDIR] [--verbose]
	extract_info_from_pileup.pl [--help]

Required arguments:
	-p | --pileup  : pileup file

Optional arguments:
	-o | --outdir  : output directory

Help:
	-v | --verbose : activate VERBOSE
	-h | --help    : print this help

Description:

	Extract information from a GFF pileup.

	1 input file is needed :
	- a pileup file

	It extracts information about each position :
	- the chromosome names,
	- the 1-based coordinates,
	- the reference base,
	- the number of reads covering the site,
	- and the read bases,
	- (base and alignment mapping qualities are not handled)

	Nucleotide, insertion(s) and deletion(s) frequencies are calculted for each position.
	Nucleotide(s) with a frequency > 20% (coverage) are kept to generate the 'IUPAC nucleotide'.
	If insertion(s) and deletion(s) frequency > 20% (coverage), the 'IUPAC nucleotide' will be in lower case, if not, it will be in UPPER case.

	1 output file is generated and it contains 7 columns :
	- the chromosome names
	- the 1-based coordinates,
	- the refenrece base,
	- the number of reads covering the site,
	- the number of A, C, G, T and indels
	- the calculated 'IUPAC nucleotide'
EOF

	exit;

}

sub complement {

	my $nuc = shift;
	$nuc =~ tr/ATGCUatgcuNnYyRrSsWwKkMmBbDdHhVv/TACGAtacgaNnRrYySsWwMmKkVvHhDdBb/;
	return $nuc;

}

sub iupac_nuc {

	# Based on http://www.bioinformatics.org/sms/iupac.html

	#my @nuc = shift; # Array must contains at least 1 element (A, C, G, T, N, or *)
my @nuc = @_;
	my @nuc_without_asterisk;
	my $asterisk = 0;
	my $iupac_nuc;

	# Create another array without the "*"
	foreach my $elem (@nuc){

		if ($elem ne "*"){
			push @nuc_without_asterisk, $elem;
		}

		if ($elem eq "*"){
			$asterisk = 1;
		}

	}

	my $nuc_number = @nuc;

	# If length array is equal to 0, it means the initial array contained only "*"
	if ( @nuc_without_asterisk == 0 ){
		return "*";

	} elsif ( @nuc_without_asterisk == 1 ){

		# if A --> A
		if    ($nuc_without_asterisk[0] eq "A"){ $iupac_nuc = "A"}

		# if C --> C
		elsif ($nuc_without_asterisk[0] eq "C"){ $iupac_nuc = "C"}

		# if G --> G
		elsif ($nuc_without_asterisk[0] eq "G"){ $iupac_nuc = "G"}

		# if T --> T
		elsif ($nuc_without_asterisk[0] eq "T"){ $iupac_nuc = "T"}

		# else N --> N
		elsif ($nuc_without_asterisk[0] eq "N"){ $iupac_nuc = "N"}

		else {
			#print STDERR "Problem with iupac_nuc subroutine : nucleotide is not A, C, G, T or N\n";
			#print STDERR $nuc_without_asterisk[0]."\n";
		}

	} elsif (@nuc_without_asterisk == 2){

		# if A or G --> R
		if    ($nuc_without_asterisk[0] eq "A" && $nuc_without_asterisk[1] eq "G" ){ $iupac_nuc = "R"}
		elsif ($nuc_without_asterisk[0] eq "G" && $nuc_without_asterisk[1] eq "A" ){ $iupac_nuc = "R"}

		# if C or T --> Y
		elsif ($nuc_without_asterisk[0] eq "C" && $nuc_without_asterisk[1] eq "T" ){ $iupac_nuc = "Y"}
		elsif ($nuc_without_asterisk[0] eq "T" && $nuc_without_asterisk[1] eq "C" ){ $iupac_nuc = "Y"}

		# if G or C --> S
		elsif ($nuc_without_asterisk[0] eq "G" && $nuc_without_asterisk[1] eq "C" ){ $iupac_nuc = "S"}
		elsif ($nuc_without_asterisk[0] eq "C" && $nuc_without_asterisk[1] eq "G" ){ $iupac_nuc = "S"}

		# if A or T --> W
		elsif ($nuc_without_asterisk[0] eq "A" && $nuc_without_asterisk[1] eq "T" ){ $iupac_nuc = "W"}
		elsif ($nuc_without_asterisk[0] eq "T" && $nuc_without_asterisk[1] eq "A" ){ $iupac_nuc = "W"}

		# if G or T --> K
		elsif ($nuc_without_asterisk[0] eq "G" && $nuc_without_asterisk[1] eq "T" ){ $iupac_nuc = "K"}
		elsif ($nuc_without_asterisk[0] eq "T" && $nuc_without_asterisk[1] eq "G" ){ $iupac_nuc = "K"}

		# if A or C --> M
		elsif ($nuc_without_asterisk[0] eq "A" && $nuc_without_asterisk[1] eq "C" ){ $iupac_nuc = "M"}
		elsif ($nuc_without_asterisk[0] eq "C" && $nuc_without_asterisk[1] eq "A" ){ $iupac_nuc = "M"}

		# else --> N
		elsif ($nuc_without_asterisk[0] eq "N" || $nuc_without_asterisk[1] eq "N" ){ $iupac_nuc = "N"}

		else { #print STDERR "Problem with iupac_nuc subroutine : nucleotide is not A, C, G, T or N\n"; print STDERR $nuc_without_asterisk[0].$nuc_without_asterisk[1]."\n";
		}

	} elsif (@nuc_without_asterisk == 3){

		# if C or G or T --> B
		if    ($nuc_without_asterisk[0] eq "C" && $nuc_without_asterisk[1] eq "G" && $nuc_without_asterisk[2] eq "T" ){ $iupac_nuc = "B"}
		elsif ($nuc_without_asterisk[0] eq "C" && $nuc_without_asterisk[1] eq "T" && $nuc_without_asterisk[2] eq "G" ){ $iupac_nuc = "B"}
		elsif ($nuc_without_asterisk[0] eq "G" && $nuc_without_asterisk[1] eq "C" && $nuc_without_asterisk[2] eq "T" ){ $iupac_nuc = "B"}
		elsif ($nuc_without_asterisk[0] eq "G" && $nuc_without_asterisk[1] eq "T" && $nuc_without_asterisk[2] eq "C" ){ $iupac_nuc = "B"}
		elsif ($nuc_without_asterisk[0] eq "T" && $nuc_without_asterisk[1] eq "C" && $nuc_without_asterisk[2] eq "G" ){ $iupac_nuc = "B"}
		elsif ($nuc_without_asterisk[0] eq "T" && $nuc_without_asterisk[1] eq "G" && $nuc_without_asterisk[2] eq "C" ){ $iupac_nuc = "B"}

		# if A or G or T --> D
		elsif ($nuc_without_asterisk[0] eq "A" && $nuc_without_asterisk[1] eq "G" && $nuc_without_asterisk[2] eq "T" ){ $iupac_nuc = "D"}
		elsif ($nuc_without_asterisk[0] eq "A" && $nuc_without_asterisk[1] eq "T" && $nuc_without_asterisk[2] eq "G" ){ $iupac_nuc = "D"}
		elsif ($nuc_without_asterisk[0] eq "G" && $nuc_without_asterisk[1] eq "A" && $nuc_without_asterisk[2] eq "T" ){ $iupac_nuc = "D"}
		elsif ($nuc_without_asterisk[0] eq "G" && $nuc_without_asterisk[1] eq "T" && $nuc_without_asterisk[2] eq "A" ){ $iupac_nuc = "D"}
		elsif ($nuc_without_asterisk[0] eq "T" && $nuc_without_asterisk[1] eq "A" && $nuc_without_asterisk[2] eq "G" ){ $iupac_nuc = "D"}
		elsif ($nuc_without_asterisk[0] eq "T" && $nuc_without_asterisk[1] eq "G" && $nuc_without_asterisk[2] eq "A" ){ $iupac_nuc = "D"}

		# if A or C or T --> H
		elsif ($nuc_without_asterisk[0] eq "A" && $nuc_without_asterisk[1] eq "C" && $nuc_without_asterisk[2] eq "T" ){ $iupac_nuc = "H"}
		elsif ($nuc_without_asterisk[0] eq "A" && $nuc_without_asterisk[1] eq "T" && $nuc_without_asterisk[2] eq "C" ){ $iupac_nuc = "H"}
		elsif ($nuc_without_asterisk[0] eq "C" && $nuc_without_asterisk[1] eq "A" && $nuc_without_asterisk[2] eq "T" ){ $iupac_nuc = "H"}
		elsif ($nuc_without_asterisk[0] eq "C" && $nuc_without_asterisk[1] eq "T" && $nuc_without_asterisk[2] eq "A" ){ $iupac_nuc = "H"}
		elsif ($nuc_without_asterisk[0] eq "T" && $nuc_without_asterisk[1] eq "A" && $nuc_without_asterisk[2] eq "C" ){ $iupac_nuc = "H"}
		elsif ($nuc_without_asterisk[0] eq "T" && $nuc_without_asterisk[1] eq "C" && $nuc_without_asterisk[2] eq "A" ){ $iupac_nuc = "H"}

		# if A or C or G --> V
		elsif ($nuc_without_asterisk[0] eq "A" && $nuc_without_asterisk[1] eq "C" && $nuc_without_asterisk[2] eq "G" ){ $iupac_nuc = "V"}
		elsif ($nuc_without_asterisk[0] eq "A" && $nuc_without_asterisk[1] eq "G" && $nuc_without_asterisk[2] eq "C" ){ $iupac_nuc = "V"}
		elsif ($nuc_without_asterisk[0] eq "C" && $nuc_without_asterisk[1] eq "A" && $nuc_without_asterisk[2] eq "G" ){ $iupac_nuc = "V"}
		elsif ($nuc_without_asterisk[0] eq "C" && $nuc_without_asterisk[1] eq "G" && $nuc_without_asterisk[2] eq "A" ){ $iupac_nuc = "V"}
		elsif ($nuc_without_asterisk[0] eq "G" && $nuc_without_asterisk[1] eq "A" && $nuc_without_asterisk[2] eq "C" ){ $iupac_nuc = "V"}
		elsif ($nuc_without_asterisk[0] eq "G" && $nuc_without_asterisk[1] eq "C" && $nuc_without_asterisk[2] eq "A" ){ $iupac_nuc = "V"}

		# else --> N
		elsif ($nuc_without_asterisk[0] eq "N" || $nuc_without_asterisk[1] eq "N" || $nuc_without_asterisk[2] eq "N" ){ $iupac_nuc = "N"}

		else { #print STDERR "Problem with iupac_nuc subroutine : nucleotide is not A, C, G, T or N\n"; print STDERR $nuc_without_asterisk[0].$nuc_without_asterisk[1].$nuc_without_asterisk[2]."\n"; 
		}

	} else {
		print STDERR "Problem with iupac_nuc subroutine\n";
		$iupac_nuc = "N";
	}

	if ($asterisk == 1){
		$iupac_nuc = lc $iupac_nuc;
	}

	return $iupac_nuc;

}

################################################################################

############
# DEFAULTS #
############

my $pileup  = "NOTDEFINED";
my $outdir  = getcwd()    ;
my $verbose = 0           ;
my $help    = 0           ;

################################################################################

###########
# OPTIONS #
###########

GetOptions (
	"p|pileup=s" => \$pileup,  # file
	"o|outdir=s" => \$outdir,  # file
	"v|verbose"  => \$verbose, # flag
	"h|help"     => \$help,    # flag
) or die ("Error in command line arguments\n");

################################################################################

if( $help || $pileup eq "NOTDEFINED" || $outdir eq "NOTDEFINED" ) {
	usage();
}

################################################################################

################
# START SCRIPT #
################

print STDERR "\n# BEGIN extract_info_from_pileup.pl\n";

################################################################################

if ( $verbose ){
	print STDERR "\t[VERBOSE] ---> pileup file     : ".$pileup."\n";
	print STDERR "\t[VERBOSE] ---> Output diretory : ".$outdir."\n";
}

################################################################################

my $seqid;

my $pos      = 0;
my $last_pos = 0;
my $outfile  = "$outdir/info_from_pileup.tsv";

open( my $FH_pileup_out, '>>', $outfile ) or die "Could not open file '$outfile' $!";

if ( open( my $FH_pileup_in, '<', $pileup ) ) {

	if ( $verbose ){ print STDERR "\t[VERBOSE] ---> Process each position\n"; }

	# Process each line
	while ( my $line = <$FH_pileup_in> ) {

		# In the pileup format (without -u or -g), each line represents a genomic position, consisting of
			# - chromosome name,
			# - 1-based coordinate,
			# - reference base,
			# - the number of reads covering the site,
			# - read bases,
			# - base qualities
			# - and alignment mapping qualities.

		#seq1 272 T 24  ,.$.....,,.,.,...,,,.,..^+. <<<+;<<<<<<<<<<<=<;<;7<&

		# Information on match, mismatch, indel, strand, mapping quality and start and end of a read are all encoded at the read base column. At this column :
			# - a dot stands for a match to the reference base on the forward strand,
			# - a comma for a match on the reverse strand,
			# - a '>' or '<' for a reference skip,
			# - 'ACGTN' for a mismatch on the forward strand
			# - 'acgtn' for a mismatch on the reverse strand.
			# - a pattern '\\+[0-9]+[ACGTNacgtn]+' indicates there is an insertion between this reference position and the next reference position. The length of the insertion is given by the integer in the pattern, followed by the inserted sequence.
			# - a pattern '-[0-9]+[ACGTNacgtn]+' represents a deletion from the reference. The deleted bases will be presented as '*' in the following lines. * (asterisk) is a placeholder for a deleted base in a multiple basepair deletion that was mentioned in a previous line by the -[0-9]+[ACGTNacgtn]+ notation (source : wikipedia)
			# - a symbol `^' marks the start of a read. The ASCII of the character following '^' minus 33 gives the mapping quality.
			# - a symbol `$' marks the end of a read segment.

		chomp $line;
		my @column  = split "\t", $line;
		my $seqid   = $column[0];
		my $pos     = $column[1];
		my $ref_nuc = $column[2];
		my $cov     = $column[3];
		my $bases   = $column[4];
		# my $quals   = $column[5];

		if( $cov > 0 ) {

			# print STDERR "Read bases at position ".$pos.": ".$bases."\n";

			my %hash_count;
			my @insertions;
			my @corrected_insertions;
			my @deletions;
			my @corrected_deletions;
			my $indel_count = 0;

			# Remove symbol which marks the start of a read (and mapping quality following the symbol)
			while ( $bases =~ m/\^./ ){
				$bases =~ s/\^.//;
				# print STDERR "Read bases at position ".$pos.": ".$bases."\n";
			}

			# Remove symbol which marks the end of a read
			while ( $bases =~ m/\$/ ){
				$bases =~ s/\$//;
				# print STDERR "Read bases at position ".$pos.": ".$bases."\n";
			}

			# Insertion(s)
			if( $bases =~ m/\+[0-9]+[ACGTNacgtn]+/ ){

				# Store all insertions in an array
				@insertions = ($bases =~ m/\+([0-9]+[ACGTNacgtn]+)/g);

				# Process each insertion
				foreach my $insertion ( @insertions ){

					#TODO: output colonne supplémentaire avec toutes les insertions séparée par des ;

					# Get its length and sequence
					my ( $length, $seq ) = ( $insertion =~ m/([0-9]+)([ACGTNacgtn]+)/ );
					# Problem with this pattern when mismatches are presents after a insertion in the 'read bases' field
					# --> wrong sequence = too long

					# If mismatch(es) after insertion == insertion sequence too long
					if( $length != length( $seq ) ){

						# Get the real insertion sequence based on length
						my $newseq = substr( $seq, 0, $length );

						# Remove the insertion in the 'read bases' field
						$bases =~ s/\+$length$newseq//;
						# print STDERR "Read bases at position ".$pos.": ".$bases."\n";

						push @corrected_insertions, $newseq;

					} else { # If insertion sequence is correct

						# Remove the insertion in the 'read bases' field
						$bases =~ s/\+$length$seq//;
						# print STDERR "Read bases at position ".$pos.": ".$bases."\n";

						push @corrected_insertions, $seq;

					}

				}

			}

			# Deletion(s)
			if( $bases =~ m/\-[0-9]+[ACGTNacgtn]+/ ){

				# Store all insertions in an array
				@deletions = ($bases =~ m/\-([0-9]+[ACGTNacgtn]+)/g);

				# Process each deletion
				foreach my $deletion ( @deletions ){

					# Get its length and sequence
					my ( $length, $seq ) = ( $deletion =~ m/([0-9]+)([ACGTNacgtn]+)/ );
					# Problem with this pattern when mismatches are presents after a deletion in the 'read bases' field
					# --> wrong sequence = too long

					# If mismatch(es) after deletion == deletion sequence too long
					if( $length != length( $seq ) ){

						# Get the real insertion sequence based on length
						my $newseq = substr( $seq, 0, $length );

						# Remove the deletion in the 'read bases' field
						$bases =~ s/\-$length$newseq//;
						# print STDERR "Read bases at position ".$pos.": ".$bases."\n";

						push @corrected_deletions, $newseq;

					} else { # If deletion sequence is correct

						# Remove the deletion in the 'read bases' field
						$bases =~ s/\-$length$seq//;
						# print STDERR "Read bases at position ".$pos.": ".$bases."\n";

						push @corrected_deletions, $seq;

					}

				}

			}

			# print STDERR "Read bases at position ".$pos.": ".$bases."\n";

			# Count MATCHES
			my $dot_number   = ($bases =~ tr/.//);
			my $comma_number = ($bases =~ tr/,//);

			# Count MISMATCHES
			my $A_number     = ($bases =~ tr/Aa//);
			my $C_number     = ($bases =~ tr/Cc//);
			my $G_number     = ($bases =~ tr/Gg//);
			my $T_number     = ($bases =~ tr/Tt//);
			my $N_number     = ($bases =~ tr/Nn//);

			# Count DELETED bases
			my $asterik_count = ($bases =~ tr/*//);

			# Count INDELS
			$indel_count = ( @insertions + $asterik_count );

			# print STDERR "A : ".$A_number."\n";
			# print STDERR "C : ".$C_number."\n";
			# print STDERR "G : ".$G_number."\n";
			# print STDERR "T : ".$T_number."\n";
			# print STDERR "N : ".$N_number."\n";

			# Store MISMATCHES count in hash_count
			$hash_count{"A"} = $A_number;
			$hash_count{"C"} = $C_number;
			$hash_count{"G"} = $G_number;
			$hash_count{"T"} = $T_number;
			$hash_count{"N"} = $N_number;
			$hash_count{"*"} = $asterik_count;

			# Store MATCHES count in hash_count
			$hash_count{$ref_nuc} = $hash_count{$ref_nuc} + $dot_number + $comma_number;
			#$hash_count{complement($ref_nuc)} = $hash_count{complement($ref_nuc)} + $comma_number;

			my $total = ( $hash_count{"A"} + $hash_count{"C"} + $hash_count{"G"} + $hash_count{"T"} + $hash_count{"N"} + $hash_count{"*"} );

			if ( $total != $cov ){
				print STDERR "Problem ! Number of read bases is different than coverage at position ".$pos."(".$total." vs. ".$cov.")\n";
				print STDERR $column[4]."\n";
				print STDERR $bases."\n\n";
			}

			# Calculate 20% threshold
			my $threshold = $cov * 20 / 100;

			my @nuc_threshold;

			#TODO: gerer cas avec uniquement *

			foreach my $nuc ( keys %hash_count ){
				if( $hash_count{$nuc} > $threshold) {
					push @nuc_threshold, $nuc;
				}
			}

			# TODO: what do we do when none nucleotide is greater than threshold ?

			my $iupac_nuc;

			# calcul IUPAC nucleotide based on base count
			if ( @nuc_threshold > 0 ){
				$iupac_nuc = iupac_nuc( @nuc_threshold );
			} else {
				$iupac_nuc = 'N';
			}

			# Output :
				# chromosome/seqid
				# position
				# nb_de_A
				# nb_de_T
				# nb_de_G
				# nb_de_T
				# nb_de_indel
				# nucleotide_format_IUPAC

			if ( $indel_count >= $threshold ){
				$iupac_nuc = lc $iupac_nuc;
			}

			# iupac_minus
			my $iupac_nuc_minus = complement($iupac_nuc);

			if (@corrected_insertions == 0){

				print $FH_pileup_out $seqid."\t".$pos."\t".$hash_count{"A"}."\t".$hash_count{"C"}."\t".$hash_count{"G"}."\t".$hash_count{"T"}."\t".$indel_count."\t".$iupac_nuc."\t".$iupac_nuc_minus."\n";


			} else {

				my $insertions = join( ";", @corrected_insertions );

				print $FH_pileup_out $seqid."\t".$pos."\t".$hash_count{"A"}."\t".$hash_count{"C"}."\t".$hash_count{"G"}."\t".$hash_count{"T"}."\t".$indel_count."\t".$iupac_nuc."\t".$iupac_nuc_minus."\t".$insertions."\n";
			}

			#print STDERR "\n";

		} else { # Cov = 0
			print $FH_pileup_out $seqid."\t".$pos."\t0\t0\t0\t0\t0\tN\tN\tcov0\n";
		}

	}

	close( $FH_pileup_in );

} else {
	warn "Could not open file '$pileup' $!";
}

close( $FH_pileup_out );

################################################################################

##############
# END SCRIPT #
##############

print STDERR "# END extract_info_from_pileup.pl\n\n";
