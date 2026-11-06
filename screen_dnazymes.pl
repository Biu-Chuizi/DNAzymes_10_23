#!/usr/bin/perl -w
use strict;

use Getopt::Long qw(:config no_ignore_case bundling);

#=begin Program Description
#
#	DNAzyme screening helper for DNAzyme_scan.pl output
#	Wraps: ranking by binding DeltaG -> logistic score -> sequence filters -> shortlist
#	Companion to DNAzyme_scan.pl (Antonio Marco, 2022) and run_analysis_DNAzymes.R
#
#	Basic usage:
#		perl screen_dnazymes.pl -i DNAzymes_HPV_16 -o shortlist
#
#	IMPORTANT - what the DeltaG column from DNAzyme_scan.pl can and cannot do:
#		DNAzyme_scan.pl reports the nearest-neighbour free energy of the arms
#		hybridised to the target, assuming a perfect, structure-free duplex. It does
#		not account for the intramolecular folding of the DNAzyme itself, for the
#		secondary structure of the target site, or for the catalytic step. On the 38
#		HPV16 DNAzymes with published cleavage data (HPV_FullLength.tab) that column
#		carries essentially no signal for actual efficiency (leave-one-out AUC 0.08),
#		whereas the 'Energy' column of that table, which was recomputed with the
#		ViennaRNA Package, reaches AUC 0.67. So use this script to enumerate and
#		prune candidates, then recompute energies with RNAup/RNAfold and rank with
#		rank_dnazymes.R. The P_ge40 column below is only meaningful once you have
#		fed in those ViennaRNA energies; computed from the scan DeltaG it is a
#		placeholder, which is why it is labelled as such in the output.
#=cut

my $infile      = 'outfile';
my $outbase     = 'shortlist';
my $min_abs_dg  = 22.0;   # keep candidates at least this stably bound (kcal/mol)
my $min_gc      = 0.20;   # per-arm GC fraction
my $max_gc      = 0.80;
my $max_homopol = 5;      # longest allowed run of one base inside an arm
my $top         = 0;      # 0 = keep everything that passes the filters
my $arm         = 9;      # arm length used by DNAzyme_scan.pl (default 9)
my $linker_len  = 15;     # ggctagctacaacga
my $b0          = -7.7642;
my $b1          = 0.3582;

GetOptions(
	'infile=s'    => \$infile,     'i=s' => \$infile,
	'outbase=s'   => \$outbase,    'o=s' => \$outbase,
	'mindg=f'     => \$min_abs_dg, 'd=f' => \$min_abs_dg,
	'mingc=f'     => \$min_gc,
	'maxgc=f'     => \$max_gc,
	'homopol=i'   => \$max_homopol,
	'armlength=i' => \$arm,        'l=i' => \$arm,
	'top=i'       => \$top,        't=i' => \$top,
	'intercept=f' => \$b0,
	'coef=f'      => \$b1,
	'help|h'      => \my $help_flag,
);

my $help = $help_flag || !-e $infile;
print helpScreen() and exit 0 if $help;

open(my $IN, '<', $infile) or die "Cannot read $infile: $!\n";
my @lines = <$IN>;
close $IN;
chomp @lines;
s/\r//g for @lines;

my $header = shift @lines;
unless ($header =~ /DNAzyme/ && $header =~ /DeltaG/) {
	die "'$infile' does not look like DNAzyme_scan.pl output (no DNAzyme/DeltaG columns).\n";
}

my @hits;
for my $line (@lines) {
	next unless $line =~ /\S/;
	my ($transcript, $type, $position, $rna, $dnazyme, $dg) = split /\t/, $line;
	next unless defined $dnazyme && defined $dg && $dg =~ /^-?\d+(\.\d+)?$/;

	# The scan output stores the DNAzyme as <left arm><linker><right arm>.
	my $core = $dnazyme;
	my $left  = substr($core, 0, $arm);
	my $right = substr($core, $arm + $linker_len);
	if (length($core) != 2 * $arm + $linker_len) {
		# Fall back to locating the linker explicitly (handles non-default arm lengths).
		my $idx = index(lc($core), 'ggctagctacaacga');
		if ($idx > 0) { $left = substr($core, 0, $idx); $right = substr($core, $idx + $linker_len); }
	}
	$left  = uc $left;
	$right = uc $right;

	my $abs_dg = abs($dg);
	my $prob   = 1 / (1 + exp(-($b0 + $b1 * $abs_dg)));

	my @flags;
	push @flags, 'weak_binding'  if $abs_dg < $min_abs_dg;
	for my $side (['left', $left], ['right', $right]) {
		my ($name, $seq) = @$side;
		my $gc = gcFraction($seq);
		push @flags, "${name}_GC_low"  if $gc < $min_gc;
		push @flags, "${name}_GC_high" if $gc > $max_gc;
		my $hp = longestHomopolymer($seq);
		push @flags, "${name}_homopolymer$hp" if $hp > $max_homopol;
		push @flags, "${name}_has_N" if $seq =~ /N/;
	}

	push @hits, {
		transcript => $transcript, type => $type, position => $position, rna => $rna,
		dnazyme => $dnazyme, dg => $dg, abs_dg => $abs_dg, prob => $prob,
		gc_left => gcFraction($left), gc_right => gcFraction($right),
		left => $left, right => $right, flags => join(',', @flags),
	};
}

# Rank: strongest binding first
@hits = sort { $b->{abs_dg} <=> $a->{abs_dg} } @hits;

my @kept = grep { $_->{flags} eq '' } @hits;
@kept = @kept[0 .. $top - 1] if $top > 0 && @kept > $top;
$_->{keep} = 1 for @kept;

my $tabfile = $outbase . '.tab';
my $fafile  = $outbase . '_dnazymes.fa';

open(my $OUT, '>', $tabfile) or die "Cannot write $tabfile: $!\n";
print $OUT join("\t", 'Rank', 'Transcript', 'Position', 'Target_site', 'DNAzyme', 'DeltaG',
	'P_ge40_unreliable', 'GC_left', 'GC_right', 'Left_arm', 'Right_arm', 'Keep', 'Filters') . "\n";
my $rank = 0;
for my $h (@hits) {
	$rank++;
	print $OUT join("\t", $rank, $h->{transcript}, $h->{position}, $h->{rna}, $h->{dnazyme},
		$h->{dg}, sprintf('%.3f', $h->{prob}), sprintf('%.2f', $h->{gc_left}),
		sprintf('%.2f', $h->{gc_right}), $h->{left}, $h->{right},
		($h->{keep} ? 'yes' : 'no'), ($h->{flags} eq '' ? 'PASS' : $h->{flags})) . "\n";
}
close $OUT;

open(my $FA, '>', $fafile) or die "Cannot write $fafile: $!\n";
my $i = 0;
for my $h (@kept) {
	$i++;
	my $id = 'DZ' . $i . '_' . $h->{position} . '_dG' . $h->{dg};
	$id =~ s/[^\w.-]/_/g;
	print $FA ">$id\n" . $h->{dnazyme} . "\n";
}
close $FA;

print "\n";
print "Scanned  : " . scalar(@hits) . " candidate DNAzymes from $infile\n";
print "Passed   : " . scalar(@kept) . " after filters (|DeltaG| >= $min_abs_dg, ${min_gc} <= GC <= $max_gc, homopolymer <= $max_homopol)\n";
print "Ranked   : $tabfile (all candidates, best binding first)\n";
print "Shortlist: $fafile ($i sequences)\n\n";
print "Next step - the DeltaG above is NOT a usable efficiency predictor on its own.\n";
print "Recompute the thermodynamic features with the ViennaRNA Package, then rank:\n\n";
print "  pip install ViennaRNA\n";
print "  python compute_features.py $fafile -o " . $outbase . "_features.tab\n\n";
print "compute_features.py uses the same energy set as the paper (Mathews 2004 DNA,\n";
print "--noLP) and reproduces the Internal column of HPV_FullLength.tab exactly.\n";
print "Equivalent command-line calls, if you prefer the ViennaRNA binaries:\n\n";
print "  RNAup -b -d2 --noLP -P dna_mathews2004.par --noconv -c 'S' $fafile > homodimer.out\n";
print "  RNAfold -p -d2 --noLP -P dna_mathews2004.par --noconv < $fafile > intramolecular.out\n\n";

exit 0;

sub gcFraction {
	my $s = $_[0];
	return 0 unless length $s;
	my $gc = ($s =~ tr/GCGC//);
	return $gc / length($s);
}

sub longestHomopolymer {
	my $s = $_[0];
	my $best = 0;
	while ($s =~ /(.)\1*/g) {
		my $len = length($1);
		$best = $len if $len > $best;
	}
	return $best;
}

sub helpScreen {
	my $o = "Usage: perl screen_dnazymes.pl -i <DNAzyme_scan_output> [-OPTIONS]\n\n";
	$o .= "Options [default]:\n";
	$o .= "  -i/--infile       output table from DNAzyme_scan.pl [outfile]\n";
	$o .= "  -o/--outbase      prefix for the ranked table and shortlist [$outbase]\n";
	$o .= "  -d/--mindg        minimum |DeltaG| to keep (kcal/mol) [$min_abs_dg]\n";
	$o .= "  --mingc/--maxgc   allowed GC fraction per arm [$min_gc - $max_gc]\n";
	$o .= "  --homopol         max run of a single base per arm [$max_homopol]\n";
	$o .= "  -l/--armlength    arm length used during the scan [$arm]\n";
	$o .= "  -t/--top          keep at most N passing candidates (0 = all) [$top]\n";
	$o .= "  --intercept/--coef  logistic coefficients for P(cleavage>=40%) [$b0, $b1]\n";
	return $o;
}
