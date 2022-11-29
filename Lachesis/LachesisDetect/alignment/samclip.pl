#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;
use File::Basename;
use Data::Dumper;

#----------------------------------------------------------------------
# globals

my $EXE = basename($0);

# SAM file TSV columns
use constant {
  SAM_RNAME => 2,
  SAM_POS   => 3,
  SAM_CIGAR => 5,
  SAM_TLEN  => 8,
  SAM_SEQ   => 9,
  SAM_OPT_FIELD => 11
};

#----------------------------------------------------------------------
# command line parameters

my $max      = 20;
my $max_read_mismatch = 2;
my $max_pair_mismatch = 3;
my $ref      = '';
my $debug    = 0;
my $progress = 100_000;

#----------------------------------------------------------------------
sub usage {
  my($exitcode) = @_;
  $exitcode=0 if !defined($exitcode) or $exitcode eq 'help';
  my $fh = $exitcode ? \*STDERR : \*STDOUT;
  print $fh
    "SYNOPSIS\n  Filter SAM file for soft & hard clipped alignments\n",
    "USAGE\n",
    "  % samclip --ref ref.fa < in.sam > out.sam\n",
    "  % minimap2 ref.fa R1.fq R2.fq | samclip --ref ref.fa | samtools sort > out.bam\n",
    "OPTIONS\n",
    "  --help         This help\n",
    "  --ref FASTA    Reference genome - needs FASTA.fai index\n",
    "  --max NUM      Maximum clip length to allow (default=$max)\n",
    "  --debug        Print verbose debug info to stderr\n",
    "  --progress N   Print progress every NUM records (default=$progress,none=0)\n",
    "";
  exit($exitcode);
}

#----------------------------------------------------------------------
# getopts

@ARGV or usage(1);

GetOptions(
  "help"       => \&usage,
  "ref=s"      => \$ref,
  "max=i"      => \$max,
  "debug"      => \$debug,
  "progress=i" => \$progress,
) or usage(1);
             
$ref or err("Please supply reference genome with --ref");
$max >= 0 or err("Please supply --max >= 0");  
$ref .= ".fai" unless $ref =~ m/\.fai$/;
-r $ref or err("Can't see '$ref' index. Run 'samtools faidx $ref' ?"); 
!@ARGV and -t STDIN and err("Please provide or pipe a SAM file to $EXE");

#----------------------------------------------------------------------
# main

# get a hash of { seqname => length }
#msg("Loading: $ref");
my $len = fai_to_dict($ref);
msg(Dumper($len)) if $debug;
#msg("Found", scalar keys %$len, "sequences in $ref");

my $total=0;
my $softclip_removed=0;
my $softclip_both_sides_removed=0;
my $mismatch_read_removed=0;
my $mismatch_pair_removed=0;
my $kept=0;
my $header=0;

# read SAM one line ar a time
SAM_LINE: while (my $line = <ARGV>) {
  # SAM header
  if ($line =~ m/^@/) {
    print $line;
    $header++;
    next;
  }
  $total++;
  #msg("Processed $total records...") if $progress and $total % $progress == 0;
  my @sam = split m/\t/, $line;
  # do a quick 'clipped?' check before heavyweight parsing
  if ($sam[SAM_CIGAR] =~ /\d[SH]/) {
    my($HL, $SL, undef, $SR, $HR) 
      = ($sam[5] =~ m/ ^ (?:(\d+)H)? (?:(\d+)S)? (.*?) (?:(\d+)S)? (?:(\d+)H)? $/x);
    $HL ||= 0; $SL ||= 0; $SR ||= 0; $HR ||= 0;
    # if either end is clipped more than --max allowed, then remove it
    # unless it is at a contig end
    my $start = $sam[SAM_POS];
    my $end = $start + length($sam[SAM_SEQ]) - 1;
    my $contiglen = $len->{$sam[SAM_RNAME]} or err("Reference", $sam[SAM_RNAME], "not in '$ref'");
    msg("CHROM=$sam[SAM_RNAME]:1-$contiglen POS=$start..$end CIGAR=$sam[SAM_CIGAR] HL=$HL SL=$SL SR=$SR HR=$HR max=$max)") if $debug;
    unless ($start == 0 or $end >= $contiglen) {
      if ($HL+$SL+ $HR+$SR > $max) {
        #print "    SOFT CLIP\n";
        $softclip_removed++;
        next;
      } elsif ($SL > 0 and $SR > 0) { # NEW: remove read if both ends are soft-clipped
        $softclip_both_sides_removed++;
        #print "    SOFT CLIP ON BOTH SIDES\n";
        next;
      }
    }

  }

  #COUNT NUMBER OF MISMATCHES
  # nM = total number of mismatches in a read pair (only substitutions, indels do not count)
  # NM = total number of mismatches and indels in a read, i.e edit distance (substitutions AND indels)
  # optional fields are in cols 12+
  # nM and NM aren't always in the same column

  foreach (@sam[SAM_OPT_FIELD .. $#sam]) {
    if ($_ =~ m/nM:i:(\d+)/) {
      my $nM_pair = $1;

      if ($nM_pair > $max_pair_mismatch) {
        #print("    nM REMOVED\n");
        $mismatch_pair_removed++; 
        next SAM_LINE;
      }
    } elsif ($_ =~ m/NM:i:(\d+)/) {
      my $NM_read = $1;

      if ($NM_read > $max_read_mismatch) {
        #print("    NM REMOVED\n");
        $mismatch_read_removed++; 
        next SAM_LINE;
      }
    }
  }
  print $line;
  $kept++;
}
# stats
msg("Total SAM records $total");
msg("clip removed $softclip_removed");
msg("clip both sides removed $softclip_both_sides_removed");
msg("mismatch read (NM) removed $mismatch_read_removed");
msg("mismatch pair (nM) removed $mismatch_pair_removed");

if ($total > 0) {
    msg(sprintf("allowed $kept (%.1f%%)", $kept / $total * 100));
} else {
    msg(sprintf("allowed $kept (NA%%)"));
}
msg("Header contained $header lines");

#----------------------------------------------------------------------
sub fai_to_dict {
  my($fname) = @_;
  open my $fai, '<', $fname or err("Can't read FAI '$fname'");
  my $len = {};
  while (<$fai>) {
    my($name, $bp) = split m/\t/;
    $len->{$name} = $bp;
  }
  close $fai;
  return $len;
}



#----------------------------------------------------------------------
sub msg {
  print STDERR "@_\n";
}

#----------------------------------------------------------------------
sub err {
  msg("ERROR:", @_);
  exit(1);
}

