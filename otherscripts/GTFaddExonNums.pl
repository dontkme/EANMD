#!/usr/bin/perl
use strict;
use warnings;

my $input_file = shift or die "Usage: $0 input.gtf\n";

open my $fh, '<', $input_file or die "Cannot open file: $input_file\n";

my %transcripts;

# First pass: Collect exons and CDS by transcript
while (<$fh>) {
    chomp;
    my @fields = split /\t/;
    next unless @fields == 9;

    my ($seqid, $source, $type, $start, $end, $score, $strand, $phase, $attributes) = @fields;
    my ($transcript_id) = $attributes =~ /transcript_id "([^"]+)"/;
    my ($gene_id) = $attributes =~ /gene_id "([^"]+)"/;

    if ($type eq 'exon' || $type eq 'CDS') {
        push @{$transcripts{$transcript_id}{$type}}, { 
            line => $_, start => $start, end => $end, strand => $strand, gene_id => $gene_id 
        };
    }
}
close $fh;

# Second pass: Add exon_id, exon_number, and gene_name
foreach my $transcript_id (keys %transcripts) {
    my $exon_number = 0;
    
    # Sort exons by strand direction
    my $strand = $transcripts{$transcript_id}{'exon'}[0]{strand};
    my @sorted_exons = $strand eq '+'
        ? sort { $a->{start} <=> $b->{start} } @{$transcripts{$transcript_id}{'exon'}}
        : sort { $b->{start} <=> $a->{start} } @{$transcripts{$transcript_id}{'exon'}};
    
    # Map exon range to exon_number and exon_id
    my %exon_map;
    foreach my $exon (@sorted_exons) {
        $exon_number++;
        my $exon_id = "${transcript_id}.exon$exon_number";
        
        # Store exon range with number and id
        $exon_map{$exon->{start}}{$exon->{end}} = { 
            number => $exon_number, 
            exon_id => $exon_id, 
            cds_id => "${transcript_id}.cds$exon_number" 
        };
        
        # Update exon annotation
        $exon->{line} =~ s/$/ gene_name \"$exon->{gene_id}\"; exon_id \"$exon_id\"; exon_number $exon_number;/;
        print $exon->{line}, "\n";
    }

    # Match CDS to correct exon_number and CDS-specific exon_id
    foreach my $cds (@{$transcripts{$transcript_id}{'CDS'}}) {
        my ($matched_exon_number, $matched_cds_id) = (0, '');
        
        foreach my $start (keys %exon_map) {
            foreach my $end (keys %{$exon_map{$start}}) {
                if ($cds->{start} >= $start && $cds->{end} <= $end) {
                    $matched_exon_number = $exon_map{$start}{$end}->{number};
                    $matched_cds_id = $exon_map{$start}{$end}->{cds_id};
                    last;
                }
            }
        }
        
        $cds->{line} =~ s/$/ gene_name \"$cds->{gene_id}\"; exon_id \"$matched_cds_id\"; exon_number $matched_exon_number;/;
        print $cds->{line}, "\n";
    }
}

__END__

# Usage:
# perl script_name.pl input.gtf > output.gff
