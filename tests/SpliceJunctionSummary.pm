package Genome::Model::Tools::Transcriptome::SpliceJunctionSummary;

use strict;
use warnings;

use Genome;
#use Bio::DB::Sam;

class Genome::Model::Tools::Transcriptome::SpliceJunctionSummary {
    is => ['Genome::Model::Tools::Transcriptome::Base'],
    has_input => [
        output_directory => {
            is => 'Text',
            doc => 'The output directory where splice junction files are output.',
        },
        observed_junctions_bed12_file => {
            is => 'Text',
            doc => 'A BED12 format file of splice junctions.  Each junction has two blocks.',
        },
        reference_fasta_file => {
            is => 'Text',
            doc => 'The FASTA format reference sequence file.  The reference FASTA must be compatible with the annotation GTF file.',
        },
        annotation_gtf_file => {
            is => 'Text',
            doc => 'The GTF format annotation file used to annotate the observed junctions and determine if observed juctions match known junctions.',
        },
        annotation_name => {
            is => 'Text',
            doc => 'The name of the annotation used.  This name is also used as a basename for the output file.',
        },
        bedtools_version => {
            is => 'Text',
            is_optional => 1,
            doc => 'The version of the BEDTools software to be used.',
            default_value => Genome::Model::Tools::BedTools->default_bedtools_version,
        },
    ],
    has_output => [
        observed_junctions_bed6_file => {
            is => 'Text',
            doc => 'The output BED6 format file of observed splice junctions.',
            is_optional => 1,
        },
        annotation_genes_bed12_file => {
            is => 'Text',
            doc => 'The output BED12 format file of Gene annotation.',
            is_optional => 1,
        },
        known_junctions_bed6_file => {
            is => 'Text',
            doc => 'The output BED6 format file of known splice junctions.',
            is_optional => 1,
        },
        known_junctions_tsv_file => {
            is => 'Text',
            doc => 'The output TSV format file of known splice junctions.',
            is_optional => 1,
        },
        annotation_exons_bed6_file => {
            is => 'Text',
            doc => 'The output BED6 format file of Exon annotation.',
            is_optional => 1,
        },
        annotation_merged_exon_bed_file => {
            is => 'Text',
            doc => 'The output BED format file of merged Exon annotation.',
            is_optional => 1,
        },
        annotated_observed_junction_tsv_file => {
            is => 'Text',
            doc => 'The output TSV format file that contains annotated information about each observed splice junction.',
            is_optional => 1,
        },
        transcript_expression_tsv_file => {
            is => 'Text',
            doc => 'The output TSV format file containing transcript level expression values using JPJM.',
            is_optional => 1,
        },
        gene_expression_tsv_file => {
            is => 'Text',
            doc => 'The output TSV format file containing gene level expression values using JPJM.',
            is_optional => 1,
        },
        _known_junctions => {is_optional => 1,},
        _known_donors => {is_optional => 1,},
        _known_acceptors => {is_optional => 1,},
        _known_ec_blocks => {is_optional => 1,},
        _observed_junctions => {is_optional => 1,},
        _faidx_index => {is_optional => 1,},
        _transcripts => {is_optional => 1,},
        _genes => {is_optional => 1,},
    ],
};

sub help_detail {
    return <<EOS
This tool will determine the splice sites for all observed and known splice junctions.
A comparison of observed versus known is performed including annotating with all observed junctions with exon, transcript and gene level information.
Expression of transcripts and genes is calculated using JPJM (see below for equation):
    JPJM = (sum of exon junction read counts for a gene / number of exon-exon junctions of that gene) => then normalized to per million junction mapped reads
Finally, an R script generates a summary spreadsheet with result R graphics written as PDF files.

EOS
}

sub execute {
    my $self = shift;
    
    $self->_resolve_output_file_paths();

    $self->_load_fasta_reference_index();

    $self->_map_transcript_and_gene_ids();

    $self->_generate_observed_junctions();

    # infer the splice site of each observed junction
    $self->infer_splice_site('observed_junctions');

    $self->_generate_known_junctions();

    # infer the splice sites of these known junctions
    $self->infer_splice_site('known_junctions');

    $self->_write_known_junction_splice_sites();

    $self->_generate_exon_content();

    #Now annotate observed exon-exon junctions against databases of known junctions by Gene Name(symbol)
    $self->annotate_observed_junctions();

    # Annotation based gene/transcript level expression calculations:
    $self->calculate_gene_expression();

    $self->generate_summary_and_plots();

    return 1;
}

sub _resolve_output_file_paths {
    my $self = shift;

    $self->observed_junctions_bed6_file($self->output_directory .'/observed.junctions.bed');
    $self->known_junctions_bed6_file($self->output_directory .'/'. $self->annotation_name.'.junctions.bed');
    $self->known_junctions_tsv_file($self->output_directory .'/'. $self->annotation_name.'.junctions.tsv');
    $self->annotation_genes_bed12_file($self->output_directory .'/'. $self->annotation_name.'.Genes.bed');
    $self->annotation_exons_bed6_file($self->output_directory .'/'. $self->annotation_name.'.Exons.bed');
    $self->annotation_merged_exon_bed_file($self->output_directory .'/'. $self->annotation_name .'.ExonContent.bed');
    $self->annotated_observed_junction_tsv_file($self->output_directory .'/observed.junctions.anno.'. $self->annotation_name .'.tsv');
    $self->transcript_expression_tsv_file($self->output_directory .'/'. $self->annotation_name .'.Junction.TranscriptExpression.tsv');
    $self->gene_expression_tsv_file($self->output_directory .'/'. $self->annotation_name .'.Junction.GeneExpression.tsv');
    return 1;
}

sub _write_known_junction_splice_sites {
    my $self = shift;

    my $headers = Genome::Utility::IO::BedReader->default_headers;
    my @known_junctions_headers = @$headers;
    push @known_junctions_headers, ('splice_site','intron_size');
    my $known_junctions_tsv_writer = Genome::Utility::IO::SeparatedValueWriter->create(
        output => $self->known_junctions_tsv_file,
        separator => "\t",
        headers => \@known_junctions_headers,
        print_headers => 1,
    );
    my $known_junctions = $self->_known_junctions;
    my %known_junctions = %{$known_junctions};
    my @sorted_known_junction_ids = sort { $known_junctions{$a}->{order} <=> $known_junctions{$b}->{order} } keys %known_junctions;
    for my $known_junction_id (@sorted_known_junction_ids) {
        my $data = $known_junctions{$known_junction_id};
        delete($data->{order});
        $known_junctions_tsv_writer->write_one($data);
    }
    $known_junctions_tsv_writer->output->close;
    return 1;
}


sub _load_fasta_reference_index {
    my $self = shift;
    my $fai = Bio::DB::Sam::Fai->load($self->reference_fasta_file);
    unless ($fai) {
        die('Could not open faidx for fasta file : '. $self->reference_fasta_file);
    }
    $self->_faidx_index($fai);
    return 1;
}

sub _map_transcript_and_gene_ids {
    my $self = shift;
    my $annotation_gtf_file = $self->annotation_gtf_file;

    $self->debug_message('Importing gene, gene name, and transcript mappings from: '. $annotation_gtf_file);
    my $gff_reader = Genome::Utility::IO::GffReader->create(
        input => $annotation_gtf_file,
    );
    unless ($gff_reader) {
        $self->error_message('Failed to load annotation GTF file: '. $annotation_gtf_file);
        die($self->error_message);
    }
    my %genes;
    my %transcripts;
    while (my $data = $gff_reader->next_with_attributes_hash_ref) {
        my $attributes = delete($data->{attributes_hash_ref});
        my $transcript_id = $attributes->{transcript_id};
        my $gene_id = $attributes->{gene_id};
        my $gene_name = $attributes->{gene_name};
        if ( $transcripts{$transcript_id} ) {
            unless ($gene_id eq $transcripts{$transcript_id}{ensg_id}) {
                die('Multiple gene_ids '. $gene_id .' and '. $transcripts{$transcript_id}{ensg_id} .' found for transcript: '. $transcript_id);
            }
            unless ($gene_name eq $transcripts{$transcript_id}{name}) {
                die('Multiple gene_names '. $gene_name .' and '. $transcripts{$transcript_id}{name} .' found for transcript: '. $transcript_id);
            }
        } else {
            $transcripts{$transcript_id}{ensg_id} = $gene_id;
            $transcripts{$transcript_id}{name} = $gene_name;
        }
        if ( $genes{$gene_id} ) {
            unless ($gene_id eq $genes{$gene_id}{ensg_id}) {
                die('Multiple gene_ids '. $gene_id .' and '. $genes{$gene_id}{ensg_id} .' found for gene: '. $gene_id);
            }
            unless ($gene_name eq $genes{$gene_id}{name}) {
                die('Multiple gene_names '. $gene_name .' and '. $genes{$gene_id}{name} .' found for gene: '. $gene_id);
            }
        } else {
            $genes{$gene_id}{ensg_id} = $gene_id;
            $genes{$gene_id}{name} = $gene_name;
        }
    }
    $gff_reader->input->close;
    $self->_transcripts(\%transcripts);
    $self->_genes(\%genes);
    return 1;
}

sub _generate_observed_junctions {
    my $self = shift;

    # Covert BED12 format to BED6 junctions
    my $tmp_junctions_bed6_file = Genome::Sys->create_temp_file_path;
    unless ( Genome::Model::Tools::Bed::ToBedJunc->execute(
        bed12_file => $self->observed_junctions_bed12_file,
        bed_file => $tmp_junctions_bed6_file,
    )->result ) {
        $self->error_message('Failed to convert BED12 file \''. $self->observed_junctions_bed12_file .'\' to junctions BED file \''. $tmp_junctions_bed6_file .'\'');
        die($self->error_message);
    }
    
    # Sort the BED6 file
    $self->debug_message('Sorting merged BED6 junctions to file: '. $self->observed_junctions_bed6_file);
    unless (Genome::Model::Tools::BedTools::Sort->execute(
        input_file => $tmp_junctions_bed6_file,
        output_file => $self->observed_junctions_bed6_file,
        use_version => $self->bedtools_version,
    )->result) {
        $self->error_message('Failed to sort tmp junctions BED file \''. $tmp_junctions_bed6_file .'\' to \''. $self->observed_junctions_bed6_file .'\'');
        die($self->error_message);
    }

    # Load the observed junctions
    $self->debug_message('Loading the sorted and merged observed junctions: '. $self->observed_junctions_bed6_file);
    my $observed_junctions_bed6_reader = Genome::Utility::IO::BedReader->create(
        input => $self->observed_junctions_bed6_file,
    );
    my %observed_junctions;
    my $observed_order = 0;
    while (my $junction_data = $observed_junctions_bed6_reader->next) {
        my $key = $junction_data->{chr} .':'. $junction_data->{start} .'-'. $junction_data->{end} .'('. $junction_data->{strand} .')';
        $junction_data->{order} = $observed_order++;
        $observed_junctions{$key} = $junction_data;
    }
    $observed_junctions_bed6_reader->input->close;
    $self->_observed_junctions(\%observed_junctions);
    return 1;
}

sub _generate_known_junctions {
    my $self = shift;
    
    # Generate BED12 format genes file
    $self->debug_message('Convert annotation GTF \''. $self->annotation_gtf_file .'\' to BED12 format genes file \''. $self->annotation_genes_bed12_file .'\'');
    unless (Genome::Model::Tools::Gtf::ToBed12->execute(
        gtf_file => $self->annotation_gtf_file,
        bed12_file => $self->annotation_genes_bed12_file,
    )->result) {
        $self->error_message('Failed to generate BED12 format genes file: '. $self->annotation_genes_bed12_file);
        die($self->error_message);
    }
    
    # Make a BED6 junctions file from the BED12 genes file
    my $tmp_known_junctions_bed_file = Genome::Sys->create_temp_file_path($self->annotation_name.'.junctions.bed');
    $self->debug_message('Convert BED12 genes file \''. $self->annotation_genes_bed12_file .'\' to temporary BED6 format file of known junctions \''. $tmp_known_junctions_bed_file .'\'');
    unless (Genome::Model::Tools::Bed::ToBedJunc->execute(
        bed12_file => $self->annotation_genes_bed12_file,
        bed_file => $tmp_known_junctions_bed_file,
    )->result ) {
        $self->error_message('Failed to convert BED12 file \''. $self->annotation_genes_bed12_file .'\' to unosrted BED6 junctions file \''. $tmp_known_junctions_bed_file .'\'');
        die($self->error_message);
    }

    # Key on chr:start-end(strand) and sum the score as a count of transcripts
    $self->debug_message('Load the BED6 known junctions file \''. $tmp_known_junctions_bed_file .'\' and count transcripts observed per junction.');
    my $known_junctions_bed_reader = Genome::Utility::IO::BedReader->create(
        input => $tmp_known_junctions_bed_file,
    );
    unless ($known_junctions_bed_reader) {
        $self->error_message('Failed to load known junctions BED6 reader for unsorted temporary known junctions: '. $tmp_known_junctions_bed_file);
        die($self->error_message);
    }
    my %unordered_known_junctions;
    while (my $data = $known_junctions_bed_reader->next) {
        # TODO: Should known junction transcript counts enforce strandedness?
        my $key = $data->{chr} .':'. $data->{start} .'-'. $data->{end} .'('. $data->{strand} .')';
        unless ($unordered_known_junctions{$key}) {
            $unordered_known_junctions{$key} = $data;
        } else {
            $unordered_known_junctions{$key}->{name} .= ','. $data->{name};
            $unordered_known_junctions{$key}->{score}++;
        }
    }

    # Write out the unique known junctions in BED6 format
    my $tmp_unsorted_known_junctions_bed_file = Genome::Sys->create_temp_file_path($self->annotation_name.'.unsorted.junctions.bed');
    $self->debug_message('Writing unique known junctions as unsorted BED6 file: '. $tmp_unsorted_known_junctions_bed_file);
    my $unsorted_known_junctions_bed_writer = Genome::Utility::IO::BedWriter->create(
        output => $tmp_unsorted_known_junctions_bed_file,
    );
    for my $key (keys %unordered_known_junctions) {
        $unsorted_known_junctions_bed_writer->write_one($unordered_known_junctions{$key});
    }
    $unsorted_known_junctions_bed_writer->output->close;

    # Sort BED6
    $self->debug_message('Sorting BED6 known junctions to file: '. $self->known_junctions_bed6_file);
    unless (Genome::Model::Tools::BedTools::Sort->execute(
        input_file => $tmp_unsorted_known_junctions_bed_file,
        output_file => $self->known_junctions_bed6_file,
        use_version => $self->bedtools_version,
    )->result) {
        $self->error_message('Failed to sort known junctions BED6 file \''. $tmp_unsorted_known_junctions_bed_file .'\' to \''. $self->known_junctions_bed6_file.'\'');
        die($self->error_message);
    }
    
    # Load the known junctions
    $self->debug_message('Loading the sorted known junctions: '. $self->known_junctions_bed6_file);
    my $sorted_known_junctions_bed6_reader = Genome::Utility::IO::BedReader->create(
        input => $self->known_junctions_bed6_file,
    );

    my $transcripts = $self->_transcripts;
    my %known_junctions;
    my %known_donors;
    my %known_acceptors;
    my $known_order = 0;
    while (my $junction_data = $sorted_known_junctions_bed6_reader->next) {
        my $key = $junction_data->{chr} .':'. $junction_data->{start} .'-'. $junction_data->{end} .'('. $junction_data->{strand} .')';
        $junction_data->{order} = $known_order++;
        $known_junctions{$key} = $junction_data;

        my $transcript_ids = $junction_data->{name};
        my $chr = $junction_data->{chr};
        my $strand = $junction_data->{strand};
        my $start = $junction_data->{start};
        my $end = $junction_data->{end};
        if ($strand eq '+'){
            $known_donors{$chr}{$strand}{$start}{transcript_ids} = $transcript_ids;
            $known_acceptors{$chr}{$strand}{$end}{transcript_ids} = $transcript_ids;
        } elsif($strand eq '-'){
            $known_donors{$chr}{$strand}{$end}{transcript_ids} = $transcript_ids;
            $known_acceptors{$chr}{$strand}{$start}{transcript_ids} = $transcript_ids;
        } else{
            $self->error_message('Unknown strand in reference junctions file.');
            die($self->error_message);
        }
    }
    $sorted_known_junctions_bed6_reader->input->close;
    $self->_known_junctions(\%known_junctions);
    $self->_known_donors(\%known_donors);
    $self->_known_acceptors(\%known_acceptors);
    return 1;
}

sub _generate_exon_content {
    my $self = shift;
    
    # Generate BED6 format Exons file
    $self->debug_message('Convert BED12 format gene file \''. $self->annotation_genes_bed12_file .'\' to exon BED6 file \''. $self->annotation_exons_bed6_file);
    # TODO: Add BedTools version
    unless (Genome::Model::Tools::BedTools::Bed12ToBed6->execute(
        input_bed12_file => $self->annotation_genes_bed12_file,
        output_bed6_file => $self->annotation_exons_bed6_file,
        use_version => $self->bedtools_version,
    )->result) {
        $self->error_message('Failed to generate BED6 format exons file: '. $self->annotation_exons_bed6_file);
        die($self->error_message);
    }
    # Generate Exon Content BED6 format file
    $self->debug_message('Merge exon BED6 file to \''. $self->annotation_merged_exon_bed_file .'\'');
    # TODO: Add BedTools version
    unless (Genome::Model::Tools::BedTools::Merge->execute(
        input_file => $self->annotation_exons_bed6_file,
        output_file => $self->annotation_merged_exon_bed_file,
        force_strandedness => 1,
        report_number => 1,
        use_version => $self->bedtools_version,
    )->result) {
        $self->error_message('Failed to generate merged exon BED6 format file: '. $self->annotation_merged_exon_bed_file);
        die($self->error_message);
    }
    $self->import_exon_content_blocks();
    
    return 1;
}


############################################################################################################################
#Infer strand from alignment to reference genome                                                                           #
############################################################################################################################
sub infer_splice_site {
    my $self = shift;
    my $junctions_type = shift;
    my $junctions_hash_ref_method =  '_'. $junctions_type;

    my $junctions_hash_ref = $self->$junctions_hash_ref_method;
    my $fai = $self->_faidx_index;

    #Determine the splice site by comparison back to the reference genome...
    $self->debug_message('Determining splice sites from reference for '. $junctions_type);
    
    #Now go through the junctions and look for donor/acceptor splice sites at the coordinates reported
    #SPLICE_SITES = ["GT-AG", "CT-AC", "GC-AG", "CT-GC", "AT-AC", "GT-AT"]
    for my $junc_key (keys %{$junctions_hash_ref}) {
        my $junc_data = $junctions_hash_ref->{$junc_key};

        my $j_chr = $junc_data->{chr};
        my $left = $junc_data->{start};
        my $right = $junc_data->{end};

        my $intron_size = ($right - $left) + 1;
        $junc_data->{intron_size} = $intron_size;

        my $left_seq_id = $j_chr .':'. ($left + 1) .'-'. ($left+2);
        my $left_dn = uc($fai->fetch($left_seq_id));
        unless (length($left_dn) == 2) {
            die('Invalid sequence length for left seq_id '. $left_seq_id);
        }
        my $right_seq_id = $j_chr .':'. ($right-2) .'-'. ($right-1);
        my $right_dn = uc($fai->fetch($right_seq_id));
        unless (length($right_dn) == 2) {
            die('Invalid sequence length for right seq_id '. $right_seq_id);
        }
        # Strand is assigned by the aligner...
        # TODO: Is the strand in the junction key ex. chr:start-end(strand)
        # TODO: Could we determine the correctness of the donor/acceptor if we know the strand already
        my $splice_site;
        if ($left_dn eq "GT" && $right_dn eq "AG"){
            $splice_site = "GT-AG";
        } elsif($left_dn eq "CT" && $right_dn eq "AC"){
            $splice_site = "GT-AG";
        } elsif($left_dn eq "GC" && $right_dn eq "AG"){
            $splice_site = "GC-AG";
        } elsif($left_dn eq "CT" && $right_dn eq "GC"){
            $splice_site = "GC-AG";
        } elsif($left_dn eq "AT" && $right_dn eq "AC"){
            $splice_site = "AT-AC";
        } elsif($left_dn eq "GT" && $right_dn eq "AT"){
            $splice_site = "AT-AC";
        } else {
            $splice_site = "Other";
        }
        $junc_data->{splice_site} = $splice_site;
    }
    $self->debug_message('Finished infering splice sites for '. $junctions_type);
    $self->$junctions_hash_ref_method($junctions_hash_ref);
    return 1;
}

sub import_exon_content_blocks {
    my $self = shift;

    $self->debug_message('Importing reference exon content blocks from file: '. $self->annotation_merged_exon_bed_file);
    my %known_ec_blocks;
    #Import the ref exon content block file.  In this file, the coordinates are always ordered regardless of strand
    #i.e. the 'left' coordinate is always a smaller number than the 'right' coordinate
    my $exon_reader = Genome::Utility::IO::BedReader->create(
        input => $self->annotation_merged_exon_bed_file,
        # After mergeBed the names are replaced by the count or score
        headers => ['chr','start','end','score','strand'],
    );
    my $exon_count = 0;
    while (my $exon_data = $exon_reader->next) {
        $exon_count++;
        if ($known_ec_blocks{$exon_data->{chr}}{$exon_data->{strand}}){
            my $ec_ref = $known_ec_blocks{$exon_data->{chr}}{$exon_data->{strand}};
            $ec_ref->{$exon_count}->{left} = $exon_data->{start};
            $ec_ref->{$exon_count}->{right} = $exon_data->{end};
        } else {
            my %tmp;
            $tmp{$exon_count}{left} = $exon_data->{start};
            $tmp{$exon_count}{right} = $exon_data->{end};
            $known_ec_blocks{$exon_data->{chr}}{$exon_data->{strand}} = \%tmp;
        }
    }
    $exon_reader->input->close;
    $self->debug_message('Imported '. $exon_count .' known exon content blocks.');
    $self->_known_ec_blocks(\%known_ec_blocks);
    return 1;
}

sub annotate_observed_junctions {
    my $self = shift;

    my $observed_junctions = $self->_observed_junctions;
    my $known_junctions = $self->_known_junctions;
    my $known_donors = $self->_known_donors;
    my $known_acceptors = $self->_known_acceptors;

    my %observed_junctions = %{$observed_junctions};
    my %known_junctions = %{$known_junctions};
    my %known_donors = %{$known_donors};
    my %known_acceptors = %{$known_acceptors};

    for my $jid (keys %observed_junctions) {
        my $chr = $observed_junctions{$jid}->{chr};
        my $strand = $observed_junctions{$jid}->{strand};
        my $left = $observed_junctions{$jid}->{start};
        my $right = $observed_junctions{$jid}->{end};
        
        $observed_junctions{$jid}->{anchored} = 'N';
        $observed_junctions{$jid}->{transcript_ids} = 'na';
        #First check for an exact match to a known junction (i.e. anchored by both Donor and Acceptor)
        if ($known_junctions{$jid}){
            $observed_junctions{$jid}->{anchored} = 'DA';
            $observed_junctions{$jid}->{transcript_ids} = $known_junctions{$jid}{name};
        } else {
            #Now check for anchoring on one side only
            #    - If the Donor matches to a gene on the positive strand, the non matching Acceptor should have a higher chromosome coordinate
            #    - If the Donor matches to a gene on the negative strand, the non matching Acceptor should have a lower chromosome coordinate
            #    - If the Acceptor matches to a gene on the positive strand, the non-matching Donor should have a lower chromosome coordinate
            #    - If the Acceptor matches to a gene on the negative strand, the non-matching Donor should have a higher chromosome coordinate
            if ($strand eq '+'){
                if ( $known_donors{$chr}{'+'}{$left} && $known_acceptors{$chr}{'+'}{$right} ) {
                    $observed_junctions{$jid}->{anchored} = "NDA";
                    $observed_junctions{$jid}->{transcript_ids} = $known_donors{$chr}{'+'}{$left}{transcript_ids};
                } elsif ( $known_donors{$chr}{'+'}{$left} ){
                    $observed_junctions{$jid}->{anchored} = "D";
                    $observed_junctions{$jid}->{transcript_ids} = $known_donors{$chr}{'+'}{$left}{transcript_ids};
                } elsif ($known_acceptors{$chr}{'+'}{$right}) {
                    $observed_junctions{$jid}->{anchored} = "A";
                    $observed_junctions{$jid}->{transcript_ids} = $known_acceptors{$chr}{'+'}{$right}{transcript_ids};
                }
            } elsif ($strand eq '-') {
                if ( $known_donors{$chr}{'-'}{$right} && $known_acceptors{$chr}{'-'}{$left} ) {
                    $observed_junctions{$jid}->{anchored} = "NDA";
                    $observed_junctions{$jid}->{transcript_ids} = $known_donors{$chr}{'-'}{$right}{transcript_ids};
                } elsif ($known_donors{$chr}{'-'}{$right}) {
                    $observed_junctions{$jid}->{anchored} = "D";
                    $observed_junctions{$jid}->{transcript_ids} = $known_donors{$chr}{'-'}{$right}{transcript_ids};
                } elsif ($known_acceptors{$chr}{'-'}{$left}){
                    $observed_junctions{$jid}->{anchored} = "A";
                    $observed_junctions{$jid}->{transcript_ids} = $known_acceptors{$chr}{'-'}{$left}{transcript_ids};
                }
            }
        }
        #Initialize the skipping values
        my $anchored = $observed_junctions{$jid}->{anchored};
        if ($anchored eq 'N'){
            $observed_junctions{$jid}->{exons_skipped} = 'na';
            $observed_junctions{$jid}->{donors_skipped} = 'na';
            $observed_junctions{$jid}->{acceptors_skipped} = 'na';
        }else{
            $observed_junctions{$jid}->{exons_skipped} = 0;
            $observed_junctions{$jid}->{donors_skipped} = 0;
            $observed_junctions{$jid}->{acceptors_skipped} = 0;
        }
    }
    $self->_observed_junctions(\%observed_junctions);
    $self->annotate_skipping();
    
    #Calculate the junction read count per million junction reads mapped (JPM)
    my $grand_count = 0;
    $observed_junctions = $self->_observed_junctions();
    %observed_junctions = %{$observed_junctions};
    foreach my $jid (keys %observed_junctions){
        $grand_count += $observed_junctions{$jid}->{score};
    }
    foreach my $jid (keys %observed_junctions){
        my $jpm = $observed_junctions{$jid}->{score} * (1000000 / $grand_count);
        $observed_junctions{$jid}->{jpm} = $jpm;
    }

    my @output_headers = qw/chr start end name score strand jpm intron_size splice_site anchored exons_skipped donors_skipped acceptors_skipped transcript_ids gene_ids gene_names/;
    my $output_tsv_writer = Genome::Utility::IO::SeparatedValueWriter->create(
        output => $self->annotated_observed_junction_tsv_file,
        headers => \@output_headers,
        separator => "\t",
        print_headers => 1,
    );
    unless ($output_tsv_writer) {
        $self->error_message('Failed to load output TSV file: '. $self->annotated_observed_junction_tsv_file);
        die($self->error_message);
    }

    my $transcripts = $self->_transcripts;
    my %transcripts = %{$transcripts};
    my @sorted_junction_ids = sort {$observed_junctions{$a}->{order} <=> $observed_junctions{$b}->{order}} keys %observed_junctions;
    foreach my $jid (@sorted_junction_ids){
        my $data = $observed_junctions{$jid};
        delete($data->{order});
        my @transcript_ids = split(',',$data->{transcript_ids});
        my @gene_ids;
        my @gene_names;
        for my $transcript_id (@transcript_ids) {
            if ($transcript_id eq 'na') {
                push @gene_ids, 'na';
                push @gene_names, 'na';
            } else {
                push @gene_ids, $transcripts{$transcript_id}{ensg_id};
                push @gene_names, $transcripts{$transcript_id}{name};
            }
        }
        $data->{gene_ids} = join(',', @gene_ids);
        $data->{gene_names} = join(',', @gene_names);
        $output_tsv_writer->write_one($data);
    }
    $output_tsv_writer->output->close;

    $self->debug_message('Printed resulting annotated junctions to: '. $self->annotated_observed_junction_tsv_file);
    return;
}



##############################################################################################################################################
#Determine the exon and splice site skipping of each junction observed (if it known, or anchored to a known splice site at one or both ends) #
##############################################################################################################################################
sub annotate_skipping {
    my $self = shift;
    
    my $observed_junctions = $self->_observed_junctions();
    my $known_ec_blocks = $self->_known_ec_blocks();
    my $known_donors = $self->_known_donors();
    my $known_acceptors = $self->_known_acceptors();

    my %observed_junctions = %{$observed_junctions};
    my %known_ec_blocks = %{$known_ec_blocks};
    my %known_donors = %{$known_donors};
    my %known_acceptors = %{$known_acceptors};
    
    $self->debug_message('Annotating skipping of each junction - exons, then acceptors, then donors - Using BEDTools');

    #Use BEDTools to determine the exons, acceptors, or donors contained within each exon-exon observed junction (i.e. intron)
    #Do this by writing a temp BED file for each pair of coordinates (of the form: chr, start, end, strand) and then parsing the results file
    
    #'intersectBed -a exons.bed -b junctions.bed -f 1.0 -s -wa -wb'  
    #The '-f 1.0' option should give the exons that are entirely overlapped by junctions
    #The '-s' option should make this overlap between things on the same strand
    #The '-wa -wb' options, write the original coordinates for exons and junctions (as opposed to the merged coordinates).  Each overlaping pair will be reported as a seperate line
    #Then process the output file and count the exon entries associated with each junction entry
    
    #Example BEDTools command
    #Run BEDTools as follows and 'cut' the column containing the junction ID value.  
    #The occurence count for each of these IDs, should be the number of exons contained within the corresponding junction (i.e. skipped)
    #Unix sort and uniq can even do this counting for us...
    #Then just parse the counts out

    my $temp_obs_junctions = Genome::Sys->create_temp_file_path('ObsJunc.tmp.bed');
    my $temp_known_exons = Genome::Sys->create_temp_file_path('KnownExonContent.tmp.bed');
    my $temp_known_donors = Genome::Sys->create_temp_file_path('KnownDonors.tmp.bed');
    my $temp_known_acceptors = Genome::Sys->create_temp_file_path('KnownAcceptors.tmp.bed');
    my $temp_result = Genome::Sys->create_temp_file_path('Result.tmp.txt');
    
    #OBSERVED JUNCTIONS
    $self->debug_message('Looking for entire exons skipped');
    my $temp_obs_junction_writer = Genome::Utility::IO::BedWriter->create(
        output => $temp_obs_junctions,
    );
    unless ($temp_obs_junction_writer) {
        $self->error_message('Failed to open BED6 writer for temporary observed junctions: '. $temp_obs_junctions);
    }
    #Print print out the observed hmmSplicer junctions (i.e. introns as a temp bed file
    foreach my $jid (sort {$observed_junctions{$a}->{order} <=> $observed_junctions{$b}->{order}} keys %observed_junctions) {
        my %tmp_data = (
            'chr' => $observed_junctions{$jid}->{'chr'},
            'start' => $observed_junctions{$jid}->{'start'},
            'end' => $observed_junctions{$jid}->{'end'},
            'name' => $jid,
            'score' => '.',
            'strand' => $observed_junctions{$jid}->{'strand'},
        );
        $temp_obs_junction_writer->write_one(\%tmp_data);
    }
    $temp_obs_junction_writer->output->close;

    #ENTIRE EXONS
    #Print out all the known exon content blocks
    my $temp_known_exons_writer = Genome::Utility::IO::BedWriter->create(
        output => $temp_known_exons,
    );
    unless ($temp_known_exons_writer) {
        $self->error_message('Failed to open BED6 writer for temporary known exons: '. $temp_known_exons);
    }
    foreach my $chr (sort keys %known_ec_blocks){
        foreach my $strand (sort keys %{$known_ec_blocks{$chr}}){
            my $ec_ref = $known_ec_blocks{$chr}{$strand};
            foreach my $ec (sort {$ec_ref->{$a}->{left} <=> $ec_ref->{$b}->{left}} keys %{$ec_ref}){
                my $ec_left = $ec_ref->{$ec}->{left};
                my $ec_right = $ec_ref->{$ec}->{right};
                my %tmp_data = (
                    'chr' => $chr,
                    'start' => $ec_left,
                    'end' => $ec_right,
                    'name' => 'EC',
                    'score' => '.',
                    'strand' => $strand,
                );
                $temp_known_exons_writer->write_one(\%tmp_data);
            }
        }
    }
    $temp_known_exons_writer->output->close;

    my $temp_exons_intersect_bed = Genome::Sys->create_temp_file_path('intersect_exons_junctions.bed');
    unless (Genome::Model::Tools::BedTools::Intersect->execute(
        force_strandedness => 1,
        intersection_type => 'write_both',
        minimum_overlap => 1,
        input_file_a => $temp_known_exons,
        input_file_a_format => 'bed',
        input_file_b => $temp_obs_junctions,
        output_file => $temp_exons_intersect_bed,
        output_file_format => 'bed',
        use_version => $self->bedtools_version,
    )->result) {
        die('Failed to intersect exons and observed junctions.');
    }
    
    # TODO: Can this be made a command or performed with Perl inline
    my $exons_cmd = "cat $temp_exons_intersect_bed | cut -f 10 | sort | uniq -c > $temp_result";
    Genome::Sys->shellcmd(cmd => $exons_cmd);
    
    open (COUNTS, "$temp_result") || die "\n\nCould not open temp results file: $temp_result\n\n";
    while(<COUNTS>){
        chomp($_);
        if ($_ =~ /(\d+)\s+(.*)/){
            my $count = $1;
            my $jid = $2;
            if ($observed_junctions{$jid}){
                my $anchored = $observed_junctions{$jid}->{anchored};
                unless ($anchored eq "N"){
                    $observed_junctions{$jid}->{exons_skipped} = $count;
                }
            }else{
                $self->error_message('Unrecognized junction id:'. $jid);
                die($self->error_message);
            }
        }else{
            $self->error_message('Entry in results file not understood: '. $_);
            die($self->error_message);
        }
    }
    close(COUNTS);


    #DONORS
    $self->debug_message('Looking for donors skipped');
    my $temp_known_donors_writer = Genome::Utility::IO::BedWriter->create(
        output => $temp_known_donors,
    );
    unless ($temp_known_donors_writer) {
        $self->error_message('Failed to open BED6 writer for temporary known donors: '. $temp_known_donors);
    }
    foreach my $chr (sort keys %known_donors){
        foreach my $strand (sort keys %{$known_donors{$chr}}){
            foreach my $donor (sort keys %{$known_donors{$chr}{$strand}}){
                my $donor_p = $donor+1;
                my %tmp_data = (
                    'chr' => $chr,
                    'start' => $donor,
                    'end' => $donor_p,
                    'name' => 'D',
                    'score' => '.',
                    'strand' => $strand,
                );
                $temp_known_donors_writer->write_one(\%tmp_data);
            }
        }
    }
    $temp_known_donors_writer->output->close;

    my $temp_donors_intersect_bed = Genome::Sys->create_temp_file_path('intersect_donors_junctions.bed');
    unless (Genome::Model::Tools::BedTools::Intersect->execute(
        force_strandedness => 1,
        intersection_type => 'write_both',
        minimum_overlap => 1,
        input_file_a => $temp_known_donors,
        input_file_a_format => 'bed',
        input_file_b => $temp_obs_junctions,
        output_file => $temp_donors_intersect_bed,
        output_file_format => 'bed',
        use_version => $self->bedtools_version,
    )->result) {
        die('Failed to intersect donors and observed junctions.');
    }
    
    # TODO: Can this be made a command or performed with Perl inline
    my $donors_cmd = "cat $temp_donors_intersect_bed | cut -f 10 | sort | uniq -c > $temp_result";
    Genome::Sys->shellcmd(cmd => $donors_cmd);
    
    open (COUNTS, "$temp_result") || die "\n\nCould not open temp results file: $temp_result\n\n";
    while(<COUNTS>){
        chomp($_);
        if ($_ =~ /(\d+)\s+(.*)/){
            my $count = $1;
            my $jid = $2;
            if ($observed_junctions{$jid}){
                my $anchored = $observed_junctions{$jid}->{anchored};
                unless ($anchored eq "N"){
                    $observed_junctions{$jid}->{donors_skipped} = $count;
                }
            }else{
                $self->error_message('Unrecognized junction id:'. $jid);
                die($self->error_message);
            }
        }else{
            $self->error_message('Entry in results file not understood:'. $_);
            die($self->error_message);
        }
    }
    close(COUNTS);


    #ACCEPTORS
    $self->debug_message('Looking for acceptors skipped');
    my $temp_known_acceptors_writer = Genome::Utility::IO::BedWriter->create(
        output => $temp_known_acceptors,
    );
    unless ($temp_known_acceptors_writer) {
        $self->error_message('Failed to open BED6 writer for temporary known acceptors: '. $temp_known_acceptors);
    }
    foreach my $chr (sort keys %known_acceptors){
        foreach my $strand (sort keys %{$known_acceptors{$chr}}){
            foreach my $acceptor (sort keys %{$known_acceptors{$chr}{$strand}}){
                my $acceptor_p = $acceptor+1;
                my %tmp_data = (
                    'chr' => $chr,
                    'start' => $acceptor,
                    'end' => $acceptor_p,
                    'name' => 'A',
                    'score' => '.',
                    'strand' => $strand,
                );
                $temp_known_acceptors_writer->write_one(\%tmp_data);
            }
        }
    }
    $temp_known_acceptors_writer->output->close;

    my $temp_acceptors_intersect_bed = Genome::Sys->create_temp_file_path('intersect_acceptors_junctions.bed');
    unless (Genome::Model::Tools::BedTools::Intersect->execute(
        force_strandedness => 1,
        intersection_type => 'write_both',
        minimum_overlap => 1,
        input_file_a => $temp_known_acceptors,
        input_file_a_format => 'bed',
        input_file_b => $temp_obs_junctions,
        output_file => $temp_acceptors_intersect_bed,
        output_file_format => 'bed',
        use_version => $self->bedtools_version,
    )->result) {
        die('Failed to intersect acceptors and observed junctions.');
    }
    
    # TODO: Can this be made a command or performed with Perl inline
    my $acceptors_cmd = "cat $temp_acceptors_intersect_bed  | cut -f 10 | sort | uniq -c > $temp_result";
    Genome::Sys->shellcmd(cmd => $acceptors_cmd);
    open (COUNTS, "$temp_result") || die "\n\nCould not open temp results file: $temp_result\n\n";
    while(<COUNTS>){
        chomp($_);
        if ($_ =~ /(\d+)\s+(.*)/){
            my $count = $1;
            my $jid = $2;
            if ($observed_junctions{$jid}){
                my $anchored = $observed_junctions{$jid}->{anchored};
                unless ($anchored eq "N"){
                    $observed_junctions{$jid}->{acceptors_skipped} = $count;
                }
            } else{
                $self->error_message('Unrecognized junction id: '. $jid);
                die($self->error_message);
            }
        } else{
            $self->error_message('Entry in results file not understood: '. $_);
            die($self->error_message);
        }
    }
    close(COUNTS);

    $self->_observed_junctions(\%observed_junctions);
    $self->_known_ec_blocks(\%known_ec_blocks);
    $self->_known_donors(\%known_donors);
    $self->_known_acceptors(\%known_acceptors);

    return 1;
}


############################################################################################################################
#Ensembl based gene/transcript level expression calculations:                                                              #
############################################################################################################################
sub calculate_gene_expression {
    my $self = shift;

    # Inputs
    my $known_junctions = $self->_known_junctions;
    my $observed_junctions = $self->_observed_junctions;
    my $transcripts = $self->_transcripts;
    my $genes = $self->_genes;
    
    # Outputs
    my $transcript_expression_file = $self->transcript_expression_tsv_file;
    my $gene_expression_file = $self->gene_expression_tsv_file;

    #my $entrez_ensembl_data = &loadEntrezEnsemblData();

    #Calculate gene/transcript level read counts and expression estimates from exon-exon junction counts
    #Gene/transcript level expression = (sum of exon junction read counts for a gene / number of exon-exon junctions of that gene) => then normalized to per million junction mapped reads (JPJM)
    #Only known 'DA' junctions among the observed junctions will be used for these calculations
    #Store the proportion of junctions of the gene that were observed at least 1X, 5X, 10X, etc.
    #Make sure the output file has both ENSG/ENST ID, Gene Symbol, and Mapped Gene Symbol
    #Also output the number of known exon-exon junctions of the gene
    #Note that one junction can correspond to multiple genes - in this implementation reads for these junctions can be counted multiple times (once for each gene they correspond to)??
    #Create a gene/transcript expression record for every gene/transcript that has at least one exon-exon junction regardless of whether it was detected or not
    
    #1.) Build a map of:
    #-   all genes to gene names and list of transcripts
    #-   all transcripts to gene names and gene ids
    #-   also store a list of known junctions keyed on jid for convenient lookups

    my %known_junctions = %{$known_junctions};
    my %transcripts = %{$transcripts};
    my %genes = %{$genes};
    my %observed_junctions = %{$observed_junctions};
    
    #2.) Build a map of all known junctions to their: gene_ids, transcript_ids.  i.e.
    #    - Add junctions to transcript hash create above -> also store known junction count
    #$self->debug_message('Importing known junction to transcript mappings from: '. $known_junction_file);
    #open (JUNC1, "$known_junction_file") || die "\n\nCould not open known junction file: $known_junction_file\n\n";
    #while(<JUNC1>){
    for my $j (keys %known_junctions) {
        my $chr = $known_junctions{$j}{chr};
        my $count = $known_junctions{$j}{score};
        my $tid_string = $known_junctions{$j}{name};
        my @tids = split(",", $tid_string);
        #Attach junction list to each transcript record
        my %gene_ids;
        foreach my $enst_id (@tids){
            if (defined($transcripts{$enst_id}{junction_list})){
                my $jlist = $transcripts{$enst_id}{junction_list};
                $jlist->{$j}->{count} = $count;
            } else {
                my %jlist;
                $jlist{$j}{count} = $count;
                $transcripts{$enst_id}{junction_list} = \%jlist;
                $transcripts{$enst_id}{chromosome} = $chr;
            }
            $gene_ids{$transcripts{$enst_id}{ensg_id}} = 1;
        }
        foreach my $ensg_id (keys %gene_ids){
            if (defined($genes{$ensg_id}{junction_list})){
                my $jlist = $genes{$ensg_id}{junction_list};
                $jlist->{$j}->{count} = $count;
            } else{
                my %jlist;
                $jlist{$j}{count} = $count;
                $genes{$ensg_id}{junction_list} = \%jlist;
                $genes{$ensg_id}{chromosome} = $chr;
            }
        }
    }

    #3.) Build a hash of all observed junctions (only store 'DA' junctions) and their read counts
    #    - Store the grand read count of the library
    #    - Store the grand number of junctions observed for the library
    my $grand_read_count = 0;
    for my $jid (keys %observed_junctions) {
        #Only store know junctions ('DA' anchored)
        unless ($known_junctions{$jid}) { next(); }
        $grand_read_count += $observed_junctions{$jid}->{score};
    }
    my $grand_observed_junction_count = keys %observed_junctions;
    $self->debug_message('Found '. $grand_observed_junction_count .' observed known junctions with a grand read count of '. $grand_read_count);

    #4.) Now go through each gene/transcript in the gene/transcript hash, and go through each junction associated with it
    #    - Count the number of junctions actually observed - remove genes/transcripts with 0 known junctions
    #    - Store the cumulative read count for this gene/transcript
    #    - Store the proportion of junctions of the gene that were observed at least 1X, 5X, 10X, etc.
    #    - Gene/transcript level expression = (sum of exon junction read counts for a gene / number of exon-exon junctions of that gene) => then normalized to per million junction mapped reads (JPJM)
    $self->debug_message('Calculating gene/transcript level expression values');
    my $genes_ref = \%genes;
    my $transcripts_ref = \%transcripts;
    my %lists;
    $lists{'transcripts'}{features} = $transcripts_ref;
    $lists{'transcripts'}{outfile} = $transcript_expression_file;
    $lists{'genes'}{features} = $genes_ref;
    $lists{'genes'}{outfile} = $gene_expression_file;
    
    foreach my $feature_type (keys %lists){
        my $features = $lists{$feature_type}{features};
        my $feature_count = keys %{$features};
        $self->debug_message('Processing '. $feature_count .' features - removing features with no junctions');

        foreach my $fid (keys %{$features}){
            my $junction_list = $features->{$fid}->{junction_list};
            my $known_junction_count = keys %{$junction_list};
      
            #print "\n\tknown_junction count: $known_junction_count";
            if ($known_junction_count > 0){
                $features->{$fid}->{known_junction_count} = $known_junction_count;

                #values to calculate
                my $feature_read_count = 0;
                my $junctions_1x = 0;
                my $junctions_2x = 0;
                my $junctions_5x = 0;
                my $junctions_10x = 0;
                my $junctions_20x = 0;
                my $junctions_50x = 0;
                my $junctions_100x = 0;
                my $junctions_500x = 0;
                my $junctions_1000x = 0;
                foreach my $jid (keys %{$junction_list}){
                    my $read_count = 0;
                    if ($observed_junctions{$jid}){
                        $read_count = $observed_junctions{$jid}->{score};
                        $feature_read_count += $read_count;
                        if ($read_count >= 1){$junctions_1x++;}
                        if ($read_count >= 2){$junctions_2x++;}
                        if ($read_count >= 5){$junctions_5x++;}
                        if ($read_count >= 10){$junctions_10x++;}
                        if ($read_count >= 20){$junctions_20x++;}
                        if ($read_count >= 50){$junctions_50x++;}
                        if ($read_count >= 100){$junctions_100x++;}
                        if ($read_count >= 500){$junctions_500x++;}
                        if ($read_count >= 1000){$junctions_1000x++;}
                    }
                }
                my $junctions_1x_p = sprintf ("%.2f", (($junctions_1x/$known_junction_count)*100));
                my $junctions_2x_p = sprintf ("%.2f", (($junctions_2x/$known_junction_count)*100));
                my $junctions_5x_p = sprintf ("%.2f", (($junctions_5x/$known_junction_count)*100));
                my $junctions_10x_p = sprintf ("%.2f", (($junctions_10x/$known_junction_count)*100));
                my $junctions_20x_p = sprintf ("%.2f", (($junctions_20x/$known_junction_count)*100));
                my $junctions_50x_p = sprintf ("%.2f", (($junctions_50x/$known_junction_count)*100));
                my $junctions_100x_p = sprintf ("%.2f", (($junctions_100x/$known_junction_count)*100));
                my $junctions_500x_p = sprintf ("%.2f", (($junctions_500x/$known_junction_count)*100));
                my $junctions_1000x_p = sprintf ("%.2f", (($junctions_1000x/$known_junction_count)*100));
                
                #calculate expression level
                #    - Gene/transcript level expression = (sum of exon junction read counts for a gene / number of exon-exon junctions of that gene) => then normalized to per million junction mapped reads (JPJM)
                my $jpj = ($feature_read_count / $known_junction_count); #Junction count normalized to number of junctions for the feature
                my $jpjm = $jpj * (1000000 / $grand_read_count); #Further normalized to be per million junction mapping reads
                
                #Store values for this feature
                $features->{$fid}->{read_count} = $feature_read_count;
                $features->{$fid}->{jpj} = $jpj;
                $features->{$fid}->{jpjm} = $jpjm;
                $features->{$fid}->{junctions_1x_p} = $junctions_1x_p;
                $features->{$fid}->{junctions_2x_p} = $junctions_2x_p;
                $features->{$fid}->{junctions_5x_p} = $junctions_5x_p;
                $features->{$fid}->{junctions_10x_p} = $junctions_10x_p;
                $features->{$fid}->{junctions_20x_p} = $junctions_20x_p;
                $features->{$fid}->{junctions_50x_p} = $junctions_50x_p;
                $features->{$fid}->{junctions_100x_p} = $junctions_100x_p;
                $features->{$fid}->{junctions_500x_p} = $junctions_500x_p;
                $features->{$fid}->{junctions_1000x_p} = $junctions_1000x_p;
            }else{
                delete $features->{$fid};
            }
        }
    }

    #5.) Create output files, one for gene-level and one for transcript-level
    $self->debug_message('Creating output files');
    foreach my $feature_type (keys %lists){
        my $features = $lists{$feature_type}{features};
        my $outfile = $lists{$feature_type}{outfile};
        my $feature_count = keys %{$features};

        $self->debug_message('Processing '. $feature_count .' features and printing to: '. $outfile);
        open (OUT, ">$outfile") || die ("Could not open output file: $outfile");
        print OUT "fid\tensg_id\tgene_name\tchromosome\tknown_junction_count\tread_count\tjpjm\tjunctions_1x_p\tjunctions_2x_p\tjunctions_5x_p\tjunctions_10x_p\tjunctions_20x_p\tjunctions_50x_p\tjunctions_100x_p\tjunctions_500x_p\tjunctions_1000x_p\n";
        foreach my $fid (sort {$features->{$b}->{jpjm} <=> $features->{$a}->{jpjm}} keys %{$features}){
            my $gene_name = $features->{$fid}->{name};
            print OUT "$fid\t$features->{$fid}->{ensg_id}\t$gene_name\t$features->{$fid}->{chromosome}\t$features->{$fid}->{known_junction_count}\t$features->{$fid}->{read_count}\t$features->{$fid}->{jpjm}\t$features->{$fid}->{junctions_1x_p}\t$features->{$fid}->{junctions_2x_p}\t$features->{$fid}->{junctions_5x_p}\t$features->{$fid}->{junctions_10x_p}\t$features->{$fid}->{junctions_20x_p}\t$features->{$fid}->{junctions_50x_p}\t$features->{$fid}->{junctions_100x_p}\t$features->{$fid}->{junctions_500x_p}\t$features->{$fid}->{junctions_1000x_p}\n";
        }
        close(OUT);

    }

    return 1;
}

