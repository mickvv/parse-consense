#!/usr/bin/env perl 

use Modern::Perl '2011';
use autodie;

use Smart::Comments '###';
use Getopt::Euclid qw( :vars );

use File::Basename;
use Path::Class 'file', 'dir';

use Bio::MUST::Core;
use aliased 'Bio::MUST::Core::SeqId';
use aliased 'Bio::MUST::Core::IdList';
use aliased 'Bio::MUST::Core::Taxonomy';
use aliased 'Bio::MUST::Core::Taxonomy::Filter';
use aliased 'Bio::MUST::Core::Taxonomy::Criterion';
use aliased 'Bio::MUST::Core::Taxonomy::Category';
use aliased 'Bio::MUST::Core::Taxonomy::Classifier';

# build taxonomy objects
my $tax = Taxonomy->new( tax_dir => $ARGV_taxdir );

### Processing OTUs...
open my $in, '<', file($ARGV_otu_file);

# build classifier from labels
my @categories;

while ( my $line = <$in> ) {    
    chomp $line;
    my ($label, $otu) = split ':', $line;

    my $list      = IdList->new( ids => [ split ',', $otu ] );
    my $filter    = $tax->tax_filter( $list );
    my $criterion = Criterion->new( tax_filter => $filter );
    my $category  = Category->new(
        label    => $label,
        criteria => [ $criterion ],
    );
    push @categories, $category;
}
close $in;

my $classifier = Classifier->new( categories => \@categories );
# ### $classifier


# declare main data structures
my %data_for;
#my @tests = ();
#my @trees = ();
my (%TP_for, %FP_for, %TN_for, %FN_for);


# list directories
my @series_dir = File::Find::Rule
    ->directory()
    ->maxdepth(1)
    ->relative()
    ->name( "*$ARGV_dir_name*" )
    ->in($ARGV_series_dir)
    ;
 ### @series_dir

for my $serie_dir (@series_dir) {
    ### Processing: $serie_dir

    # Process consense files
    # ######################

    my $in_strip = '-' . $serie_dir;
    ### $in_strip

    my @consfiles = File::Find::Rule
        ->file()
        ->name( qr{ \. consense \z}xmsi )
        ->in($serie_dir)
    ;
    ### @consfiles


    # Parse output files from consense...
    # ...data structures:
    my @results;
    my $spec_hash;
    my $spec_n = 0;
    my $factor = 0.0;

    INFILE:
    for my $infile (@consfiles) {


        #-------------------------------------------------------------------------------
        # custom filter to remove when final run
        #-------------------------------------------------------------------------------
        next unless $infile =~ m/copy\-20/xmsg;
#        next unless $infile =~ m/OG\d+599/xmsg;
        ### Processing: $infile
     
        # counter for bipartitions
        my $bip;
        my %deja_vu;
     
        my ($OG, $serie, $ref_mul) = $infile =~ m/ (OG\d+) ([\w\-]+) - (\d \.? \d?)/xmsg;
        my ($pres, $copies) = get_batch($infile);
        $serie  =~ tr/-//d;
        ### $OG
        ### $serie
        ### $pres
        ### $copies
        ### $ref_mul
     
        
        # build hashes of species and bootstrap values using consense output
        my %boot_hash;
        my %delta_for;
     
        open my $in, '<', $infile;
        LINE:
        while (my $line = <$in>) {
            chomp ($line);
            next LINE if $line =~ m/^$/xmsg;
            next LINE if $line =~ m/^ Species | Consensus | Set /xmsg;
            
            # compute BP normalization factor
            if ($line =~ /How many times out of\s+([0-9\.]+)/) {
                $factor = 100.0 / $1;
            }
            
            # process species
            if ($line =~ /^\s*(\d+)\.\s+(\S+.*)$/) {
                my $index = $1;
                my $spec = $2;
                $spec =~ s/_{2,}//g;
                $spec =~ s/ /_/g;			# replace spaces by underscores
                $spec_hash->{$index} = $spec;
                next LINE;
            }
            
            # process bipartitions
            my @indexes;
            if ($line =~ /^((\.|\*).*?)(\d+\.\d+)$/) {
                my $boot = $3;
                # skip bipartitions with a bootstrap <= 50
                next LINE unless $boot > 50; 
                
                my $bips = $1;
                   $bips =~ s/ //g;			# remove spaces
                
                my $dots  = $bips =~ tr/.//; # compute bipartition size
                my $stars = $bips =~ tr/*//; # compute bipartition size
                my $delta = abs $dots - $stars;
                $delta_for{$bips} = $delta;
                
                next LINE;
            }
            ### %delta_for
            last LINE if $line =~ m/^Extended/xmsg;
        }
        close $in;
     
        # sort bipartitions from the most balanced to the less
        my @bips = sort { $delta_for{$a} <=> $delta_for{$b} } keys %delta_for;
        # my @vals = @delta_for{@keys};
        BIPS:
        for my $bip ( @bips ) {
            ### Computing: $bip
            # compute indexes for both partitions
            my $dots_idxs  = char_idxs($bip, '.');
            my $stars_idxs = char_idxs($bip, '*');
            
            MONO:
            for my $idxs ( ($dots_idxs, $stars_idxs) ) {
                ### Computing monopartition
                # identify species in %spec_hash and create corresponding SeqId objects  
                ### $idxs
                my $size = @$idxs;
                ### $size
                my @seq_ids = map { SeqId->new( full_id => $_ ) } 
                              map { $spec_hash->{$_} } @$idxs
                            ;
#                ### @seq_ids
                ### seq_ids: map { $_->full_id } @seq_ids
                
                my @groups = map { $classifier->classify($_) } @seq_ids;
                # skip partitions not containing any org from the considered 
                # group (i.e. orgs belonging to removed group)
                next MONO unless grep { $_ eq $serie } @groups; 
                
                # process partitions of size 2
                if ( $size == 2 ) {
                    # identify scenario # test if one, both or no seqs are from 
                    # genomic data and if groups are identical or not 
               
                    # skip partition with already processed organisms
                    next MONO unless grep { ! $deja_vu{$_->full_id} } @seq_ids;

                    my $is_genomic = _is_genomic(@seq_ids);
                    my $groups_st  = _groups_st(@groups);
                    ### scenario: [$is_genomic, $groups_st]
                    
                    if ($is_genomic eq 'none') {
                     
                        my $inc = $groups_st eq 'same' ? 2 : 1;
                        $data_for{$OG}{$serie}{$pres}{$copies}{$ref_mul}{FN} += $inc;  
                        ### $inc
                        
                        $deja_vu{$_->full_id}++ for @seq_ids; 
                        
                        ### %data_for
                        next MONO;
                    }
                 
                    $data_for{$OG}{$serie}{$pres}{$copies}{$ref_mul}{TP}++    if $is_genomic eq 'one'  && $groups_st eq 'same';
                    $data_for{$OG}{$serie}{$pres}{$copies}{$ref_mul}{FP}++    if $is_genomic eq 'one'  && $groups_st eq 'diff';
                    $data_for{$OG}{$serie}{$pres}{$copies}{$ref_mul}{FP} += 2 if $is_genomic eq 'both' && $groups_st eq 'same';
                    
                    $deja_vu{$_->full_id}++ for @seq_ids; 
                  
                    ### %data_for
                    next MONO;
                }
               
                # skip partition with already processed organisms
                next MONO if grep { $deja_vu{$_->full_id} } @seq_ids;
                
                # skip partition if heterogeneous in term of phylum
                next MONO if grep { $_ ne $serie } @groups; 
                ### $serie
                ### @groups
                
                my @orgs = map { $_->org } @seq_ids;
                my $ref_org = $orgs[0];
                ### @orgs
                ### $ref_org
                # skip partition if not unique species
                next MONO if grep { $_ ne $ref_org } @orgs;
                say 'YOU WIN ! \o/';
                my @isnt_new;
                @isnt_new = grep { ! $_->is_new } @seq_ids;
                ### isnt_new: map { $_->full_id } @isnt_new
                my $inc = @isnt_new;
                ### $inc
                $data_for{$OG}{$serie}{$pres}{$copies}{$ref_mul}{TP} += $inc;
                ### %data_for
                
                $deja_vu{$_->full_id}++ for @seq_ids; 
            }
        }
        ### %data_for
        ### %deja_vu
    }
}

#            my @dots_seq_ids  = map { SeqId->new( full_id => $_ ) } 
#                                map { $spec_hash->{$_} } @indexes
#                              ;
#            my @stars_seq_ids = map { SeqId->new( full_id => $_ ) } 
#                                map { $spec_hash->{$_} } @indexes
#                              ;
#            ### $dots_idxs
#            ### $stars_idxs
##        ### @seq_ids
##        ### seq_ids: map { $_->full_id } @seq_ids
#
#            my @dots_groups = map { $classifier->classify($_) } @dots_seq_ids;
#            next PART unless grep { $_ eq $serie } @dots_groups;
#            ### @groups
#            ### $serie




                


            # process a bipartition line
#            my @indexes = ();
#            if ($line =~ /^((\.|\*).*?)(\d+\.\d+)$/) {
#                my $boot = $3;
#                # skip bipartitions with a bootstrap <= 50
#                next LINE unless $boot > 50; 
#
#                my $bips = $1;
#                   $bips =~ s/ //g;			# remove spaces
#                my $size = $bips =~ tr/*//; # compute bipartition size
#                next unless $size == 2;     # consider only bipartitions of size = 2  
##                ### $bips
##                ### $size
#
#                # compute indexes for bipartition to identify species in %spec_hash
#                my $offset = 0;
#                my $star   = '*';
#                my $index  = index($bips, $star, $offset);
#                push @indexes, ++$index;
#                $offset = $index;
#                $index  = index($bips, $star, $offset);
#                push @indexes, ++$index;
##               ### @indexes
#
#                # identify species thanks to indexes and create a SeqId object on 
#                # the fly
#                my @seq_ids = map { SeqId->new( full_id => $_ )} 
#                              map { $spec_hash->{$_} } @indexes
#                              ;
##               ### @seq_ids
#               ### seq_ids: map { $_->full_id } @seq_ids
#                
#                # identify taxonomic groups and continue only if at least one is the 
#                # same as the removed seqs
#                my @groups = map { $classifier->classify($_) } @seq_ids;
#                next unless grep { $_ eq $serie } @groups;
##               ### @groups
##               ### $serie
#                $bip++;
#
#                # identify scenario
#                # test if one, both or no seqs are from genomic data...
#                # .. and if groups are identical or not
#                my $is_genomic = _is_genomic(@seq_ids);
#                my $groups_st = _groups_st(@groups);
#                ### scenario: [$is_genomic, $groups_st]
#                if ($is_genomic eq 'none') {
#
#                    my $rep = $groups_st eq 'same' ? 2 : 1;
#                    $data_for{$OG}{$serie}{$pres}{$copies}{$ref_mul}{FN}++ for 1 .. $rep;  
##                    ### FN: $data_for{$OG}{$serie}{$pres}{$copies}{$ref_mul}{FN}
#                    next LINE;
#                }
#             
#                $data_for{$OG}{$serie}{$pres}{$copies}{$ref_mul}{TP}++ if $is_genomic eq 'one'  && $groups_st eq 'same';
#                $data_for{$OG}{$serie}{$pres}{$copies}{$ref_mul}{FP}++ if $is_genomic eq 'one'  && $groups_st eq 'diff';
#                $data_for{$OG}{$serie}{$pres}{$copies}{$ref_mul}{FP} += 2 if $is_genomic eq 'both' && $groups_st eq 'same';
#            }
#    ### Processed bipartitions: $bip

        #pour l'arbre X, for each org à ajouter
        #for each bipar 
        #    if une des deux faces contient les org d'intérêts
        #        next if limite = bipar de taille 2 
        #        next unless org intérêt
        #        next if bp > 50
        #        if seq n°2= homologue génomique m(ême org)? -> TP
        #                  = autre génomique (pas le même)-> FP
        #                  = autre prot -> FN
        #        FP -> si j'ai rajouté une génomique au mauvais endroit
        #
        #        TN -> corriger les calculs (limiter le pool aux séqs homologues)


    ### Writing output...
    open my $out, '>', 'outfile.test';

    say {$out} join "\t", 'OG', 'serie', 'pres', 'copies', 'ref_mul', 'TPs', 'FPs', 'FNs', 'precision', 'recall';

    for my $OG      ( sort keys %data_for ) {
    for my $serie   ( sort keys %{ $data_for{$OG} } ) {
    for my $pres    ( sort keys %{ $data_for{$OG}{$serie} } ) {
    for my $copies  ( sort keys %{ $data_for{$OG}{$serie}{$pres} } ) {
    for my $ref_mul ( sort keys %{ $data_for{$OG}{$serie}{$pres}{$copies} } ) {

#        my @seqs   = map { $data_for{$OG}{$serie}{$pres}{$copies}{$ref_mul}{$_}{seqs}   // '0'  } @all_orgs;
#        my @size   = map { $seqs_for{$_}                                                // '0'  } @all_orgs;
#        my @p_seqs = map { $data_for{$OG}{$serie}{$pres}{$copies}{$ref_mul}{$_}{p_seqs} // '0'  } @all_orgs;
#        my @mask   = map { $data_for{$OG}{$serie}{$pres}{$copies}{$ref_mul}{$_}{mask}   // '0'  } @all_orgs;
#        my @TPs    = map { $data_for{$OG}{$serie}{$pres}{$copies}{$ref_mul}{$_}{TP}     // '0'  } @all_orgs;
#        my @FPs    = map { $data_for{$OG}{$serie}{$pres}{$copies}{$ref_mul}{$_}{FP}     // '0'  } @all_orgs;
#        my @TNs    = map { $data_for{$OG}{$serie}{$pres}{$copies}{$ref_mul}{$_}{TN}     // '0'  } @all_orgs;
#        my @FNs    = map { $data_for{$OG}{$serie}{$pres}{$copies}{$ref_mul}{$_}{FN}     // '0'  } @all_orgs;
        my $TPs = $data_for{$OG}{$serie}{$pres}{$copies}{$ref_mul}{TP} // '0';
        my $FPs = $data_for{$OG}{$serie}{$pres}{$copies}{$ref_mul}{FP} // '0';
        my $FNs = $data_for{$OG}{$serie}{$pres}{$copies}{$ref_mul}{FN} // '0';

#        my $seqs     = $e_data_for{$serie}{$OG}{'seqs'};
#        my $warnings = $e_data_for{$serie}{$OG}{'warnings'};
#        my $len      = $ali_len_for{$OG};
        
        #    # PR-curve : y=precision x=recall
        #    # ROC-curve : y=TPR x=FPR
        #    my @precision = # = TP/(TP+FP)
        #    my @recall    = # sensitivity = TPR = TP/P = TP/(TP+FN)

        my $precision = eval { $TPs/($TPs+$FPs) } // '0'; # use eval function to 
        my $recall    = eval { $TPs/($TPs+$FNs) } // '0'; # counter annoying undef
                                                           # when dividing by zero
#        say {$out} join "\t", $OG, $serie, $pres, $copies, $ref_mul, (sum @p_seqs), (sum @TPs), (sum @FPs), (sum @TNs), (sum @FNs), $precision, $recall, $fallout, $warnings, $seqs, $len;
        say {$out} join "\t", $OG, $serie, $pres, $copies, $ref_mul, $TPs, $FPs, $FNs, $precision, $recall;
     
                    }
                }
            } 
        }
    } 





sub _is_genomic {
    my ($seq_id_a, $seq_id_b) = @_;
    return 'both' if $seq_id_a->is_new && $seq_id_b->is_new;
    return 'one'  if $seq_id_a->is_new || $seq_id_b->is_new;
    return 'none';  
}
sub _groups_st {
    my ($group_a, $group_b) = @_;
    return 'same' if $group_a eq $group_b;
    return 'diff' if $group_a ne $group_b;
}
sub get_batch {
    my $infile = shift;

    my ($batch) = join "-", split "/", change_path($infile, undef, 2);
    return  (split "-", $batch)[1,5]; 
}
sub change_path {
    my $infile          = shift;
    my $new_directories = shift;
    my $level           = shift;

    my ($filename, $directories) = fileparse($infile);
    # path::class::file
    my @wanted_directories = File::Spec->catdir(splice [split "/", $directories], -$level)
        if $level;

    my $outfile = file($new_directories ? $new_directories : '.', 
                       $level ? @wanted_directories : '',
                       $filename
                       );
    return $outfile->stringify;
}
sub char_idxs {
    my @indexes;
    my $offset = 0;
    my ($bips, $char) = @_;

    my $index = index($bips, $char, $offset);
    push @indexes, ++$index;

    while ($index != 0) {

        $offset = $index;
        $index  = index($bips, $char, $offset);

        return \@indexes if $index == -1;

        push @indexes, ++$index;
    }

    return \@indexes;
}

#die ("ERROR! Empty or unrecognized input file\n") if ($spec_n == 0);
#
#print "# Species list [$spec_n]\n";
#foreach my $spec (
#	sort { $spec_hash->{$a} <=> $spec_hash->{$b} } keys %$spec_hash) {
#	printf "# %2d. %s\n", $spec_hash->{$spec}, $spec;
#}
#
#if ($mode ne "puzzle") {
#	my $boot_hash = $results[0];
#	print "# Bipartition list [" . scalar (keys %$boot_hash) . "]\n";
## 	foreach my $key (sort { $boot_hash->{$b} <=> $boot_hash->{$a} }
## 		keys %$boot_hash) {
## 		printf "# %s\t\%d\n", $key, $boot_hash->{$key};
## 	}
#}
#
## build complementary bipartitions
#foreach my $boot_hash (@results) {
#	my @keys = keys %$boot_hash;
#	foreach my $key (@keys) {
#		my $boot = $boot_hash->{$key};
#		$key =~ tr/\.\*/\*\./;
#		$boot_hash->{$key} = $boot;
#	}
#
## normalize BPs for all results
#	foreach my $key (keys %$boot_hash) {
#		my $boot = $boot_hash->{$key};
#		$boot *= $factor;
#		$boot_hash->{$key} = $boot;
#	}
#}
#
#### PART 3: Compute statistics and output results
#
#if ($mode ne "puzzle") {
#	my $boot_hash = $results[0];
#
#	# loop trough all groups and output bootstrap values
#	print "# Group statistics [" . scalar (@look_array) . "]\n";
#	foreach my $otu (@look_array) {
#		print "$otu\t";
## 		print $look_hash{$otu} . "\t";
#		# handle empty otus
#		if (defined ($empt_hash{$otu})) {
#			printf ("n.a.%s", $empt_hash{$otu});
#		}
#
#		else {
#			# handle combinations not found
#			if (!defined ($boot_hash->{$look_hash{$otu}})
#				|| $boot_hash->{$look_hash{$otu}} < $threshold) {
#				print "-";
#			}
#			# handle combinations found
#			else {
#				# print bootstrap
#				printf ("%.1f", $boot_hash->{$look_hash{$otu}});
#			}
#			# add one '#' for each missing group from current otu
#			print ("#" x $miss_hash{$otu}) if (defined ($miss_hash{$otu}));
#		}
#		print "\n";
#	}
#}
#
#else {
#
#	# setup topology filter
#	my @candidates = ();
#	for (my $i = 0; $i < @results; $i++) {
#		push @candidates, $i;
#	}
#
#	# loop trough all groups and filter topologies
#	foreach my $otu (@look_array) {
#
#		# loop through all topologies
#		my @newcomers = ();
#		for (my $i = 0; $i < @results; $i++) {
#			my $boot_hash = $results[$i];
#
#			# keep topology if group exists or is uninformative
#			if (defined ($boot_hash->{$look_hash{$otu}}) or
#				(defined ($empt_hash{$otu}) && $empt_hash{$otu} eq 'U')) {
#				$tests[$i] =~ /^\s*(\d+)/;
#				push @newcomers, ($1 - 1);
#			}
#		}
#
#		# keep topologies compatible with all groups
#		@candidates = @{get_intersection (\@candidates, \@newcomers)};
#	}
#
#	# extract gene and group names from infile name
#	my @subs = split (/\./, $infile);
#	@subs = split (/\-/, $subs[0]);
#	my $group = pop @subs;
#	my $gene = join ("-", @subs);
#
#	# check number of retained topologies
#	die ("# No topology retained for $gene-$group!\n") if (@candidates == 0);
#	die ("ERROR! More than one topology retained!\n") if (@candidates > 1);
#
#	# print results
#	my $id = shift @candidates;
#	printf "# gene   spec/group      tree   log L   difference    S.E."
#		. "      p-1sKH     p-SH       c-ELW      2sKH     topology\n";
#	printf "%-8s %-15s %-80s %s\n", $gene, $group, $tests[$id], $trees[$id];
#}
#
#
#
#sub get_intersection  {
#	my ($array1, $array2) = @_;
#
#    my @intersection = ();
#    my @difference = ();
#    my %count = ();
#    foreach my $element (@$array1, @$array2) {
#    	$count{$element}++;
#    }
#    foreach my $element (keys %count) {
#		push @{$count{$element} > 1 ? \@intersection : \@difference }, $element;
#    }
#
#	return \@intersection;
#}
#
#
#
#sub get_topology {
#	my ($tree) = @_;
#
#	$tree =~ s/_{2,}//g;				# delete series of underscores
#	$tree =~ s/\)[0-9Ee\.\-]+/)/g;		# delete BPs
#	$tree =~ s/:[0-9Ee\.\-]+//g;		# delete branch lengths
#
#	if ($tree =~ /.*?(\(.*?\);)/) {		# not fool-proof for bad trees
#		$tree = $1;
#	} else {
#		die ($badtree);
#	}
#
#	return $tree;
#}
#
#
#
#sub get_species_hash {
#	my ($topology) = @_;
#
#	# extract species list
#	$topology =~ s/[\(\)\;]+//g;
#	my @specs = split (',', $topology);
#
#	# build species hash (numbered from 1 to spec_n)
#	my $spec_n = 0;
#	my $hash = {};
#	foreach my $spec (@specs) {
#		$spec =~ s/^\s+|\s+$//g;	# chop leading and trailing spaces
#		$hash->{$spec} = ++$spec_n;
#	}
#
#	return $hash;
#}



=head1 NAME

parse_consense.pl

=head1 VERSION

This documentation refers to parse_consense.pl version 0.0.1

=head1 USAGE

parse_consense.pl --report-dir=<dir> --otu[-file]=<file> --taxdir=<dir>

=head1 REQUIRED ARGUMENTS

=over

=item --series-dir=<dir>

Path to input ALI files [repeatable argument].

=for Euclid: 
    dir.type: str

#=item <infiles>
#
#Path to input consense files [repeatable argument].
#
#=for Euclid: 
#    infiles.type: readable
#    repeatable

=item --taxdir=<dir>

Path to local NCBI taxonomy DB.

=for Euclid:
    dir.type: string

=item --otu[-file]=<file>

Path to artificial groups' file (user defined).

=for Euclid:
    file.type: string

#=item --score=<str>
#
#Choose either percent_identity or bit_score to compute delta
#
#=for Euclid:
#    str.type: str

=back

=head1 OPTIONAL ARGUMENTS

=over

=item --dir-name=<dir>

Specify dirs name.

=for Euclid:
    dir.type: str
    dir.default: ''

=item --hsp[-length]=<int>

Filter hits according to a min hsp length.

=for Euclid:
    int.type: int
    int.default: 0

=item --evalue=<num>

Filter hits according to an e-value threshold.

=for Euclid:
    num.type: number
    num.default: 1e-10

=item --perc-id=<int>

Filter hits according to a percent_identity threshold.
 
=for Euclid:
    int.type: int
    int.default: 0
    
=back
