""" Helper functions for processing the t7 sequencing
"""

import sys
import re
import itertools # Used to do pairwise comparisons of umis for similarity

from collections import defaultdict

import numpy as np
import pandas as pd

from pysam.libcalignmentfile import AlignmentFile

import warnings
from Bio import BiopythonDeprecationWarning

# Suppress only BiopythonDeprecationWarning
warnings.simplefilter("ignore", BiopythonDeprecationWarning)
from Bio import pairwise2  # pip install biopython

def get_t7_count_matrix(cell_barcodes_dict, canonical_to_edits, canonical_edit_bc_cell_umis_filtered):
    """ Performs the UMI counting per edit site, either for the barcoded or non-barcoded reads!
    """

    t7_barcoded_counts = np.zeros((len(cell_barcodes_dict), len(canonical_to_edits)),
                                  dtype=np.uint16)

    #edit_bc_cell_umis_filtered # edit to cell barcodes to umis
    #canonical_to_edits # canonical-edits to edits
    
    # REQUIRE DICT SINCE IT MAINTAINS ORDER
    # cell_barcodes_dict = {key: value for value, key in enumerate(cell_barcodes_list)}
    
    edit_site_names = []
    for edit_site_index, (edit_site, edits) in enumerate(canonical_to_edits.items()):

        edit_site_names.append( f'{edit_site.chrom}:{edit_site.ref_pos}' )

        # for edit in edits:
        #
        # NOTE edit_bc_cell_umis_filtered[edit] MUST NOW BE A DICTIONARY THAT STORES SETS AS VALUES
        t7_barcoded_counts[:, edit_site_index] += get_cell_counts_from_umi_dict(
                                                    canonical_edit_bc_cell_umis_filtered[edit_site], cell_barcodes_dict)

    t7_barcoded_counts = pd.DataFrame(t7_barcoded_counts, index=list(cell_barcodes_dict.keys()), columns=edit_site_names)

    return t7_barcoded_counts


def get_cell_counts_from_umi_dict(cell_bc_to_umis, cell_barcodes_dict):
    """Gets the counts across cells for a given barcode to umi dict
    """
    cell_umi_counts = np.zeros((len(cell_barcodes_dict)), dtype=np.uint16)

    for cell_barcode, cell_index in cell_barcodes_dict.items():

        umi_set = cell_bc_to_umis[cell_barcode]

        # umi_set = set(umis) # stored as set automatically now
        num_unique_umi = len(umi_set)

        if num_unique_umi == 1:
            cell_umi_counts[cell_index] += 1
        elif num_unique_umi == 2:
            umi_1, umi_2 = umi_set
            
            if within_single_mismatch(umi_1, umi_2):
                cell_umi_counts[cell_index] += 1
            else:
                cell_umi_counts[cell_index] += 2
        else:

            # Kinda complex, fix if it becomes an issue later
            umis_to_match_umis = {umi: {umi} for umi in umi_set}
            
            for umi_1, umi_2 in itertools.combinations(umi_set, 2):

                if within_single_mismatch(umi_1, umi_2):
                    umis_to_match_umis[umi_1] = umis_to_match_umis[umi_1].union(umis_to_match_umis[umi_2])
                    umis_to_match_umis[umi_2] = umis_to_match_umis[umi_2].union(umis_to_match_umis[umi_1])

                    # syncing!
                    # Need to also pull in the umis with edit distance 1 that are also edit distance 1 from
                    # this umi!
                    for umi in umis_to_match_umis[umi_1]:
                        umis_to_match_umis[umi_1] = umis_to_match_umis[umi_1].union(umis_to_match_umis[umi])
                        umis_to_match_umis[umi_2] = umis_to_match_umis[umi_2].union(umis_to_match_umis[umi])

                    # Now need to update all the other umis with these umis edit dist 1
                    for umi in umis_to_match_umis[umi_1]:
                        umis_to_match_umis[umi] = umis_to_match_umis[umi].union(umis_to_match_umis[umi_1])

            unique_umi_sets = []
            for umi_set_ in umis_to_match_umis.values():
                if umi_set_ not in unique_umi_sets:
                    unique_umi_sets.append(umi_set_)
            
            cell_umi_counts[cell_index] += len(unique_umi_sets)

    return cell_umi_counts

def bam_count_gene_umis(bam_file, cell_barcodes_dict, gene_names, allt7_reads=None,
                        max_reads=None, # Just for testing purposes...
                        ):
    """Get's gene UMI counts per cell, returning as pandas dataframe of counts!
    """
    
    # I SHOULD REWRITE THIS WHOLE THING WHEN I GET A CHANCE
    
    gene_set = set(gene_names)

    # genes_to_bcs_to_umis = defaultdict(dict)
    genes_to_bcs_to_umis = defaultdict(lambda: defaultdict(set))
    
    # If the bam has not reads for a given region, just doesn't iterate
    with AlignmentFile(bam_file, "rb") as bam:
        nreads = 0
        for i, read in enumerate( bam ):
            nreads += 1

            cell_barcode = read.get_tag('CB')
            
            # Using filt bam means I shouldnt have to check against allt7
            if (cell_barcode not in cell_barcodes_dict):
                continue

            gene_name_tag = read.get_tag('GN')
            gene_id_tag = read.get_tag('GX')

            if (gene_name_tag != '') and (gene_name_tag in gene_set):
                gene_or_id = gene_name_tag
            elif (gene_id_tag != '') and (gene_id_tag in gene_set):
                gene_or_id = gene_id_tag
            else:
                continue

            # Should be cleaner and faster
            genes_to_bcs_to_umis[gene_or_id][cell_barcode].add(read.get_tag('pN'))

    gene_names_list = list(gene_names)
    gene_names_dict = {key: value for value, key in enumerate(gene_names_list)}
    cell_by_gene_umi_counts = np.zeros((len(cell_barcodes_dict), len(gene_names)), dtype=np.uint16)
    for gene in genes_to_bcs_to_umis:
        gene_index = gene_names_dict[gene]
        cell_by_gene_umi_counts[:, gene_index] = get_cell_counts_from_umi_dict(genes_to_bcs_to_umis[gene], cell_barcodes_dict)

    cell_by_gene_umi_counts = pd.DataFrame(cell_by_gene_umi_counts,
                                           index=list(cell_barcodes_dict.keys()),
                                           columns=gene_names_list)

    return cell_by_gene_umi_counts


def within_single_mismatch(seq1, seq2):  
    diff_score = 0
    for i1, i2 in zip(seq1, seq2):
        if i1 != i2:
            if diff_score > 0:
                return False
            diff_score+=1
    return True

def bio_edit_distance(seqA, seqB, start_from_first_smallest_seq_aln=True, alns_to_compare=None):
    """More biologically relevant edit distance, which does a global alignment between two strings, then takes the
        edit distance as the number of gaps introduced minus the difference in sequence length. Will result
        in a 0 score if one sequence is a subset of the other.
    """
    ### Puts the sequences together better, at the 5' end
    aln = pairwise2.align.localms(seqA, seqB,
                                  1, # score for match
                                  -1, # mismatch penalty
                                  -.5, # gap-open penalty
                                  -.5, # gap-extension penalty
                                  one_alignment_only=True)[0]

    seq_lens = [len(seqA), len(seqB)]
    shorter_seq_index = np.argmin(seq_lens)
    #### Counting mismatches till all of shorter sequence accounted for
    aln_seqs = [aln.seqA, aln.seqB]
    n_seen = 0
    edit_dist = 0
    start_count = start_from_first_smallest_seq_aln == False # Only start counting once the alignment of the shorter sequence starts...
    if type(alns_to_compare)==type(None):
        alns_to_compare = len(aln.seqA) # Compare all of the global alignment, not just from one-side
    else: # Min between the max alignment length and the requested number to compare.
        alns_to_compare = min([len(aln.seqA), alns_to_compare])

    for i in range( alns_to_compare ):

        if start_count and aln.seqA[i] != aln.seqB[i]:
            edit_dist += 1

        if aln_seqs[shorter_seq_index][i] != '-': # Aligned character
            n_seen += 1
            start_count = True

            if n_seen == seq_lens[shorter_seq_index]: # seen all characters
                break

    return edit_dist

def get_edit_sets(edit_set):
    """Gets set of edits, version which did not considered subsequence problem, where partial overlap of subsequences
        between larger t7 edits results in very bad over-estimation of allelic number.
    """
    # Keep track of which edits look like each other, so only have a count of 1 for these edits.
    # This de-duplicates for edits that might represent slightly offset t7 priming, which can be quite common.
    # Basically the criteria is that the end of the read (forward read) or start of the read (reverse read)
    # is different between the two edit sites, indicating extra sequence content not present in the other.
    # Otherwise, could indicate different priming / read truncation due to the library fragmentation.

    edits_to_edit_matches = {edit: {edit} for edit in edit_set}  # initialise with just a set of itself

    for edit_1, edit_2 in itertools.combinations(edit_set, 2):

        # If they map in opposite directions, definitely different alleles edited
        # (donor-insert went in different direction at both alleles)
        # Also different alleles if the genomic proportion of the read is different, indicates different insert site
        reflen = len(edit_1.ref_seq)
        if (edit_1.forward != edit_2.forward) or (edit_1.ref_pos != edit_2.ref_pos):
            continue

        # Both forwards, so could be from the same allele with different variation in start site
        elif edit_1.forward and edit_2.forward:  # Reverse the direction to get the end of transcription
            edit_1_seq = edit_1.alt_seq[::-1][0:-reflen]  # Cutting off the added reference sequence
            edit_2_seq = edit_2.alt_seq[::-1][0:-reflen]

        elif edit_2.forward == False and edit_2.forward == False:
            edit_1_seq = edit_1.alt_seq[0:-reflen]  # Already in opposite direction of transcription
            edit_2_seq = edit_2.alt_seq[0:-reflen]

        #### The same insert at a different position is clear case of different insert
        if edit_1.alt_seq == edit_2.alt_seq and edit_1.ref_pos != edit_2.ref_pos:
            continue

        # Getting the common sequence at the end of the read in the direction of transcription
        # edit_seqs = [edit_1_seq, edit_2_seq]

        ##### Original method, does not consider sequencing error in either edit, below version does, even
        ####    for cases where the sequences are of different length!
        # longest_common_seq = os.path.commonprefix( edit_seqs )
        # # if the longest common sequence between the two is equal to one of the sequences, and ref
        # # position is the same, indicating the genomic part of the read is the same, the one sequence
        # # is likely a truncated read generated of the other, but generated from the same allele
        # if longest_common_seq in edit_seqs:
        #     edits_to_edit_matches[edit_1] = edits_to_edit_matches[edit_1].union({edit_2})
        #     edits_to_edit_matches[edit_2] = edits_to_edit_matches[edit_2].union({edit_1})
        #     continue

        # dist_between_seqs = nltk.edit_distance(edit_1_seq, edit_2_seq) - abs(len(edit_1_seq) - len(edit_2_seq))

        # The levenshtein distance above can under-estimate the difference between biological sequences, because
        # mutates them in ways that don't quite make biological sense which resulted in edit distances that
        # didn't make sense in some cases.
        # Below does a global alignment and gets the gap count in both aligned sequences as the edit distance,
        # subtracting the expected number of gaps due to differences in sequence length.
        dist_between_seqs = bio_edit_distance(edit_1_seq, edit_2_seq)

        if dist_between_seqs > 1:
            # HOMO-OLIGO sequencing error correction.
            # Need to correct for homo-oligo sequencing error. Can see cases where it OVERESTIMATES the number
            # of alleles, where the only different BETWEEN the edits is the NUMBER of a particular base IF there are
            # consecutive versions of that base. Illumina sequences is prone to this error. SO if that is the source
            # of difference, will remove this consecuative run of bases and re-calculate edit distance.
            # Example seen during debugging:
            # ATAATACTCTCCCTATTCACTCTGCGT
            # ATAATACTCTCCCCCCATTCACTC
            # -> Aside from the consecutive runs of 'CC', there is only an edit distance of 1 between these two sequences.
            #  and the 'CC' is not accurate, so would remove these consecutive runs and re-calculate the edit distance.
            # Find subsequences with more than 2 identical characters in a row, and replace them with one-copy of the
            # character BEFORE the edit-distance calculation between the sequences, to handle the homo-polymer
            # illumina sequencing error
            edit_1_consecs = re.findall(r"((.)\2{2,})", edit_1_seq)
            edit_2_consecs = re.findall(r"((.)\2{2,})", edit_2_seq)
            # Only a problem if one of them has a consecuative run of a particular base, the 'GGG' or 'CCC' of the
            # superb-seq barcode in particular seems to sometimes cause this issue, if is a non-canonical donor
            # sequence.
            if len(edit_1_consecs) > 0 or len(edit_2_consecs) > 0:
                # Example case am trying to catch here:
                # edit_1_seq: ATAATACTCTCCCCTTCACTCTGCGTTGATAC
                # edit_2_seq: ATAATACTCTCCCCCCATTCACT
                # edit_1_consecs: [('CCCC', 'C')]
                # edit_2_consecs: [('CCCCCC', 'C')]
                edit_1_seq_homo_fix = edit_1_seq
                for consec, base in edit_1_consecs:
                    #### I noticed that where there is a homopolymer, the next base
                    edit_1_seq_homo_fix = edit_1_seq_homo_fix.replace(consec, base)  # make it not a polymer

                edit_2_seq_homo_fix = edit_2_seq
                for consec, base in edit_2_consecs:
                    edit_2_seq_homo_fix = edit_2_seq_homo_fix.replace(consec, base)  # make it not a polymer
                # edit_1_seq_homo_fix: ATAATACTCTCTTCACTCTGCGTTGATAC
                # edit_2_seq_homo_fix: ATAATACTCTCATTCACT
                # the edit dist here is now:
                # nltk.edit_distance(edit_1_seq_homo_fix, edit_2_seq_homo_fix) - abs( len(edit_1_seq_homo_fix)-len(edit_2_seq_homo_fix) )
                # 2
                # due to the substitution cost, it prefers to flip two characters rather than create a gap,
                # which would yield a 1 edit distance. To enable this, will lower the substitution cost
                # dist_between_seqs = (nltk.edit_distance(edit_1_seq_homo_fix, edit_2_seq_homo_fix, substitution_cost=0.5) -
                #                      abs(len(edit_1_seq_homo_fix) - len(edit_2_seq_homo_fix)))
                dist_between_seqs = bio_edit_distance(edit_1_seq_homo_fix,
                                                      edit_2_seq_homo_fix) - 1  # -1 because the homopolymer extra error

        # Base-pair sequencing error correction
        # Some cases where the edits are VERY similar, off by 1bp, likely sequencing error!
        # Below is using the levenstein distance, so number of edits used to make one sequence into the other.
        # Taking the difference with the sequence lengths, to account for extra bp which would be added
        # to make the sequences the same if they are of different lengths. So would give 0 if 1 seq is a subset
        # of the other !
        if dist_between_seqs <= 1:  # Only one edit distance from each other, accounts for short indel at start also
            edits_to_edit_matches[edit_1] = edits_to_edit_matches[edit_1].union({edit_2})
            edits_to_edit_matches[edit_2] = edits_to_edit_matches[edit_2].union({edit_1})

    # Unique edits will be the number of unique sets of edits generated from the above dictionary!
    edit_sets = [set(edits_) for edits_ in edits_to_edit_matches.values()]
    # These unique edit sets will be inflated in terms of the unique edits, IF there are two edits that have
    # shared subsequence, with sub-edits shared between them. It causes them to be added to shared sets,
    # but they will be non-overlapping.
    unique_edit_sets = []
    [unique_edit_sets.append(edit_set_) for edit_set_ in edit_sets if edit_set_ not in unique_edit_sets]

    if len(unique_edit_sets) == 1:
        return unique_edit_sets
        # cell_allelic_edits[celli, edit_sitei] = len(unique_edit_sets)
        # continue

    # Need to de-duplicate these edit sets, since shared subsequences between longer edits can cause them to
    # create unique sets. So to de-duplicate, need to go through and REMOVE edit sets which represent
    # combinations of other edit sets!
    unique_edit_sets_deoverlapped = unique_edit_sets.copy()
    for edit_set_1, edit_set_2 in itertools.combinations(unique_edit_sets, 2):

        union_set = edit_set_1.union(edit_set_2)
        if union_set in unique_edit_sets_deoverlapped:
            unique_edit_sets_deoverlapped.pop(unique_edit_sets_deoverlapped.index(union_set))

    return unique_edit_sets_deoverlapped

def get_longest_edits(edit_set):
    """This method tries to get the longest t7 edits that are distinct from one another, overcomes subsequence problem
        of previous method.
    """
    # In this version, just for looking
    edits_to_subedits = {edit: set() for edit in edit_set}  # initialise with list

    ##### Will keep track of the longer edits in a set, where if two edits are matched, only store the longer
    ##### edit, and if the shorter edit is currently called as a longer edit from previous comparisons, pop
    ##### it from the long-edit list.
    edit_set = list(edit_set)
    edit_lens = np.array([len(edit.alt_seq) for edit in edit_set])
    edit_order = np.argsort(edit_lens)  # Will go from smallest to largest

    # By starting with the shortest edits, comparing to everything larger than it, the shorter sequences will initially
    # be called the longer edits, but if a longer match is found will downgrade to a sub-edit, with the final
    # set of long edits being the edits that could not be matched to a longer edit.

    longest_edits = []  # These will be our final edit set!
    sub_edits = [] # These are the edits that have already been called subedits, so don't add them as longer edit!
    for i, edit_1_i in enumerate(edit_order):
        edit_1 = edit_set[edit_1_i]

        for edit_2_i in edit_order[i+1:]: # Compare against all longer edits not yet seen
            edit_2 = edit_set[edit_2_i]

            # If they map in opposite directions, definitely different alleles edited
            # (donor-insert went in different direction at both alleles)
            # Also different alleles if the genomic proportion of the read is different, indicates different insert site
            reflen = len(edit_1.ref_seq)
            if (edit_1.forward != edit_2.forward) or (edit_1.ref_pos != edit_2.ref_pos):
                if edit_1 not in sub_edits: # Hasn't been logged as a subedit.
                    longest_edits.append( edit_1 )
                if edit_2 not in sub_edits:
                    longest_edits.append( edit_2 )

                continue

            # Both forwards, so could be from the same allele with different variation in start site
            elif edit_1.forward == False and edit_2.forward == False:  # Reverse the direction to get the end of transcription
                edit_1_seq = edit_1.alt_seq[reflen:]  # Cutting off the added reference sequence
                edit_2_seq = edit_2.alt_seq[reflen:]
                # Need to reverse the direciton, since the 'bio_edit_distance' function will ONLY count the edit distance
                # alignement correctly if the shorter subsequence maps to the FRONT of the longer subsequence.
                edit_1_seq = edit_1_seq[::-1]
                edit_2_seq = edit_2_seq[::-1]

            elif edit_2.forward and edit_2.forward:
                edit_1_seq = edit_1.alt_seq[0:-reflen]  # Already in opposite direction of transcription
                edit_2_seq = edit_2.alt_seq[0:-reflen]

            # Below does a global alignment and gets the gap count in both aligned sequences as the edit distance,
            # subtracting the expected number of gaps due to differences in sequence length.
            dist_between_seqs = bio_edit_distance(edit_1_seq, edit_2_seq)

            if dist_between_seqs > 1:
                # HOMO-OLIGO sequencing error correction.
                # Need to correct for homo-oligo sequencing error. Can see cases where it OVERESTIMATES the number
                # of alleles, where the only different BETWEEN the edits is the NUMBER of a particular base IF there are
                # consecutive versions of that base. Illumina sequences is prone to this error. SO if that is the source
                # of difference, will remove this consecuative run of bases and re-calculate edit distance.
                # Example seen during debugging:
                # ATAATACTCTCCCTATTCACTCTGCGT
                # ATAATACTCTCCCCCCATTCACTC
                # -> Aside from the consecutive runs of 'CC', there is only an edit distance of 1 between these two sequences.
                #  and the 'CC' is not accurate, so would remove these consecutive runs and re-calculate the edit distance.
                # Find subsequences with more than 2 identical characters in a row, and replace them with one-copy of the
                # character BEFORE the edit-distance calculation between the sequences, to handle the homo-polymer
                # illumina sequencing error
                edit_1_consecs = re.findall(r"((.)\2{2,})", edit_1_seq)
                edit_2_consecs = re.findall(r"((.)\2{2,})", edit_2_seq)
                # Only a problem if one of them has a consecuative run of a particular base, the 'GGG' or 'CCC' of the
                # superb-seq barcode in particular seems to sometimes cause this issue, if is a non-canonical donor
                # sequence.
                if len(edit_1_consecs) > 0 or len(edit_2_consecs) > 0:
                    # Example case am trying to catch here:
                    # edit_1_seq: ATAATACTCTCCCCTTCACTCTGCGTTGATAC
                    # edit_2_seq: ATAATACTCTCCCCCCATTCACT
                    # edit_1_consecs: [('CCCC', 'C')]
                    # edit_2_consecs: [('CCCCCC', 'C')]
                    edit_1_seq_homo_fix = edit_1_seq
                    for consec, base in edit_1_consecs:
                        #### I noticed that where there is a homopolymer, the next base
                        edit_1_seq_homo_fix = edit_1_seq_homo_fix.replace(consec, base)  # make it not a polymer

                    edit_2_seq_homo_fix = edit_2_seq
                    for consec, base in edit_2_consecs:
                        edit_2_seq_homo_fix = edit_2_seq_homo_fix.replace(consec, base)  # make it not a polymer
                    # edit_1_seq_homo_fix: ATAATACTCTCTTCACTCTGCGTTGATAC
                    # edit_2_seq_homo_fix: ATAATACTCTCATTCACT
                    # the edit dist here is now:
                    # nltk.edit_distance(edit_1_seq_homo_fix, edit_2_seq_homo_fix) - abs( len(edit_1_seq_homo_fix)-len(edit_2_seq_homo_fix) )
                    # 2
                    # due to the substitution cost, it prefers to flip two characters rather than create a gap,
                    # which would yield a 1 edit distance. To enable this, will lower the substitution cost
                    dist_between_seqs = bio_edit_distance(edit_1_seq_homo_fix, edit_2_seq_homo_fix) - 1  # -1 because the homopolymer adds extra error

            # They are different from one another in the sense that one is not a subset of the other (i.e. not
            # different due to different t7 primering).
            # Need to now consider if they only differ at the 5' end of the read, if so likely only different
            # due to sequencing adapters / TSO stil present at the end of the read!
            # Getting the edit distance from the first 5 base-pairs of the 3' end of the read
            if dist_between_seqs > 2:
                # They are different from one another in the sense that one is not a subset of the other (i.e. not
                # different due to different t7 primering).
                # Need to now consider if they only differ at the 5' end of the read, if so likely only different
                # due to sequencing adapters / TSO stil present at the end of the read!
                # Getting the edit distance from the first 5 base-pairs of the 3' end of the read
                three_prime_edit_dist = bio_edit_distance(edit_1_seq[::-1], edit_2_seq[::-1],
                                                          start_from_first_smallest_seq_aln=False,
                                                          alns_to_compare=10)

                # Only 1 mismatch at the 3' end! so let's consider them the same sequence, allocating shorter as sub...
                if three_prime_edit_dist <= 1:
                    dist_between_seqs = three_prime_edit_dist

            # Base-pair sequencing error correction
            # Some cases where the edits are VERY similar, off by 1bp, likely sequencing error!
            # Below is using the levenstein distance, so number of edits used to make one sequence into the other.
            # Taking the difference with the sequence lengths, to account for extra bp which would be added
            # to make the sequences the same if they are of different lengths. So would give 0 if 1 seq is a subset
            # of the other !
            if dist_between_seqs <= 2:  # edit distance cutoff from each other, accounts for short indel at start also

                # Recording the longer edit, assigning the shorter edit to the sub-edit list
                edit_pair = [edit_1, edit_2]
                edit_pair_lens = np.array([len(edit_1.alt_seq), len(edit_2.alt_seq)])
                if edit_pair_lens[0] != edit_pair_lens[1]: # different lengths, so clear which is the subedit
                    long_index = np.argmax([len(edit_1.alt_seq), len(edit_2.alt_seq)])
                else: # same length, choose one with more kmer matches to barcode
                    n_kmers = [len(edit_1.kmer_matches), len(edit_2.kmer_matches)]
                    if n_kmers[0] != n_kmers[1]: # Different kmer matches
                        long_index = np.argmax(n_kmers)
                    else: # Match from this criteria to, just choose the first then
                        long_index = 0

                long_edit = edit_pair.pop(long_index)
                long_edits = [long_edit] # So can process as a list below!

                subedit = edit_pair[0] # remaining edit is the subedit

                # Add the subedit!
                edits_to_subedits[long_edit] = edits_to_subedits[long_edit].union({subedit})

                while subedit in longest_edits: # Need to remove the subedit from list of long edits, found a longer match
                    longest_edits.pop( longest_edits.index(subedit) )

                if subedit not in sub_edits: # Need to add to subedit list
                    sub_edits.append( subedit )

            else: # They are different by all our criteria! including at 3' end. So let's consider them both long edits
                long_edits = [edit_1, edit_2]

            # Now adding to the longest_edits list if any long_edits found, if they were not already categorised as
            # subedits!!!!
            for long_edit in long_edits:
                # not previously found as a subedit, and found here is a longer edit, so add it to longest edits list!
                if long_edit not in longest_edits and long_edit not in sub_edits:
                    longest_edits.append( long_edit )

    # The surviving longest edits did not have a match with a longer t7 insert, so they are the set of unique edit sites
    longest_edits = list( set(longest_edits) )

    return longest_edits

def bed_file_flag_edits(bed_file, canonical_edit_sites, keep_sites, whitelist,
                          edit_dist, # Extra error around the edit site specification to allow for overlap with the bed regions
                          verbosity=1,
                        ):
    """ whitelist = True flag to flag any of the canonical edit sites for REMOVAL if intersect with the bed regions.
        whitelist = False to flag any canonical edit sites that intersect with the bed regions to KEEP.
    """
    if whitelist:
        filter_type = "whitelist"
        keep_indices = np.array( list(range(len(canonical_edit_sites))) ) # Since need to iterate through every site.
    else:
        filter_type = "blacklist"
        # Updating the keep_sites to remove those that intersect blacklisted regions...
        keep_indices = np.where(keep_sites)[0]

    blacklist_regions = pd.read_csv(bed_file, sep='\t', header=None)
    blacklist_chroms = blacklist_regions.values[:, 0].astype(str)
    blacklist_starts = blacklist_regions.values[:, 1].astype(int)
    blacklist_ends = blacklist_regions.values[:, 2].astype(int)

    black_chrom_set = set(list(blacklist_chroms))
    blacklist_chrom_to_indices = {chrom: np.where(blacklist_chroms == chrom)[0]
                                  for chrom in black_chrom_set}

    flagged_indices = []
    for index in keep_indices:
        edit_site = canonical_edit_sites[index]

        # Checking if intersects blacklist region.
        edit_chrom, edit_pos = edit_site.chrom, edit_site.ref_pos
        if edit_chrom in black_chrom_set:
            black_chrom_indices = blacklist_chrom_to_indices[edit_chrom]
            starts_, ends_ = blacklist_starts[black_chrom_indices], blacklist_ends[black_chrom_indices]
            # +/- edit_dist to the edit site, so can account for a level of error in edit site determination.
            overlap_bool = np.logical_and(starts_ <= (edit_pos + edit_dist), ends_ >= (edit_pos - edit_dist))

            if np.any(overlap_bool):
                keep_sites[index] = whitelist # Setting the keep_site to false, so now we filter this edit site !!!
                flagged_indices.append( index )
                print(f"Overlap found: {edit_site}: ", file=sys.stdout, flush=True) if verbosity >= 2 else None
                for overlap_index in np.where(overlap_bool)[0]:
                    print(f"{filter_type} region {edit_chrom}: {starts_[overlap_index]}-{ends_[overlap_index]}",
                                               file=sys.stdout, flush=True) if verbosity >= 2 else None

    return flagged_indices