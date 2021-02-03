# cython: language_level=3
# Author: Avraham Kahan (v3 author) and Ruslana Frazer (v1/v2 author)
import pysam, csv, random
from pysam.libcalignedsegment import AlignedSegment  # for type delcarations
from pysam.libcalignmentfile import IteratorRowRegion  # for type delcarations
from collections import namedtuple, defaultdict
from typing import List, Tuple

# depending on version of python being used may be better to use normal classes
CachedEndRead = namedtuple("CachedEndRead", "read real_end")
Locus = namedtuple("Locus", "chrom start end pattern num_repeats mapped_reads")  # type declaration Locus
Result = namedtuple("Result", "locus MS_occurence_count") # MS_occurence_count = (length between locus start and end)/(length of MS)


class LatestLocus:
    # class to hold position and mapped reads of current latest locus
    def __init__(self):
        self.mapped_reads = []
        self.start = -1


class CIGAR_OPTIONS:
    ALG_MATCH = 0
    INSERTION = 1
    DELETION = 2
    REF_SKIP = 3
    SOFT_CLIP = 4
    HARD_CLIP = 5
    PADDING = 6
    SEQ_MATCH = 7
    SEQ_MISMATCH = 8


class FLAG_OPTIONS:
    MUTLIPLE_SEGMENTS_EXIST = 0x1
    ALL_PROPERLY_ALIGNED = 0x2
    SEG_UNMAPPED = 0x4
    NEXT_SEG_UNMAPPED = 0x8
    SEQ_REV_COMP = 0x10
    NEXT_SEQ_REV_COMP = 0x20
    FIRST_SEQ_IN_TEMP = 0x40
    LAST_SEQ_IN_TEMP = 0x80
    SECONDARY_ALG = 0x100
    POOR_QUALITY = 0x200
    DUPLICATE_READ = 0x400
    SUPPLEMENTARY_ALG = 0x800


def get_read_repeat_length(read: AlignedSegment, locus: Locus) -> float:
    #print("*************\n start: {0}\n end:{1} \n Cigar {2}".format(read.reference_start, read.reference_end, read.cigartuples))
    read_position = read.reference_start+1
    indel_bases = 0 # number of added/deleted bases in MS locus
    for cigar_op in read.cigartuples:
        if cigar_op[0] in [CIGAR_OPTIONS.ALG_MATCH, CIGAR_OPTIONS.SEQ_MATCH, CIGAR_OPTIONS.SEQ_MISMATCH]:
            read_position+=cigar_op[1]
        elif cigar_op[0] == CIGAR_OPTIONS.INSERTION:
            if locus.start <= read_position <= locus.end:
                indel_bases+=cigar_op[1]
        elif cigar_op[0] == CIGAR_OPTIONS.DELETION:
            if locus.start <= read_position <= locus.end:
                indel_bases-=min(cigar_op[1], locus.end-read_position) # in case deletion goes past end of MS locus
            elif read_position < locus.start:
                indel_bases -= max(0, read_position + cigar_op[1] - locus.start) # case where deletion extends into locus from flanking
            read_position+=cigar_op[1]
    return max(locus.num_repeats + indel_bases/len(locus.pattern), 0) # so is never negative


def strip_chromosome(chromosome: str) -> str:
    """ return chromosome as only its integer (ex. chr16 -> 16) """
    if chromosome.upper().find("X") != -1:
        return "X"
    elif chromosome.upper().find("Y") != -1:
        return "Y"
    else:
        numeric_filter = filter(str.isdigit, chromosome)
        return "".join(numeric_filter)


def format_result(result: Result) -> str:
    """
    :return: formatted string result for given histogram
    """
    results_string = ", ".join(
        ["{0}_{1}".format(MS_repeats, read_count) for MS_repeats, read_count  in result.MS_occurence_count.items()])
    return "{0}:{1}:{2}:{3}:{4}, {5}".format(strip_chromosome(result.locus.chrom), result.locus.start, result.locus.end,
                                                            result.locus.pattern, result.locus.num_repeats,
                                                            results_string)


def probability_filter(p_exclude) -> bool:
    return p_exclude < random.random()


def process_MS_locus(cur_locus: Locus, p_exclude: float) -> str:
    MS_histogram = defaultdict(lambda: 0)  # default value for unadded key is 0; ex. res_dict[<unadded_key>]+=1 will create key value pair of 1
    for cached_end_read in cur_locus.mapped_reads:
        read = cached_end_read.read  # initial read, without end data
        if probability_filter(p_exclude):
            # the way of calculating repeat lenghts is one that is up for discussion
            MS_occurences = round(get_read_repeat_length(read, cur_locus) + .00001)  # so x.5 rounds to x+1
            MS_histogram[MS_occurences] += 1
    return format_result(Result(locus=cur_locus, MS_occurence_count=MS_histogram))


def simple_filter(read: AlignedSegment) -> bool:
    return not (read.flag & FLAG_OPTIONS.SECONDARY_ALG or read.flag & FLAG_OPTIONS.POOR_QUALITY
                or read.flag & FLAG_OPTIONS.DUPLICATE_READ or read.flag & FLAG_OPTIONS.SUPPLEMENTARY_ALG or not read.cigartuples)


def get_next_mapped_read(reads_iterator: IteratorRowRegion) -> AlignedSegment:
    cur_read = next(reads_iterator, None)
    while cur_read is not None:
        if cur_read.cigartuples is None:  # unaligned
            cur_read = next(reads_iterator, None)
        else:
            return cur_read
    return None


def backtrack_reads(latest_locus_reads: List[CachedEndRead], locus_start, locus_end: int, flanking) -> List[
    CachedEndRead]:
    # get reads from last locus' reads
    return [adjusted_read for adjusted_read in latest_locus_reads if
            adjusted_read.read.reference_start + 1 <= locus_start - flanking and locus_end + flanking <= adjusted_read.real_end]


def find_mapped_reads(first_read: AlignedSegment, reads_iterator: IteratorRowRegion, latest_locus_reads:
List[CachedEndRead], locus_start: int, locus_end: int, flanking: int, removed_bases: int) \
        -> Tuple[List[CachedEndRead], AlignedSegment]:

    loci_mapped_reads = backtrack_reads(latest_locus_reads, locus_start, locus_end,
                                        flanking)  # maps read from previous latest start
    cur_read = first_read  # first read is used so that the first read that didn't map to the previous locus can still be used to map to this one; elsewise, we would have lost it in the iterator
    while cur_read is not None:
        if locus_start - flanking < cur_read.reference_start + 1:
            return loci_mapped_reads, cur_read
        read_end = cur_read.reference_end+1
        if cur_read.reference_start + 1 <= locus_start - flanking and locus_end + flanking <= read_end:
            while cur_read.reference_start + 1 <= locus_start - flanking:
                # not part of while loop since this way if a read with the same start location but shorter length appears the while loop does not break
                if locus_end + flanking <= read_end and simple_filter(cur_read):
                    loci_mapped_reads.append(CachedEndRead(cur_read, read_end))
                cur_read = get_next_mapped_read(reads_iterator)
                if cur_read is None:  # reads iterator is exhausted
                    break
                read_end = cur_read.reference_end+1 - removed_bases
            return loci_mapped_reads, cur_read

        else:
            cur_read = get_next_mapped_read(reads_iterator)

    return loci_mapped_reads, None


def find_prefix(bam_file: pysam.AlignmentFile) -> str:
    """ returns the prefix of contigs in the BAM (ex. Chr, nothing [''], or chr)"""
    prefixes = ['', 'chr', 'Chr']
    for prefix in prefixes:
        try:
            test_fetch = bam_file.fetch("{0}1".format(prefix), 10_000)
            print("PREFIX FOR CHROMOSOMES IS: {0}".format(prefix))
            return prefix
        except ValueError: # different prefix
            continue


def main(input_bam: str, output_file: str, csv_loci: str, flanking: int, removed_bases: int, p_exclude: float):
    loci_iterator = csv.reader(open(csv_loci), dialect="excel-tab")
    bam_file = pysam.AlignmentFile(input_bam, "rb")  # load BAM index file
    prefix = find_prefix(bam_file)
    results = []
    current_chrom = '0'
    prev_start = 0
    for locus in loci_iterator:
        chrom = locus[0]
        locus_start = int(locus[3])
        locus_end = int(locus[4])
        # get new iterator for new chromosome or if locus is 8,000 away, since a fetch will go through 8,000 irrelevant alignments on average
        if chrom != current_chrom or abs(prev_start - locus_start) > 8000:
            try:
                reads_iterator = bam_file.fetch("{0}{1}".format(prefix, chrom), start=locus_start)  # hop to part of chromosome with relevant reads
                first_read = get_next_mapped_read(reads_iterator)
            except ValueError:
                print("NO READS FOR CHROMOSOME \'{0}\'  IN BAM".format(chrom))
                first_read = None
                reads_iterator = None

            current_chrom = chrom
            latest_locus = LatestLocus()



        # first read is used so that the first read that didn't map to the previous locus can still be used to map to this one; elsewise, we would have lost it in the iterator
        current_reads, first_read = find_mapped_reads(first_read, reads_iterator, latest_locus.mapped_reads,
                                                      locus_start, locus_end,
                                                      flanking, removed_bases)
        #print("********************************")
        #for read in current_reads:
            #print(read.read.reference_start, read.real_end)
        results.append(process_MS_locus(Locus(chrom=chrom, start=locus_start, end=locus_end, pattern=locus[12],
                                              num_repeats=float(locus[6]), mapped_reads=current_reads), p_exclude))
        prev_start = locus_start
        if latest_locus.start < locus_start:
            latest_locus.start = locus_start
            latest_locus.mapped_reads = current_reads

    with open(output_file, "w") as results_file:
        for formatted_result in results:
            results_file.write(formatted_result + "\n")

