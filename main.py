import HistogramGenerator 
import argparse, time


def create_parser() -> argparse.ArgumentParser:
    """
    :return: creates parser with all command line arguments arguments
    """
    parser = argparse.ArgumentParser(description='create INDEL information from a sample on a specific set of loci.')
    parser.add_argument("-I", "--input_BAM", help="Input BAM file.")
    parser.add_argument("-O", "--output_file", help="Output file.", default="output.txt")
    parser.add_argument("-l", "--loci_list", help="List of loci to be processed and included in the output.")
    parser.add_argument("-f", "--flanking", help="Length of flanking on both sides of an accepted read.", default=10)
    parser.add_argument("-r", "--removed_bases", help="Number of bases to be removed at the end of the read.", default=0)
    parser.add_argument("-e", "--exclude", help="The probability that a read will be randomly removed.", default=0)
    return parser


if __name__=='__main__':
    print(time.localtime())
    parser = create_parser()
    args = parser.parse_args()
    HistogramGenerator.main(args.input_BAM, args.output_file, args.loci_list, int(args.flanking), int(args.removed_bases), float(args.exclude))
    print(time.localtime())
