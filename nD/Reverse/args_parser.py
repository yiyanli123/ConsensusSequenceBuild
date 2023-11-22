import argparse

def set_parser():
    parser = argparse.ArgumentParser()

    parser.add_argument( '--fq_file', type=str, help='.fq file path.')
    parser.add_argument( '--res_dir', type=str, help='output directory.')
    parser.add_argument( '--hairpin_seq', type=str, help='hairpin seq.')
    parser.add_argument( '--err_threshold', type=float, default=0.3, help='err threshold.')
    parser.add_argument( '--step_len', type=int, default=1, help='step len')
    parser.add_argument( '--adjusted_len', type=int, default=200, help='adjust len.')
    parser.add_argument( '--subread_len_filter_thre', type=float, default=0.8, help='subread_len_filter_thre.')
    return  parser.parse_known_args()[0]