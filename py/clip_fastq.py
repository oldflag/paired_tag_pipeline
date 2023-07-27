from argparse import ArgumentParser
import os 

def get_args():
    parser = ArgumentParser()
    parser.add_argument('in_fastq')
    parser.add_argument('out_fastq')
    parser.add_argument('--bp', type=int, default=25)

    return parser.parse_args()


def clip(fastq, out, nbp):
    out_hdl = open(out, 'wt')
    in_hdl = open(fastq, 'rt')

    for line in in_hdl:
        if line[0] in {'@', '+'}:
            out_hdl.write(line)
        else:
            out_hdl.write(line[:nbp].strip() + '\n') 

    in_hdl.close()
    out_hdl.close()


if __name__ == '__main__':
    args = get_args()
    clip(args.in_fastq, args.out_fastq, args.bp)
