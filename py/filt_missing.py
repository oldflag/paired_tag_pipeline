import gzip
from argparse import ArgumentParser

def get_args():
    parser = ArgumentParser()
    parser.add_argument('r1_in')
    parser.add_argument('r2_in')
    parser.add_argument('r1_out')
    parser.add_argument('r2_out')

    return parser.parse_args()


def main(r1_in, r2_in, r1_out, r2_out):
    r1_hdl = gzip.open(r1_in, 'rt') if r1_in[-2:] == 'gz' else open(r1_in, 'rt')
    r2_hdl = gzip.open(r2_in, 'rt') if r2_in[-2:] == 'gz' else open(r2_in, 'rt')

    r1_o = open(r1_out, 'wt')
    r2_o = open(r2_out, 'wt')

    r1_lines, r2_lines = list(), list()
    for i, (r1l, r2l) in enumerate(zip(r1_hdl, r2_hdl)):
        r1_lines.append(r1l)
        r2_lines.append(r2l)
        if i % 4 == 3:
            # check the sequence length
            if len(r1_lines[1]) <= 1:
                r1_lines[1] = 'N\n'
                r1_lines[3] = '!\n'
            elif len(r1_lines[1]) < len(r1_lines[3]):
                r1_lines[1] = r1_lines[1].strip() + 'N' * (len(r1_lines[3]) - len(r1_lines[1])) + '\n'
            elif len(r1_lines[1]) > len(r1_lines[3]):
                r1_lines[3] = r1_lines[3].strip() + '!' * (len(r1_lines[1]) - len(r1_lines[3])) + '\n'
            # write out
            r1_o.write(''.join(r1_lines))
            r2_o.write(''.join(r2_lines))
            r1_lines, r2_lines = list(), list()
    for g in (r1_hdl, r2_hdl, r1_o, r2_o):
        g.close()


if __name__ == '__main__':
    args = get_args()
    main(args.r1_in, args.r2_in, args.r1_out, args.r2_out)
