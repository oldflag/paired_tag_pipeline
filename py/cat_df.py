from argparse import ArgumentParser


parser = ArgumentParser()
parser.add_argument('infiles', help='input CSV files', nargs='+')
parser.add_argument('--out', help='the output file')

args = parser.parse_args()

out_hdl = open(args.out, 'wt')
infiles = [x for x in args.infiles if not x.startswith('input')]
for i, if_ in enumerate(infiles):
    with open(if_, 'rt') as inf:
        for j, line in enumerate(inf):
            if j == 0:
                if i == 0:
                    out_hdl.write(line)
            else:
                out_hdl.write(line)


out_hdl.close()
