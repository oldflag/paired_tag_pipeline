"""
Test for functions within `add_tags.py`
"""

import pysam
import add_tags as atg


def test_transform_read():
    """
    Tests for read transformation
    """
    READ_NAMES = [
      'read1',
      'FQZZDDAAXX whatever',
      'DFAR!#@#$TA',
      'this:one:has:colons and other stuff',
      'HNGAFED#@!::',
      'agg83128y4fDFEA agh'
    ]
    UMI = ['AACCGGTT', 'TTAACCGGAAAA', 'GGGGCCCCGG',
           'NNNNNNNNNNNNNNNN']
    BC = [('GAGAGAG', 'BC4'), ('GGGGGGGG', 'BC6'),
          ('AAAAAAAAA', 'BC2'), ('ATGACGTAC', 'BC92'),
          ('TACCGATAC', 'BC32')]
    SM = [('AAAA', 'A0'), ('TTTT', 'C4'), ('GGGG', 'my_name_is'),
          ('TATA', 'some sample'), ('TCGC', 'another name')]
    seq = 'A' * 50
    
    for rname in READ_NAMES:
        for umi in UMI:
            for bc1 in BC:
                for bc2 in BC:
                    for sam in SM:
                        rg = f'{sam[1]}_RG'
                        rgm = {sam[1]: rg}
                        rec = pysam.AlignedSegment()
                        rec.query_sequence = seq
                        fullname = rname + '|' + umi + ':' + bc1[1] + ':' + bc2[1] + ':' + sam[1] + ':rna'
                        fullname += '|' + umi + ':' + bc1[0] + ':' + bc2[0] + ':' + sam[0] + ':GCCGCCG'
                        rec.query_name = fullname
                        newrec = atg.transform_read(rec, 1, rgm)
                        assert newrec.get_tag('MI') == umi
                        assert newrec.get_tag('BC') == sam[0]
                        assert newrec.get_tag('CR') == umi + ':' + bc1[0] + ':' + bc2[0] + ':' + sam[0] + ':GCCGCCG'
                        assert newrec.get_tag('CB') == bc1[1] + '.' + bc2[1] + '.' + sam[1]
                        assert newrec.get_tag('XX') == '000000001'
                        
    
