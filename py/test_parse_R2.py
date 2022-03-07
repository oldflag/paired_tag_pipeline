"""
Test the read parser by building up R2 reads and putting them through the parsers
"""
import os
from skbio.alignment import StripedSmithWaterman
import parse_R2 as pr2

def test_r2_parsing():
    LINKERS_R1 = ['ACTGACTAGATCAAACTATC', 'GGATTACTATCGCGACAT',
                  'CCATACGATATTACGATACCGAT', 'GGGATACTTCACGTATTACGAGGC']
    LINKERS_R2 = ['CTAGTGACGAACTACCCTGA', 'CCATTTGGTACCTTATAC',
                  'ACTGAGCCTAAACGTACCGGCTG', 'TACGGACATGAGTTCAGCGGCATC']
    
    BARCODES = ['AACGATAC', 'CCGACAGT', 'GGATACAG', 'TACATTAC', 'GTATATTG',
                'CTGATCTA', 'ATTCAGTA', 'AGAGTAGA', 'CATGTTAT', 'TACCTCAC',
                'CCAGATC', 'TAGGAT']
    
    SAM_BARCODES = ['AAAA', 'TTTT', 'CAGC', 'GATT', 'TTAC', 'GACA', 'CAG', 'GCG']
    
    UMI = ['ACTGACTAGA', 'GGATATGCAT', 'ATACGGATAC', 'TTACGGGATA',
           'AAAAGGATAC', 'TTACGATATA', 'GGGGATATGA', 'AAGAGTCCAC',
           'GGATCATC', 'GGCATTAT', 'TTATTC']
   
    GAPS = ['', 'C', 'T', 'AC', 'AG', 'TCA', 'ACT', 'TTC', 'GGG']
    PRIMERS = ['AAC', 'AA', 'A', '']
    

    for l1 in LINKERS_R1:
      l1_sw = StripedSmithWaterman(l1)
      for l2 in LINKERS_R2:
        if len(l2) != len(l1):
          continue
        l2_sw = StripedSmithWaterman(l2)
        for bc1 in BARCODES:
          for bc2 in BARCODES:
            if len(bc1) != len(bc2):
              continue
            for sbc in SAM_BARCODES:
              for umi in UMI:
                for gp in GAPS:
                  for pm in PRIMERS:
                    r2_seq = f'{umi}{bc1}{gp}{l1}{bc2}{l2}T{sbc}GCGGCCGC{pm}'
                    parsed = pr2.extract_barcodes(r2_seq, l1_sw, l2_sw, len(umi),
                                                  len(bc1), len(sbc), len(l1))
                    print(r2_seq)
                    assert umi == parsed[0], (umi, parsed[0])
                    assert bc1 == parsed[1], (bc1, parsed[1])
                    assert bc2 == parsed[2], (bc2, parsed[2])
                    assert sbc == parsed[3], (sbc, parsed[3])
                    assert 'TGCGGCCGC' == parsed[4], ('TGCGGCCGC', parsed[4])
                   
 
