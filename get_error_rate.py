import sys
import re

def main(sam_file):
    ncorrect_reads = 0
    ncorrect_mm = 0
    nincorrect_reads = 0
    nincorrect_mm = 0
    with open(sam_file) as infile:
        for row in infile:
            if row.startswith('READ'):
                cols = row.split('\t')
                readid, chromo, poso = cols[0].split('_')
                chromo = chromo.strip().upper()
                poso = int(poso)
                chrom = cols[2].strip().upper()
                pos = int(cols[3])-1
                
                if pos != -1:
                    if chromo == chrom or poso == pos:
                        ncorrect_reads += 1
                        m = re.match('NM:i:(\d+)', cols[12])
                        mm = int(m.groups()[0]) if m else 0
                        ncorrect_mm += mm
                    else:
                        m = re.match('NM:i:(\d+)', cols[12])
                        mm = int(m.groups()[0]) if m else 0
                        nincorrect_mm += mm
                        nincorrect_reads += 1

    print('Error rate, correctly mapped reads (%): ',
          round(ncorrect_mm*100/(ncorrect_reads*50), 2))
    print('Error rate, incorrectly mapped reads (%): ',
          round(nincorrect_mm*100/(nincorrect_reads*50), 2))
                

if __name__ == '__main__':
    main(sys.argv[1])