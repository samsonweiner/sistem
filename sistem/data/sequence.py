import subprocess

from sistem.genome import SNVChromosome
from sistem.utilities.utilities import get_reg_id

def build_cell_ref(chromosome, ref, num_regions, region_len, out_path):
    nucs = ['A', 'C', 'G', 'T']
    with open(out_path, 'w+') as f:
        f.write(f'>{chromosome.name}\n')
        for r in chromosome.seq:
            q = get_reg_id(r)
            if q == num_regions - 1:
                ref_seq = list(ref[chromosome.name][int(q*region_len):].seq.upper())
            else:
                ref_seq = list(ref[chromosome.name][int(q*region_len):int((q+1)*region_len)].seq.upper())
            if isinstance(chromosome, SNVChromosome):
                if r in chromosome.SNVs:
                    positions = chromosome.SNVs[r]
                    for p,idx in positions.items():
                        if p < len(ref_seq) - 1:
                            if ref_seq[p] != 'N':
                                alt_nucs = nucs[:]
                                alt_nucs.remove(ref_seq[p])
                                alt = alt_nucs[idx]
                                ref_seq[p] = alt
            f.write("".join(ref_seq))
    call = subprocess.run(['samtools', 'faidx', out_path])
