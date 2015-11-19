import numpy as np
import cyvcf2
print cyvcf2
from cyvcf2 import VCF

def is_phased(gt_array, ploidy):
    return gt_array[1::ploidy]&1

def gt_allele(gt_array, ploidy):
    v = (gt_array>>1) - 1
    v[v < 0] = -1
    return v.reshape((len(v)/ploidy), ploidy)


def is_missing(gt_array, ploidy):
    vals = []
    for i in range(ploidy):
        vals.append(((gt_array[i::ploidy]>>1) == 0) | (gt_array[i::ploidy] < 0))
    return np.array(vals).T



for v in VCF("test-hemi.vcf"):
    print v.ploidy, v.gt_array
    alleles = (gt_allele(v.gt_array, v.ploidy))
    print alleles
    print "phase:", is_phased(v.gt_array, v.ploidy)
    print is_missing(v.gt_array, v.ploidy)
