import pysam

from utils import GenotypeData


def calculate_pr_d_given_g(reads, genotype, ref_base):
        """
        计算 Pr[D|G]：给定基因型的情况下观察到的reads的概率
        :param reads: 对应位置的reads
        :param genotype: 假设的基因型
        :param ref_base: 参考基因组碱基
        :return: 概率值
        """
        prob = 1.0
        if genotype == "AA" or genotype == "TT" or genotype == "GG" or genotype == "CC":
            allele = genotype[0]
            for read in reads:
                prob *= (1 - self.epsilon) if read == allele else self.epsilon
        else:
            # For heterozygous genotypes (e.g., AC, AG, etc.)
            for read in reads:
                if read == ref_base:
                    prob *= (1 - self.epsilon) / 2 + self.epsilon / 3
                else:
                    prob *= self.epsilon / 3 + (1 - self.epsilon) / 2
        return prob

def calculate_pr_g_given_r(ref_base):
    """计算 Pr[G|R]：给定参考碱基的基因型先验概率"""
    if ref_base == "A":
        return {"AA": 0.999, "AC": 0.0008, "AG": 0.0002}
    elif ref_base == "C":
        return {"CC": 0.999, "CG": 0.0008, "CT": 0.0002}
    elif ref_base == "G":
        return {"GG": 0.999, "GT": 0.0008, "GC": 0.0002}
    elif ref_base == "T":
        return {"TT": 0.999, "TG": 0.0008, "TC": 0.0002}

def call_genotype(contig, position):
    """调用基因型"""
    ref_base = self.get_reference_at_position(contig, position)
    reads = self.get_reads_at_position(contig, position)

    possible_genotypes = ["AA", "AC", "AG", "AT", "CC", "CG", "CT", "GG", "GT", "TT"]
    genotype_posteriors = {}

    for genotype in possible_genotypes:
        pr_d_given_g = self.calculate_pr_d_given_g(reads, genotype, ref_base)
        pr_g_given_r = self.calculate_pr_g_given_r(ref_base).get(genotype, 0)
        genotype_posteriors[genotype] = pr_d_given_g * pr_g_given_r

    most_likely_genotype = max(genotype_posteriors, key=genotype_posteriors.get)
    alt_base = most_likely_genotype[1]
    genotype_value = most_likely_genotype.count(alt_base)

    return contig, position, ref_base, alt_base, genotype_value


