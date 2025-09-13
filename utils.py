import pysam

class GenotypeData:
    def __init__(self, fasta_path, bam_path):
        """
        初始化 GenotypeData 对象，读取并解析 FASTA 和 BAM 文件
        :param fasta_path: 参考基因组 FASTA 文件路径
        :param bam_path: BAM 文件路径
        """
        self.fasta_path = fasta_path
        self.bam_path = bam_path
        self.reference = self._load_reference_genome(fasta_path)
        self.bam_data = self._load_bam_data(bam_path)

    def _load_reference_genome(self, fasta_path):
        """
        读取 FASTA 文件，加载参考基因组序列
        :param fasta_path: 参考基因组 FASTA 文件路径
        :return: 一个字典，key 为 contig 名称，value 为对应的基因组序列
        """
        reference = {}
        try:
            with pysam.FastaFile(fasta_path) as ref_file:
                for contig in ref_file.references:
                    reference[contig] = ref_file.fetch(contig)
            return reference
        except Exception as e:
            raise ValueError(f"无法加载 FASTA 文件 {fasta_path}: {str(e)}")

    def _load_bam_data(self, bam_path):
        """
        读取 BAM 文件，获取其内容
        :param bam_path: BAM 文件路径
        :return: pysam.AlignmentFile 对象，包含 BAM 文件内容
        """
        try:
            bam_file = pysam.AlignmentFile(bam_path, "rb")
            return bam_file
        except Exception as e:
            raise ValueError(f"无法加载 BAM 文件 {bam_path}: {str(e)}")

    def get_reference_at_position(self, contig, position):
        """
        获取参考基因组在某个位置的碱基
        :param contig: 染色体/基因组片段名（例如 'chr17'）
        :param position: 位置（1-based）
        :return: 参考基因组的碱基
        """
        try:
            return self.reference[contig][position - 1]  # 1-based -> 0-based
        except KeyError:
            raise ValueError(f"无法找到 contig {contig} 的参考序列")
        except IndexError:
            raise ValueError(f"位置 {position} 超出了 {contig} 的序列范围")

    def get_reads_at_position(self, contig, position):
        """
        获取在指定位置上所有比对到该位置的 reads
        :param contig: 染色体/基因组片段名（例如 'chr17'）
        :param position: 位置（1-based）
        :return: 所有在该位置上比对的 reads 列表
        """
        reads = []
        try:
            for pileupcolumn in self.bam_data.pileup(contig, position - 1, position, truncate=True):
                if pileupcolumn.reference_pos == position - 1:  # 0-based
                    for p in pileupcolumn.pileups:
                        if not p.is_del and not p.is_refskip:
                            # 使用 p.alignment.query_sequence 获取 reads 序列
                            reads.append(p.alignment.query_sequence[p.query_position])
            return reads
        except ValueError as e:
            raise ValueError(f"无法获取 {contig} 在位置 {position} 的 reads: {str(e)}")


# 示例用法
def main():
    fasta_path = '/root/DSAA6000-Assignment1/chr17-1.fa.gz'
    bam_path = '/root/DSAA6000-Assignment1/chr17_41000000_42000000-1-1.bam'

    genotype_data = GenotypeData(fasta_path, bam_path)

    # 获取参考基因组在某个位置的碱基
    ref_base = genotype_data.get_reference_at_position("17", 41000001)
    print(f"参考基因组在位置 41000001 的碱基是: {ref_base}")

    # 获取比对到某个位置的 reads
    reads = genotype_data.get_reads_at_position("17", 41000001)
    print(f"在位置 41000001 上比对到的 reads 数量: {len(reads)}")

if __name__ == "__main__":
    main()
