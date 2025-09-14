import math
import pysam
import argparse

from tqdm import tqdm
from collections import Counter
from utils import GenotypeData


def calculate_pr_g_given_r(ref_base):
    all_genos = ["AA","AC","AG","AT","CC","CG","CT","GG","GT","TT"]
    pri = {g: 0.0 for g in all_genos}

    if ref_base not in ("A","C","G","T"):
        return pri  # N 等非标准碱基直接置 0

    # 同型参考：xx（x==R）
    pri[ref_base * 2] = 0.999

    # 杂合参考：x==R, y!=R（3种，平分0.0008）
    others = [b for b in ("A","C","G","T") if b != ref_base]
    for y in others:
        g = "".join(sorted(ref_base + y))  # 统一成字典序，如 AC 而非 CA
        pri[g] = 0.0008 / 3.0

    # 同型非参考：yy（y!=R，3种，平分0.0002）
    for y in others:
        pri[y * 2] = 0.0002 / 3.0

    # 双非参考杂合不赋值（作业设定为 0）
    return pri


# =========================
# 固定 ε 似然（快）
# =========================
def calculate_pr_d_given_g(reads, genotype, epsilon=0.01):
    """Pr[D | G]，对 reads 做计数乘积（快；区域不深时足够稳定）"""
    reads = [b for b in reads if b in ("A","C","G","T")]
    if not reads:
        return 0.0

    cnt = Counter(reads)
    bases = ("A","C","G","T")

    if genotype[0] == genotype[1]:  # homozygous xx
        x = genotype[0]
        probs = {b: (1 - epsilon) if b == x else (epsilon / 3.0) for b in bases}
    else:  # heterozygous xy
        x, y = genotype[0], genotype[1]
        probs = {b: (0.5 - epsilon/3.0) if b in (x, y) else (epsilon / 3.0) for b in bases}

    p = 1.0
    for b, c in cnt.items():
        p *= probs[b] ** c
    return p


# =========================
# 逐碱基质量版对数似然（更准，稍慢）
# =========================
def loglik_reads_given_g_with_quals(reads, quals, genotype):
    """sum log Pr[D_i | G]，用 Phred 质量 Q_i 转 ε_i=10^{-Q/10}"""
    logp = 0.0
    for b, q in zip(reads, quals):
        if b not in ("A","C","G","T") or q is None:
            continue
        eps = max(1e-6, min(0.25, 10 ** (-q / 10)))  # clamp 防极端
        if genotype[0] == genotype[1]:
            x = genotype[0]
            p = (1 - eps) if b == x else (eps / 3.0)
        else:
            x, y = genotype[0], genotype[1]
            p = (0.5 - eps / 3.0) if b in (x, y) else (eps / 3.0)
        logp += math.log(p)
    return logp


# =========================
# 用“已有 reads + ref”直接判型（不再 pileup）
# =========================
def call_genotype_from_reads(reads, ref_base, epsilon=0.01, use_quals=False, quals=None):
    genos = ["AA","AC","AG","AT","CC","CG","CT","GG","GT","TT"]
    pri_map = calculate_pr_g_given_r(ref_base)

    best_g, best_score = None, -float("inf")

    if use_quals and quals is not None:
        # 对数后验：loglike + logprior
        for g in genos:
            prior = pri_map.get(g, 0.0)
            if prior <= 0.0:
                continue
            ll = loglik_reads_given_g_with_quals(reads, quals, g)
            score = ll + math.log(prior)
            if score > best_score:
                best_score = score
                best_g = g
    else:
        # 直接乘积后验：like * prior（用 log 比较更稳，但为速度保留乘积）
        for g in genos:
            prior = pri_map.get(g, 0.0)
            if prior <= 0.0:
                continue
            like = calculate_pr_d_given_g(reads, g, epsilon)
            score = like * prior
            if score > best_score:
                best_score = score
                best_g = g

    if not best_g:
        return None

    # 映射到 Alt 与 0/1/2
    if best_g == ref_base * 2:
        alt, geno_count = ref_base, 0
    elif ref_base in best_g and best_g[0] != best_g[1]:
        alt = best_g.replace(ref_base, "")
        geno_count = 1
    elif best_g[0] == best_g[1] and best_g[0] != ref_base:
        alt, geno_count = best_g[0], 2
    else:
        # 双非参考杂合（理论上先验为 0）
        return None

    return alt, geno_count, best_g, best_score


# =========================
# 主流程：单次 pileup + 进度条 + 多线程
# =========================
def run_region(fasta_path, bam_path,
               region_start=41000000, region_end=42000000,
               epsilon=0.01, min_depth=3, out_path=None,
               show_progress=True, threads=4, min_baseq=13, min_mapq=0,
               max_depth=10000, use_quals=False):
    gd = GenotypeData(fasta_path, bam_path)

    # BAM 解压多线程（若 HTSlib 支持）
    try:
        gd.bam_data.set_threads(threads)
    except Exception:
        pass

    # 选择 contig 并规范输出名
    bam_refs = set(gd.bam_data.references)
    fasta_refs = set(gd.reference.keys())

    contig, chrom_out = None, None
    for cand in ("chr17", "17"):
        if cand in bam_refs and cand in fasta_refs:
            contig = cand
            break
    if contig is None:
        # 回退：任取交集里一个
        common = list(bam_refs & fasta_refs)
        if not common:
            raise RuntimeError("BAM 与 FASTA 的 contig 不匹配！")
        contig = common[0]

    chrom_out = "chr17" if "17" in contig else contig

    results = []
    total_bp = region_end - region_start + 1
    last_pos_to_update = region_start
    pbar = tqdm(total=total_bp, disable=not show_progress,
                desc=f"{contig}:{region_start}-{region_end}",
                unit="bp", dynamic_ncols=True, miniters=1)

    # 只做一次 pileup，Python 层不再重复
    for col in gd.bam_data.pileup(
        contig, region_start, region_end,
        truncate=True,
        stepper="samtools",
        min_base_quality=min_baseq,
        ignore_overlaps=True,
        max_depth=max_depth
    ):
        pos1 = col.reference_pos + 1
        if pos1 < region_start or pos1 > region_end:
            continue

        # 进度推进（包括未覆盖区间的一次性跳过）
        if pos1 >= last_pos_to_update:
            pbar.update(pos1 - last_pos_to_update)
            last_pos_to_update = pos1 + 1

        # 收集该位点碱基 / 质量（过滤 indel/跳跃；可加映射质量过滤）
        reads_here, quals_here = [], []
        for p in col.pileups:
            if p.is_del or p.is_refskip:
                continue
            if min_mapq and p.alignment.mapping_quality < min_mapq:
                continue
            base = p.alignment.query_sequence[p.query_position]
            if base in ("A","C","G","T"):
                reads_here.append(base)
                if use_quals:
                    q = p.alignment.query_qualities[p.query_position]
                    quals_here.append(q)

        if len(reads_here) < min_depth:
            continue

        ref_base = gd.get_reference_at_position(contig, pos1).upper()
        res = (call_genotype_from_reads(reads_here, ref_base,
                                        epsilon=epsilon,
                                        use_quals=use_quals,
                                        quals=quals_here if use_quals else None))
        if not res:
            continue
        alt, geno_count, best_g, _ = res

        if geno_count >= 1:
            results.append((chrom_out, pos1, ref_base, alt, geno_count))
            pbar.set_postfix(variants=len(results))

    if last_pos_to_update <= region_end:
        pbar.update(region_end - last_pos_to_update + 1)
    pbar.close()

    # 输出
    if out_path:
        with open(out_path, "w") as f:
            for r in results:
                f.write(f"{r[0]} {r[1]} {r[2]} {r[3]} {r[4]}\n")
    else:
        for r in results:
            print(f"{r[0]} {r[1]} {r[2]} {r[3]} {r[4]}")

    return results


# =========================
# CLI
# =========================
def parse_args():
    ap = argparse.ArgumentParser(
        description="Bayesian genotyper (fast, single-pileup, tqdm, HTSlib threads).")
    ap.add_argument("--fasta", default='./chr17-1.fa.gz', help="Reference FASTA (indexed)")
    ap.add_argument("--bam", default='./chr17_41000000_42000000-1-1.bam', help="BAM (indexed)")
    ap.add_argument("--start", type=int, default=41000000, help="Region start (1-based)")
    ap.add_argument("--end", type=int, default=42000000, help="Region end (1-based, inclusive)")
    ap.add_argument("--epsilon", type=float, default=0.01, help="Base error rate (fixed ε, if not using quals)")
    ap.add_argument("--min-depth", type=int, default=3, help="Min effective depth")
    ap.add_argument("--min-baseq", type=int, default=13, help="Min base quality (pileup filter)")
    ap.add_argument("--min-mapq", type=int, default=0, help="Min mapping quality (read-level filter)")
    ap.add_argument("--threads", type=int, default=4, help="HTSlib threads for BAM decompression")
    ap.add_argument("--max-depth", type=int, default=10000, help="Pileup max depth cap")
    ap.add_argument("--use-quals", action="store_true",
                    help="Use per-base Phred qualities (log-likelihood) for ε_i")
    ap.add_argument("--out", default='./output.txt', help="Output path (default: print to stdout)")
    ap.add_argument("--no-progress", action="store_true", help="Disable tqdm progress bar")
    return ap.parse_args()


def main():
    args = parse_args()
    run_region(
        fasta_path=args.fasta,
        bam_path=args.bam,
        region_start=args.start,
        region_end=args.end,
        epsilon=args.epsilon,
        min_depth=args.min_depth,
        out_path=args.out,
        show_progress=not args.no_progress,
        threads=args.threads,
        min_baseq=args.min_baseq,
        min_mapq=args.min_mapq,
        max_depth=args.max_depth,
        use_quals=args.use_quals
    )


if __name__ == "__main__":
    main()