#!/bin/env python3
import sys
from collections import defaultdict
import gzip
import pysam
import re
import itertools
import misc
import time
import pandas as pd
import math
from scipy.stats import mannwhitneyu
from cliffs_delta import cliffs_delta
from statistics import mean, stdev
from statistics import mean, variance
import scipy.special as sp
import os
import copy

logger = misc.create_logger("bridge")
  
class PhasedBlock(object):
    def __init__(self, chrom, ps):
        self.chrom = chrom
        self.ps = ps
        self.variants = []
    def add(self, pos):
        self.variants.append(pos)
    def start(self):
        return self.variants[0]
    def end(self):
        return self.variants[-1]

class PhasedBlockSet(object):

    @staticmethod

    def from_vcf(fname, overlapping=False):
        blk_dict = {}
        CHROM, FORMAT, SAMPLE = 0, -2, -1
        logger.info("PhasedBlockSet.form_vcf")
        for rec in pysam.VariantFile(fname):
             vid = (rec.chrom, rec.pos)
             assert(len(rec.samples) == 1)
             if "PS" in rec.format and type(rec.samples[0]["PS"]) == int:
                 bid = (rec.chrom, rec.samples[0]["PS"])
                 if bid not in blk_dict:
                     blk_dict[bid] = PhasedBlock(bid[0], bid[1])
                 blk_dict[bid].add(rec.pos)
        logger.info("Load %d phased blocks from %s" % (len(blk_dict), fname))
    
        if not overlapping:
            prev = None
            trivial = set()
            for k, v in sorted(blk_dict.items()):
                if prev != None and prev[0][0] == k[0] and prev[1].end() > k[1]:
                    trivial.add(k)
                    logger.info("detected overlapping blocks: %s:%d" % k)
                else:
                    prev = (k, v)
            [blk_dict.pop(k) for k in trivial]
        blk_list = sorted(blk_dict.values(), key=lambda x: (x.chrom, x.ps))
        
        return blk_list 


class SNP(object):
    # SNP对象：保存某个位点在各read上的碱基支持记录
    def __init__(self, chrom, pos, ref, alt):
        self.chrom = chrom
        self.pos = pos
        self.ref = ref.upper()
        self.alt = alt.upper()
        self.reads = []

    def add(self, qname, is_forward, read_base, support_type, ps, hp):
        read_base = read_base.upper()
        mutation = f"{self.ref}>{read_base}"
        self.reads.append((qname, is_forward, read_base, support_type, mutation, ps, hp))

    def get_phased_distribution_beta(self, pseudocount=0.5):
        # 基于已定相read构建每个(ps,hp)下ALT概率的Beta posterior参数
        # 编码规则：ALT=1, REF=0（OTHER不参与）
        by_hap = defaultdict(list)
        for qname, is_fwd, read_base, support_type, mutation, ps, hp in self.reads:
            if hp not in (1, 2) or ps == 0:
                continue
            if support_type == "ALT":
                by_hap[(ps, hp)].append(1.0)
            elif support_type == "REF":
                by_hap[(ps, hp)].append(0.0)

        result = []
        for (ps, hp), values in by_hap.items():
            if len(values) >= 2:
                alt_count = sum(values)
                ref_count = len(values) - alt_count
                alpha = alt_count + pseudocount
                beta = ref_count + pseudocount
                result.append((ps, hp, alpha, beta))
        return result

    def get_mean_per_hap(self, pseudocount=0.5):
        """Return per-haplotype beta means as {hp: (ps, mean)}."""
        means = {}
        for ps, hp, alpha, beta in self.get_phased_distribution_beta(pseudocount=pseudocount):
            means[hp] = (ps, alpha / (alpha + beta))
        return means

    def valid_phased(self):
        means = self.get_mean_per_hap()
        return 1 in means and 2 in means

    def test_phased(self, sval, pseudocount=0.5):
        """Assign REF/ALT-coded value (0/1) to haplotype using beta means."""
        means = self.get_mean_per_hap(pseudocount=pseudocount)
        if 1 not in means or 2 not in means:
            raise ValueError("SNP site lacks phased support for both haplotypes")
        ps1, mean1 = means[1]
        ps2, mean2 = means[2]
        ps = ps1 if ps1 != 0 else ps2
        midpoint = (mean1 + mean2) / 2.0
        hp = 1 if abs(sval - mean1) <= abs(sval - mean2) else 2
        return (ps, hp, abs(sval - midpoint))
    
    def get_consensus(self):
        counts = {1: [0, 0], 2: [0, 0]}
        for tup in self.reads:
            if len(tup) == 7:  # 兼容旧格式
                _, _, _, st, _, ps, hp = tup
            else:
                _, _, _, st, ps, hp = tup
            if hp not in (1, 2) or ps == 0: continue
            idx = 1 if st == "ALT" else 0
            counts[hp][idx] += 1

        consensus = {}
        for hp in (1, 2):
            ref_c, alt_c = counts[hp]
            if ref_c + alt_c >= 2:
                consensus[hp] = 1 if alt_c > ref_c else 0
        return consensus


mean = lambda x: sum(x) / len(x)

def snp_beta_binomial_log_predictive(support_type, alpha, beta):
    """Beta-Binomial posterior predictive log probability for one REF/ALT SNP observation."""
    total = alpha + beta
    if total <= 0:
        return float("-inf")
    if support_type == "ALT":
        prob = alpha / total
    elif support_type == "REF":
        prob = beta / total
    else:
        return float("-inf")
    return math.log(max(prob, 1e-12))

def kimean(scores):

    center = lambda x: sum(x) / len(x)
    split_set = lambda s, m:[[i for i in s if i < m],[i for i in s if i >= m]]

    split_set2 = lambda s, m: [[i for i in s if abs(i - m[0]) < abs(i-m[1])], [i for i in s if abs(i-m[0]) >= abs(i-m[1])]]

    sz = [0, len(scores)]
    sets = split_set(scores, center(scores))
    while sz[0] != len(sets[0]):
        sz = len(sets[0]), len(sets[1])
        sets = split_set2(scores, [center(sets[0]), center(sets[1])])

    return sets




class Bridger(object):
    def __init__(self):
        self.options = {}
        self.options["margin"] = 100000

    def load_read_haplotype(self, fname):
        bf = pysam.AlignmentFile(fname)

        rtypes = {}
        for al in bf:
            if al.has_tag("HP") and al.has_tag("PS"):
                rtypes[al.qname] = (al.get_tag("PS"), al.get_tag("HP"))

        return rtypes

 
    def get_gap_between_blocks(self, blk0, blk1):
        if(blk0.chrom == blk1.chrom):
            return max(blk0.start(), blk0.end() - self.options["margin"]), min(blk1.start() + self.options["margin"], blk1.end())        
 

    def extend_block_in_vcf(self, ifname, ofname, links, switchs = {}):
        '''根据links标记得两个相邻得blocks扩展vcf'''
        # rewrite VCF
        #print(links) 
        haplotype = {} # block: block, hp
    
        for blk1, (blk0, flip) in sorted(links.items(), key = lambda x: int(x[0][1])):
            #print(blk0, blk1)
            assert int(blk1[1]) > int(blk0[1]) and blk1[0] == blk0[0]
            swt = 0 if blk0 not in switchs else len(switchs[blk0])
            #swt = 0
            if blk0 in haplotype:
                 blk2, flip1 = haplotype[blk0]
                 haplotype[blk1] = (blk2, flip1 ^ flip ^ (swt % 2 != 0))
            else:
                 haplotype[blk1] = (blk0, flip ^ (swt % 2 != 0))
   
        ivcf = pysam.VariantFile(ifname)
        ovcf = pysam.VariantFile(ofname, "w", header=ivcf.header)
        for rec in ivcf.fetch():
            #assert len(rec.samples) == 1
            if "PS" in rec.format and type(rec.samples[0]["PS"]) == int:
                 blk = (rec.chrom, rec.samples[0]["PS"])
                 swt = [] if blk not in switchs else switchs[blk]
                 gt = 0
                 if swt != None:
                    for s in swt:
                        if rec.pos > s:
                            gt += 1

                 if blk in haplotype:
                    nblk, flip = haplotype[blk]
                    assert blk[0] == nblk[0]  # chrom

                    rec.samples[0]["PS"] = nblk[1]
            
                    flip = flip ^ (gt % 2 == 1)
                    if flip:
                        gt = rec.samples[0]["GT"]
                        rec.samples[0]["GT"] = (gt[1], gt[0])
                        rec.samples[0].phased = True
                 else:
                    flip = (gt % 2 == 1)
                    if flip:
                        gt = rec.samples[0]["GT"]
                        rec.samples[0]["GT"] = (gt[1], gt[0])
                        rec.samples[0].phased = True
                                              
            ovcf.write(rec)
                 
class MethBridger(Bridger):
    def __init__(self, threads):
        # super().__init__()
        # self.options = {}
        # self.options["threads"] = threads
        # self.options["margin"] = 100000
        
        # # read 有多少meth就可以判断单倍型
        # self.options["read_count"] = 5
        # self.options["read_rate"] = 0.6
        # self.find_switch_time = 0
        # self.bridge_two_blocks_time = 0
        super().__init__()
        self.options = {
            "threads": threads,
            "margin": 100000,
            "read_count": 5,
            "read_rate": 0.6,
            "switch_min_error_reduction": 3,      
            "switch_min_support_reads": 5,       
            "switch_min_distance_bp": 2000,      
            "switch_min_length_weight": 1.0,
            "joint_phase": False
        }

   

    def bridge2(self, vcf_fname, bam_fname, ofname, detect_switch=True, joint_phase=False):
        self.options["vcf"] = vcf_fname
        self.options["joint_phase"] = joint_phase
        blocks = PhasedBlockSet.from_vcf(vcf_fname)
        print("start bridge2")  

        bf = pysam.AlignmentFile(bam_fname)
        jobs = [] # (job_type, job_data)
        if detect_switch:
            print("detect_switch")
            for blk in blocks:
                block_alignment = bf.fetch(blk.chrom, blk.start(), blk.end(), multiple_iterators=True)
                jobs.append([0, blk, len(list(block_alignment))])

        for blk0, blk1 in zip(blocks[:-1], blocks[1:]):

            if(blk0.chrom == blk1.chrom):
                start, end = max(blk0.start(), blk0.end() - self.options["margin"]), min(blk1.start() + self.options["margin"], blk1.end()) 
                block_alignment = bf.fetch(blk0.chrom, start, end, multiple_iterators=True)

                jobs.append([1, blk0, blk1, len(list(block_alignment))])

        jobs.sort(key = lambda x: -x[-1] )

        # import multiprocessing as mp
        # with mp.Pool(self.options["threads"]) as pool:
        #     result = pool.starmap(self.run_job, zip(jobs, [bam_fname]*len(jobs)), chunksize=1)
        
        import time 
        start_time_pool = time.time()  
        import multiprocessing as mp
        with mp.Pool(self.options["threads"]) as pool:
            result = pool.starmap(self.run_job, zip(jobs, [bam_fname]*len(jobs)), chunksize=1)
        end_time_pool = time.time()
        run_time_pool = end_time_pool - start_time_pool
        print("----------------------------------------------------------")
        print("Thread pool run time: {:.2f} seconds(the bridge2 )".format(run_time_pool))
        print("----------------------------------------------------------")

        links, switchs = {}, {}
        all_switches = []
        for job, rz in zip(jobs, result):
            if job[0] == 0:
                swt, blk = rz, job[1]
                switchs[(blk.chrom, blk.ps)] = swt
                for pos in swt:
                    all_switches.append((blk.chrom, blk.ps, pos))
            elif job[0] == 1:
                lnk, blk0, blk1 = rz, job[1], job[2]
                if lnk == 1:
                    links[(blk1.chrom, blk1.ps)] = ((blk0.chrom, blk0.ps), False)
                elif lnk == -1:
                    links[(blk1.chrom, blk1.ps)] = ((blk0.chrom, blk0.ps), True)
        # for job, rz in zip(jobs, result):
        #     if job[0] == 1:
        #         lnk, blk0, blk1 = rz, job[1], job[2]
        #         if lnk == 1:
        #             links[(blk1.chrom, blk1.ps)] = ((blk0.chrom, blk0.ps), False)
        #         elif lnk == -1:
        #             links[(blk1.chrom, blk1.ps)] = ((blk0.chrom, blk0.ps), True)


        self.extend_block_in_vcf(vcf_fname, ofname, links, switchs)


    # def run_job(self, job, bam_fname):
    #     logger.info("XXX start %d %d" % (job[0], job[1].ps))
    #     if job[0] == 0:
    #         result = self.find_switch_in_block(job[1], bam_fname)
    #     elif job[0] == 1:
    #         result =  self.bridge_two_blocks(job[1], job[2], bam_fname)

    #     logger.info("XXX end %d %d" % (job[0], job[1].ps))
    #     return result
    
    #统计job的运行时间

    def run_job(self, job, bam_fname):
        start_time = time.time()
        # if job[0]==0:
        #     logger.info("find_switch_in_block start %d %d %s" % (job[0], job[1].ps,job[1].chrom))
        # elif job[0]==1:
        #     logger.info("bridge_two_blocks start %d %d %s" % (job[0], job[1].ps,job[1].chrom))
        if job[0] == 0:
            # print("不执行find_switch_in_block")
            # result = None
            result = self.find_switch_in_block(job[1], bam_fname)
        elif job[0] == 1:
            result = self.bridge_two_blocks(job[1], job[2], bam_fname)
        end_time = time.time()
        # if job[0]==0:
        #     logger.info("find_switch_in_block end %d %d %s" % (job[0], job[1].ps,job[1].chrom))
        # elif job[0]==1:
        #     logger.info("bridge_two_blocks end %d %d %s" % (job[0], job[1].ps,job[1].chrom))
        elapsed_time = end_time - start_time
        logger.info("Job %d %d %s took %d seconds to complete" % (job[0], job[1].ps,job[1].chrom, elapsed_time))
        return result

    @staticmethod
    def are_two_blocks_adjacent(blk0, blk1):
        if blk0.chrom != blk1.chrom:
            return False
        # if not (blk0.start() <= blk0.end() and blk0.end() < blk1.start() and blk1.start() <= blk1.end()):
        #     print(blk0.start(), blk0.end(), blk0.end(), blk1.start(), blk1.start(), blk1.end())
        assert blk0.start() <= blk0.end() and blk0.end() < blk1.start() and blk1.start() <= blk1.end()
        return True

    def bridge_two_blocks(self, blk0, blk1, bam_fname):
        bf = pysam.AlignmentFile(bam_fname)
        gap = self.Gap(bf, blk0, blk1, self.options)
        return gap.bridge() 

    def find_switch_in_block(self, blk, bam_fname):
        bf = pysam.AlignmentFile(bam_fname)
        block = self.Block(bf, blk, self.options)
        return block.find_switch()

    class Gap(object):
        """Gap between two adjacent blocks"""
        def __init__(self, bf, blk0, blk1, opts):
            assert MethBridger.are_two_blocks_adjacent(blk0, blk1)
            self.bf = bf
            self.chrom = blk0.chrom
            self.blk0 = blk0
            self.blk1 = blk1
            self.opts = opts
            self.start, self.end = max(blk0.start(), blk0.end() - self.opts["margin"] ), min(blk1.start() + self.opts["margin"], blk1.end())        
            #self.info("position=%s:%d-%d" % (self.blk0.chrom, self.start, self.end))

        # def info(self, msg):
        #     logger.info("Gap(%d-%d): %s" %(self.blk0.ps, self.blk1.ps, msg))      
        def get_meths_linkscores(self, meths):
            #logger.info("get_meths_linkscores(Gap)")
            chrom, start, end, blk0, blk1 = self.chrom, self.start, self.end, self.blk0, self.blk1
            links = defaultdict(lambda: [0, 0, 0, 0])
            block_alignment = self.bf.fetch(self.chrom, start, end, multiple_iterators=True)
            for read in block_alignment:
               valid_meth = []
               for rloc, m in MethBridger.collect_meths_in_read(read):
                   if rloc in meths and meths[rloc].valid():
                       ps, hp, d = meths[rloc].test(m/256)
                       valid_meth.append((rloc, hp))
               for (p0, h0), (p1, h1) in itertools.combinations(valid_meth, 2):
                    assert p0 < p1 and (h0 == 1 or h0 == 2) and (h1 == 1 or h1 == 2)
                    links[(p0, p1)][h0-1 + (h1-1)*2] += 1
            return links
        
        def extract_consistent_meths(self, meths, links):
            print("extract_consistent_meths")
            def extract_top(links, cands, rate) :
                count = defaultdict(lambda: [0,0])
                for k, v in links.items():
                    if k[0] not in cands or k[1] not in cands: continue
                    c = v[0] + v[3], v[1] + v[2]
                    r = c[0] / sum(c)
                    if abs(r - (1-r)) >= 0.6 and sum(c) >= 10:
                        count[k[0]][0] += 1
                        count[k[1]][0] += 1
                    count[k[0]][1] += 1
                    count[k[1]][1] += 1
                scount = sorted(count.items(), key=lambda x: -x[1][0]/(x[1][1]+1))
            
                ext = set()
                for x in scount:
                    ext.add(x[0])
                    if len(ext) >= len(scount)*rate:
                            break
                return ext

            valid_locs = set([k[0] for k in links]) | set([k[1] for k in links])
            total = len(valid_locs)
            # valid_locs = extract_top(links, valid_locs, 0.5)
            # i = 1 
            # while i < 3:
            #     valid_locs = extract_top(links, valid_locs, 0.5)
            #     i = i + 1

            while len(valid_locs) > total*0.15 and len(valid_locs) >= 10:
                valid_locs = extract_top(links, valid_locs, 0.5)
            return valid_locs

        def meth_options(self):
            return 8, 0.3, 0.3
        
        def is_meth_cliffs_delta(self, meth):
            hap1_values = []
            hap2_values = []
            for qname, is_forward, meth_val, ps, hp in meth.infos:
                if ps == 0 or hp == 0:
                    continue
                if hp == 1:
                    hap1_values.append(meth_val)
                elif hp == 2:
                    hap2_values.append(meth_val)

            if len(hap1_values) >= 2 and len(hap2_values) >= 2:
                delta, size = cliffs_delta(hap1_values, hap2_values)
                return abs(delta) > 0.23
            else:
                return False
            
        # Keep the methylation-only path behaviorally identical to cliff_bayes.py.
        def bridge_methyl_only(self):
            meths0, meths1, unphased = self.collect_meths()
            meths3 = self.collect_middle_block_meths()
            cliffs_delta_meths1 = {}
            cliffs_delta_meths2 = {}
            extend_meths1 = {}
            extend_meths2 = {}
            meths_ASM = {}
            meths0_ASM = {}
            meths1_ASM = {}

            phased = []

            print("the length of the meths0, meths1, meths3", len(meths0), len(meths1), len(meths3))
            
            for key in meths0:
                if self.is_meth_cliffs_delta(meths0[key]):
                    cliffs_delta_meths1[key] = meths0[key]
            for key in meths1:
                if self.is_meth_cliffs_delta(meths1[key]):
                    cliffs_delta_meths2[key] = meths1[key]

            
       
            [v.verify_phased(*self.meth_options()) for v in itertools.chain(meths0.values(), meths1.values())]
            [v2.verify(*self.meth_options()) for v2 in itertools.chain(meths3.values())]

            print("cliffs_delta_meths1", len(cliffs_delta_meths1))
            print("cliffs_delta_meths2", len(cliffs_delta_meths2))

            extend_meths1 = self.extend_phased_reads_bayes_beta(cliffs_delta_meths1, unphased)
            extend_meths2 = self.extend_phased_reads_bayes_beta(cliffs_delta_meths2, unphased)
            [v.verify_cliffs_phased() for v in itertools.chain(extend_meths1.values(), extend_meths2.values())]
            # phased0 = self.get_phased_reads_bayes_beta(extend_meths1, unphased)
            # phased1 = self.get_phased_reads_bayes_beta(extend_meths2, unphased)
            self.extend_phased_reads(meths0, unphased)
            self.extend_phased_reads(meths1, unphased)

            print("extend_meths1", len(extend_meths1))
            print("extend_meths2", len(extend_meths2))
            # print("the length of the meths0, meths1, meths3", len(meths0), len(meths1), len(meths3))
            # links = self.get_meths_linkscores(meths3)
            # valid_locs = self.extract_consistent_meths(meths3, links)
            # for key in valid_locs:
            #     if key in meths3:
            #         meths_ASM[key] = meths3[key]

            # for key0 in extend_meths1:
            #     if key0 in meths_ASM:
            #         meths0_ASM[key0] = extend_meths1[key0]
            # print("length of the meths0_ASM:", len(meths0_ASM))
            # for key1 in extend_meths2:
            #     if key1 in meths_ASM:
            #         meths1_ASM[key1] = extend_meths2[key1]
            # print("length of the meths1_ASM:", len(meths1_ASM))               
            logger.info("------------------------------------------------------------------")
            logger.info("get_link")
            #return self.get_link(meths0, meths1)
            return self.get_link_ASM(extend_meths1, extend_meths2)
            #return self.get_link_ASM(extend_meths1, extend_meths2)
            #return self.get_link_by_spanning_reads(extend_meths1, extend_meths2)
            #return self.get_link_by_intersection(meths0_ASM, meths1_ASM, unphased)
            
        def bridge(self):
            if not self.opts.get("joint_phase", False):
                return self.bridge_methyl_only()

            meths0, meths1, unphased = self.collect_meths()
            snps0, snps1 = self.collect_snps()
            meths3 = self.collect_middle_block_meths()
            cliffs_delta_meths1 = {}
            cliffs_delta_meths2 = {}
            extend_meths1 = {}
            extend_meths2 = {}

            print("the length of the meths0, meths1, meths3", len(meths0), len(meths1), len(meths3))
            print("the length of snps0, snps1", len(snps0), len(snps1))
            print("the length of unphased(meth)", len(unphased))
            
            for key in meths0:
                if self.is_meth_cliffs_delta(meths0[key]):
                    cliffs_delta_meths1[key] = meths0[key]
            for key in meths1:
                if self.is_meth_cliffs_delta(meths1[key]):
                    cliffs_delta_meths2[key] = meths1[key]

            
       
            [v.verify_phased(*self.meth_options()) for v in itertools.chain(meths0.values(), meths1.values())]
            [v2.verify(*self.meth_options()) for v2 in itertools.chain(meths3.values())]

            print("cliffs_delta_meths1", len(cliffs_delta_meths1))
            print("cliffs_delta_meths2", len(cliffs_delta_meths2))

            # joint 方法从同一批初始可区分甲基化位点开始，不能复用 methyl-only 扩展后的模型
            joint_seed_meths1 = copy.deepcopy(cliffs_delta_meths1)
            joint_seed_meths2 = copy.deepcopy(cliffs_delta_meths2)

            extend_meths1 = self.extend_phased_reads_bayes_beta(cliffs_delta_meths1, unphased)
            extend_meths2 = self.extend_phased_reads_bayes_beta(cliffs_delta_meths2, unphased)
            [v.verify_cliffs_phased() for v in itertools.chain(extend_meths1.values(), extend_meths2.values())]
            [v.verify_phased(*self.meth_options()) for v in itertools.chain(extend_meths1.values(), extend_meths2.values())]

            # 第一阶段：保留现有甲基化独立分型流程
            phased_meth0 = self.get_phased_reads_bayes_beta(extend_meths1, unphased)
            phased_meth1 = self.get_phased_reads_bayes_beta(extend_meths2, unphased)
            meth_qnames0 = set([qname for qname, _ in phased_meth0])
            meth_qnames1 = set([qname for qname, _ in phased_meth1])

            remaining0 = set([read for read in unphased if read.qname not in meth_qnames0])
            remaining1 = set([read for read in unphased if read.qname not in meth_qnames1])

            print("extend_meths1", len(extend_meths1))
            print("extend_meths2", len(extend_meths2))
            print("meth_only_phased0", len(phased_meth0))
            print("meth_only_phased1", len(phased_meth1))
            print("meth_failed0", len(remaining0))
            print("meth_failed1", len(remaining1))

            self.meth_only0 = phased_meth0
            self.meth_only1 = phased_meth1

            print("joint_phase", "enabled")
            joint_unphased = self.collect_unphased_reads_for_joint(unphased)
            print("the length of unphased(joint_all)", len(joint_unphased))

            joint_meths1 = copy.deepcopy(joint_seed_meths1)
            joint_meths2 = copy.deepcopy(joint_seed_meths2)
            joint_snps0 = copy.deepcopy(snps0)
            joint_snps1 = copy.deepcopy(snps1)

            phased_joint0_raw, joint_meths1, joint_snps0, joint_debug0_raw = self.extend_phased_reads_bayes_beta_joint(
                joint_meths1, joint_snps0, joint_unphased,
                meth_weight=1.0, snp_weight=1.5,
                min_meth_sites=1, min_snp_sites=1,
                min_total_sites=2, min_margin=1.5
            )
            joint_debug0 = [("block0",) + row for row in joint_debug0_raw]

            phased_joint1_raw, joint_meths2, joint_snps1, joint_debug1_raw = self.extend_phased_reads_bayes_beta_joint(
                joint_meths2, joint_snps1, joint_unphased,
                meth_weight=1.0, snp_weight=1.5,
                min_meth_sites=1, min_snp_sites=1,
                min_total_sites=2, min_margin=1.5
            )
            joint_debug1 = [("block1",) + row for row in joint_debug1_raw]

            phased_joint0 = [(read.qname, sel) for read, sel in phased_joint0_raw]
            phased_joint1 = [(read.qname, sel) for read, sel in phased_joint1_raw]
            print("joint_iterative_phased0", len(phased_joint0))
            print("joint_iterative_phased1", len(phased_joint1))
            self.print_phasing_comparison("block0", phased_meth0, phased_joint0)
            self.print_phasing_comparison("block1", phased_meth1, phased_joint1)
            self.print_joint_evidence_summary("block0", joint_debug0_raw)
            self.print_joint_evidence_summary("block1", joint_debug1_raw)

            # 可选：导出联合rescue成功read的调试信息
            # 输出列包含 meth_n/snp_n/best_hap/margin，便于快速复盘
            joint_debug_rows = joint_debug0 + joint_debug1
            joint_debug_tsv = self.opts.get("joint_rescue_tsv")
            if joint_debug_tsv:
                if type(joint_debug_tsv) != str:
                    joint_debug_tsv = "joint_rescue_reads.tsv"
                with open(joint_debug_tsv, "w", encoding="utf-8") as fout:
                    fout.write("block\tround\tread_id\tmeth_n\tsnp_n\tbest_hap\tbest_hp\ttop1_score\ttop2_score\tmargin\n")
                    for block_name, round_idx, qname, meth_n, snp_n, best_hap, top1_score, top2_score, margin in joint_debug_rows:
                        fout.write(
                            "%s\t%d\t%s\t%d\t%d\t%s\t%s\t%f\t%f\t%f\n" %
                            (block_name, round_idx, qname, meth_n, snp_n, best_hap[0], best_hap[1], top1_score, top2_score, margin)
                        )
                print("joint_rescue_file", joint_debug_tsv, "rows", len(joint_debug_rows))

            # 保存两套独立分型结果，便于外部调试或后续流程比较
            self.joint0 = phased_joint0
            self.joint1 = phased_joint1
            self.joint_debug0 = joint_debug0
            self.joint_debug1 = joint_debug1
            self.final0 = phased_joint0
            self.final1 = phased_joint1

            logger.info("------------------------------------------------------------------")
            logger.info("get_link")
            return self.get_link_ASM(joint_meths1, joint_meths2)

        def get_link_ASM(self, meths0_ASM, meths1_ASM):
            # unzip self params
            chrom, start, end, blk0, blk1 = self.chrom, self.start, self.end, self.blk0, self.blk1
            print("start, end:",  blk0.end(), blk1.start())
    
            sup_count, sup_rate = self.opts["read_count"], self.opts["read_rate"]
    
            scores = defaultdict(int)
            block_alignment = self.bf.fetch(chrom, start, end, multiple_iterators=True)
            for read in block_alignment:
                sup0, sup1 = defaultdict(int), defaultdict(int)
                for rloc, m in MethBridger.collect_meths_in_read(read):

                    for meths, sup in ((meths0_ASM, sup0), (meths1_ASM, sup1)):
                        if rloc in meths and meths[rloc].valid_phased():
                            (ps, hp, d) = meths[rloc].test_phased(m/256)
                            if d >= 0.3:
                                sup[(ps, hp)] += 1

                    # for meths, sup in ((meths0_ASM, sup0), (meths1_ASM, sup1)):
                    #     if rloc in meths and meths[rloc].valid_cliffs():
                    #         (ps, hp, d) = meths[rloc].test_phased_cliff(m/256)
                    #         if d >= 0.3:
                    #             sup[(ps, hp)] += 1                            
                        
                # if sup0 is not None:
                #     print("sup0", sup0)
                # if sup1 is not None:
                #     print("sup1", sup1)
                sel0 = MethBridger.get_best_support(sup0, sup_count, sup_rate)
                # if sel0 is not None:
                #     print("sel0: ", sel0)
                sel1 = MethBridger.get_best_support(sup1, sup_count, sup_rate)
                # if sel1 is not None:
                #     print("sel1: ", sel1)
                if sel0 != None and sel1 != None: 
                    consist = sel0[1] == sel1[1]
                    scores[(sel0[0], sel1[0], consist)] += 1
            print("scores", scores)
            sel_scores = MethBridger.get_best_support(scores, 2, 0.5)
            print("sel_scores", sel_scores)
            if sel_scores != None:
                if sel_scores[2]: # consist
                    print("----------1")
                    return 1
                else:
                    print("---------- -1")
                    return -1
                    
            else:
                print("没有连接block")
                return 0     
            
        
        def extend_phased_reads(self, meths, unphased):
            phased = set()
            while True:
                phasing = self.phase_reads(meths, unphased - phased)
                if len(phasing) > 0:
                    self.add_phasing_reads(meths, phasing)
                    [v.verify_phased(*self.meth_options()) for v in meths.values()]
                    [phased.add(r[0]) for r in phasing]
                else:
                    break

        def extend_phased_reads_bayes_beta(self, meths, unphased):#返回一个meth字典就可以了
            phased = set()
            extend_meths = {}
            while True:
                phasing = self.phase_reads_bayes_beta(meths, unphased - phased)
                if len(phasing) > 0:
                    self.add_phasing_reads(meths, phasing)
                    for key in meths:
                        if self.is_meth_cliffs_delta(meths[key]):
                            extend_meths[key] = meths[key]                        
                    [phased.add(r[0]) for r in phasing]
                    print("extend_phased_reads_bayes_beta, phasing reads:", len(phasing))
                else:
                    break
            return extend_meths
           
        def get_phased_reads_bayes_beta(self, meths, unphased):
            logger.info("Gap_get_phased_reads")
            phased_tag = set()
            phased = []
            while True:
                phasing = self.phase_reads_bayes_beta(meths, unphased - phased_tag)

                if len(phasing) > 0:
                    for read, (ps ,hp)in phasing:
                        sel = (ps ,hp) 
                        phased_tag.add(read)
                        phased.append((read.qname, sel))
                else:
                    break
            return phased
        
        def collect_meths(self):
            # unzip self params
            chrom, start, end, blk0, blk1 = self.chrom, self.start, self.end, self.blk0, self.blk1

            meths0, meths1 = {}, {}
            unphased = set()
    
            block_alignment = self.bf.fetch(chrom, start, end, multiple_iterators=True)
            for read in block_alignment:
                ps, hp = 0, 0
                if read.has_tag("PS") and read.has_tag("HP"):
                    ps = read.get_tag("PS")
                    hp = read.get_tag("HP")
                    
                    for rloc, m in MethBridger.collect_meths_in_read(read):
                        for meths, blk in ((meths0, blk0), (meths1, blk1)):
                            if blk.ps != ps: continue
    
                            if rloc not in meths:
                                 meths[rloc] = MethBridger.Meth(blk.chrom, rloc)
                            meths[rloc].add(read.qname, read.is_forward, m/256, ps, hp)
        
                else:
                    if read.modified_bases != None: 
                        unphased.add(read)
            #self.info("collect meths, left=%d, right=%d, unphased=%d" % (len(meths0), len(meths1), len(unphased))) 
            return meths0, meths1, unphased

        def collect_snps(self):
            # 按block分开收集SNP：snps0对应blk0，snps1对应blk1
            # 位点范围使用start~end（覆盖两个block之间区域），再按read的PS归属到对应block
            chrom, start, end, blk0, blk1 = self.chrom, self.start, self.end, self.blk0, self.blk1
            snps0, snps1 = {}, {}
            self.unphased_snp_sites = {}
            vcf_fname = self.opts.get("vcf")
            if not vcf_fname:
                return snps0, snps1

            # 统一加载窗口内所有未定相SNP位点（包含中间区域）
            snp_sites = MethBridger.load_unphased_snp_sites_from_vcf(vcf_fname, chrom, start, end)
            self.unphased_snp_sites = snp_sites
            if not snp_sites:
                return snps0, snps1

            block_alignment = self.bf.fetch(chrom, start, end, multiple_iterators=True)
            for read in block_alignment:
                ps = read.get_tag("PS") if read.has_tag("PS") else 0
                hp = read.get_tag("HP") if read.has_tag("HP") else 0

                for snps, blk in ((snps0, blk0), (snps1, blk1)):
                    if blk.ps != ps:
                        continue
                    for snp, ref, alt, read_base, support_type in MethBridger.collect_snps_in_read(read, snp_sites):
                        if snp not in snps:
                            snps[snp] = SNP(chrom, snp, ref, alt)
                        snps[snp].add(read.qname, read.is_forward, read_base, support_type, ps, hp)

            return snps0, snps1

        def collect_unphased_reads_for_snp(self):
            # SNP使用独立的未定相read集合：
            # 只要read没有PS/HP标签（或标签为0）就纳入，不要求有甲基化信息
            chrom, start, end = self.chrom, self.start, self.end
            unphased_snp = set()
            block_alignment = self.bf.fetch(chrom, start, end, multiple_iterators=True)
            for read in block_alignment:
                ps = read.get_tag("PS") if read.has_tag("PS") else 0
                hp = read.get_tag("HP") if read.has_tag("HP") else 0
                if ps == 0 or hp == 0:
                    unphased_snp.add(read)
            return unphased_snp

        def collect_unphased_reads_for_joint(self, methyl_unphased=None):
            # joint 方法使用窗口内全部未定相read；不要求read带甲基化标记，SNP-only read也纳入
            chrom, start, end = self.chrom, self.start, self.end
            reads_by_qname = {}

            if methyl_unphased is not None:
                for read in methyl_unphased:
                    reads_by_qname.setdefault(read.qname, read)

            block_alignment = self.bf.fetch(chrom, start, end, multiple_iterators=True)
            for read in block_alignment:
                ps = read.get_tag("PS") if read.has_tag("PS") else 0
                hp = read.get_tag("HP") if read.has_tag("HP") else 0
                if ps == 0 or hp == 0:
                    reads_by_qname.setdefault(read.qname, read)

            return set(reads_by_qname.values())

        def unique_reads_by_qname(self, reads):
            reads_by_qname = {}
            for read in reads:
                reads_by_qname.setdefault(read.qname, read)
            return set(reads_by_qname.values())


            
        def collect_middle_block_meths(self):
            # unzip self params
            chrom, start, end, blk0, blk1 = self.chrom, self.start, self.end, self.blk0, self.blk1
            meths3 = {}           
            blk1_start = blk1.start()
            blk0_end = blk0.end()
            block_alignment = self.bf.fetch(chrom, blk0_end, blk1_start, multiple_iterators=True)
            # block_alignment = self.bf.fetch(chrom, start, end, multiple_iterators=True)
            for read in block_alignment:
                ps, hp = 0, 0
                if read.has_tag("PS") and read.has_tag("HP"):
                    ps = read.get_tag("PS")
                    hp = read.get_tag("HP")
                    for rloc, m in MethBridger.collect_meths_in_read(read):
                        if rloc not in meths3:
                            meths3[rloc] = MethBridger.Meth(blk0.chrom, rloc)
                        meths3[rloc].add(read.qname, read.is_forward, m/256, ps, hp)
                else:         
                    for rloc, m in MethBridger.collect_meths_in_read(read):
                        if rloc not in meths3:
                            meths3[rloc] = MethBridger.Meth(blk0.chrom, rloc)
                        meths3[rloc].add(read.qname, read.is_forward, m/256, ps, hp)                    
            return meths3                    

        def beta_log_likelihood(self, x, alpha, beta):
            import scipy.special as sp
            if x <= 0 or x >= 1:
                x = min(max(x, 1e-4), 1 - 1e-4)  # 避免 log(0)
            return (alpha - 1) * math.log(x) + (beta - 1) * math.log(1 - x) - math.lgamma(alpha) - math.lgamma(beta) + math.lgamma(alpha + beta)
        
        def phase_reads_bayes_beta(self, meths, unphased):
            # 甲基化独立判别：仅使用甲基化位点进行Beta-Bayes后验累加
            phasing = []

            for read in unphased:
                posterior = defaultdict(float)

                for rloc, m in MethBridger.collect_meths_in_read(read):
                    if rloc in meths and self.is_meth_cliffs_delta(meths[rloc]):
                        meth_obj = meths[rloc]
                        mval = m / 256
                        haps = meth_obj.get_phased_distribution_beta()
                        for ps, hp, alpha, beta in haps:
                            logp = self.beta_log_likelihood(mval, alpha, beta)
                            posterior[(ps, hp)] += logp

                if len(posterior) == 0:
                    continue

                best_hap, logprob = max(posterior.items(), key=lambda x: x[1])

                # 过滤低置信度结果（可调整阈值）
                if len(posterior) >= 3 and logprob < -5:
                    continue

                phasing.append((read, best_hap))

            return phasing

        def phase_reads_bayes_beta_snp(self, snps, unphased):
            # SNP独立判别：使用Beta-Binomial posterior predictive累加REF/ALT证据
            phasing = []
            if not snps:
                return phasing

            snp_sites = {pos: (snp.ref, snp.alt) for pos, snp in snps.items()}
            for read in unphased:
                posterior = defaultdict(float)

                for rloc, ref, alt, read_base, support_type in MethBridger.collect_snps_in_read(read, snp_sites):
                    # SNP编码规则：ALT=1, REF=0, OTHER直接丢弃
                    if support_type not in ("REF", "ALT"):
                        continue
                    snp_obj = snps.get(rloc)
                    if snp_obj is None:
                        continue
                    haps = snp_obj.get_phased_distribution_beta()
                    for ps, hp, alpha, beta in haps:
                        logp = snp_beta_binomial_log_predictive(support_type, alpha, beta)
                        posterior[(ps, hp)] += logp

                if len(posterior) == 0:
                    continue

                best_hap, logprob = max(posterior.items(), key=lambda x: x[1])
                if len(posterior) >= 3 and logprob < -5:
                    continue
                phasing.append((read, best_hap))

            return phasing

        def filter_joint_meths(self, meths):
            # joint迭代中每轮重新筛选可区分单倍型的甲基化位点
            filtered = {pos: meth for pos, meth in meths.items() if self.is_meth_cliffs_delta(meth)}
            for meth in filtered.values():
                meth.verify_cliffs_phased()
                meth.verify_phased(*self.meth_options())
            return filtered

        def is_snp_informative(self, snp, min_delta=0.3):
            # SNP需要两个单倍型都有支持，且REF/ALT编码均值存在足够差异
            means = snp.get_mean_per_hap()
            if 1 not in means or 2 not in means:
                return False
            return abs(means[1][1] - means[2][1]) >= min_delta

        def filter_joint_snps(self, snps):
            # joint迭代中每轮重新筛选可区分单倍型的SNP位点
            if not snps:
                return {}
            return {pos: snp for pos, snp in snps.items() if self.is_snp_informative(snp)}

        def extend_phased_reads_bayes_beta_joint(self, meths, snps, unphased,
                                                 meth_weight=1.0, snp_weight=1.5,
                                                 min_meth_sites=1, min_snp_sites=1,
                                                 min_total_sites=2, min_margin=1.5):
            # joint迭代分型：单轮joint判别 -> 回填甲基化/SNP模型 -> 重新筛选位点 -> 下一轮
            phased_qnames = set()
            phased_results = []
            debug_rows = []
            round_idx = 0
            all_unphased = self.unique_reads_by_qname(unphased)

            current_meths = self.filter_joint_meths(meths)
            current_snps = self.filter_joint_snps(snps)

            while True:
                round_idx += 1
                candidates = set([read for read in all_unphased if read.qname not in phased_qnames])
                if len(candidates) == 0:
                    break

                phasing = self.phase_reads_bayes_beta_joint(
                    current_meths, current_snps, candidates,
                    meth_weight=meth_weight, snp_weight=snp_weight,
                    min_meth_sites=min_meth_sites, min_snp_sites=min_snp_sites,
                    min_total_sites=min_total_sites, min_margin=min_margin
                )
                if len(phasing) == 0:
                    break

                self.add_phasing_reads(meths, phasing)
                self.add_phasing_reads_snp(snps, phasing, getattr(self, "unphased_snp_sites", None))

                for read, sel in phasing:
                    phased_qnames.add(read.qname)
                    phased_results.append((read, sel))

                for row in getattr(self, "joint_rescue_debug_rows", []):
                    debug_rows.append((round_idx,) + row)

                current_meths = self.filter_joint_meths(meths)
                current_snps = self.filter_joint_snps(snps)
                print(
                    "extend_phased_reads_bayes_beta_joint, round:",
                    round_idx,
                    "phasing reads:",
                    len(phasing),
                    "informative_meths:",
                    len(current_meths),
                    "informative_snps:",
                    len(current_snps)
                )

            self.joint_rescue_debug_rows = debug_rows
            return phased_results, current_meths, current_snps, debug_rows

        def print_phasing_comparison(self, block_name, meth_phased, joint_phased):
            meth_map = {qname: sel for qname, sel in meth_phased}
            joint_map = {qname: sel for qname, sel in joint_phased}
            meth_qnames = set(meth_map.keys())
            joint_qnames = set(joint_map.keys())
            both = meth_qnames & joint_qnames
            same = sum(1 for qname in both if meth_map[qname] == joint_map[qname])
            different = len(both) - same
            print(
                "phasing_compare",
                block_name,
                "meth_only",
                len(meth_qnames),
                "joint",
                len(joint_qnames),
                "both",
                len(both),
                "same",
                same,
                "different",
                different,
                "meth_only_unique",
                len(meth_qnames - joint_qnames),
                "joint_unique",
                len(joint_qnames - meth_qnames)
            )

        def print_joint_evidence_summary(self, block_name, joint_debug_rows):
            counts = defaultdict(int)
            for round_idx, qname, meth_n, snp_n, best_hap, top1_score, top2_score, margin in joint_debug_rows:
                if meth_n > 0 and snp_n > 0:
                    counts["methyl+snp"] += 1
                elif meth_n > 0:
                    counts["methyl_only"] += 1
                elif snp_n > 0:
                    counts["snp_only"] += 1
                else:
                    counts["no_evidence"] += 1

            print(
                "joint_evidence_types",
                block_name,
                "total",
                sum(counts.values()),
                "methyl+snp",
                counts["methyl+snp"],
                "methyl_only",
                counts["methyl_only"],
                "snp_only",
                counts["snp_only"],
                "no_evidence",
                counts["no_evidence"]
            )

        def phase_reads_bayes_beta_joint(self, meths, snps, unphased,
                                         meth_weight=1.0, snp_weight=1.5,
                                         min_meth_sites=1, min_snp_sites=1,
                                         min_total_sites=2, min_margin=1.5):
            # 联合判别：先分别累计甲基化与SNP证据，再按权重合并到同一个后验
            phasing = []
            debug_rows = []
            snp_sites = {pos: (snp.ref, snp.alt) for pos, snp in snps.items()} if snps else {}

            for read in unphased:
                meth_posterior = defaultdict(float)
                snp_posterior = defaultdict(float)
                meth_informative_locs = set()
                snp_informative_locs = set()

                # 甲基化证据：仅使用可判别的位点，累加每个(ps, hp)的log-likelihood
                for rloc, m in MethBridger.collect_meths_in_read(read):
                    if rloc not in meths:
                        continue
                    meth_obj = meths[rloc]
                    if not self.is_meth_cliffs_delta(meth_obj):
                        continue

                    haps = meth_obj.get_phased_distribution_beta()
                    if len(haps) == 0:
                        continue

                    meth_informative_locs.add(rloc)
                    mval = m / 256
                    for ps, hp, alpha, beta in haps:
                        meth_posterior[(ps, hp)] += self.beta_log_likelihood(mval, alpha, beta)

                # SNP证据：仅使用REF/ALT，按Beta-Binomial posterior predictive打分
                if len(snp_sites) > 0:
                    for rloc, ref, alt, read_base, support_type in MethBridger.collect_snps_in_read(read, snp_sites):
                        if support_type not in ("REF", "ALT"):
                            continue
                        snp_obj = snps.get(rloc)
                        if snp_obj is None:
                            continue
                        if not self.is_snp_informative(snp_obj):
                            continue

                        haps = snp_obj.get_phased_distribution_beta()
                        if len(haps) == 0:
                            continue

                        snp_informative_locs.add(rloc)
                        for ps, hp, alpha, beta in haps:
                            snp_posterior[(ps, hp)] += snp_beta_binomial_log_predictive(support_type, alpha, beta)

                # 合并后验：posterior[(ps, hp)] = meth_weight * meth + snp_weight * snp
                posterior = defaultdict(float)
                all_haps = set(meth_posterior.keys()) | set(snp_posterior.keys())
                for hap in all_haps:
                    posterior[hap] = meth_weight * meth_posterior.get(hap, 0.0) + snp_weight * snp_posterior.get(hap, 0.0)

                meth_n = len(meth_informative_locs)
                snp_n = len(snp_informative_locs)

                # 过滤：无后验/证据总数不足/两类证据都不足
                if len(posterior) == 0:
                    continue
                if meth_n + snp_n < min_total_sites:
                    continue
                if meth_n < min_meth_sites and snp_n < min_snp_sites:
                    continue

                ranked = sorted(posterior.items(), key=lambda x: -x[1])
                best_hap, top1_score = ranked[0]
                if len(ranked) >= 2:
                    top2_score = ranked[1][1]
                    margin = top1_score - top2_score
                else:
                    top2_score = float("-inf")
                    margin = top1_score

                # 过滤：top1-top2间隔过小，不输出低置信度结果
                if margin < min_margin:
                    continue

                phasing.append((read, best_hap))
                debug_rows.append((read.qname, meth_n, snp_n, best_hap, top1_score, top2_score, margin))
                print("joint_rescue_read", read.qname, "meth_n", meth_n, "snp_n", snp_n, "best_hap", best_hap, "margin", margin)

            self.joint_rescue_debug_rows = debug_rows
            return phasing

        def add_phasing_reads_snp(self, snps, phasing, snp_sites=None):
            # 将本轮新定相read回填到SNP模型，供下一轮继续迭代
            if len(phasing) == 0:
                return

            if snp_sites is None:
                snp_sites = {pos: (snp.ref, snp.alt) for pos, snp in snps.items()}
            if not snp_sites:
                return

            for read, (ps, hp) in phasing:
                for rloc, ref, alt, read_base, support_type in MethBridger.collect_snps_in_read(read, snp_sites):
                    # 仅使用REF/ALT更新模型，OTHER直接丢弃
                    if support_type == "OTHER":
                        continue
                    if rloc not in snps:
                        snps[rloc] = SNP(self.chrom, rloc, ref, alt)
                    snps[rloc].add(read.qname, read.is_forward, read_base, support_type, ps, hp)

        def extend_phased_reads_bayes_beta_snp(self, snps, unphased):
            # SNP迭代定相：单轮判别 -> 回填模型 -> 下一轮，直到没有新增read
            phased = set()
            phased_results = []
            while True:
                phasing = self.phase_reads_bayes_beta_snp(snps, unphased - phased)
                if len(phasing) == 0:
                    break

                self.add_phasing_reads_snp(snps, phasing)
                for read, sel in phasing:
                    phased.add(read)
                    phased_results.append((read.qname, sel))

                print("extend_phased_reads_bayes_beta_snp, phasing reads:", len(phasing))

            return phased_results

        def get_phased_reads_bayes_beta_snp(self, snps, unphased):
            logger.info("Gap_get_phased_reads_snp")
            # 输出格式与甲基化判别保持一致：[(qname, (ps, hp)), ...]
            # 采用SNP迭代定相，而不是单轮定相
            return self.extend_phased_reads_bayes_beta_snp(snps, unphased)

        def get_meth_posterior_detail_in_read(self, read, meths):
            # 计算单条read在甲基化模型下的后验及有效位点数
            posterior = defaultdict(float)
            informative_locs = set()
            for rloc, m in MethBridger.collect_meths_in_read(read):
                if rloc in meths and self.is_meth_cliffs_delta(meths[rloc]):
                    haps = meths[rloc].get_phased_distribution_beta()
                    if len(haps) > 0:
                        informative_locs.add(rloc)
                    mval = m / 256
                    for ps, hp, alpha, beta in haps:
                        posterior[(ps, hp)] += self.beta_log_likelihood(mval, alpha, beta)
            return posterior, len(informative_locs)

        def get_snp_posterior_detail_in_read(self, read, snps):
            # 计算单条read在SNP模型下的后验及统计信息
            posterior = defaultdict(float)
            informative_locs = set()
            ref_count, alt_count, other_count = 0, 0, 0
            if not snps:
                return posterior, 0, ref_count, alt_count, other_count

            snp_sites = {pos: (snp.ref, snp.alt) for pos, snp in snps.items()}
            for rloc, ref, alt, read_base, support_type in MethBridger.collect_snps_in_read(read, snp_sites):
                if support_type == "REF":
                    ref_count += 1
                elif support_type == "ALT":
                    alt_count += 1
                else:
                    other_count += 1
                    continue

                snp_obj = snps.get(rloc)
                if snp_obj is None:
                    continue
                haps = snp_obj.get_phased_distribution_beta()
                if len(haps) > 0:
                    informative_locs.add(rloc)
                for ps, hp, alpha, beta in haps:
                    posterior[(ps, hp)] += snp_beta_binomial_log_predictive(support_type, alpha, beta)

            return posterior, len(informative_locs), ref_count, alt_count, other_count

        def summarize_posterior_top2(self, posterior):
            # 返回后验前两名及margin，便于判断冲突是否由低置信度导致
            if len(posterior) == 0:
                return ("NA", "NA", 0.0, "NA", "NA", 0.0, 0.0)

            ranked = sorted(posterior.items(), key=lambda x: -x[1])
            (top_ps, top_hp), top_score = ranked[0]
            if len(ranked) >= 2:
                (second_ps, second_hp), second_score = ranked[1]
                margin = top_score - second_score
            else:
                second_ps, second_hp, second_score = "NA", "NA", float("-inf")
                margin = top_score

            return top_ps, top_hp, top_score, second_ps, second_hp, second_score, margin
        
        def phase_reads(self, meths, unphased):
            phasing = []
    
            remove_abnormal = lambda x: x if type(x) ==int else 0
            # params
            sup_count, sup_rate = self.opts["read_count"], self.opts["read_rate"]
            for read in sorted(unphased, key=lambda x: remove_abnormal(x.reference_start)):
                checks = defaultdict(int)
                for rloc, m in MethBridger.collect_meths_in_read(read):
                    if rloc in meths and meths[rloc].valid_phased():
                         (ps, hp, w) = meths[rloc].test_phased(m/256) 
                         if w >= 0.3: # TODO 检查  
                             checks[(ps, hp)] += 1
        
                sel = MethBridger.get_best_support(checks, sup_count, sup_rate)
                if sel != None:
                    phasing.append((read, sel)) #read和他read上每个甲基化位点的12值
            return phasing
        
        def add_phasing_reads(self, meths, phasing):
            for read, (ps, hp) in phasing:
           
                for rloc, m in MethBridger.collect_meths_in_read(read):
                    if rloc not in meths:
                        meths[rloc] = MethBridger.Meth(self.chrom, rloc) 
                    meths[rloc].add(read.qname, read.is_forward, m/256, ps, hp)

        def add_ASM(self, meths):
            chrom, start, end, blk0, blk1 = self.chrom, self.start, self.end, self.blk0, self.blk1
            meths3 = {}           
            blk1_start = blk1.start()
            blk0_end = blk0.end()
            block_alignment = self.bf.fetch(chrom, blk0_end, blk1_start, multiple_iterators=True)
            for read in block_alignment:
                valid_meth = []
                for rloc, m in MethBridger.collect_meths_in_read(read):
                    if rloc in meths and meths[rloc].valid():
                        ps, hp, d = meths[rloc].test(m/256) 
                        valid_meth.append((rloc, hp, ps, d))
                if len(valid_meth) != 0:
                    print("valid_meth", valid_meth)           

        def statistics_meths_ASM(self, meths0, meths1, meths0_ASM, meths1_ASM):
            chrom, start, end, blk0, blk1 = self.chrom, self.start, self.end, self.blk0, self.blk1    
            block_alignment = self.bf.fetch(chrom, start, end, multiple_iterators=True)
            num_meth0_ASM = len(meths0_ASM)
            num_meth1_ASM = len(meths1_ASM)
            for key, value  in meths0.items():
                if key in meths0_ASM:
                    if meths0[key].valid_phased():
                        num_meth0_ASM -= 1
                        meths0_ASM.pop(key)

            for key , value in meths1.items():
                if key in meths1_ASM:
                    if meths1[key].valid_phased():
                        num_meth1_ASM -= 1
                        meths1_ASM.pop(key)
            print("num_meth0_ASM, num_meth1_ASM", num_meth0_ASM, num_meth1_ASM)
            return meths0_ASM, meths1_ASM
    

        def get_link(self, meths0, meths1):
            # unzip self params
            chrom, start, end, blk0, blk1 = self.chrom, self.start, self.end, self.blk0, self.blk1
    
            sup_count, sup_rate = self.opts["read_count"], self.opts["read_rate"]
    
            scores = defaultdict(int)
            block_alignment = self.bf.fetch(chrom, start, end, multiple_iterators=True)
            for read in block_alignment:
                sup0, sup1 = defaultdict(int), defaultdict(int)
                for rloc, m in MethBridger.collect_meths_in_read(read):
                    for meths, sup in ((meths0, sup0), (meths1, sup1)):
                        if rloc in meths and meths[rloc].valid_phased:
                            (ps, hp, d) = meths[rloc].test_phased(m/256)
                            if d >= 0.3:
                                sup[(ps, hp)] += 1
                sel0 = MethBridger.get_best_support(sup0, sup_count, sup_rate)
                sel1 = MethBridger.get_best_support(sup1, sup_count, sup_rate)
                if sel0 != None and sel1 != None:
                    consist = sel0[1] == sel1[1]
                    scores[(sel0[0], sel1[0], consist)] += 1
            print(scores)
            sel_scores = MethBridger.get_best_support(scores, 2, 0.5)
            if sel_scores != None:
                if sel_scores[2]: # consist
                    return 1
                else:
                    return -1
            else:
                return 0          
    


    @staticmethod
    def collect_meths_in_read(read):
        meths = []
        if read.modified_bases != None:
            ref_loc = read.get_reference_positions(full_length=True)
            for m, locs in read.modified_bases.items():
                if m[0] == 'C' and m[2] == 'm':
                    for lc in locs:
                        if ref_loc[lc[0]] != None:  
                            rloc = ref_loc[lc[0]] + 1 if read.is_forward else ref_loc[lc[0]]
                            meths.append((rloc, lc[1]))
        return meths

    @staticmethod
    def is_biallelic_snp(rec):
        # 仅保留单碱基替换且只有一个ALT的SNP位点
        if rec.alts is None:
            return False
        if len(rec.ref) != 1:
            return False
        if len(rec.alts) != 1:
            return False
        if len(rec.alts[0]) != 1:
            return False
        return True

    @staticmethod
    def is_heterozygous_snp(rec):
        # 只有杂合位点才携带单倍型区分信息；0/0 和 1/1 等纯合位点会干扰分型
        if len(rec.samples) == 0:
            return False

        sample = list(rec.samples.values())[0]
        gt = sample.get("GT")
        if gt is None or len(gt) != 2:
            return False
        if any(allele is None for allele in gt):
            return False
        return set(gt) == {0, 1}

    @staticmethod
    def is_unphased_snp(rec):
        # 仅保留未定相且杂合的SNP；纯合位点不提供单倍型信息
        if not MethBridger.is_biallelic_snp(rec):
            return False
        if not MethBridger.is_heterozygous_snp(rec):
            return False

        sample = list(rec.samples.values())[0]
        phased = getattr(sample, "phased", False)
        return not phased

    @staticmethod
    def load_unphased_snp_sites_from_vcf(vcf_file, chrom, start, end):
        # 从指定区间加载未定相SNP，返回 {pos: (REF, ALT)}
        snp_sites = {}
        vf = pysam.VariantFile(vcf_file)
        for rec in vf.fetch(chrom, start - 1, end):
            if rec.pos < start or rec.pos > end:
                continue
            if MethBridger.is_unphased_snp(rec):
                snp_sites[rec.pos] = (rec.ref.upper(), rec.alts[0].upper())
        return snp_sites

    @staticmethod

    
    def collect_snps_in_read(read, snp_sites, only_alt=False):
        # 在read比对结果中提取SNP支持类型（REF/ALT/OTHER）
        snps = []
        if read.query_sequence is None:
            return snps

        aligned_pairs = read.get_aligned_pairs(matches_only=False)
        for qpos, rpos in aligned_pairs:
            if qpos is None or rpos is None:
                continue

            rloc = rpos + 1
            if rloc not in snp_sites:
                continue

            ref, alt = snp_sites[rloc]
            read_base = read.query_sequence[qpos].upper()

            if read_base == ref:
                support_type = "REF"
            elif read_base == alt:
                support_type = "ALT"
            else:
                support_type = "OTHER"

            if only_alt and support_type != "ALT":
                continue
            snps.append((rloc, ref, alt, read_base, support_type))

        return snps


    @staticmethod
    def get_best_support(sup, count, rate):
        '''support must be >= max(count, total*rate)
        '''
        
        sorted_sup = sorted(sup.items(), key = lambda x: -x[1])
        sum_sup = sum([i[1] for i in sorted_sup])
 
        if sum_sup > 0 and sorted_sup[0][1] >= max(sum_sup*rate, count):
            return sorted_sup[0][0]
        else:
           return None

    @staticmethod
    def get_phased_distribution(self):


        by_hap = defaultdict(list)
        for qname, is_fwd, mval, ps, hp in self.infos:
            if hp in (1, 2):
                by_hap[(ps, hp)].append(mval)

        result = []
        for (ps, hp), values in by_hap.items():
            if len(values) >= 2:
                mu = mean(values)
                sigma = stdev(values) if len(values) > 1 else 0.01
                result.append((ps, hp, mu, sigma))
        return result   
      
    @staticmethod
    def is_biallelic_snp(rec):
        if rec.alts is None: return False
        if len(rec.ref) != 1: return False
        if len(rec.alts) != 1 or len(rec.alts[0]) != 1: return False
        return True

    @staticmethod
    def is_phased_snp(rec):
        """关键修改：只检查 PS，不再要求 HP"""
        if not MethBridger.is_biallelic_snp(rec): return False
        sample = list(rec.samples.values())[0]
        return "PS" in rec.format and isinstance(sample.get("PS"), int)

    @staticmethod
    def load_phased_snp_sites_from_vcf(vcf_file, chrom, start, end):
        """只要求有 PS（已定相）即可"""
        snp_sites = {}
        vf = pysam.VariantFile(vcf_file)
        for rec in vf.fetch(chrom, start - 1, end):
            if rec.pos < start or rec.pos > end: continue
            if MethBridger.is_phased_snp(rec):
                snp_sites[rec.pos] = (rec.ref.upper(), rec.alts[0].upper())
        return snp_sites


    class Meth(object):
        def __init__(self, chrom, position):
            self.chrom = chrom
            self.position = position
            self.infos = []
            self.centroids = None
            self.centroids_phased = None
            self.centroids_cliffs = None
            pass
         
        @staticmethod
        def mean(x):

            ava = sum(x) / len(x)
            #return ava
            xlen = len(x)
            off = xlen // 5
            x_sorted = sorted(x)
            assert x_sorted[-1] >= x_sorted[0] 
            if ava < 0.5:
                x_core = x_sorted[0:xlen-off]
                #print(x, ava)
                #print(x_core, sum(x_core)/len(x_core))
                assert ava >= sum(x_core) / len(x_core)
                return sum(x_core) / len(x_core)
            else:
                x_core = x_sorted[off:]
                #print(x, ava)
                #print(x_core, sum(x_core)/len(x_core))
                assert ava <= sum(x_core) / len(x_core)
     
                return sum(x_core) / len(x_core)
 
        def add(self, qname, is_forward, mvalue, ps=0, hp=0):
            assert len(qname) > 10 and mvalue <= 1
            self.infos.append((qname, is_forward, mvalue, ps, hp))

            # TODO
            self.ps = ps
            
        def verify(self, sup_count, sup_rate, distance):
            signals = [i[2] for i in self.infos]
            splited_signals = kimean(signals)
            #print("verify", self.position, splited_signals[0], splited_signals[1])
            
            self.centroids = self._verify_two_signal_sets(splited_signals[0], splited_signals[1], sup_count, sup_rate, distance)

        def verify_phased(self, sup_count, sup_rate, distance):
            signals0 = [i[2] for i in self.infos if i[4] == 1]
            signals1 = [i[2] for i in self.infos if i[4] == 2]
            #print("verify_phase", self.position, signals0, signals1)
            self.centroids_phased = self._verify_two_signal_sets(signals0, signals1, sup_count, sup_rate, distance)

        def verify_cliffs_phased(self):
            signals0 = [i[2] for i in self.infos if i[4] == 1]
            signals1 = [i[2] for i in self.infos if i[4] == 2]
            #print("verify_phase", self.position, signals0, signals1)
            self.centroids_cliffs = self._verify_two_signal_sets_cliffs(signals0, signals1)

        def is_meth_cliffs_delta(self):
            hap1_values = []
            hap2_values = []
            for qname, is_forward, meth_val, ps, hp in self.infos:
                if ps == 0 or hp == 0:
                    continue
                if hp == 1:
                    hap1_values.append(meth_val)
                elif hp == 2:
                    hap2_values.append(meth_val)

            if len(hap1_values) >= 2 and len(hap2_values) >= 2:
                delta, size = cliffs_delta(hap1_values, hap2_values)
                return abs(delta) > 0.33
            else:
                return False

        def get_phased_distribution_beta(self):
            from collections import defaultdict

            by_hap = defaultdict(list)
            for qname, is_fwd, mval, ps, hp in self.infos:
                if hp in (1, 2):
                    by_hap[(ps, hp)].append(mval)

            result = []
            for (ps, hp), values in by_hap.items():
                if len(values) >= 2:
                    alpha, beta = self.get_beta_parameters(values)
                    result.append((ps, hp, alpha, beta))
            return result
        
        def _verify_two_signal_sets(self, signals0, signals1, sup_count, sup_rate, distance):
            '''验证两个集合是否是有效的两类，如果是则返回质心，否则返回None
               sup_count: 每一类的最小支持数目
               sup_rate:  每一类的最小支持比例
               distance: 中心的最小距离 [0,1]
            '''

            if len(signals0) > 0 and len(signals1) > 0:
                mean0, mean1 = self.mean(signals0), self.mean(signals1)
                size_total = len(signals0) + len(signals1)
                size_threshold = max(sup_count, sup_rate*size_total)
                #print("vvv", size_threshold, mean0, mean1, distance)
                if len(signals0) > size_threshold and len(signals1) >= size_threshold and abs(mean0-mean1) >= distance:
                    # 返回集合的质心
                    return [mean0, mean1]
            
            return None
        
        def _verify_two_signal_sets_cliffs(self, signals0, signals1):
            '''验证两个集合是否是有效的两类，如果是则返回质心，否则返回None
            '''

            if len(signals0) > 0 and len(signals1) > 0:
                mean0, mean1 = self.mean(signals0), self.mean(signals1)
                if len(signals0) >= 2 and len(signals1) >= 2:
                    return [mean0, mean1]
            
            return None
              

            

        def valid(self):
            return self.centroids != None

        def valid_phased(self):
            return self.centroids_phased != None
        
        def valid_cliffs(self):
            return self.centroids_cliffs != None

        def test(self, m):
            assert self.valid() and m <= 1
            return self._test(m, self.centroids)
      
        def test_phased(self, m):
            assert self.valid_phased() and m <= 1
            return self._test(m, self.centroids_phased)
        
        def test_phased_cliff(self, m):
            assert self.valid_cliffs() and m <= 1
            return self._test(m, self.centroids_cliffs)
        
      
        def _test(self, m, centroids):
            hp1_ave, hp2_ave = centroids
            t = m - (hp1_ave + hp2_ave) / 2

            if abs(hp1_ave - m) < abs(m - hp2_ave):
                return (self.ps, 1, abs(t))
            else:
                return (self.ps, 2, abs(t))

        def get_beta_parameters(self, values, pseudocount=0.5):
            n = len(values)
            if n < 2:
                return (1.0, 1.0)  # fallback to uniform
            
            m = mean(values)
            v = variance(values)

            # 修复边界情况
            if v == 0 or m == 0 or m == 1:
                # 完全甲基化或非甲基化或无变异
                if m <= 0.01:
                    return (0.5 + pseudocount, 5.0 + pseudocount)  # 非甲基化
                elif m >= 0.99:
                    return (5.0 + pseudocount, 0.5 + pseudocount)  # 甲基化
                else:
                    return (1.0 + pseudocount, 1.0 + pseudocount)  # 中性

            try:
                alpha = ((1 - m) / v - 1 / m) * m ** 2 + pseudocount
                beta = alpha * (1 / m - 1) + pseudocount
                return max(alpha, 0.1), max(beta, 0.1)
            except ZeroDivisionError:
                # 极端情况下 fallback
                return (1.0 + pseudocount, 1.0 + pseudocount)
            

    class Block(object):
        def __init__(self, bf, blk, opts):
            self.bf = bf
            self.blk = blk
            self.opts = opts

        def meth_options(self):
            return 10, 0.3, 0.4

        def collect_meths(self):
            #logger.info("collect_meths(Block)")
            bf, blk = self.bf, self.blk
            meths = {} # { Meth } 
            unphased = set()
            #print(blk.start(), blk.end())
            block_alignment = bf.fetch(blk.chrom, blk.start(), blk.end(), multiple_iterators=True)

            for read in block_alignment:
               hp = read.get_tag("HP") if read.has_tag("HP") else 0
               ps = read.get_tag("PS") if read.has_tag("PS") else 0
               
               for rloc, m in MethBridger.collect_meths_in_read(read):
                   if rloc not in meths:
                       meths[rloc] = MethBridger.Meth(blk.chrom, rloc)
                   meths[rloc].add(read.qname, read.is_forward, m/256, ps, hp)
               if ps == 0:
                   unphased.add(read)
            return meths, unphased
        
        def collect_snps_block(self):
            snps = {}
            vcf_fname = self.opts.get("vcf")
            if not vcf_fname: return snps

            snp_sites = MethBridger.load_phased_snp_sites_from_vcf(vcf_fname, self.blk.chrom, self.blk.start(), self.blk.end())
            if not snp_sites: return snps

            for read in self.bf.fetch(self.blk.chrom, self.blk.start(), self.blk.end(), multiple_iterators=True):
                ps = read.get_tag("PS") if read.has_tag("PS") else 0
                hp = read.get_tag("HP") if read.has_tag("HP") else 0
                for rloc, ref, alt, read_base, support_type in MethBridger.collect_snps_in_read(read, snp_sites):
                    if rloc not in snps:
                        snps[rloc] = SNP(self.blk.chrom, rloc, ref, alt)
                    snps[rloc].add(read.qname, read.is_forward, read_base, support_type, ps, hp)
            return snps
        
        def beta_log_likelihood(self, x, alpha, beta):
                if x <= 0 or x >= 1:
                    x = min(max(x, 1e-4), 1 - 1e-4)
                return (alpha - 1) * math.log(x) + (beta - 1) * math.log(1 - x) \
                    - math.lgamma(alpha) - math.lgamma(beta) + math.lgamma(alpha + beta)

        def test_switch_probabilistic(self, meths, snps=None):
            """SNP-only 模式：完全去掉甲基化位点，只使用 SNP 进行全局概率 switch 检测
               保留 read-length 加权 + 全局 HMM-like 打分 + 多重过滤"""
            if not snps:
                logger.warning("No SNP sites found in block, switch detection skipped.")
                return {}

            site_betas = {}
            for rloc, snp_obj in snps.items():
                if not snp_obj.valid_phased(): continue
                haps = snp_obj.get_phased_distribution_beta()
                betas = {hp: (alpha, beta) for ps, hp, alpha, beta in haps}
                if 1 in betas and 2 in betas:
                    site_betas[rloc] = betas

            if not site_betas: 
                logger.warning("No valid phased SNP sites after filtering.")
                return {}

            positions = sorted(site_betas.keys())
            improvement = defaultdict(lambda: [0.0, 0, 0.0])   # [total_delta, num_reads, total_length_weight]

            for read in self.bf.fetch(self.blk.chrom, self.blk.start(), self.blk.end(), multiple_iterators=True):
                if not (read.has_tag("HP") and read.has_tag("PS")): continue
                hp_current = read.get_tag("HP")
                if hp_current not in (1, 2): continue

                read_len = read.query_length or 1000
                length_weight = read_len / 1000.0                     # kb 为单位加权（长 read 贡献更大）

                # ==================== 只收集 SNP REF/ALT 观测 ====================
                obs = []
                if snps:
                    snp_sites = {pos: (snp.ref, snp.alt) for pos, snp in snps.items()}
                    for rloc, _, _, _, st in MethBridger.collect_snps_in_read(read, snp_sites):
                        if st not in ("REF", "ALT"): continue
                        if rloc in site_betas:
                            obs.append((rloc, st))

                if len(obs) < 5: continue
                obs.sort(key=lambda x: x[0])

                # 无 switch 全局 log-likelihood
                logp_no = sum(snp_beta_binomial_log_predictive(st, *site_betas[rloc][hp_current])
                              for rloc, st in obs if rloc in site_betas)

                # 对每个候选 switch 位置计算翻转后似然
                for p in positions:
                    if p not in [rloc for rloc, _ in obs]: continue
                    logp_switch = 0.0
                    hp_flipped = 3 - hp_current
                    for rloc, st in obs:
                        hp_use = hp_current if rloc <= p else hp_flipped
                        logp_switch += snp_beta_binomial_log_predictive(st, *site_betas[rloc][hp_use])

                    delta = logp_switch - logp_no
                    if delta > 0:
                        weighted_delta = delta * length_weight
                        improvement[p][0] += weighted_delta
                        improvement[p][1] += 1
                        improvement[p][2] += length_weight

            return improvement

        def test_switch_joint(self, meths, snps, valid_locs=None, meth_conf=0.3, snp_conf=0.3, min_sites=4):
            """Joint methylation+SNP switch test inside one phased block."""
            bf, blk = self.bf, self.blk
            improvement = defaultdict(lambda: [0, 0, 0, 0])
            snp_sites = {pos: (snp.ref, snp.alt) for pos, snp in snps.items()} if snps else {}

            block_alignment = bf.fetch(blk.chrom, blk.start(), blk.end(), multiple_iterators=True)
            for read in block_alignment:
                hap = {}

                for rloc, m in MethBridger.collect_meths_in_read(read):
                    if rloc in meths and meths[rloc].valid_phased() and (valid_locs is None or rloc in valid_locs):
                        ps, hp, d = meths[rloc].test_phased(m / 256)
                        if d >= meth_conf:
                            hap[rloc] = hp

                for rloc, ref, alt, read_base, support_type in MethBridger.collect_snps_in_read(read, snp_sites):
                    if support_type not in ("REF", "ALT"):
                        continue
                    snp = snps.get(rloc)
                    if snp is None or not snp.valid_phased():
                        continue
                    sval = 1.0 if support_type == "ALT" else 0.0
                    ps, hp, d = snp.test_phased(sval)
                    if d >= snp_conf:
                        hap[rloc] = hp

                if len(hap) < min_sites:
                    continue

                c = (
                    sum(1 for i in hap.values() if i == 1),
                    sum(1 for i in hap.values() if i == 2),
                )
                sorted_hap = sorted(hap.items())
                left_c = [0, 0]
                for p, v in sorted_hap[:-1]:
                    if v == 1:
                        left_c[0] += 1
                    elif v == 2:
                        left_c[1] += 1

                    right_c = [c[0] - left_c[0], c[1] - left_c[1]]
                    if sum(left_c) < 2 or sum(right_c) < 2:
                        continue
                    new_c = [left_c[0] + right_c[1], left_c[1] + right_c[0]]
                    improvement[p][0] += min(new_c) - min(c)
                    improvement[p][1] += sum(new_c)
                    improvement[p][2] += sum(left_c)
                    improvement[p][3] += sum(right_c)
            return improvement

        def extract_consistent_meths(self, meths, links):
            #logger.info("extract_consistent_meths(Block)")
            def extract_top(links, cands, rate) :
                count = defaultdict(lambda: [0,0])
                for k, v in links.items():
                    if k[0] not in cands or k[1] not in cands: continue
                    c = v[0] + v[3], v[1] + v[2]
                    r = c[0] / sum(c)
                    if abs(r - (1-r)) >= 0.6 and sum(c) >= 10:
                        count[k[0]][0] += 1
                        count[k[1]][0] += 1
                    count[k[0]][1] += 1
                    count[k[1]][1] += 1
                scount = sorted(count.items(), key=lambda x: -x[1][0]/(x[1][1]+1))
            
                ext = set()
                for x in scount:
                    ext.add(x[0])
                    if len(ext) >= len(scount)*rate:
                          break
                return ext
    
            valid_locs = set([k[0] for k in links]) | set([k[1] for k in links])
            total = len(valid_locs)
    
            while len(valid_locs) > total*0.15 and len(valid_locs) >= 10:
                valid_locs = extract_top(links, valid_locs, 0.5)
            
            return valid_locs

        def get_meths_linkscores(self, bf, blk, meths):
            logger.info("get_meths_linkscores(Block)")
            links = defaultdict(lambda: [0, 0, 0, 0])
            block_alignment = bf.fetch(blk.chrom, blk.start(), blk.end(), multiple_iterators=True)
            start_time = time.time()
            for read in block_alignment:
               valid_meth = []
               #counter_3 += 1
               for rloc, m in MethBridger.collect_meths_in_read(read):
                   if rloc in meths and meths[rloc].valid():
                       ps, hp, d = meths[rloc].test(m/256)
                       valid_meth.append((rloc, hp))
                       #counter_1 += 1
               for (p0, h0), (p1, h1) in itertools.combinations(valid_meth, 2):
                    assert p0 < p1 and (h0 == 1 or h0 == 2) and (h1 == 1 or h1 == 2)
                    links[(p0, p1)][h0-1 + (h1-1)*2] += 1
            end_time = time.time()
            run_time = end_time - start_time
            print("get_meths_linkscores run time: {:.2f} seconds".format(run_time))
            return links

        def compare_phased_meths(self, meths, valid_locs):
            #logger.info("compare_phased_meths(Block)")
            count = [0, 0, 0, 0]
            consist = {} 
            for pos, m in meths.items():
                if m.valid_phased():
                    count[0] += 1
                
                if m.valid(): 
                    count[1] += 1
    
                if m.valid_phased() and m.valid():
                    count[2] += 1
    
                if m.valid_phased() and pos in valid_locs:
                    count[3] += 1
            
                if pos in valid_locs or m.valid_phased():
                    consist[pos] = (pos in valid_locs, m.valid_phased())
            
            for p, c in sorted(consist.items()):
                print("CC", p, c)
            print(count, count[3] / len(valid_locs))

        def find_switch(self):
            logger.info("find_switch(Block)")
            bf, blk = self.bf, self.blk
            #block_alignment = bf.fetch(blk.chrom, blk.start(), blk.end(), multiple_iterators=True)
            meths = {} # Meth
            meths, unphased = self.collect_meths()
            snps = {}
            if self.opts.get("joint_phase", False):
                snps = self.collect_snps_block()
    
            for p, v in sorted(meths.items()):
                v.verify(*self.meth_options())
                v.verify_phased(*self.meth_options())

            #self.extend_phased_reads(meths, unphased)
            improvement = self.test_switch(meths) 


            for p, c in sorted(improvement.items()):
                imp = c[0]/c[1] if c[1] > 0 else 0
                # if True or imp < 0:
                #     print(p, c, imp)    
                    
            # links = self.get_meths_linkscores(bf, blk, meths)
            # valid_locs = self.extract_consistent_meths(meths, links)
            # sorted(valid_locs)
            # for v in sorted(valid_locs):
            #     print("valid", v)
            #self.compare_phased_meths(meths, valid_locs)
            if self.opts.get("joint_phase", False):
                improvement_1 = self.test_switch_joint(meths, snps)
                self.detect_switch(improvement_1)
    
            #improvement = self.test_switch(meths, valid_locs)     
            return self.detect_switch(improvement)
        
        def detect_switch(self, improvement):
            ranges = []
            rr = []
            for p, c in sorted(improvement.items()):
                imp = c[0]/c[1] if c[1] > 0 else 0
                if imp < 0:
                    rr.append((p, c)) 
                else:
                    if len(rr) > 0:
                        ranges.append(rr)
                        rr = []

            switchs = []            
            for rr in ranges:
                if len(rr) >= 3:
                    sorted_rr = sorted(rr, key = lambda x: -x[1][0])
                    for ir in sorted_rr:
                        print(ir)
                        print("ir")
                        if ir[1][0] <= -10 and ir[1][0] / ir[1][1] <= -0.03:
                            switchs.append(ir[0])
                            break
            print("detect_switch", switchs)
            return switchs

        def test_switch_fraction(self, meths, snps_beta, valid_locs=None, mismatch_threshold=0.50, beta_factor=0.25):
            """
            Fraction‑based mismatch method with stronger Beta effect (v5).

            Parameters
            ----------
            meths : dict
                Methylation sites.
            snps_beta : dict
                SNP Beta summaries.
            valid_locs : set, optional
                Positions to include.
            mismatch_threshold : float
                Base fraction threshold (default 0.50).
            beta_factor : float
                Beta adjustment factor (default 0.25).

            Returns
            -------
            defaultdict
                Mapping from position to [diff, weight].
            """
            bf, blk = self.bf, self.blk
            positions = set(meths.keys()) | set(snps_beta.keys())
            if not positions:
                return defaultdict(lambda: [0.0, 0.0])
            sorted_positions = sorted(positions)
            improvement_beta = self.test_switch(meths)
            snp_sites = None
            if snps_beta:
                snp_sites = {pos: (None, None) for pos in snps_beta.keys()}
            pair_counts = defaultdict(lambda: [0, 0])
            for read in bf.fetch(blk.chrom, blk.start(), blk.end(), multiple_iterators=True):
                assignments = []
                # methylation assignments
                for rloc, m in MethBridger.collect_meths_in_read(read):
                    if rloc in meths and meths[rloc].valid_phased() and (valid_locs is None or rloc in valid_locs):
                        ps_val, hp_val, d_val = meths[rloc].test_phased(m / 256)
                        if d_val >= 0.3:
                            assignments.append((rloc, hp_val))
                if snp_sites:
                    for rloc, ref, alt, read_base, support_type in MethBridger.collect_snps_in_read(read, snp_sites):
                        if support_type not in ("REF", "ALT"):
                            continue
                        if rloc not in snps_beta:
                            continue
                        sval = 1.0 if support_type == "ALT" else 0.0
                        means = snps_beta[rloc]["means"]
                        if 1 not in means or 2 not in means:
                            continue
                        mean1, mean2 = means[1], means[2]
                        mid = (mean1 + mean2) / 2.0
                        d_val = abs(sval - mid)
                        if d_val < 0.3:
                            continue
                        hp_val = 1 if abs(mean1 - sval) < abs(sval - mean2) else 2
                        assignments.append((rloc, hp_val))
                if len(assignments) < 2:
                    continue
                assignments.sort(key=lambda x: x[0])
                assign_dict = dict(assignments)
                for i in range(len(sorted_positions) - 1):
                    pos_i = sorted_positions[i]
                    pos_j = sorted_positions[i + 1]
                    if pos_i in assign_dict and pos_j in assign_dict:
                        hp_i = assign_dict[pos_i]
                        hp_j = assign_dict[pos_j]
                        if hp_i == hp_j:
                            pair_counts[(pos_i, pos_j)][0] += 1
                        else:
                            pair_counts[(pos_i, pos_j)][1] += 1
            improvement = defaultdict(lambda: [0.0, 0.0])
            for (pos_i, pos_j), (matches, mismatches) in pair_counts.items():
                weight = matches + mismatches
                if weight == 0:
                    continue
                mismatch_fraction = mismatches / weight
                base_diff = matches - mismatches
                # determine sign based on mismatch fraction
                if mismatch_fraction > mismatch_threshold:
                    diff = -base_diff
                else:
                    diff = base_diff
                # incorporate Beta effect from pos_j
                val_beta = improvement_beta.get(pos_j, [0.0, 0.0, 0.0, 0.0])
                diff_beta = val_beta[0]
                weight_beta = val_beta[1]
                beta_effect = (diff_beta / weight_beta) if weight_beta > 0 else 0.0
                diff += beta_factor * beta_effect * weight
                improvement[pos_j][0] += diff
                improvement[pos_j][1] += weight
            return improvement
                    

        def test_switch(self, meths, valid_locs=None):
            #logger.info("test_switch(Block)")
            bf, blk = self.bf, self.blk

            improvement = defaultdict(lambda: [0, 0, 0, 0])
            block_alignment = bf.fetch(blk.chrom, blk.start(), blk.end(), multiple_iterators=True)
            for read in block_alignment:
                hap = {}
                for rloc, m in MethBridger.collect_meths_in_read(read):
                    if rloc in meths and meths[rloc].valid_phased() and (valid_locs==None or rloc in valid_locs):
                        (ps, hp, d) = meths[rloc].test_phased(m/256)
                        if d >= 0.3: # TODO
                            hap[rloc] = hp
                if len(hap) < 4: continue
    
                c = sum([1 for i in hap.values() if i == 1]), sum([1 for i in hap.values() if i == 2])
                sorted_hap = sorted(hap.items()) # by position
                left_c = [0, 0]
                for p, v in sorted_hap[0:-1]:
                    if v == 1:
                        left_c[0] += 1
                    elif v == 2:
                        left_c[1] += 1
    
                    right_c = [c[0] - left_c[0], c[1] - left_c[1]]
                    if sum(left_c) < 2 or sum(right_c) < 2: continue 
                    new_c = [left_c[0] + right_c[1], left_c[1] + right_c[0]] # switch
                    improvement[p][0] += min(new_c) - min(c)
                    improvement[p][1] += sum(new_c)
                    improvement[p][2] += sum(left_c)
                    improvement[p][3] += sum(right_c)
            
            return improvement


        def extend_phased_reads(self, meths, unphased):
            #logger.info("extend_phased_reads(Block)")  
            phased = set()
            while True:
                phasing = self.phase_reads(meths, unphased - phased)
                if len(phasing) > 0:
                    #logger.info(str(len(phasing)))
                    self.add_phasing_reads(meths, phasing)
                    [v.verify_phased(*self.meth_options()) for v in meths.values()]
                    [phased.add(r[0]) for r in phasing]
                else:
                    break
           

        def phase_reads(self, meths, unphased):
            #logger.info("phase_reads(Block)")
            phasing = []
            remove_abnormal = lambda x: x if type(x) ==int else 0
            # params
            sup_count, sup_rate = self.opts["read_count"], self.opts["read_rate"]
            for read in sorted(unphased, key=lambda x: remove_abnormal(x.reference_start)):
                checks = defaultdict(int)
                for rloc, m in MethBridger.collect_meths_in_read(read):
                    if rloc in meths and meths[rloc].valid_phased():
                         (ps, hp, w) = meths[rloc].test_phased(m/256) 
                         if w >= 0.3: # TODO 检查  
                             checks[(ps, hp)] += 1
        
                sel = MethBridger.get_best_support(checks, sup_count, sup_rate)
                if sel != None:
                    phasing.append((read, sel))
            return phasing
        
        def add_phasing_reads(self, meths, phasing):
            for read, (ps, hp) in phasing:
           
                for rloc, m in MethBridger.collect_meths_in_read(read):
                    if rloc not in meths:
                        meths[rloc] = MethBridger.Meth(self.chrom, rloc) 
                    meths[rloc].add(read.qname, read.is_forward, m/256, ps, hp)
        

def main():
    import argparse
    parse = argparse.ArgumentParser(description='BayesPhase.')
    parse.add_argument('vcf', help='The VCF file.')
    parse.add_argument('bam', help='The BAM file.')
    parse.add_argument("out_vcf", help="The output VCF file")
    parse.add_argument("-t", "--threads", help="number of threads", default=1, nargs="?", const=1, type=int)
    parse.add_argument("-jointPhase", "--jointPhase", help="use methylation+SNP joint read phasing", action="store_true")
    parse.add_argument("--str", help="")
    args = parse.parse_args()

    logger.info("Input VCF: %s" % args.vcf)
    logger.info("Input BAM: %s" % args.bam)
    logger.info("Output VCF: %s" % args.out_vcf)

    bridger = MethBridger(args.threads)
    #bridger = SVBridger()
    #bridger.bridge_by_str(args.vcf, args.bam, args.str, args.out_vcf) 
    #bridger.check_blocks(args.vcf, args.bam, args.out_vcf) 
    bridger.bridge2(args.vcf, args.bam, args.out_vcf, joint_phase=args.jointPhase)
    #bridger.bridge2_in_memory(args.vcf, args.bam, args.out_vcf)
    print("---------------------------------------------------")
    #bridger.find_switch(args.vcf, args.bam)


if __name__ == '__main__':
    main()
