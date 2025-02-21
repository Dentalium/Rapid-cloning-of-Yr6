#!/usr/bin/env/ python3
#coding=utf-8

import argparse
import os.path
import re

import threading
import subprocess
import time


parser=argparse.ArgumentParser(prog="MutCandidator",
                               description="This script is to identify candidate contigs with vcf files produced by minipileup.")
parser.add_argument('-w', '--work_path', type=str, default='.', help='Work path')
parser.add_argument('-i', '--input', required=True, help='Input vcf, sorting is not a must.')
parser.add_argument('-n', '--min_mutants', type=int, 
                    default=2, dest="paran",
                    help='Minimum number of mutants to report a contig. Default is 2.')
parser.add_argument('-c', '--min_coverage', type=int, 
                    default=10, dest="parac",
                    help='Minimum coverage for position to be regarded. Default is 10.')
parser.add_argument('-a', '--max_ref_allele_freq', type=float, 
                    default=0.01, dest="paraa",
                    help='Maximum reference allele frequency of a SNP to be reported. Default is 0.01.')
parser.add_argument('-z', '--max_mutant_lines', type=int, 
                    default=2, dest="paraz",
                    help='Number mutant lines that are allowed to have SNV in same position. Default is 2.')
parser.add_argument('-m', '--param', type=int, 
                    default=0, dest="param",
                    help='Maximum SNP number per contig per mutant be able to hold after filter. Default is 0 (unlimited).')
parser.add_argument('-s', '--strict', action="store_true", help="Strict mode, with only C->T/G->A type permitted.")

args = parser.parse_args()

outputf=f"{args.input.split(sep='/')[-1]}.a{args.paraa}c{args.parac}z{args.paraz}.filter.vcf"
if args.strict:
    output_report=f"{args.input.split(sep='/')[-1]}.n{args.paran}a{args.paraa}c{args.parac}z{args.paraz}m{args.param}.s.report.txt"
else:
    output_report=f"{args.input.split(sep='/')[-1]}.n{args.paran}a{args.paraa}c{args.parac}z{args.paraz}m{args.param}.report.txt"


contigs={}
report=[]

def split_vcf_records(line) -> list:
    # fields: [0]: chr; [1]: pos; [3]: ref; [4]: alt; [7]: info; [8]: format; [9+]: samples
    fields = line.strip().split("\t")
    return fields

def read_vcf_records(sample) -> list:
    # infos: [0]:allele1(num); [1]:allele2(num); [2]:dp0(ref); [3]:dp1(alt1); ...
    infos=list(map(int, re.split(r"[/:,]",sample)))
    return infos

class SNP:
    def __init__(self, sample, pos, ref, alt, infos:list):
        self.sample=sample
        self.pos=pos
        self.ref=ref
        self.alt=alt
        self.infos=infos

    def report(self):
        gt1=self.infos[0]
        gt2=self.infos[1]
        dp=self.infos[2:]
        alleles=[self.ref]+self.alt

        result=f"\t\tPos: {self.pos}\t{self.ref}->"
        if gt1 == gt2:
            result+=f"{alleles[gt1]}\tdepth: {dp[gt1]}\n"
        else:
            result+=f"{alleles[gt1]}/{alleles[gt2]}\tdepth: {dp[gt1]}/{dp[gt2]}\n"

        return result


class DisplayProgress:
    def __init__(self, inputfile) -> None:
        self.inputfile=inputfile
        self.current_line=0
        self.thread=None
        # 终止进程的事件
        self.end=threading.Event()

    def stat_lines(self):
        line_num=int(subprocess.getoutput("wc -l %s" % self.inputfile).split()[0])
        print("stat finished!")
        return line_num
    
    def add_current(self):
        self.current_line+=1

    def display_progress(self, all_line:int, end_event, fresh_time=2):
        time.sleep(fresh_time)
        while self.current_line < all_line:
            print(f"File proccessed: {self.current_line/all_line*100:.2f}%",end="\r")
            if end_event.is_set():
                break
            time.sleep(fresh_time)

    def gen_thread(self,all_line="",fresh_time=2):
        if all_line == "":
            all_line=self.stat_lines()
        self.thread=threading.Thread(target=self.display_progress, args=(all_line, self.end, fresh_time))
        self.thread.start()
        
    def ter_thread(self):
        self.end.set()
        self.thread.join()


def filter_snp(obj):
    with open(args.input,"r") as src, open(outputf,"w") as dst:
        headerflag=False

        samples = {}
        for i in src:
            if i.startswith("#") and not i.startswith("##"):
                headerflag=True
                print("Header proccessed!")

                header = split_vcf_records(i)
                for j in range(9, len(header)):
                    header[j]=header[j].split(sep="/")[-1]
                
                dst.write("\t".join(header)+"\n")

                samplenum=len(header)-9

            if headerflag and not i.startswith("#"):
                fields = split_vcf_records(i)

                if fields[3]=="N":
                    continue

                bflag=0
                aflag=samplenum

                for j in range(9, len(header)):
                    sample_name=header[j]
                    this_sample=fields[j]

                    if this_sample.startswith("."):
                        samples[sample_name]={"field":"."}
                        bflag+=1

                    else:
                        samples[sample_name]={"field":this_sample, "infos":read_vcf_records(this_sample)}

                        if samples[sample_name]["infos"][0] != 0:
                            bflag+=1
                        
                        dp=sum(samples[sample_name]["infos"][2:])

                        if dp <= args.parac or samples[sample_name]["infos"][2]/dp >= args.paraa:
                            samples[sample_name]["field"]="."

                    if samples[sample_name]["field"] == ".":
                        aflag-=1

                if bflag <= 0.5*samplenum and aflag > 0 and aflag <= args.paraz:
                    filtered_line="\t".join(fields[0:9])
                    for j in samples:
                        filtered_line+=f"\t{samples[j]['field']}"
                    dst.write(f"{filtered_line}\n")

            obj.add_current()
        
        print("\nFiltering finished!")


# main body

myDP=DisplayProgress(args.input)

if os.path.exists(outputf):
    print(f"{outputf} already exitsts!")
else:
    print("filtering SNP...")
    myDP.gen_thread(fresh_time=10)
    filter_snp(myDP)
    myDP.ter_thread()

# debug
"""
print("filtering SNP...")
filter_snp()
print("exit")
exit()
"""

with open(outputf,"r") as file:
    for i in file:
        # header
        if i.startswith("#"):
            header = split_vcf_records(i)
            samplenum=len(header)-9
        # body
        else:
            fields = split_vcf_records(i)
            fields[4]=fields[4].split(sep=",")
            this_contig=fields[0]
            
            if this_contig not in contigs.keys():
                contigs[this_contig]={}
            
            for j in range(9, samplenum+9):
                this_sample=fields[j]
                sample_name=header[j]
                if this_sample != ".":
                    # strict mode
                    if args.strict:
                        this_sample_infos=read_vcf_records(this_sample)
                        alt_allele=this_sample_infos[0]-1
                        if not ((fields[3]=="C" and fields[4][alt_allele]=="T") or (fields[3]=="G" and fields[4][alt_allele]=="A")):
                            continue

                    if sample_name not in contigs[this_contig]:
                        contigs[this_contig][sample_name]=[]

                    contigs[this_contig][sample_name].append(SNP(sample=sample_name, pos=fields[1], ref=fields[3], alt=fields[4], infos=read_vcf_records(fields[j])))


print("selecting candidates...")

for contig_name,this_contig in contigs.items():
    flag=True
    for i in this_contig.values():
        if args.param > 0 and len(i) >= args.param:
            flag=False
            break

    # paran
    if flag and len(this_contig) >= args.paran:
        #print(this_contig)
        report.append(f"Candidate contig: {contig_name}\n{len(this_contig)} mutants\n")
        for sample in this_contig.keys():
            report.append(f"\t{sample}:\n")
            for i in this_contig[sample]:
                report.append(i.report())
            report.append("\n")

with open(output_report,"w") as file:
    for i in report:
        file.write(i)

print(f"Done! Report is {output_report}")
