#!/usr/bin/env/ python3
#coding=utf-8
'''
vcfChopperMultiplerv2
python=3.12.2
pysam=0.22.1
'''
import pysam
import multiprocessing
import argparse
import sys
import os

parser=argparse.ArgumentParser(description="This script splits the vcf file into overlapping intervals. \
                               Output is unsorted vcf file, which will be written to stdout.")
parser.add_argument("-i", required=True, help="Input bgziped vcf file, must end with an extension \".vcf.gz\". Must be sorted and indexed.")
parser.add_argument("-l", type=int, default=50000, help="Length of interval, default is 50000.")
parser.add_argument("-p", type=int, default=5000, help="Length of overlap, default is 5000.")
parser.add_argument("-t", type=int, default=0, help="Offset of the first interval, default is 0.")
#parser.add_argument("-o", default="", help="Prefix of output chopped vcf file, default is same as input file.")
parser.add_argument("-f", default="", help="File of a list of chromosomes to process parallelly. Choromesomes not mentioned will be deprecated.")
parser.add_argument("-j", type=int, default=1, help="Number of blocks to be processed parallely. Default is 1 (with out parallely processing).")

args=parser.parse_args()

if args.f and args.j !=1:
    parser.error("You can't define a chromosome list while using parameter j!")

#if not args.o:
#    args.o=args.i[0:-4]

paral=args.l # 50000 default
parap=args.p # 5000 default
offset=args.t # 0 default

paraf=args.f
paraj=args.j

input_file=args.i
#output_file=args.o+".choped.vcf"

def chop_contigs(inputdict:dict):
    contigs_chopped={}
    for i in inputdict:
        contigs_chopped[i]={}
        contig_length=inputdict[i]
        cumu_length=0
        
        if contig_length <= 2*paral-parap+offset:
            contigs_chopped[i][f"{i}_00"]=[1, contig_length]
        else:
            while contig_length-cumu_length >= paral-parap:
                if cumu_length==0:
                    contigs_chopped[i][f"{i}_0"]=[1 ,paral+offset]
                    cumu_length=cumu_length+paral+offset-parap
                else:
                    contigs_chopped[i][f"{i}_{cumu_length}"]=[cumu_length + 1 ,cumu_length + paral]
                    cumu_length=cumu_length+paral-parap 

            contigs_chopped[i][list(contigs_chopped[i].keys())[-1]][1]=contig_length

    return contigs_chopped

def contigLen(inputlist:list):
    return inputlist[1]-inputlist[0]+1

def modify_header(vcffile):
    contigs={}
    with pysam.VariantFile(vcffile) as vcf_in:
        for i in vcf_in.header.contigs:
            contigs[i]=vcf_in.header.contigs[i].length
        new_header=pysam.VariantHeader()
        for i in vcf_in.header.records:
            if i.key != 'contig':
                #print(i)
                new_header.add_record(i)
        new_header.add_samples(vcf_in.header.samples)

    chr_list=[]
    if paraf:
        with open(paraf, "r") as file:
            for i in file:
                chr_list.append(i.rstrip())

        contigs={k: v for k, v in contigs.items() if k in chr_list}

    contigs_chopped=chop_contigs(contigs)
    
    print("Chopping calculation finished!", file=sys.stderr)

    for i in contigs_chopped:
        for j in contigs_chopped[i]:
            #print(j, contigs_chopped[i][j])
            new_header.contigs.add(j, contigLen(contigs_chopped[i][j]))
    metatext=f"paral={paral}, parap={parap}, offset={offset}"
    if not paraf:
        metatext=metatext+f", parallel={paraj}"
    new_header.add_meta("ChoppingParameters", metatext)

    print("Header processed!", file=sys.stderr)

    return (contigs, chr_list, new_header, contigs_chopped)

def modify_records(vcf:pysam.VariantFile, oldchr, start, end, newchr):
    recs=[]
    for rec in vcf.fetch(oldchr, start-1, end):
        rec.contig=newchr
        try:
            rec.pos=rec.pos-start+1
        except ValueError:
            #print("ValueError!", rec, "\n\trange:", start, end, file=sys.stderr)
            pass
        recs.append(rec)
    return recs

def split_chr_group(chrs:list,number):
    grouped_chrs=[]
    chrs_length=len(chrs)
    group_size=chrs_length//number
    group_remain=chrs_length%number

    if group_size==0:
        print(f"Error!\n\tParameter j shouldn't be larger than {chrs_length}", file=sys.stderr)
        exit(1)

    for _ in range(number):
        grouped_chrs.append([])
    
    while group_size > 0:
        for i in range(number):
            grouped_chrs[i].append(chrs.pop(0))
        group_size-=1
    
    for i in range(group_remain):
        grouped_chrs[i].append(chrs.pop(0))

    return grouped_chrs

def mtp_vcf(tmpfile, vcffile, this_chrs:list, contigs_chopped:dict, new_header:pysam.VariantHeader):
    with pysam.VariantFile(vcffile) as vcf_in:
        with open(tmpfile, "w") as vcf_tmp:
            vcf_tmp.write(str(new_header))
            for i in this_chrs:
                for j in contigs_chopped[i]:
                    vcf_in.header.contigs.add(j,paral*2)#;print(tmpfile,i,j,file=sys.stderr)
                    this_range=contigs_chopped[i][j]
                    modified_recs=modify_records(vcf=vcf_in,oldchr=i,start=this_range[0],end=this_range[1],newchr=j)
                    for k in modified_recs:
                        vcf_tmp.write(str(k))

def main(vcffile):
    contigs, chr_list, new_header, contigs_chopped=modify_header(vcffile)

    if (not args.f) and paraj==1:
        print("Single process", file=sys.stderr)
        with pysam.VariantFile(vcffile) as vcf_in:
            with open(sys.stdout.fileno(), "w") as vcf_out:
                vcf_out.write(str(new_header))
                for i in contigs_chopped:
                    for j in contigs_chopped[i]:
                        vcf_in.header.contigs.add(j,paral*2)
                        this_range=contigs_chopped[i][j]
                        modified_recs=modify_records(vcf=vcf_in,oldchr=i,start=this_range[0],end=this_range[1],newchr=j)
                        
                        for rec in modified_recs:
                            vcf_out.write(str(rec))
    else:
        if paraj!=1:
            chr_list=split_chr_group(chrs=list(contigs.keys()),number=paraj)
        else:
            chr_list=[[x] for x in chr_list]

        #print(chr_list, file=sys.stderr)

        processes=[]
        tmpfiles=[]
        for i in range(paraj):
            tmpfiles.append(f"{i}.tmp.vcf")
            processes.append(multiprocessing.Process(target=mtp_vcf, args=(tmpfiles[i], input_file, chr_list[i], contigs_chopped, new_header)))

        for p in processes:
            p.start()
        for p in processes:
            p.join()

        print("Writing output...", file=sys.stderr)        
        with pysam.VariantFile("-","w",header=new_header) as vcf_out:
            for i in tmpfiles:
                with pysam.VariantFile(i,"r") as tmp_vcf:
                    for j in tmp_vcf.fetch():
                        vcf_out.write(j)

        print("Removing tmp...", file=sys.stderr)
        for i in tmpfiles:
            os.remove(i)

if __name__ == '__main__':
    print(f"parameters:\n\tinterval length: {paral}\n\tinterval overlap: {parap}\n\toffset of first interval: {offset}", file=sys.stderr)
    print("processing...", file=sys.stderr)
    main(vcffile=input_file)
    print("Done!", file=sys.stderr)
