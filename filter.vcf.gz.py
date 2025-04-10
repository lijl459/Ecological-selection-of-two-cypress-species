#!/usr/bin/env python
import re, os, gzip, sys
if len(sys.argv) != 3:
    print("Usage: <Script> <in.vcf or vcf.gz> <out.vcf or vcf.gz>")
    exit()
k = re.compile('\s+')
inroot = str(sys.argv[1])
outroot = str(sys.argv[2])

if inroot.endswith("vcf.gz"):
    infile = gzip.open(inroot, 'r')
elif inroot.endswith("vcf"):
    infile = open(inroot, 'r')
else:
    print("Plese import vcf or vcf.gz file")
    exit()

if outroot.endswith("vcf.gz"):
    outfile = gzip.open(outroot, 'w')
elif outroot.endswith("vcf"):
    outfile = open(outroot, 'w')
else:
    print("Plese export vcf or vcf.gz file")
    exit()
    

dic = {}
st = set()

for line in infile:
    if line[0] == "#":
        pass
    else:
        l = line.split("\t")
        if  "INDEL" in l[7]:
            if l[0] in dic.keys():
                pass
            else:
                dic[l[0]] = []
            for a in range(int(l[1])-5, int(l[1])+6):
                dic[l[0]].append(a)
        else:
            pass

infile.seek(0, 0)

for line in infile:
    if line[0] == "#":
        outfile.write(line)
    else:
        line = line.strip()
        l = line.split("\t")
        if len(l[3]) >1:
            pass
        elif len(l[4]) > 1:
            pass
        elif l[0] in dic.keys() and int(l[1]) in dic[l[0]]:
            pass
        elif float(l[5]) <30:
            pass
        else:
            d = l[8].split(":").index("DP")
            for i in range(9, len(l)):
                dp = l[i].split(":")[d]
                if dp.isdigit():
                    if int(dp) <10:
                        l[i] = "./."
                    else:
                        pass
                else:
                    pass
            outline = '\t'.join(l) + '\n'
            outfile.write(outline)

infile.close()
outfile.close()
            
            
