#!/usr/bin/env python

#Find cluster based on HMM profile

import argparse
import glob
import os
import subprocess as sp

hmm_path = "hmmsearch"


def run_command(cmd,ommit):
    """This function will execute a command in bash"""
    if ommit:
        try: process = sp.Popen(cmd,shell=True)
        except: pass
        process.communicate("Y\n")
        if process.wait() != 0: print("Error ocurred, but you chose to ommit it")
    else:
        try: process = sp.Popen(cmd,shell=True)
        except OSError: sys.exit("Error: Execution cmd failed")
        process.communicate("Y\n")
        if process.wait() != 0: sys.exit("ERROR: Execution cmd failed")

		
def create_folder(name):
    """Checks if a folder exists and if not creates it"""
    if not os.path.exists(name):
        cmd = "mkdir "+name
        try:
            run_command(cmd,False)
        except:
            print("Unable to create directory ",name)

def get_clusters(info,outfile,families):
    """ List potential clusters based on the HMM results """
    genes = list(info.keys())
    if len(genes) != 0:
        genes = sorted(genes,key=lambda x:int(x.split("_")[2]))

        groups = []
        new_group = [genes[0]]
        for gene in genes[1:]:
            c1 = gene.split("_")[1]
            c2 = new_group[-1].split("_")[1]
            p1 = int(gene.split("_")[2])
            p2 = int(new_group[-1].split("_")[2])
            diff = p1 - p2
            if c1 == c2:
                if diff <= 3:
                    new_group.append(gene)
                    to_add = False
                else:
                    to_add = True
            else:
                to_add = True
            if to_add:
                if len(new_group) >= 4:
                    groups.append(new_group)
                new_group = [gene]

        for group in groups:
            found_fams = set([])
            for gene in group:
                found_fams.add(info[gene])
            p = round(len(found_fams) / len(families) * 100)
            if p >= args.percentage:
                print(";".join(group),file=outfile)

def filter_HMM(HMM_file,evalthr1,evalthr2,limit_hits):
    """ obtain list of valid hits in the HMM search """
    codes = {}
    stats = {}
    for line in open(HMM_file):
        line = line.strip()
        if "#" not in line and "--" not in line:
            dades = line.split()
            if float(dades[4]) < evalthr1 and float(dades[7]) < evalthr2:
                if dades[0] not in codes:
                    codes[dades[0]] = dades[2]
                    stats[dades[0]] = float(dades[5])
                else:
                    if float(dades[5]) > stats[dades[0]]:
                        codes[dades[0]] = dades[2]
                        stats[dades[0]] = float(dades[5])
    if limit_hits != 0:
        codes2 = {}
        for code in codes:
            f = codes[code]
            if f not in codes2:
                codes2[f] = set([])
            codes2[f].add(code)

        codes = {}
        for f in codes2:
            genes = list(codes2[f])
            genes = sorted(genes,key=lambda x: stats[x],reverse=True)
            for gene in genes[:limit_hits]:
                codes[gene] = f

    return codes    

def get_families(HMMprofile):
    """ Obtain potential families from original profile file"""
    fams = set([])
    for line in open(HMMprofile):
        line = line.strip()
        if "NAME" in line:
            dades = line.split()
            fams.add(dades[1])
    return fams

parser = argparse.ArgumentParser(description="Search the presence of a cluster based on a HMM profile")
parser.add_argument("-i",dest="HMMprofile",action="store",required=True,help="HMM profile of the cluster you're searching for")
parser.add_argument("-p",dest="proteome",action="store",required=True,help="Proteome of the species you want to search the cluster in")
parser.add_argument("-o",dest="outfolder",action="store",default=".",help="Name for the output folder")
parser.add_argument("-e",dest="evalue",action="store",default=1e-05,type=float,help="Evalue filter")
parser.add_argument("-n",dest="numHits",action="store",default=0,type=int,help="Limit the number of hits taken, set to 0 if no limit is considered")
parser.add_argument("-c",dest="percentage",action="store",default=70,type=int,help="Percentage of protein families that need to be present for a cluster to be considered")
args = parser.parse_args()

create_folder(args.outfolder)

fams = get_families(args.HMMprofile)
HMM = args.HMMprofile.split("/")[-1].split(".")[0]
outFolder = args.outfolder
outfile = open(outFolder+"/results_search."+HMM+".txt","w")
cmd = hmm_path+" --tblout "+outFolder+"/search"+HMM+".tbl "+args.HMMprofile+" "+args.proteome+" >log"
run_command(cmd,False)
info = filter_HMM(outFolder+"/search"+HMM+".tbl",args.evalue,args.evalue,args.numHits)
get_clusters(info,outfile,fams)
outfile.close()

detected_clusters = [x.split()[0] for x in open(outFolder+"/results_search."+HMM+".txt")]
spe = args.proteome.split("/")[-1].split(".")[0]
cl = args.HMMprofile.split("/")[-1].split(".")[0]
if len(detected_clusters) != 0:
    print("Check file "+outFolder+"/results_search."+HMM+".txt for the resulting cluster")
else:
    print("No clusters detected")
 

