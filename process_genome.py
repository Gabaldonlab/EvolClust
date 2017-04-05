#!/usr/bin/env python 

"""
  Evolclust - automated pipeline to detect regions of conserved gene
  order using comparative genomics.
  Copyright (C) 2017 - Marina Marcet-Houben, Toni Gabaldon
  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

#Process genome and obtain the data needed from it to make the following analysis

import argparse
import glob
import os
import subprocess as sp
from Bio import SeqIO

#Checks if a folder exists and if not creates it		
def create_folder(name):
	if not os.path.exists(name):
		cmd = "mkdir "+name
		try:
			run_command(cmd,False)
		except:
			print "Unable to create directory ",name

#Run a command in linux
def run_command(cmd,ommit):
	if ommit:
		try: process = sp.Popen(cmd,shell=True)
		except: pass
		process.communicate("Y\n")
		if process.wait() != 0: print "Error ocurred, but you chose to ommit it"
	else:
		try: process = sp.Popen(cmd,shell=True)
		except OSError,e: sys.exit("Error: Execution cmd failed")
		process.communicate("Y\n")
		if process.wait() != 0: sys.exit("ERROR: Execution cmd failed")

#Prints a sequence into a file
def print_sequence(code,sequence,outfile):
	print >>outfile,">"+code
	i = 0
	if sequence[-1] == "*":
		sequence = sequence[:-1]
	while i < len(sequence):
		print >>outfile,sequence[i:i+60]
		i += 60

#Parses a .gb file and obtains the gff and the proteome re-coded in a gene-order analysis format
def parse_gb(gbFileName,outfileGFF,outfileProt,outfileConversion,tag):
	outfile = open(outfileGFF,"w")
	outfileP = open(outfileProt,"w")
	numContig = 1
	numProt = 1
	conversion = {"contigs":{},"proteins":{}}
	for record in SeqIO.parse(gbFileName,"gb"):
		contig = record.name
		species = record.annotations["organism"]
		sequence = str(record.seq)
		baseName = "%s_%.4d" % (tag,numContig)
		conversion["contigs"][contig] = str(numContig)
		numContig += 1
		source = "NCBI"
		for feature in record.features:
			if feature.type == "CDS":
				geneID = None
				transID = None
				if "gene" in feature.qualifiers:
					geneID = feature.qualifiers["gene"][0]
				if "locus_tag" in feature.qualifiers:
					transID = feature.qualifiers["locus_tag"][0]
				if not transID:
					pass
				elif "translation" in feature.qualifiers:
					protName = "%s_%.5d" %(baseName,numProt)
					numProt += 1
					comment = 'gene_id "%s"; transcript_id "%s"' % (geneID,transID)
					if "product" in feature.qualifiers:
						comment += '; protein_id "%s"' % feature.qualifiers["product"][0]
					elif "protein_id" in feature.qualifiers:
						comment += '; protein_id "%s"' % feature.qualifiers["protein_id"][0]
					if int(feature.strand) > 0: 
						frame = "+"
					else:
						frame = "-"
					positions = []
					for position in feature.location.parts:
						start = position.start + 1
						stop = position.end
						p = "%.10d-%.10d" %(start,stop)
						positions.append(p)
					positions = sorted(positions,key=lambda x:int(x.split("-")[0]))
					for p in positions:
						print >>outfile,contig+"\t"+source+"\t"+feature.type+"\t"+str(int(p.split("-")[0]))+"\t"+str(int(p.split("-")[1]))+"\t.\t"+frame+"\t.\t"+comment+"; conversion "+protName
					print_sequence(protName,feature.qualifiers["translation"][0],outfileP)
					conversion["proteins"][protName] = comment					
	outfile.close()
	outfileP.close()
	outfile = open(outfileConversion,"w")
	print >>outfile,tag+"\t"+species
	for code in conversion["contigs"]:
		print >>outfile,"Contigs\t"+conversion["contigs"][code]+"\t"+code
	proteins = conversion["proteins"].keys()
	proteins = sorted(proteins,key=lambda x:int(x.split("_")[2]))
	for code in proteins:
		print >>outfile,"Proteins\t"+code+"\t"+conversion["proteins"][code]
	outfile.close()
					
parser = argparse.ArgumentParser(description="Process a genome to the needed format to obtain the known conserved clusters")
parser.add_argument("-g","--gb",dest="gbFileName",action="store",required=True,help="Genbank fileName")
parser.add_argument("-t","--tag",dest="tag",action="store",required=True,help="TAG for the species in order to name the proteins and outfiles")
parser.add_argument("-d","--dirName",dest="dirProt",action="store",required=True,help="Directory where the proteome and conversion files will be stored")
args = parser.parse_args()

##Creates needed directories if not present
tag = args.tag
dirName = args.dirProt+"/"
create_folder(dirName)
dirGFF = dirName+"gff/"
create_folder(dirGFF)
dirProteome = dirName+"proteomeFiles/"
create_folder(dirProteome)
dirConversion = dirName+"conversionFiles/"
create_folder(dirConversion)

#Parse .gb to obtain the relevant files
parse_gb(args.gbFileName,dirGFF+"/"+tag+".gff",dirProteome+"/"+tag+".aa.txt",dirConversion+"/"+tag+".txt",tag)

			

