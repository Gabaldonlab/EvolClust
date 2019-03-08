#!/usr/bin/env python

"""
  Evolclust v1.0 - automated pipeline to detect regions of conserved gene
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

import argparse
import glob
import numpy
import os
from random import randint
import subprocess as sp
import sys

#For non parametric tests uncomment this part
#import rpy2
#import rpy2.robjects as robjects
#from rpy2.robjects.packages import importr


########################################################################
# Operational modules
########################################################################

#Checks whether a folder exists or else creates it
def create_folder(name):
	if not os.path.exists(name):
		cmd = "mkdir "+name
		try:
			run_command(cmd,False)
			# ~ print "Created Folder"
		except:
			print "Unable to create directory ",name

#Run a bash command
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



########################################################################
# Load information modules
########################################################################

#Load clusters into memory
def load_clusters(fileName):
	refClusters = {}
	species = set([])
	for line in open(fileName):
		line = line.strip()
		dades = line.split("\t")
		spe = dades[0].split("_")[1]
		species.add(spe)
		refClusters[dades[0]] = dades[1].split(";")
	return refClusters,species

#Load predicted clusters into memory
def load_cluster_families(fileName):
	clusters = {}
	num = 1
	for line in open(fileName):
		line = line.strip()
		dades = line.split()
		if len(dades) > 1:
			code = "CF_%.6d" %(num)
			num += 1
			clusters[code] = {}
			for name in dades:
				spe = name.split("_")[1]
				if spe not in clusters[code]:
					clusters[code][spe] = []
				clusters[code][spe].append(name)
	return clusters

#Load conversion files into memory	
def load_conversion(dirName,species):
	conversion = {}
	for code in species:
		fileName = dirName+"/"+code+".txt"
		#~ code = fileName.split("/")[-1].split(".")[0]
		conversion[code] = {}
		for line in open(fileName):
			line = line.strip()
			dades = line.split()
			contig = dades[0].split("_")[1]
			if contig not in conversion[code]:
				conversion[code][contig] = {}
			conversion[code][contig][dades[0]] = dades[1]
	return conversion
	
#Load homologous pairs of proteins into memory
def load_pairs(spe1,spe2,outDir):
	pairs = []
	for line in open(outDir+"/"+spe1+"/"+spe2+".txt"):
		line = line.strip()
		dades = line.split()
		p = dades[2]+"-"+dades[3]
		pairs.append(p)
	return pairs

#~ #Load homologous pairs of proteins into memory
#~ def load_pairs_by_protein(spe1,spe2,outDir):
	#~ pairs = {}
	#~ for line in open(outDir+"/"+spe1+"/"+spe2+".txt"):
		#~ line = line.strip()
		#~ dades = line.split()
		#~ p1 = dades[2]
		#~ p2 = dades[3]
		#~ if p1 not in pairs:
			#~ pairs[p1] = []
		#~ pairs[p1].append(p2)
		#~ if p2 not in pairs:
			#~ pairs[p2] = []
		#~ pairs[p2].append(p1)
	#~ return pairs

########################################################################
#Build files
########################################################################

#Create the conversion files between the protein names and the number of the family they belong in
def build_conversion_files(fileName,outDir):
	info = {}
	num = 1
	for line in open(fileName):
		line = line.strip()
		dades = line.split()
		for d in dades:
			spe = d.split("_")[0]
			if d.count("_") != 2:
				pass
			if spe not in info:
				info[spe] = {}
			info[spe][d] = str(num)
		num += 1
	for code in info:
		outfile = open(outDir+"/"+code+".txt","w")
		genes = info[code].keys()
		genes = sorted(genes,key=lambda x: int(x.split("_")[2]))
		for gene in genes:
			print >>outfile,gene+"\t"+info[code][gene]
		outfile.close()

#Creates a list of all the proteins
def build_complete_protein_list(fileName,outfileName):
	outfile = open(outfileName,"w")
	for line in open(fileName):
		line = line.strip()
		if ">" in line:
			print >>outfile,line.replace(">","")
	outfile.close()

#Creates a full list of pairwise homologous genes according to the mcl
def build_pairs(fileName,outDir):
	species = {}
	info = {}
	num = 1
	for line in open(fileName):
		line = line.strip()
		codes = line.split("\t")                      
		if len(codes) > 1:         
			info[num] = codes
			for c in codes:      
				s = c.split("_")[0]
				if s not in species:
					species[s] = {}   
				species[s][c] = num    
			num += 1
	for spe in species:                     
		outfolder = outDir+"/"+spe
		create_folder(outfolder)
		pairs = {}           
		for code in species[spe]:
			num = species[spe][code]
			for name in info[num]:  
				s = name.split("_")[0]
				if s not in pairs:     
					pairs[s] = set([])
				p = code+"-"+name
				pairs[s].add(p)
		for s in pairs:
			outfile = open(outfolder+"/"+s+".txt","w")
			for p in pairs[s]:
				print >>outfile,spe,s,p.split("-")[0],p.split("-")[1]
			outfile.close()


#Creates jobs list that will automatize the process
def create_jobs(outDir,tagName,thrMode):
	#Create folder
	outJob = outDir+"/jobs/"
	create_folder(outJob)
	if tagName == "initial":
		#Obtain a list of all the species
		species = []
		for fileName in glob.glob(outDir+"/pairs_files/*"):
			code = fileName.split("/")[-1]
			species.append(code)
		#For each pair of species
		outfile = open(outJob+"jobs.step1.txt","w")
		for i in range(0,len(species)):
			spe1 = species[i]
			for j in range(0,len(species)):
				if i != j:
					spe2 = species[j]
					if i > j:
						print >>outfile,"python "+args.pathEvolClust+" -s1 "+spe1+" -s2 "+spe2+" --get_pairwise_clusters -d "+args.outDir+" --threshold "+thrMode
		outfile.close()
		outfileName = outJob+"jobs.step1.txt"
	elif tagName == "comparison":
		outfileName = outJob+"jobs.step2.txt"
		outfile = open(outJob+"jobs.step2.txt","w")
		for fileName in glob.glob(outDir+"/clusters_by_spe/*"):
			print >>outfile,"python "+args.pathEvolClust+" -i "+fileName+" -d "+outDir+" --cluster_comparison"
		outfile.close()
	return outfileName

########################################################################
# Check file presence
########################################################################

def check_files(folderName):
	jobs = {}
	for line in open(folderName+"/jobs/jobs.step1.txt"):
		line = line.strip()
		s1 = line.split("-s1 ")[1].split(" ")[0]
		s2 = line.split("-s2 ")[1].split(" ")[0]
		if s1 not in jobs:
			jobs[s1] = {}
		jobs[s1][s2] = line
	species = [x.split("/")[-1] for x in glob.glob(folderName+"/pairs_files/*")]
	ok = True
	outfile = open(folderName+"/jobs/jobs.step1.missing.txt","w")
	for a in range(0,len(species)):
		fileName = folderName+"/clusters_from_pairs/"+species[a]+"/"
		for b in range(0,len(species)):
			if a != b:
				fileName2 = fileName+species[b]+".txt"
				if not os.path.exists(fileName2):
					if species[a] in jobs:
						print >>outfile,jobs[species[a]][species[b]]
					elif species[b] in jobs:
						print >>outfile,jobs[species[b]][species[a]]
					else:
						print "File "+fileName+" not found but jobs command also not found"
					ok = False
	if not ok:
		exit("Some files are missing from the previous step, please check out the file "+folderName+"/jobs/jobs.step1.missing.txt")
	

########################################################################
# Genome walking - Main script
########################################################################

#Obtain conserved cluster given a pair of homologous proteins
def get_clusters(p1,p2,conversion,set_difference):
	spe1 = p1.split("_")[0]
	spe2 = p2.split("_")[0]
	if spe1 != spe2:
		cluster1,cluster2,score = calculate_cluster(p1,p2,conversion,set_difference)
	return cluster1,cluster2,score

def calculate_cluster(s1,s2,conversion,set_difference):
	#Obtain contigs where the seed is located
	spe1,contig1 = s1.split("_")[0],s1.split("_")[1]
	#Get all the proteins in the contig
	cluster1 = conversion[spe1][contig1].keys()
	#Sort them
	cluster1 = sorted(cluster1,key=lambda x: int(x.split("_")[2]))
	#Convert them into their protein family codes
	cluster1Trans = [conversion[spe1][contig1][x] for x in cluster1 if x in conversion[spe1][contig1]]
	#Repeat for second contig
	spe2,contig2 = s2.split("_")[0],s2.split("_")[1]
	cluster2 = conversion[spe2][contig2].keys()
	cluster2 = sorted(cluster2,key=lambda x: int(x.split("_")[2]))
	cluster2Trans = [conversion[spe2][contig2][x] for x in cluster2 if x in conversion[spe2][contig2]]
	#Trim edges according to presence of homologs in both clusters
	cluster1 = trim_clusters(cluster1,cluster1Trans,cluster2Trans)
	cluster1Trans = [conversion[spe1][contig1][x] for x in cluster1 if x in conversion[spe1][contig1]]
	cluster2 = trim_clusters(cluster2,cluster2Trans,cluster1Trans)
	cluster2Trans = [conversion[spe2][contig2][x] for x in cluster2 if x in conversion[spe2][contig2]]
	#Save current cluster length
	size1 = len(cluster1)
	size2 = len(cluster2)
	continua = True
	while continua:
		#Delete proteins that are not present in the opposite cluster
		cluster1 = delete_non_homologs(cluster1,cluster1Trans,cluster2Trans)
		#Cut the cluster if there are stretches of more than three proteins in between, will return the piece where the seed is present as long as it's still larger or equal than the minSize
		cluster1 = cut_cluster(cluster1,s1,set_difference)
		cluster1Trans = [conversion[spe1][contig1][x] for x in cluster1 if x in conversion[spe1][contig1]]
		#If the cluster still exists then the same is done for the opposite cluster
		if len(cluster1) != 0:
			cluster2 = delete_non_homologs(cluster2,cluster2Trans,cluster1Trans)
			cluster2 = cut_cluster(cluster2,s2,set_difference)
			cluster2Trans = [conversion[spe2][contig2][x] for x in cluster2 if x in conversion[spe2][contig2]]
		else:
			#If no piece remains then the clusters will be returned empty
			cluster2 = []
			cluster2Trans = []
		#We check whether we still have two clusters
		if len(cluster1) == 0 or len(cluster2) == 0:
			continua = False
		#We check if the size of the clusters has changed during the last iteration, if not this will finish
		elif len(cluster1) == size1 and len(cluster2) == size2:
			continua = False
		#We continue
		if continua:
			size1 = len(cluster1)
			size2 = len(cluster2)
	#If, after the whole process we still have two clusters left then we complete the possible gaps and calculate their score
	if len(cluster1) != 0 and len(cluster2) != 0:
		cluster1 = complete_cluster(cluster1)
		cluster2 = complete_cluster(cluster2)
		#print "Complete cluster found, calculating score"
		score = calculate_score(cluster1,cluster2,conversion)
	else:
		cluster1 = []
		cluster2 = []
		score = 0
	return cluster1,cluster2,score

#Fill in non-homologous proteins
def complete_cluster(cluster):
	start = int(cluster[0].split("_")[2])
	stop = int(cluster[-1].split("_")[2])
	base = "_".join(cluster[0].split("_")[:2])
	new_cluster = []
	for i in range(start,stop+1):
		n = "%s_%.5d" % (base,i)
		new_cluster.append(n)
	return new_cluster

#Delete proteins that don't have homologs in the opposite cluster	
def delete_non_homologs(cluster1,cluster1Trans,cluster2Trans):
	new_cluster = []
	for i in range(0,len(cluster1Trans)):
		c = cluster1Trans[i]
		if c in cluster2Trans:
			new_cluster.append(cluster1[i])
	return new_cluster

#Cut the cluster in pieces by parts where more than three non-homologous proteins are found
def cut_cluster(cluster,seed,set_difference):
	new_group = [cluster[0]]
	group_found = []
	for code in cluster[1:]:
		p1 = int(new_group[-1].split("_")[2])
		p2 = int(code.split("_")[2])
		diff = p2 - p1
		if diff <= set_difference:
			new_group.append(code)
		else:
			if seed in new_group and len(new_group) >= minSize:
				group_found = new_group
			new_group = [code]
	if seed in new_group and len(new_group) >= minSize:
		group_found = new_group
	return group_found

#Delete non-homologous edges
def trim_clusters(cluster1,cluster1Trans,cluster2Trans):
	trimmed = []
	cl1 = set(cluster1Trans)
	cl2 = set(cluster2Trans)
	common = cl1.intersection(cl2)
	if len(common) != 0:
		nums = []
		for i in range(0,len(cluster1Trans)):
			c = cluster1Trans[i]
			if c in common:
				nums.append(i)
		minim = min(nums)
		maxim = max(nums)+1
		trimmed = cluster1[minim:maxim]
	else:
		trimmed = []
	return trimmed	

########################################################################	
# Calculate thresholds and clusters
########################################################################

#Calculates all clusters and from there it obtains the thresholds and the list of clusters that pass the filters.
def get_clusters_and_thresholds(pairs,minSize,maxSize,conversion,set_difference,thr_mode):
	allClusters = {}
	total = len(pairs)
	num = 0
	thresholds1 = {}
	thresholds2 = {}
	for a in range(minSize,maxSize+1):
		thresholds1[a] = {}
		thresholds2[a] = {}
	allClusters = {}
	for p in pairs:
		if num % 1000 == 0:
			print "Processed: ",num,"out of",total
		num += 1
		p1,p2 = p.split("-")
		if args.species1 in p1:
			pass
		else:
			p1,p2 = p2,p1
		cl1,cl2,score = get_clusters(p1,p2,conversion,set_difference)
		allClusters[p] = [cl1,cl2,score]
		size1,size2 = len(cl1),len(cl2)
		if size1 >= minSize and size2 >= minSize:
			thresholds1 = get_threshold_scores(cl1,cl2,conversion,thresholds1,p1)
			thresholds2 = get_threshold_scores(cl2,cl1,conversion,thresholds2,p2)
	thresholds1 = print_thresholds(args.species1,args.species2,thresholdsPath,thresholds1,all_proteins1,thr_mode)
	thresholds2 = print_thresholds(args.species2,args.species1,thresholdsPath,thresholds2,all_proteins2,thr_mode)
	#Filter clusters based on thresholds
	path = args.outDir+"/clusters_from_pairs/"
	pathAll = args.outDir+"/all_cluster_predictions/"
	create_folder(path)
	create_folder(pathAll)
	clusters1 = set([])
	clusters2 = set([])
	pathAll1 = pathAll+"/"+args.species1+"/"
	create_folder(pathAll1)
	pathAll2 = pathAll+"/"+args.species2+"/"
	create_folder(pathAll2)
	outfileAll1 = open(pathAll1+"/"+args.species2+".txt","w")
	print >>outfileAll1,"#Code1\tCode2\tClusterSize1\tClusterSize2\tThreshold1\tThreshold2\tScore\tCluster1\tCluster2"
	outfileAll2 = open(pathAll2+"/"+args.species1+".txt","w")
	print >>outfileAll2,"#Code1\tCode2\tClusterSize1\tClusterSize2\tThreshold1\tThreshold2\tScore\tCluster1\tCluster2"
	for p in allClusters:
		cl1 = allClusters[p][0]
		cl2 = allClusters[p][1]
		score = allClusters[p][2]
		size1 = len(cl1)
		size2 = len(cl2)
		p1,p2 = p.split("-")
		if len(cl1) != 0:
			print >>outfileAll1,p1+"\t"+p2+"\t"+str(size1)+"\t"+str(size2)+"\t"+str(score)+"\t"+";".join(cl1)+"\t"+";".join(cl2)
		else:
			print >>outfileAll1,p1+"\t"+p2+"\tNone"
		if len(cl2) != 0:
			print >>outfileAll2,p2+"\t"+p1+"\t"+str(size2)+"\t"+str(size1)+"\t"+str(score)+"\t"+";".join(cl2)+"\t"+";".join(cl1)
		else:
			print >>outfileAll2,p2+"\t"+p1+"\tNone"
		if size1 >= minSize and size2 >= minSize:
			if size1 > maxSize or size2 > maxSize:
				pass
			else:
				if score > thresholds1[size1] and score > thresholds2[size2]:
					cl1 = ";".join(cl1)
					cl2 = ";".join(cl2)
					clusters1.add(cl1)
					clusters2.add(cl2)
	outfileAll1.close()
	outfileAll2.close()
	path1 = path+"/"+args.species1
	create_folder(path1)
	outfile1 = open(path1+"/"+args.species2+".txt","w")
	for cl1 in clusters1:
		print >>outfile1,cl1
	outfile1.close()
	path2 = path+"/"+args.species2
	create_folder(path2)
	outfile2 = open(path2+"/"+args.species1+".txt","w")
	for cl2 in clusters2:
		print >>outfile2,cl2
	outfile1.close()

#Given two found clusters it obtains the values of the thresholds we're going to use
#It will allways save the best score possible for a given protein, no matter how many homologs
#it has.
def get_threshold_scores(cl1,cl2,conversion,thresholds,prot):
	size = len(cl1)
	if size > maxSize:
		size = maxSize
	for s in range(minSize,size+1):
		if prot not in thresholds[s]:
			vesFent = True
		else:
			if thresholds[s][prot] == 1.0:
				vesFent = False
			else:
				vesFent = True
		if vesFent:
			spe1 = cl1[0].split("_")[0]
			contig1 = cl1[0].split("_")[1]
			spe2 = cl2[0].split("_")[0]
			contig2 = cl2[0].split("_")[1]
			cl1T = [conversion[spe1][contig1][x] for x in cl1 if x in conversion[spe1][contig1]]
			cl2T = [conversion[spe2][contig2][x] for x in cl2 if x in conversion[spe2][contig2]]
			current_score = 0.0
			for a in range(0,len(cl1)-s+1):
				if current_score != 1.0:
					if prot in cl1[a:a+s]:
						cluster2 = trim_clusters(cl2,cl2T,cl1T[a:a+s])
						if len(cluster2) != 0:
							score = calculate_score(cl1[a:a+s],cluster2,conversion)
						else:
							score = 0.0
						if score > current_score:
							current_score = score
			if prot in thresholds[s]:
				if current_score > thresholds[s][prot]:
					thresholds[s][prot] = current_score
			else:
				thresholds[s][prot] = current_score
	return thresholds

#Summarize and print thresholds in a file for record keeping
def print_thresholds(spe1,spe2,outDir,thresholds,all_proteins,thr_mode):
	outDir = outDir+"/"+spe1
	create_folder(outDir)
	outfile = open(outDir+"/"+spe2+".txt","w")
	print >>outfile,"Cluster_size\tThreshold"
	sizes = thresholds.keys()
	sizes.sort()
	thresholds2 = {}
	for size in sizes:
		values = []
		for code in all_proteins:
			if code not in thresholds[size]:
				values.append(0.0)
			else:
				values.append(thresholds[size][code])
		if thr_mode == "2stdv":
			thr = calculate_threshold1(values)
		elif thr_mode == "3stdv":
			thr = calculate_threshold2(values)
		elif thr_mode == "1stdv":
			thr = calculate_threshold3(values)
		elif thr_mode == "90percent":
			thr = calculate_threshold4(values)
		elif thr_mode == "75percent":
			thr = calculate_threshold5(values)
		elif thr_mode == "non_parametric":
			thr = calculate_threshold6(values)
		print >>outfile,"%d\t%.3f" %(size,thr)
		thresholds2[size] = thr
	outfile.close()
	return thresholds2

#Threshold 1: average + 2stdv
def calculate_threshold1(values):
	average = numpy.average(values)
	std = numpy.std(values)
	thr = average + 2.0*std
	if thr > 1.0:
		thr = 1.0
	return thr

#Threshold 2: average + 3stdv
def calculate_threshold2(values):
	average = numpy.average(values)
	std = numpy.std(values)
	thr = average + 3.0*std
	if thr > 1.0:
		thr = 1.0
	return thr

#Threshold 3: average + 1stdv
def calculate_threshold3(values):
	average = numpy.average(values)
	std = numpy.std(values)
	thr = average + 1*std
	if thr > 1.0:
		thr = 1.0
	return thr

#Threshold 4: takes value below which 90% of the data can be found
def calculate_threshold4(values):
	values.sort()
	num = int(float(len(values)) * 0.9)
	print len(values),num
	thr = values[num]
	return thr

#Threshold 5: takes value below which 75% of the data can be found
def calculate_threshold5(values):
	values.sort()
	num = int(float(len(values)) * 0.75)
	thr = values[num]
	return thr

#Threshold 6: calculates non-parametric tolerance intervals using the wilks method
def calculate_threshold6(values):
	vals = robjects.FloatVector(values)
	results = tolerance.nptol_int(vals,alpha=0.05,P=0.99,side=1,method="WILKS")
	return results[3][0]
	



########################################################################	
# Calculate scores
########################################################################

def calculate_score(cluster1,cluster2,conversion):
	cl1 = []
	cl2 = []
	spe1 = cluster1[0].split("_")[0]
	contig1 = cluster1[0].split("_")[1]
	spe2 = cluster2[0].split("_")[0]
	contig2 = cluster2[0].split("_")[1]
	missing = 0
	if len(conversion) == 0:
		cl1 = list(cluster1)
		cl2 = list(cluster2)
		print "MISSING conversion"
	else:
		#print "Conversion ok"
		for cl in cluster1:
			if cl in conversion[spe1][contig1]:
				cl1.append(conversion[spe1][contig1][cl])
			else:
				name = "m_"+str(missing)
				missing += 1
				cl1.append(name)
		for cl in cluster2:
			if cl in conversion[spe2][contig2]:
				cl2.append(conversion[spe2][contig2][cl])
			else:
				name = "m_"+str(missing)
				missing += 1
				cl2.append(name)
	common1 = len([x for x in cl1 if x in cl2])
	common2 = len([x for x in cl2 if x in cl1])
	#print "Number of common proteins:",common1,common2
	if common1 < 2 or common2 < 2:
		#At least the cluster has to have two homologous proteins for it to be considered
		score = 0.0
	else:
		size = len(cluster1)*len(cluster2)*2
		missing1 = len(cluster1)-common1
		missing2 = len(cluster2)-common2
		missing = missing1*missing1+missing2*missing2
		common = common1*common1+common2*common2
		score = float(common - missing)/float(size)
		if score < 0.0:
			score = 0.0
		elif score > 1.0:
			score = 1.0
	return score

#Calculates score between pairs of clusters
def calculate_partial_scores(clustersSpe,clustersAll,conversion,outfileName):
	outfile = open(outfileName,"w")
	codesSpe = clustersSpe.keys()
	codesAll = clustersAll.keys()
	for i in range(0,len(codesSpe)):
		cluster1 = clustersSpe[codesSpe[i]]
		print i,"out of",len(codesSpe)
		for j in range(0,len(codesAll)):
			cluster2 = clustersAll[codesAll[j]]
			score = calculate_score(cluster1,cluster2,conversion)
			if score > 0.0:
				print >>outfile,codesSpe[i]+"\t"+codesAll[j]+"\t"+str(score)
	outfile.close()

########################################################################	
# Filter families
########################################################################

def divide_overlaps(families,clusters):
	dividedFamilies = {}
	for spe in families:
		if len(families[spe]) == 1:
			dividedFamilies[spe] = families[spe]
		else:
			fam = families[spe]
			fam = sorted(fam,key=lambda x: int(clusters[x][0].split("_")[2]))
			all_groups = []
			new_group = [fam[0]]
			current_proteins = set(clusters[fam[0]])
			for f in fam[1:]:
				prots = set(clusters[f])
				common = current_proteins.intersection(prots)
				if len(common) == 0:
					all_groups.append(new_group)
					new_group = [f]
					current_proteins = prots
				else:
					new_group.append(f)
					current_proteins = current_proteins.union(prots)
			all_groups.append(new_group)
			if len(all_groups) == 1:
				dividedFamilies[spe] = new_group
			else:
				num = 1
				for n in all_groups:
					name = "%s_%.3d" % (spe,num)
					num += 1
					dividedFamilies[name] = n
	return dividedFamilies

def delete_outliers(family,clusters):
	info = {}
	for spe in family:
		for code in family[spe]:
			info[code] = len(clusters[code])
			#~ print code,len(clusters[code])
	delete = []
	#~ print "INITIAL:"
	#~ for code in info:
		#~ print code,info[code]
	average,std,thrUp,thrDown = get_stats(info.values())
	size = 0
	saved = []
	#~ print info.keys()
	saved,info = delete_unique(info,saved)
	num = 1
	if len(info) != 0:
		while len(info) != size:
			size = len(info)
			#~ print num,delete,info.keys()
			info,delete = delete_outlier(info,average,thrUp,thrDown,std,delete)
			if len(info) != size:
				average,std,thrUp,thrDown = get_stats(info.values())
				saved,info = delete_unique(info,saved)
		#~ print delete
	#~ print "TO DELETE:",delete
	family2 = {}
	for spe in family:
		family2[spe] = []
		for code in family[spe]:
			if code not in delete:
				family2[spe].append(code)
	return family2

def get_stats(nums):
	average = numpy.average(nums)
	std = numpy.std(nums)
	thrUp = average+(std*2.0)
	thrDown = average-(std*2.0)
	return average,std,thrUp,thrDown

def delete_unique(info,saved):
	counter = {}
	info2 = {}
	for code in info:
		spe = code.split("_")[1]
		if spe not in counter:
			counter[spe] = []
		counter[spe].append(code)
	for code in counter:
		if len(counter[code]) == 1:
			saved.append(counter[code][0])
		else:
			for name in counter[code]:
				info2[name] = info[name]
	return saved,info2

def delete_outlier(info,average,thrUp,thrDown,std,delete):
	diffUp = 0
	chosenUp = None
	diffDown = 0
	chosenDown = None
	for code in info:
		val = info[code]
		if val > average:
			if val > thrUp:
				d = val - thrUp
				if d > diffUp:
					diffUp = d
					chosenUp = code
		else:
			if val < thrDown:
				d = thrDown - val
				if d > diffDown:
					diffDown = d
					chosenDown = code
	#This is added in case that two codes of the same species are chosen and there were only two codes for this species to begin with
	if chosenUp != None and chosenDown != None:
		sUp = chosenUp.split("_")[1]
		sDown = chosenDown.split("_")[1]
		if sUp == sDown:
			if diffUp < diffDown:
				chosenUp = None
			else:
				chosenDown = None
	if chosenUp != None:
		#~ print "CHOSEN OUTLIER UP:",chosenUp,diffUp,info[chosenUp],average,thrUp,std
		delete.append(chosenUp)
		del info[chosenUp]
	if chosenDown != None:
		#~ print "CHOSEN OUTLIER DOWN:",chosenDown,diffDown,info[chosenDown],average,thrDown,std
		delete.append(chosenDown)
		del info[chosenDown]
	return info,delete

#Checks that the cluster is not only formed by duplications of the same gene
def delete_multiple_duplications(spe,cl,conversion):
	counter = set([])
	for code in cl:
		contig = code.split("_")[1]
		if code in conversion[contig]:
			g = conversion[contig][code]
			counter.add(g)
	if len(counter) >= (minSize-1):
		PASS = True
	else:
		PASS = False
	#print spe,PASS,counter
	return PASS

########################################################################
# Fuse clusters from the same species
########################################################################

def run_mcl(clusters,clDir,tmpDir,outfileMCL):
	clustList = clusters.keys()
	outfile = open(tmpDir+"/mcl_list.txt","w")
	for a in range(0,len(clustList)):
		name1 = clustList[a]
		codes1 = set(clusters[name1])
		for b in range(0,len(clustList)):
			name2 = clustList[b]
			codes2 = set(clusters[name2])
			common = codes1.intersection(codes2)
			overlap = float(len(common)) / (float(len(codes1)+len(codes2))/2.0)
			if overlap >= 0.33:
				print >>outfile,name1+"\t"+name2+"\t"+str(overlap)
			else:
				pass
	outfile.close()
	cmd = "mcl "+tmpDir+"/mcl_list.txt --abc"
	run_command(cmd,False)
	cmd = "mv out.mcl_list.txt.I20 "+outfileMCL
	run_command(cmd,False)

def list_all_genes(clusters):
	prots = set([])
	for cl in clusters:
		for p in clusters[cl]:
			prots.add(p)
	return prots

def split_in_groups(genes,set_difference):
	groups = []
	new_group = [genes[0]]
	for gene in genes[1:]:
		contig1 = new_group[-1].split("_")[1]
		contig2 = gene.split("_")[1]
		if contig1 != contig2:
			groups.append(new_group)
			new_group = [gene]
		else:
			pos1 = int(new_group[-1].split("_")[2])
			pos2 = int(gene.split("_")[2])
			diff = pos2 - pos1
			if diff <= set_difference:
				new_group.append(gene)
			else:
				groups.append(new_group)
				new_group = [gene]
	groups.append(new_group)
	groups = fill_gaps(groups)
	return groups

def fill_gaps(groups):
	groups2 = []
	for group in groups:
		start = int(group[0].split("_")[2])
		stop = int(group[-1].split("_")[2]) + 1
		base = "_".join(group[0].split("_")[:2])
		new_group = []
		for a in range(start,stop):
			name = "%s_%.5d" % (base,a)
			new_group.append(name)
		groups2.append(new_group)
	return groups2

def get_most_representative_clusters(groups,clusters,spe,outfile,set_difference):
	numT = 1
	taken_genes = set([])
	j = 0
	for group in groups:
		if len(group) == 1:
			g = group[0]
			if len(clusters[g]) >= minSize and len(clusters[g]) <= maxSize:
				name = "FCL_"+spe+"_"+str(numT)
				numT += 1
				print >>outfile,name+"\t"+";".join(clusters[group[0]])
				for a in clusters[g]:
					taken_genes.add(a)
		else:
			#Gather all proteins and count how many times each protein appears in the different clusters
			proteins = set([])
			counter = {}
			for g in group:
				for prot in clusters[g]:
					proteins.add(prot)
					if prot not in counter:
						counter[prot] = set([])
					counter[prot].add(g)
			#Calculate the presence threshold
			presence = []
			for c in counter:
				s = len(counter[c])
				presence.append(s)
			threshold = int(numpy.average(presence))
			if threshold < 1:
				threshold = 1
			#Gather those proteins that have a presence above the threshold
			proteins2 = []
			for p in proteins:
				a = len(counter[p])
				if a >= threshold:
					proteins2.append(p)
			proteins2 = list(proteins2)
			proteins2 = sorted(proteins2,key=lambda x:int(x.split("_")[2]))
			groups2 = split_in_groups(proteins2,set_difference)
			for g in groups2:
				if len(g) >= minSize and len(g) <= maxSize:
					name = "FCL_"+spe+"_"+str(numT)
					numT += 1
					print >>outfile,name+"\t"+";".join(g)
					for a in g:
						taken_genes.add(a)
	return taken_genes


def fishing(lost,clusters,outfileREPR,clDir,tmpDir,set_difference):
	a = 0
	outfile = open(clDir+"/"+spe+".fishing","w")
	clusters2 = {}
	for cl in clusters:
		genes = set(clusters[cl])
		common = genes.intersection(lost)
		diff = len(genes) - len(common)
		if diff == 0 or diff == 1 or diff == 2:
			print >>outfile,cl+"\t"+";".join(clusters[cl])
			clusters2[cl] = clusters[cl]
	outfile.close()
	outfileMCL = clDir+"/"+spe+".out.2.mcl"
	run_mcl(clusters2,clDir,tmpDir,outfileMCL)
	groups = []
	for line in open(outfileMCL):
		line = line.strip()
		dades = line.split()
		groups.append(dades)
	taken_genes = get_most_representative_clusters(groups,clusters2,spe,outfileREPR,set_difference)
	return taken_genes

########################################################################
# Add clusters that were discarded
########################################################################

#Gather information of all the precomputed clusters
def compile_info(inDirCL,inDirThr):
	#Create a single threshold file
	species = set([])
	outfileTHR = open(inDirThr+"/all_thresholds.txt","w")
	species = [x.split("/")[-1] for x in glob.glob(inDirCL+"/*") if ".txt" not in x]
	counter = 0
	for spe1 in species:
		#print counter,len(species)
		counter += 1
		dirName = inDirCL+"/"+spe1+"/"
		#Join all clusters predicted for this species into a single file
		outfileCL = open(dirName+"/all_clusters.txt","w")
		for spe2 in species:
			if spe1 != spe2:
				fileName = inDirCL+"/"+spe1+"/"+spe2+".txt"
				for line in open(fileName):
					line = line.strip()
					if "None" not in line and "#" not in line:
						print >>outfileCL,line
		outfileCL.close()
		#Put all thresholds in a file
		for spe2 in species:
			if spe1 != spe2:
				fileName = inDirThr+"/"+spe1+"/"+spe2+".txt"
				for line in open(fileName):
					line = line.strip()
					dades = line.split()
					if "Cluster_size" not in line:
						print >>outfileTHR,spe1+"\t"+spe2+"\t"+dades[0]+"\t"+dades[-1]
	outfileTHR.close()
	return species

#Load predicted clusters into memory

def load_cluster_families2(fileName):
	clFM = {}
	protsCL = set([])
	for line in open(fileName):
		line = line.strip()
		dades = line.split("\t")
		if "##" in line:
			pass
		elif "# " in line:
			name = line.split()[1]
			clFM[name] = {}
		else:
			dades = line.split("\t")
			clFM[name][dades[0]] = dades[1].split(";")
			codes = dades[1].split(";")
			for c in codes:
				protsCL.add(c)
	return clFM,protsCL
	
#Load thresholds into memory

def load_thresholds(fileName):
	thresholds = {}
	for line in open(fileName):
		line = line.strip()
		dades = line.split("\t")
		c1 = dades[0]
		c2 = dades[1]
		if c1 not in thresholds:
			thresholds[c1] = {}
		if c2 not in thresholds[c1]:
			thresholds[c1][c2] = {}
		thresholds[c1][c2][int(dades[2])] = float(dades[-1])
	return thresholds

#Load clusters into memory discarding those clusters that are already into memory
def get_all_clusters(inDir,species,protsCL):
	clusters = {}
	a = 1
	for spe1 in species:
		#print a,spe1,len(species),len(clusters)
		a += 1
		for line in open(inDir+"/"+spe1+"/all_clusters.txt"):
			line = line.strip()
			dades = line.split("\t")
			spe2 = dades[1].split("_")[0]
			if dades[1] not in protsCL:
				codes = dades[-1].split(";")
				common = set(codes).intersection(protsCL)
				overlap = float(len(common)) / float(len(codes))
				if overlap < 0.1:
					if dades[0] not in clusters:
						clusters[dades[0]] = {}
					if spe2 not in clusters[dades[0]]:
						clusters[dades[0]][spe2] = set([])
					n = [codes[0].split("_")[0],codes[0].split("_")[1],codes[0].split("_")[2],codes[-1].split("_")[2]]
					n = "\t".join(n)
					clusters[dades[0]][spe2].add(n)
	return clusters

#Select candidates to add to the cluster family
def get_candidates(cf,all_clusters,thresholds,maxSize,set_difference):
	species_included = set([x.split("_")[0] for x in cf])
	candidates = {}
	for spe in cf:
		cl1 = set(cf[spe])
		#For each protein included in the cluster 
		for prot in cf[spe]:
			spe1 = spe.split("_")[0]
			if prot in all_clusters:
				#For each species that has an homologous cluster predicted from this protein
				for spe2 in all_clusters[prot]:
					#If the species is not already included in the family
					if spe2 not in species_included:
						for info in all_clusters[prot][spe2]:
							#For each alternative cluster obtain the information and rebuild the cluster
							info = info.split("\t")
							cl2 = set([])
							for a in range(int(info[2]),int(info[3])+1):
								b = "%s_%s_%.5d" %(info[0],info[1],a)
								cl2.add(b)
							#See whether the two clusters are homologs and obtain the score
							cluster2,score = compare_clusters(cl1,cl2,conversion,set_difference)
							if score:
								size = len(cluster2)
								if score >= 0.95 and size <= maxSize:
									thr = thresholds[spe2][spe1][size]
									if score >= thr:
										if spe1 not in candidates:
											candidates[spe1] = {}
										if spe2 not in candidates[spe1]:
											candidates[spe1][spe2] = []
										if len(candidates[spe1][spe2]) == 0:
											candidates[spe1][spe2] = [cluster2,score]
										else:
											if score > candidates[spe1][spe2][1]:
												candidates[spe1][spe2] = [cluster2,score]
	return candidates


def compare_clusters(cl1,cl2,conversion,set_difference):
	spe1,contig1 = list(cl1)[0].split("_")[0],list(cl1)[0].split("_")[1]
	#Sort them
	cluster1 = sorted(list(cl1),key=lambda x: int(x.split("_")[2]))
	#Convert them into their protein family codes
	cluster1Trans = [conversion[spe1][contig1][x] for x in cluster1 if x in conversion[spe1][contig1]]
	#Repeat for second contig
	spe2,contig2 = list(cl2)[0].split("_")[0],list(cl2)[0].split("_")[1]
	#Sort them
	cluster2 = sorted(list(cl2),key=lambda x: int(x.split("_")[2]))
	#Convert them into their protein family codes
	cluster2Trans = [conversion[spe2][contig2][x] for x in cluster2 if x in conversion[spe2][contig2]]	
	#Trim edges
	cluster2 = trim_clusters(cluster2,cluster2Trans,cluster1Trans)
	cluster2Trans = [conversion[spe2][contig2][x] for x in cluster2 if x in conversion[spe2][contig2]]
	if len(cluster2) !=0:
		groups = cut_cluster_without_seed(cluster2,cluster1Trans,set_difference,conversion[spe2][contig2])
		min_common = 0
		chosen = None
		if len(groups) != 0:
			for group in groups:
				cl2Trans = set([conversion[spe2][contig2][x] for x in group if x in conversion[spe2][contig2]])
				cl1Trans = set(cluster1Trans)
				common = cl1Trans.intersection(cl2Trans)
				if len(common) >= 3:
					if len(common) > min_common:
						chosen = group
						min_common = len(common)
		if chosen != None:
			cluster2 = chosen
			score = calculate_score(cluster1,cluster2,conversion)
		else:
			score = None
	return cluster2,score

#Cut the cluster in pieces by parts where more than three non-homologous proteins are found
def cut_cluster_without_seed(cluster,cluster1Trans,set_difference,conversion):
	minSize = 5
	new_cluster = []
	for cl in cluster:
		if cl in conversion:
			t = conversion[cl]
			if t in cluster1Trans:
				new_cluster.append(cl)
	if len(new_cluster) != 0:
		new_group = [new_cluster[0]]
		groups = []
		for code in new_cluster[1:]:
			p1 = int(new_group[-1].split("_")[2])
			p2 = int(code.split("_")[2])
			diff = p2 - p1
			if diff <= set_difference:
				new_group.append(code)
			else:
				if len(new_group) >= minSize:
					groups.append(new_group)
				new_group = [code]
		if len(new_group) >= minSize:
			groups.append(new_group)
	else:
		groups = []
	return groups

def fuse_clusters(group,protsCL):
	all_codes = set([])
	for g in group:
		codes = g.split(";")
		for a in codes:
			all_codes.add(a)
	common = all_codes.intersection(protsCL)
	overlap = float(len(common)) / float(len(all_codes))
	if overlap < 0.1:
		all_codes = list(all_codes)
		all_codes = sorted(all_codes,key=lambda x: int(x.split("_")[2]))
		cl = ";".join(all_codes)
	else:
		cl = None
	return cl

parser = argparse.ArgumentParser(description="Will perform the genome walking")
parser.add_argument("-i","--infile",dest="inFile",action="store",default=None,help="mcl or pairs file, can hold multiple or single families.")
parser.add_argument("-f","--fastafile",dest="fastaFile",action="store",default=None,help="Fasta file that contains the complete proteome database.")
parser.add_argument("-s1","--species1",dest="species1",action="store",default=None,help="Species tag needed to run the threshold calculation")
parser.add_argument("-s2","--species2",dest="species2",action="store",default=None,help="Species tag needed to run the threshold calculation")
parser.add_argument("-d","--outdir",dest="outDir",action="store",default="./genome_walking/",help="basepath folder where the results will be stored")
parser.add_argument("--minSize",dest="minSize",action="store",default=5,help="Minimum size a cluster needs to have to be considered. Default is set to 5")
parser.add_argument("--maxSize",dest="maxSize",action="store",default=35,help="Maximum size a cluster needs to have to be considered. Default is set to 35")
parser.add_argument("--non_homologs",dest="non_homologs",action="store",default=3,help="Number of non-homologous genes needed to split a cluster")
parser.add_argument("--threshold",dest="thr",action="store",choices=["2stdv","3stdv","1stdv","90percent","75percent","non_parametric"],default="2stdv",help="Way to calculate the thresholds to accept a conserved region as cluster")
parser.add_argument("--build_jobs",dest="buildJobs",action="store_true",help="Will prepare the files so that trees can be build. Needs to have the -i and -d options as full paths")
parser.add_argument("--initial_files",dest="conv",action="store_true",help="Will prompt the program to create the initial files")
parser.add_argument("--get_pairwise_clusters",dest="calc_scores_pairs",action="store_true",help="Main program in the genome walking - starts from pairs files")
parser.add_argument("--calculate_thresholds",dest="calc_thr",action="store_true",help="Will calculate the thresholds for each pair of species")
parser.add_argument("--cluster_comparison",dest="cluster_comparison",action="store_true",help="Will compare clusters")
parser.add_argument("--filter_clusters",dest="filter_clusters",action="store_true",help="Will filter out identical clusters and name them all; will also create the files needed to run the final cluster comparison")
parser.add_argument("--filter_cluster_families",dest="clusterFam_filter",action="store_true",help="Will delete redundancy among clusters within cluster families")
parser.add_argument("--complete_families",dest="clusterFam_complete",action="store_true",help="Will complete cluster families adding clusters for closely related species that were discarded.")
parser.add_argument("--path_evolclust",dest="pathEvolClust",action="store",default="./evolclust.py",help="Path to the python program evolclust.py")
parser.add_argument("--path_mcl",dest="pathMCL",action="store",default="mcl",help="Path to mcl")
parser.add_argument("--local",dest="local",action="store_true",help="This will run the whole evolclust pipeline without splitting it in steps. It is only recommended for small datasets")
args = parser.parse_args()

#Create output folder if it doesn't exist
create_folder(args.outDir)
minSize = int(args.minSize)
maxSize = int(args.maxSize)
set_difference = args.non_homologs
if args.local:
	print "STEP1: Create initial files"
	cmd = "python "+args.pathEvolClust+" -i "+args.inFile+" -d "+args.outDir+" -f "+args.fastaFile+" --initial_files --path_evolclust "+args.pathEvolClust+" --threshold "+args.thr
	run_command(cmd,False)
	print "STEP2: Calculate thresholds and obtain pairwise clusters (This can take a long time)"
	for line in open(args.outDir+"/jobs/jobs.step1.txt"):
		line = line.strip()
		run_command(line,False)
	print "STEP3: Filter paiwise clusters"
	cmd = "python "+args.pathEvolClust+" -d "+args.outDir+" --filter_clusters"
	run_command(cmd,False)
	print "STEP4: All to all comparison of clusters"
	cmd = "python "+args.pathEvolClust+" -d "+args.outDir+" -i "+args.outDir+"/complete_cluster_list.txt --cluster_comparison"
	run_command(cmd,False)
	print "STEP5: Put clusters into families (uses mcl)"
	cmd = args.pathMCL+" "+args.outDir+"/cluster_comparison/complete_cluster_list.txt --abc"
	run_command(cmd,False)
	cmd = "mv out.complete_cluster_list.txt.I20 "+args.outDir+"/complete_comparison.mcl"
	run_command(cmd,False)
	print "STEP6: Clean cluster families"
	cmd = "python "+args.pathEvolClust+" -d "+args.outDir+" -i "+args.outDir+"/complete_comparison.mcl --filter_cluster_families"
	run_command(cmd,False)
	print "STEP7: Complete cluster families"
	cmd = "python "+args.pathEvolClust+" -d "+args.outDir+" -i "+args.outDir+"/final_clusters.txt --complete_families"
	run_command(cmd,False)

#Creates conversion files that will correlate each protein to the corresponding protein family
#RUNS FOR ALL DATA - Checked
if args.conv:
	outDir = args.outDir+"/conversion_files/"
	create_folder(outDir)
	build_conversion_files(args.inFile,outDir)
	outfileName = args.outDir+"/complete_protein_list.txt"
	build_complete_protein_list(args.fastaFile,outfileName)
	#Creates files where all the pairs of homologous proteins between two species according to the mcl can be found
	outDir = args.outDir+"/pairs_files/"
	create_folder(outDir)
	build_pairs(args.inFile,outDir)
	jobsOutfileName = create_jobs(args.outDir,"initial",args.thr)

#This script will obtain a list of clusters that have a higher conservation score than the cluster
#RUN IN CLUSTER FOR PAIRS OF SPECIES
if args.calc_scores_pairs:
	#Load previously calculated data
	species = [args.species1,args.species2]
	pairsDir = args.outDir+"/pairs_files/"
	pairs = load_pairs(species[0],species[1],pairsDir)
	all_proteins1 = set([line.strip() for line in open(args.outDir+"/complete_protein_list.txt") if line.split("_")[0] == args.species1])
	all_proteins2 = set([line.strip() for line in open(args.outDir+"/complete_protein_list.txt") if line.split("_")[0] == args.species2])
	conversion = load_conversion(args.outDir+"/conversion_files/",species)
	#Get clusters and thresholds
	thresholdsPath = args.outDir+"/thresholds/"
	create_folder(thresholdsPath)
	if args.thr == "non_parametric":
		try:
			tolerance = importr("tolerance")
		except:
			exit("For this threshold mode you need to install the R-package tolerance: https://cran.r-project.org/web/packages/tolerance/tolerance.pdf")
	get_clusters_and_thresholds(pairs,minSize,maxSize,conversion,set_difference,args.thr)

#This will create a single list with all the clusters found. For each species it performs a mcl to remove redundancy and then creates a unique list
#NEEDS TO BE RUN AFTER ALL CLUSTERS HAVE BEEN CALCULATED
if args.filter_clusters:
	path = args.outDir+"/clusters_from_pairs/"
	#Check that all files are present, else terminate the program and send a report of the missing files
	if not args.local:
		check_files(args.outDir)
	#Load clusters
	species = set([])
	allClusters = {}
	totalNumCl = 0
	for dirName in glob.glob(path+"/*"):
		spe = dirName.split("/")[-1]
		clusters = set([])
		for fileName in glob.glob(path+"/"+spe+"/*txt"):
			for line in open(fileName):
				line = line.strip()
				clusters.add(line)
		clDir = dirName+"/fused_clusters/"
		create_folder(clDir)
		tmpDir = clDir+"/tmp/"
		create_folder(tmpDir)
		clFile = clDir+"/cluster_list.txt"
		outfile = open(clFile,"w")
		clusters2 = {}
		num = 1
		name = "CL_"+spe
		print "Processing "+spe
		for cl in clusters:
			code = "%s_%.5d" % (name,num)
			num += 1
			clusters2[code] = cl.split(";")
			print >>outfile,code+"\t"+cl
		outfile.close()
		#Run the mcl analysis on the clusters predicted within each species
		outfileMCL = clDir+"/"+spe+".out.mcl"
		run_mcl(clusters2,clDir,tmpDir,outfileMCL)
		groups = []
		for line in open(outfileMCL):
			line = line.strip()
			dades = line.split()
			groups.append(dades)
		outfileREPR = open(clDir+"/"+spe+".representative.txt","w")
		all_genes = list_all_genes(clusters2)
		taken_genes = get_most_representative_clusters(groups,clusters2,spe,outfileREPR,set_difference)
		lost_genes = all_genes.difference(taken_genes)
		#Get back whole clusters that were discarded
		taken_genes = fishing(lost_genes,clusters2,outfileREPR,clDir,tmpDir,set_difference)
		outfileREPR.close()
		#Load the newly fused clusters into memory
		allClusters[spe] = {}
		for line in open(clDir+"/"+spe+".representative.txt"):
			line = line.strip()
			dades = line.split("\t")
			allClusters[spe][dades[0]] = dades[1]
			totalNumCl += 1
	if totalNumCl == 0:
		exit("No putative clusters were found in this dataset")
	else:
		#Split clusters into files to process them in a cluster
		pathAllClusters = args.outDir+"/complete_cluster_list.txt"
		pathSpeClusters = args.outDir+"/clusters_by_spe/"
		create_folder(pathSpeClusters)
		outfile = open(pathAllClusters,"w")
		for spe in allClusters:
			numFile = 1
			num = 0
			for name in allClusters[spe]:
				if num % 250 == 0:
					if num != 0:
						outfileSpe.close()
					pathName = "%s/%s_%.5d" % (pathSpeClusters,spe,numFile)
					numFile += 1
					outfileSpe = open(pathName,"w")
				num += 1
				print >>outfile,name+"\t"+allClusters[spe][name]
				print >>outfileSpe,name+"\t"+allClusters[spe][name]
			outfileSpe.close()
		outfile.close()
		create_jobs(args.outDir,"comparison",args.thr)

#Make all the comparisons between the detected clusters	
if args.cluster_comparison:
	pathAllClusters = args.outDir+"/complete_cluster_list.txt"
	pathSpeClusters = args.inFile
	name = args.inFile.split("/")[-1].split(".")[0]
	refSpeClusters,species = load_clusters(pathSpeClusters)
	refAllClusters,species = load_clusters(pathAllClusters)
	conversion = load_conversion(args.outDir+"/conversion_files/",species)
	path1 = args.outDir+"/cluster_comparison/"
	create_folder(path1)
	pathOutFile = path1+"/"+name+".txt"
	calculate_partial_scores(refSpeClusters,refAllClusters,conversion,pathOutFile)

#After this is complete a mcl needs to be run in order to cluster the clusters.

#Filter results and create the final list of cluster families
if args.clusterFam_filter:
	#Load data into memory
	clusterFams = load_cluster_families(args.inFile)
	pathAllClusters = args.outDir+"/complete_cluster_list.txt"
	clusters,species = load_clusters(pathAllClusters)
	conversion = load_conversion(args.outDir+"/conversion_files/",species)
	total = len(clusterFams)
	num = 1
	#Open outfile
	outfile = open(args.outDir+"/final_clusters.txt","w")
	for fam in clusterFams:
		num +=1
		if num % 1000 == 0:
			print num,total
		#Delete outliers
		filteredFam = delete_outliers(clusterFams[fam],clusters)
		families = divide_overlaps(filteredFam,clusters)
		counter = {}
		#Counts how many species have a given protein family
		for spe in families:
			for cl in families[spe]:
				proteins = clusters[cl]
				for pr in proteins:
					contig = pr.split("_")[1]
					speName = pr.split("_")[0]
					if pr in conversion[speName][contig]:
						f = conversion[speName][contig][pr]
						if f not in counter:
							counter[f] = {}
						if spe not in counter[f]:
							counter[f][spe] = set([])
						counter[f][spe].add(pr)
		#The threshold determines that proteins need to be in at least 1/3 of the species in order to be kept
		thr = round(float(len(clusterFams[fam]))*0.33,0)
		new_clusters = {}
		for f in counter:
			if float(len(counter[f])) > thr:
				for spe in counter[f]:
					if spe not in new_clusters:
						new_clusters[spe] = []
					for code in counter[f][spe]:
						new_clusters[spe].append(code)
		#Complete cluster and delete clusters that are formed only by duplicated genes
		nearly_last = []
		for spe in new_clusters:
			cl = new_clusters[spe]
			cl = sorted(cl,key=lambda x: int(x.split("_")[2]))
			PASS = delete_multiple_duplications(spe,cl,conversion[spe.split("_")[0]])
			if PASS:
				cl = complete_cluster(cl)
				string = spe+"\t"+";".join(cl)
				nearly_last.append(string)
		#Print clusters that are found in at least two species
		if len(nearly_last) > 1:
			print >>outfile,"# "+fam
			for string in nearly_last:
				print >>outfile,string
	outfile.close()

#Completes cluster families with clusters that have been previously discarded

if args.clusterFam_complete:
	#First step is to compile all the information from the clusters previously predicted
	print "Compiling information, this can take a few hours depending on the amount of data you have"
	inDirThr = args.outDir+"/thresholds/"
	inDirCL = args.outDir+"/all_cluster_predictions/"
	species = [x.split("/")[-1] for x in glob.glob(inDirCL+"/*") if ".txt" not in x]
	species = compile_info(inDirCL,inDirThr)

	#Load the clusters into memory
	cluster_families,protsCL = load_cluster_families2(args.inFile)
	#Load the thresholds into memory
	thresholds = load_thresholds(inDirThr+"/all_thresholds.txt")
	#Open the outfile
	outfile = open(args.outDir+"/final_clusters.complemented.txt","w")
	#Load conversion information
	conversion = load_conversion(args.outDir+"/conversion_files/",species)
	#Load all previously predicted clusters into memory as long as they overlap less than 10% with the current cluster collection
	all_clusters = get_all_clusters(inDirCL,species,protsCL)
	for cf in cluster_families.keys():
		print >>outfile,"# "+cf
		species = set([x.split("_")[0] for x in cluster_families[cf]])
		#Search for candidates suitable to add to the family.
		candidates = get_candidates(cluster_families[cf],all_clusters,thresholds,maxSize,set_difference)
		#print candidates
		#Join all the predictions from the different clusters
		single_candidates = {}
		for spe1 in candidates:
			for spe2 in candidates[spe1]:
				cl = ";".join(candidates[spe1][spe2][0])
				if spe2 not in single_candidates:
					single_candidates[spe2] = set([])
				single_candidates[spe2].add(cl)
		#Fuse them together
		for spe in single_candidates:
			group = single_candidates[spe]
			group = fuse_clusters(group,protsCL)
			if group:
				print >>outfile,spe+"\t"+group
		for spe in cluster_families[cf]:
			print >>outfile,spe+"\t"+";".join(cluster_families[cf][spe])
	outfile.close()

		

