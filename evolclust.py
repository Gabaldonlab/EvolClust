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
#~ import ete3

########################################################################
# Operational modules
########################################################################

#Checks whether a folder exists or else creates it
def create_folder(name):
	if not os.path.exists(name):
		cmd = "mkdir "+name
		try:
			run_command(cmd,False)
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
		pairs.append(line)
	return pairs

#Load homologous pairs of proteins into memory
def load_pairs_by_protein(spe1,spe2,outDir):
	pairs = {}
	for line in open(outDir+"/"+spe1+"/"+spe2+".txt"):
		line = line.strip()
		p1,p2 = line.split("-")
		if p1 not in pairs:
			pairs[p1] = []
		pairs[p1].append(p2)
		if p2 not in pairs:
			pairs[p2] = []
		pairs[p2].append(p1)
	return pairs

#Load pre-calculated thresholds into memory
def load_thresholds(dirName,species):
	thresholds = {}
	fileName = dirName+"/"+species[0]+"/"+species[1]+".txt"
	if os.path.exists(fileName):
		for line in open(fileName):
			line = line.strip()
			if "Spe1" in line:
				pass
			else:
				dades = line.split()
				if dades[0] not in thresholds:
					thresholds[dades[0]] = {}
				if dades[1] not in thresholds[dades[0]]:
					thresholds[dades[0]][dades[1]] = {}
				thresholds[dades[0]][dades[1]][dades[2]] = float(dades[-1])
	fileName = dirName+"/"+species[1]+"/"+species[0]+".txt"
	if os.path.exists(fileName):
		for line in open(fileName):
			line = line.strip()
			if "Spe1" in line:
				pass
			else:
				dades = line.split()
				if dades[0] not in thresholds:
					thresholds[dades[0]] = {}
				if dades[1] not in thresholds[dades[0]]:
					thresholds[dades[0]][dades[1]] = {}
				thresholds[dades[0]][dades[1]][dades[2]] = float(dades[-1])
	return thresholds

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
			elif len(spe) != 5:
				pass
			else:
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
def build_pairs(fileName,outfileName):
	outfile = open(outfileName,"w")
	for line in open(fileName):
		line = line.strip()
		#~ print line
		dades = line.split()
		if len(dades) != 1:
			codes = []
			for d in dades:
				if "_" in d:
					spe = d.split("_")[0]
					if len(spe) == 5:
						codes.append(d)
			for i in range(0,len(codes)):
				c1 = codes[i]
				for j in range(i+1,len(codes)):
					c2 = codes[j]
					spe1 = c1.split("_")[0]
					spe2 = c2.split("_")[0]
					if spe1 != spe2:
						print >>outfile,spe1,spe2,c1,c2
	outfile.close()

#Divides the full list of pairs into pairwise species files
def split_file(outfileName,outDir):
	num = 1
	info = {}
	for line in open(outfileName):
		line = line.strip()
		dades = line.split()
		spe1,spe2 = dades[0],dades[1]
		if dades[2] > dades[3]:
			pair = dades[2]+"-"+dades[3]
		else:
			pair = dades[3]+"-"+dades[2]
		if spe1 not in info:
			info[spe1] = {}
		if spe2 not in info[spe1]:
			info[spe1][spe2] = set([])
		if spe2 not in info:
			info[spe2] = {}
		if spe1 not in info[spe2]:
			info[spe2][spe1] = set([])
		info[spe1][spe2].add(pair)
		info[spe2][spe1].add(pair)
		if num % 100000000 == 0:
			for spe1 in info:
				folderName = outDir+"/"+spe1
				create_folder(folderName)
				for spe2 in info[spe1]:
					outfile = open(folderName+"/"+spe2+".txt","a")
					for p in info[spe1][spe2]:
						print >>outfile,p
					outfile.close()
			info = {}
		num += 1
	for spe1 in info:
		folderName = outDir+"/"+spe1
		create_folder(folderName)
		for spe2 in info[spe1]:
			outfile = open(folderName+"/"+spe2+".txt","a")
			for p in info[spe1][spe2]:
				print >>outfile,p
			outfile.close()

#Creates jobs list that will automatize the process
def create_jobs(outDir,tagName):
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
		outfile1 = open(outJob+"jobs.step2.txt","w")
		for i in range(0,len(species)):
			spe1 = species[i]
			for j in range(0,len(species)):
				if i != j:
					spe2 = species[j]
					if i > j:
						print >>outfile,"python "+args.pathEvolClust+" -i "+args.inFile+" -s1 "+spe1+" -s2 "+spe2+" --calculate_thresholds -d "+args.outDir
						print >>outfile1,"python "+args.pathEvolClust+" -s1 "+spe1+" -s2 "+spe2+" --get_pairwise_clusters -d "+args.outDir
		outfile.close()
		outfile1.close()
	elif tagName == "comparison":
		outfile = open(outJob+"jobs.step3.txt","w")
		for fileName in glob.glob(outDir+"/clusters_by_spe/*"):
			print >>outfile,"python "+args.pathEvolClust+" -i "+fileName+" -d "+outDir+" --cluster_comparison"
		outfile.close()

########################################################################
#Genome walking - Main script
########################################################################

#Obtain conserved cluster given a pair of homologous proteins
def get_clusters(group,conversion,thresholds,outfile,printOut):
	allClustersA = {}
	allClustersB = {}
	allClustersC = {}
	seed1,seed2 = p.split("-")
	spe1 = seed1.split("_")[0]
	spe2 = seed2.split("_")[0]
	if spe1 != spe2:
		cluster1,cluster2,score = calculate_cluster(seed1,seed2,conversion)
		if len(cluster1) != 0 and len(cluster2) != 0:
			if printOut:
				print >>outfile,";".join(cluster1)+"\t"+";".join(cluster2)+"\t"+str(score)
			allClustersA,allClustersB,allClustersC = save_data(spe1,spe2,score,cluster1,thresholds,allClustersA,allClustersB,allClustersC)
			allClustersA,allClustersB,allClustersC = save_data(spe2,spe1,score,cluster2,thresholds,allClustersA,allClustersB,allClustersC)
	return allClustersA,allClustersB,allClustersC

def calculate_cluster(s1,s2,conversion):
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
		#Cut the cluster if there are stretches of more than three proteins in between, will return the piece where the seed is present as long as it's still larger or equal than 5 proteins
		cluster1 = cut_cluster(cluster1,s1)
		cluster1Trans = [conversion[spe1][contig1][x] for x in cluster1 if x in conversion[spe1][contig1]]
		#If the cluster still exists then the same is done for the opposite cluster
		if len(cluster1) != 0:
			cluster2 = delete_non_homologs(cluster2,cluster2Trans,cluster1Trans)
			cluster2 = cut_cluster(cluster2,s2)
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

#Checks the suitability of the data (size, passes threshold)
def save_data(spe1,spe2,score,cluster,thresholds,allClustersA,allClustersB,allClustersC):
	size = len(cluster)
	cluster = sorted(cluster,key=lambda x:int(x.split("_")[2]))
	cluster = ";".join(cluster)
	if int(size) > 35:
		if spe1 not in allClustersC:
			allClustersC[spe1] = set([])
		allClustersC[spe1].add(cluster)
		thr = None
	else:
		size = str(size)
		if size in thresholds[spe1][spe2]:
			thr = thresholds[spe1][spe2][size]
		else:
			thr = 0.0
		if score > thr:
			if spe1 not in allClustersA:
				allClustersA[spe1] = set([])
			allClustersA[spe1].add(cluster)
		else:
			if spe1 not in allClustersB:
				allClustersB[spe1] = set([])
			allClustersB[spe1].add(cluster)
	return allClustersA,allClustersB,allClustersC

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
def cut_cluster(cluster,seed):
	new_group = [cluster[0]]
	group_found = []
	for code in cluster[1:]:
		p1 = int(new_group[-1].split("_")[2])
		p2 = int(code.split("_")[2])
		diff = p2 - p1
		if diff <= 3:
			new_group.append(code)
		else:
			if seed in new_group and len(new_group) >= 5:
				group_found = new_group
			new_group = [code]
	if seed in new_group and len(new_group) >= 5:
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
# Calculate thresholds
########################################################################

def calculate_thresholds(outfileName,species,outDir,conversion,protein_list):
	#Calculates all thresholds between each pair of species. Thresholds are not bidirectional
	outfile = open(outfileName,"w")
	print >>outfile,"#Spe1 Spe2 Size Num_clusters Average Threshold"
	#Load data into memory
	listProteins = [x for x in protein_list if x.split("_")[0] == species[0]]
	listPairs = load_pairs_by_protein(species[0],species[1],outDir)
	randomProteins = get_random_proteins(listProteins,listPairs,1000)
	allScores = {}
	a = 0
	#For each of the selected proteins
	for s1 in randomProteins:
		if a % 250 == 0:
			p = int(float(a) / float(len(randomProteins)) * 100.00)
			print "THR - Processed: "+str(a)+" out of "+str(len(randomProteins))+"( "+str(p)+"%)"
		a += 1
		#Get the list of clusters sized 0 to 35
		clusters1 = create_clusters_different_sizes_forward(s1,protein_list)
		scores = {}
		#For each homologous protein to s1
		for s2 in listPairs[s1]:
			#Get a full list of clusters sized 5 to 35 that include this protein
			clusters2 = create_clusters_different_sizes_all(s2,protein_list)
			#For each size
			for num in clusters2:
				#Check that size is present in both cluster lines
				if num in clusters1 and num in clusters2:
					cl1 = clusters1[num].split(",")
					#For each of the pre-calculated cluster2
					for cluster2 in clusters2[num]:
						cl2 = cluster2.split(",")
						#Calculate the score for each pair of clusters
						score = calculate_score(cl1,cl2,conversion)
						#print cl1,cl2,score
						#Keep the score if it's the best one
						if num not in scores:
							scores[num] = score
						elif score > scores[num]:
							scores[num] = score
		#Save all scores
		for num in range(5,36):
			if num not in allScores:
				allScores[num] = []
			if num not in scores:
				allScores[num].append(0.0)
			else:
				allScores[num].append(scores[num])
	#Calculate all averages and thresholds
	for num in allScores:
		average = numpy.average(allScores[num])
		std = numpy.std(allScores[num])
		threshold = float(average) + float((std*2.0))
		if threshold > 1.0:
			threshold = 1.0
		print >>outfile,species[0],species[1],num,len(allScores[num]),average,threshold
	outfile.flush()
	outfile.close()

#Get the list of clusters starting from the seed and advancing a given number of proteins
def create_clusters_different_sizes_forward(seed,listProts):
	spe = seed.split("_")[0]
	contig = seed.split("_")[1]
	start = int(seed.split("_")[2])
	growing_cluster = []
	allClusters = {}
	for i in range(0,35):
		position = start + i
		code = "%s_%s_%.5d" % (spe,contig,position)
		if code in listProts:
			growing_cluster.append(code)
			if len(growing_cluster) >= 5:
				allClusters[len(growing_cluster)] = ",".join(growing_cluster)
	return allClusters

#Create all possible clusters that contain a seed protein given different cluster sizes
def create_clusters_different_sizes_all(seed,listProts):
	spe = seed.split("_")[0]
	contig = seed.split("_")[1]
	base = spe+"_"+contig
	allClusters = {}
	#For all cluster sizes
	for a in range(5,36):
		allClusters[a] = []
		#Obtain the starting point
		start = int(seed.split("_")[2]) - a + 1
		if start < 1:
			start = 1
		#For each position
		for b in range(start,start+a):
			#Obtain the cluster
			cluster = [("%s_%.5d" % (base,c)) for c in range(b,b+a)]
			#Check that all proteins exist (i.e. they are in the contig)
			cluster2 = [x for x in cluster if x in listProts]
			if len(cluster2) == a:
				cluster2 = ",".join(cluster2)
				allClusters[a].append(cluster2)
	return allClusters

#Obtain a list of random pairs
def get_random_proteins(listProteins,listPairs,numProts):
	tried = set([])
	if len(listProteins) <= numProts:
		randomProteins = listProteins
	else:
		randomProteins = set([])
		while len(randomProteins) < numProts and len(tried) < len(listProteins):
			num = randint(0,len(listProteins)-1)
			p = listProteins[num]
			if p in listPairs:
				randomProteins.add(listProteins[num])
				tried.add(p)
	return randomProteins

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
	if common1 < 5 or common2 < 5:
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
	#~ if score > 0:
	#~ print ",".join(list(cluster1)),"|",",".join(list(cluster2)),common,missing,size,score
	return score

#Calculates score between pairs of clusters
def calculate_partial_scores(clustersSpe,clustersAll,conversion,outfileName):
	outfile = open(outfileName,"w")
	codesSpe = clustersSpe.keys()
	codesAll = clustersAll.keys()
	for i in range(0,len(codesSpe)):
		cluster1 = clustersSpe[codesSpe[i]]
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
	if len(counter) >= 4:
		PASS = True
	else:
		PASS = False
	#print spe,PASS,counter
	return PASS

parser = argparse.ArgumentParser(description="Will perform the genome walking")
parser.add_argument("-i","--infile",dest="inFile",action="store",default=None,help="mcl or pairs file, can hold multiple or single families.")
parser.add_argument("-f","--fastafile",dest="fastaFile",action="store",default=None,help="Fasta file that contains the complete proteome database.")
parser.add_argument("-s1","--species1",dest="species1",action="store",default=None,help="Species tag needed to run the threshold calculation")
parser.add_argument("-s2","--species2",dest="species2",action="store",default=None,help="Species tag needed to run the threshold calculation")
parser.add_argument("-d","--outdir",dest="outDir",action="store",default="/users/tg/mmarcet/HGT3_proteomes/genome_walking3/",help="basepath folder where the results will be stored")
parser.add_argument("--build_jobs",dest="buildJobs",action="store_true",help="Will prepare the files so that trees can be build. Needs to have the -i and -d options as full paths")
parser.add_argument("--initial_files",dest="conv",action="store_true",help="Will prompt the program to create the initial files")
parser.add_argument("--get_pairwise_clusters",dest="calc_scores_pairs",action="store_true",help="Main program in the genome walking - starts from pairs files")
parser.add_argument("--calculate_thresholds",dest="calc_thr",action="store_true",help="Will calculate the thresholds for each pair of species")
parser.add_argument("--cluster_comparison",dest="cluster_comparison",action="store_true",help="Will compares clusters")
parser.add_argument("--filter_clusters",dest="filter_clusters",action="store_true",help="Will filter out identical clusters and name them all; will also create the files needed to run the final cluster comparison")
parser.add_argument("--filter_cluster_families",dest="clusterFam_filter",action="store_true",help="Will delete redundancy among clusters within cluster families")
parser.add_argument("--path_evol_cluster",dest="pathEvolClust",action="store",default="./evolclust.py",help="Path to the python program evolclust.py")
args = parser.parse_args()

#Create output folder if it doesn't exist
create_folder(args.outDir)

#Creates conversion files that will correlate each protein to the corresponding protein family
#RUNS FOR ALL DATA - Checked
if args.conv:
	outDir = args.outDir+"/conversion_files/"
	create_folder(outDir)
	build_conversion_files(args.inFile,outDir)
	outfileName = args.outDir+"/complete_protein_list.txt"
	build_complete_protein_list(args.fastaFile,outfileName)
	#Creates files where all the pairs of homologous proteins between two species according to the mcl can be found
	outfileName = args.outDir+"/complete_pair_file.txt"
	outDir = args.outDir+"/pairs_files/"
	create_folder(outDir)
	if not os.path.exists(outfileName):
		build_pairs(args.inFile,outfileName)
		split_file(outfileName,outDir)
	create_jobs(args.outDir,"initial")

#This script will calculate the thresholds for a pair of species and clusters of different sizes.
#RUN IN CLUSTER FOR PAIRS OF SPECIES - checked
if args.calc_thr:
	outDir = args.outDir+"/pairs_files/"
	if not os.path.exists(outDir):
		exit("Before calculating the thresholds the --conversion function needs to be run")
	if args.species1 == None or args.species2 == None:
		exit("Species tags need to be established for calculating their threshold")
	#Create folders and file names
	outThr = args.outDir+"/thresholds/"
	create_folder(outThr)
	outThr1 = outThr+"/"+args.species1+"/"
	create_folder(outThr1)
	outfileName1 = outThr1+"/"+args.species2+".txt"
	outThr2 = outThr+"/"+args.species2+"/"
	create_folder(outThr2)
	outfileName2 = outThr2+"/"+args.species1+".txt"
	#Loads the list of all proteins encoded in the proteomes in use
	protein_list = set([line.strip() for line in open(args.outDir+"/complete_protein_list.txt")])
	#Loads the conversion files for the species to be treated
	species1 = [args.species1,args.species2]
	species2 = [args.species2,args.species1]
	conversion = load_conversion(args.outDir+"/conversion_files/",species1)
	#Calculate thresholds from spe1 to spe2
	calculate_thresholds(outfileName1,species1,outDir,conversion,protein_list)
	#Calculate thresholds from spe2 to spe1
	calculate_thresholds(outfileName2,species2,outDir,conversion,protein_list)

#This script will obtain a list of clusters that have a higher conservation score than the cluster
#RUN IN CLUSTER FOR PAIRS OF SPECIES
if args.calc_scores_pairs:
	#Load previously calculated data
	species = [args.species1,args.species2]
	pairsDir = args.outDir+"/pairs_files/"
	pairs = load_pairs(species[0],species[1],pairsDir)
	conversion = load_conversion(args.outDir+"/conversion_files/",species)
	thresholds = load_thresholds(args.outDir+"/thresholds/",species)
	#Create folders
	path = args.outDir+"/clusters_from_pairs/"
	create_folder(path)
	#This will print extra information about the threshold
	detailsOutfile = args.outDir+"/all_clusters/"
	create_folder(detailsOutfile)
	detailsOutfile = detailsOutfile+"/"+species[0]
	create_folder(detailsOutfile)
	detailsOutfile = detailsOutfile+"/"+species[1]+".txt"
	outfile = open(detailsOutfile,"w")
	allClusters = {}
	total = len(pairs)
	num = 0
	for p in pairs:
		if num % 1000 == 0:
			print "Processed: ",num,"out of",total
		num += 1
		printOut = True
		clustersPASS,clustersFAIL,clustersLONG = get_clusters(p,conversion,thresholds,outfile,printOut)
		for spe in clustersPASS:
			if spe not in allClusters:
				allClusters[spe] = set([])
			for cl in clustersPASS[spe]:
				allClusters[spe].add(cl)
	outfile.close()
	for spe in allClusters:
		if spe == species[1]:
			path1 = path+"/"+species[1]
			create_folder(path1)
			outfile = open(path1+"/"+species[0]+".txt","w")
		else:
			path1 = path+"/"+species[0]
			create_folder(path1)
			outfile = open(path1+"/"+species[1]+".txt","w")
		for cl in allClusters[spe]:
			print >>outfile,cl
		outfile.close()

#This will create a single list with all the clusters found, filtering out identical clusters.
#NEEDS TO BE RUN AFTER ALL CLUSTERS HAVE BEEN CALCULATED
if args.filter_clusters:
	path = args.outDir+"/clusters_from_pairs/"
	#Load clusters
	species = set([])
	clusters = {}
	for fileName in glob.glob(path+"/*/*"):
		spe = fileName.split("/")[-2]
		if spe not in clusters:
			clusters[spe] = set([])
		for line in open(fileName):
			line = line.strip()
			clusters[spe].add(line)
	#Assign name to clusters
	pathAllClusters = args.outDir+"/complete_cluster_list.txt"
	pathSpeClusters = args.outDir+"/clusters_by_spe/"
	create_folder(pathSpeClusters)
	outfile = open(pathAllClusters,"w")
	for spe in clusters:
		name = "CL_"+spe
		num = 0
		numFile = 1
		for cl in clusters[spe]:
			if num % 250 == 0:
				if num != 0:
					outfileSpe.close()
				pathName = "%s/%s_%.5d" % (pathSpeClusters,spe,numFile)
				numFile += 1
				outfileSpe = open(pathName,"w")
			code = "%s_%.5d" % (name,num)
			num += 1
			print >>outfile,code+"\t"+cl
			print >>outfileSpe,code+"\t"+cl
		outfileSpe.close()
	outfile.close()
	create_jobs(args.outDir,"comparison")

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
			print >>outfile,"###############################################"
			print >>outfile,"# "+fam
			print >>outfile,"###############################################"
			for string in nearly_last:
				print >>outfile,string
	outfile.close()

		

