#!/usr/bin/env python

#Calculate support for clusters based on creating 100 random genomes and seeing if the clusters are found there.

import argparse
import glob
import random
import os,sys
import subprocess as sp

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

#Load the original conversion information for a given species
def load_conversion(fileName):
	species = set([])
	conversion = {}
	proteins = {}
	spe = fileName.split("/")[-1].split(".")[0]
	species.add(spe)
	conversion[spe] = {}
	proteins[spe] = []
	for line in open(fileName):
		line = line.strip()
		dades = line.split("\t")
		conversion[spe][dades[0]] = dades[1]
		proteins[spe].append(dades[0])
	return species,proteins,conversion

#Create a set of random genomes and print them into files. They will serve as conversion files for these random genomes
def create_random_genomes(proteins,conversion,num_random_genomes,outFolder):
	prot_families = conversion.values()
	genome_num = 0
	while genome_num < num_random_genomes:
		random.shuffle(prot_families)
		outfile = open(outFolder+"/genome_"+str(genome_num)+".txt","w")
		for a in range(0,len(proteins)):
			print >>outfile,"r"+proteins[a]+"\t"+prot_families[a]
		outfile.close()
		genome_num += 1

#Create the matching pair files for the random genomes. The pairs contain set of homologous proteins
def create_pairs_files(conversionFolder,conversion,pairsFolder):
	protFams1 = {}
	for code in conversion:
		val = conversion[code]
		if val not in protFams1:
			protFams1[val] = []
		protFams1[val].append(code)
	for fileName in glob.glob(conversionFolder+"/*"):
		tag = fileName.split("/")[-1].split("_")[1].split(".")[0]
		protFams2 = {}
		for line in open(fileName):
			line = line.strip()
			dades = line.split("\t")
			if dades[1] not in protFams2:
				protFams2[dades[1]] = []
			protFams2[dades[1]].append(dades[0])
		outfile = open(pairsFolder+"/pairs_"+tag+".txt","w")
		for fam in protFams1:
			for c1 in protFams1[fam]:
				for c2 in protFams2[fam]:
					print >>outfile,c1.split("_")[0]+"\tr"+c2.split("_")[0]+"\t"+c1+"\t"+c2
		outfile.close()

#Load clusters into memory		
def load_clusters(fileName):
	clusters = {}
	for line in open(fileName):
		line = line.strip()
		if "# " in line:
			CLname = line.split(" ")[1]
		elif "##" in line:
			pass
		else:
			dades = line.split("\t")
			spe = dades[0].split("_")[0]
			if spe not in clusters:
				clusters[spe] = {}
			if CLname not in clusters[spe]:
				clusters[spe][CLname] = {}
			clusters[spe][CLname][dades[0]] = dades[1]
	return clusters 

#Load the conversion file; either the original or the random genome. They are both loaded into the same dictionary
def load_conversion2(fileName,conversion):
	for line in open(fileName):
		line = line.strip()
		dades = line.split("\t")
		spe = dades[0].split("_")[0]
		contig = dades[0].split("_")[1]
		if spe not in conversion:
			conversion[spe] = {}
		if contig not in conversion[spe]:
			conversion[spe][contig] = {}
		conversion[spe][contig][dades[0]] = dades[1]
	return conversion

#Load pairs of homologous proteins into memory
def load_pairs(fileName,list_proteins):
	pairs = {}
	for line in open(fileName):
		line = line.strip()
		dades = line.split("\t")
		if dades[2] in list_proteins:
			p = dades[2]+"-"+dades[3]
			if dades[2] not in pairs:
				pairs[dades[2]] = set([])
			pairs[dades[2]].add(p)
	return pairs

#Return a list of all proteins found in the clusters	
def get_prots_in_clusters(clusters):
	prots = set([])
	for code in clusters:
		for tag in clusters[code]:
			ps = clusters[code][tag].split(";")
			prots = prots.union(ps)
	return prots

#Search for clusters in the random genome	
def search_for_clusters(clusters,genomeRandomFolder,pairsRandomFolder,original_conversion_file):
	#Get the list of proteins in the clusters
	list_proteins = get_prots_in_clusters(clusters)
	found_clusters = {}
	for fileName in glob.glob(genomeRandomFolder+"/*"):
		#Load conversion file
		conversion1 = {}
		conversion1 = load_conversion2(original_conversion_file,conversion1)
		tag = fileName.split("/")[-1].split("_")[1].split(".")[0]
		conversion1 = load_conversion2(fileName,conversion1)
		#Load pairs
		pairs = load_pairs(pairsRandomFolder+"/pairs_"+tag+".txt",list_proteins)
		#For each family
		for CF in clusters:
			#For each cluster in the family
			for tag in clusters[CF]:
				prots = clusters[CF][tag].split(";")
				#For each protein in the cluster
				for p in prots:
					if p in pairs:
						#For each pair of homologous proteins
						for pair in pairs[p]:
							p1,p2 = pair.split("-")
							if p1 in prots:
								#Get a potential cluster in the random genome and its score
								cl1,cl2,score = get_clusters(p1,p2,conversion1,set_difference,prots)
								#If the score is above 0.1 we consider that a cluster has been found and this is saved into memory
								if score >= 0.1:
									clA = ";".join(cl1)
									clB = ";".join(cl2)
									if clA not in found_clusters:
										found_clusters[clA] = {}
									if fileName not in found_clusters[clA]:
										found_clusters[clA][fileName] = [clB,score]
	return found_clusters
									

#Get matching clusters and their conservation score		
def get_clusters(p1,p2,conversion,set_difference,cluster1):
	spe1 = p1.split("_")[0]
	spe2 = p2.split("_")[0]
	if spe1 != spe2:
		cluster1,cluster2,score = calculate_cluster(p1,p2,conversion,set_difference,cluster1)
	return cluster1,cluster2,score

#Given a cluster and a pair of homologous proteins, find a cluster in the other genome.
def calculate_cluster(s1,s2,conversion,set_difference,cluster1):
	#Obtain contigs where the seed is located
	spe1,contig1 = s1.split("_")[0],s1.split("_")[1]
	#Sort proteins in the cluster
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
	# ~ print s1,"Enter loop ",size1,size2
	continua = True
	while continua:
		#Delete proteins that are not present in the opposite cluster
		cluster1 = delete_non_homologs(cluster1,cluster1Trans,cluster2Trans)
		#Cut the cluster if there are stretches of more than three proteins in between, will return the piece where the seed is present as long as it's still larger or equal than the minSize
		#If the cluster still exists then the same is done for the opposite cluster
		if len(cluster1) != 0:
			cluster1 = cut_cluster(cluster1,s1,set_difference)
			cluster1Trans = [conversion[spe1][contig1][x] for x in cluster1 if x in conversion[spe1][contig1]]
			cluster2 = delete_non_homologs(cluster2,cluster2Trans,cluster1Trans)
			if len(cluster2) != 0:
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

#Calculate the CS between the two clusters
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

parser = argparse.ArgumentParser(description="Will perform the genome walking")
parser.add_argument("-s","--species_tag",dest="speTAG",action="store",required=True,help="Species tag")
parser.add_argument("-i","--clustersFile",dest="clustersFile",action="store",required=True,help="File with final clusters")
parser.add_argument("-d","--inFolder",dest="inFolder",action="store",required=True,help="Folder where the analysis was run")
parser.add_argument("--non_homologs",dest="non_homologs",action="store",default=3,help="Number of non-homologous genes needed to split a cluster")
parser.add_argument("--minSize",dest="minSize",action="store",default=5,help="Minimum size a cluster needs to have to be considered. Default is set to 5")
parser.add_argument("--numGenomes",dest="numGenomes",action="store",default=100,help="Number of random genomes you want to consider")
args = parser.parse_args()

#Load the conversion between proteins and their protein families for the chosen species			
species,proteins,conversion = load_conversion(args.inFolder+"/conversion_files/"+args.speTAG+".txt")

#Create the random folders where analyses and results will be stored
randomFolder = args.inFolder+"/random"
create_folder(randomFolder)
randomFolderResults = args.inFolder+"/random/results/"
create_folder(randomFolderResults)

#Load clusters into memory
clusters = load_clusters(args.clustersFile)


species = list(species)
minSize = args.minSize
set_difference = args.non_homologs
num_genomes = args.numGenomes
for spe in list(species):
	if spe in clusters:
		#Create random folders for the species
		speRandomFolder = randomFolder+"/"+spe
		create_folder(speRandomFolder)
		genomeRandomFolder = speRandomFolder+"/random_genomes/"
		create_folder(genomeRandomFolder)
		#Create a certain number of random genomes
		create_random_genomes(proteins[spe],conversion[spe],num_genomes,genomeRandomFolder)
		#Obtain the list of homologous proteins for the newly created genomes
		pairsRandomFolder = speRandomFolder+"/random_pairs/"
		create_folder(pairsRandomFolder)
		create_pairs_files(genomeRandomFolder,conversion[spe],pairsRandomFolder)
		#Search for the clusters in the newly create genomes
		found_clusters = search_for_clusters(clusters[spe],genomeRandomFolder,pairsRandomFolder,args.inFolder+"/conversion_files/"+spe+".txt")
		#If any cluster has been found, print it into a results file
		if len(found_clusters) != 0:
			outfile = open(randomFolderResults+"/"+spe+".txt","w")
			for clA in found_clusters:
				for tag in found_clusters[clA]:
					print >>outfile,clA+"\t"+tag+"\t"+found_clusters[clA][tag][0]+"\t"+str(found_clusters[clA][tag][1])
			outfile.close()
		else:
			outfile = open(randomFolderResults+"/"+spe+".txt","w")
			print >>outfile,"No problematic clusters were found"
			outfile.close()
