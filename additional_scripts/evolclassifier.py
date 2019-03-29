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

import argparse
import ete3
import glob
import os,sys
import subprocess as sp
import numpy

########################################################################
# Operational modules
########################################################################

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

########################################################################
# Loading data modules
########################################################################

#Load clusters into memory - general method
def load_clusters_general(fileName):
	clusters = {}
	for line in open(fileName):
		line = line.strip()
		if "##" in line:
			pass
		elif "#" in line:
			name = line.split()[1]
			clusters[name] = {}
		else:
			codes = line.split()[1].split(";")
			clusters[name][line.split()[0]] = codes
	return clusters

#Load one single cluster family into memory and create pairwise combinations - detailed_method
def load_one_cluster(fileName,famName,outfileName):
	clFam = {}
	proteins = {}
	outfile = open(outfileName,"w")
	for line in open(fileName):
		line = line.strip()
		if "##" in line:
			pass
		elif "#" in line:
			name = line.split()[1]
			if name == famName:
				remember = True
				print >>outfile,line
			else:
				remember = False
		else:
			if remember:
				print >>outfile,line
				codes = line.split()[1].split(";")
				clFam[line.split()[0]] = codes
				spe = line.split()[0].split("_")[0]
				if spe not in proteins:
					proteins[spe] = set([])
				for c in codes:
					proteins[spe].add(c)
	outfile.close()
	clusters = {}
	species = clFam.keys()
	for i in range(0,len(species)):
		spe1 = species[i].split("_")[0]
		code1 = species[i]
		for j in range(i+1,len(species)):
			spe2 = species[j].split("_")[0]
			code2 = species[j]
			if spe1 != spe2:
				name = code1+"*"+code2
				clusters[name] = {}
				clusters[name][spe1] = clFam[code1]
				clusters[name][spe2] = clFam[code2]
	return clusters,proteins

#Combines clusters grouping those that evolve vertically
def combine_clusters(clusterFile,groups,nameFam):
	clFam = {}
	for line in open(clusterFile):
		line = line.strip()
		if "##" in line:
			pass
		elif "#" in line:
			name = line.split()[1]
			if name == nameFam:
				take = True
			else:
				take = False
		else:
			if take:
				code = line.split("\t")[0]
				taken = None
				for group in groups:
					if code in group:
						taken = "|".join(group)
				if not taken:
					clFam[code] = {}
					clFam[code][code] = line.split()[1].split(";")
				else:
					if taken not in clFam:
						clFam[taken] = {}
					clFam[taken][code] = line.split()[1].split(";")
	clusters = {}
	species = clFam.keys()
	for i in range(0,len(species)):
		code1 = species[i]
		for j in range(i+1,len(species)):
			code2 = species[j]
			name = code1+"*"+code2
			clusters[name] = {}
			for p in clFam[code1]:
				clusters[name][p] = clFam[code1][p]
			for p in clFam[code2]:
				clusters[name][p] = clFam[code2][p]
	return clusters	


#Load trees into memory and correlate them to their family number
def speName(node):
	return node.split("_")[0]

def whole_name(node):
	return node

def load_all_trees(fileName):
	trees = {}
	for line in open(fileName):
		line = line.strip()
		dades = line.split()
		t = ete3.PhyloTree(dades[1],sp_naming_function=speName)
		code = dades[0]
		trees[code] = t
	return trees

def load_detailed_trees(fileName,families):
	trees = {}
	for line in open(fileName):
		line = line.strip()
		dades = line.split()
		if dades[0] in families:
			t = ete3.PhyloTree(dades[1],sp_naming_function=speName)
			code = dades[0]
			trees[code] = t
	return trees

#Loads the species tree
def load_spTree(fileName):
	spTree = ete3.PhyloTree(fileName,sp_naming_function=whole_name)
	for leaf in spTree.iter_leaves():
		if len(leaf.name) != 5:
			leaf.delete()
	return spTree

#Loads conversion into memory
def load_all_conversion(conversionFolder):
	conversion = {}
	for fileName in glob.glob(conversionFolder+"/*txt"):
		spe = fileName.split("/")[-1].split(".")[0]
		conversion[spe] = {}
		for line in open(fileName):
			line = line.strip()
			dades = line.split()
			conversion[spe][dades[0]] = dades[1]
	return conversion

#Loads conversion into memory for a set of proteins
def load_detailed_conversion(conversionFolder,proteins):
	conversion = {}
	families = set([])
	for spe in proteins:
		fileName = conversionFolder+"/"+spe+".txt"
		conversion[spe] = {}
		for line in open(fileName):
			line = line.strip()
			dades = line.split()
			if dades[0] in proteins[spe]:
				conversion[spe][dades[0]] = dades[1]
				families.add(dades[1])
	return conversion,families

#Load taxonomic information into memory
def load_taxonomy(fileName,THR):
	taxonomy = {}
	tax_thr = set([])
	for line in open(fileName):
		line = line.strip()
		dades = line.split("\t")
		taxonomy[dades[1]] = dades[3:]
		if THR == "order":
			for t in dades[7:]:
				tax_thr.add(t)
		elif THR == "class":
			for t in dades[6:]:
				tax_thr.add(dades[6])
	return taxonomy,tax_thr

########################################################################
# Obtain main statistics
########################################################################

#Count the number of proteins in the cluster belong to each protein family and send back the list
def count_protein_presence(clusters,conversion):
	protCounter = {}
	for code in clusters:
		spe = code.split("_")[0]
		for name in clusters[code]:
			if name in conversion[spe]:
				nameTrans = conversion[spe][name]
				if nameTrans not in protCounter:
					protCounter[nameTrans] = []
				protCounter[nameTrans].append(name)
	return protCounter

#Check which proteins appear in a given tree
def check_protein_presence(t,proteins):
	proteins2 = []
	counter = {}
	for p in proteins:
		if p in t.get_leaf_names():
			proteins2.append(p)
			spe = p.split("_")[0]
			if spe not in counter:
				counter[spe] = 0
			counter[spe] += 1
	found = False
	for code in counter:
		if counter[spe] != 1:
			found = True
	if found:
		proteins2 = []
	return proteins2

#Check which proteins appear in a given tree and see whether both groups we are comparing have at least one representative
def check_protein_presence1(t,proteins,clusterName):
	speciesA = [x.split("_")[0] for x in clusterName.split("*")[0].split("|")]
	speciesB = [x.split("_")[0] for x in clusterName.split("*")[1].split("|")]
	foundA = False
	foundB = False
	proteins2 = []
	for p in proteins:
		if p in t.get_leaf_names():
			spe = p.split("_")[0]
			if spe in speciesA:
				foundA = True
			if spe in speciesB:
				foundB = True
			proteins2.append(p)
	return proteins2,foundA,foundB

#Checks whether a group of proteins / species are monophyletic in a given tree	
def check_monophyly(t,proteins):
	#~ print proteins
	anc = t.get_common_ancestor(proteins)
	if len(anc.get_leaves()) == len(proteins):
		monophyly = True
	else:
		monophyly = False
	return monophyly,anc

#Counts the number of losses inferred in a species tree given a list of presences
def count_losses(spTree,speciesList):
	anc = spTree.get_common_ancestor(speciesList)
	taken = []
	for node in anc.traverse():
		leaves = set(node.get_leaf_names())
		present_leaves = leaves.intersection(speciesList)
		if len(present_leaves) == 0:
			found = False
			for species in taken:
				present_leaves = leaves.intersection(species)
				if len(present_leaves) != 0:
					found = True
			if not found:
				taken.append(leaves)
	return len(taken)

#Check the taxonomic common ancestor for a group of species
def check_taxonomic_scope(additionalGenes,taxonomy,group):
	species = [x.split("_")[0] for x in additionalGenes]
	base_tax = taxonomy[species[0]]
	minim = len(base_tax)
	for spe in species[1:]:
		for i in range(0,len(base_tax)):
			b = base_tax[i]
			if b not in taxonomy[spe]:
				if i < minim:
					minim = i
	base_tax = base_tax[:minim]
	found = False
	for b in base_tax:
		if b in group:
			found = True
	return found,base_tax
	
#Obtain statistics of a protein family within a cluster		
def analyse_protein_family(protCounter,fprots,t,proteins,prot):
	mainInfo = {}
	mainInfo["prot"] = len(protCounter[prot])
	mainInfo["protInTree"] = len(fprots)
	#Check whether gene tree is monophyletic
	monophyly,anc = check_monophyly(t,fprots)
	mainInfo["GTmonophyly"] = monophyly
	#Check whether the common ancestor of the species involved is a speciation node or a duplication node
	species1 = set([x for x in anc.get_children()[0].get_species()])
	species2 = set([x for x in anc.get_children()[1].get_species()])
	common = species1.intersection(species2)
	if len(common) > 0:
		mainInfo["EvolNode"] = "Duplication"
	else:
		mainInfo["EvolNode"] = "Speciation"
	#Check number of species
	species = set([x.split("_")[0] for x in protCounter[prot]])
	mainInfo["species"] = len([x for x in species if x in spTree.get_leaf_names()])
	#Get the number of additional proteins in the tree
	allGenes = [x for x in anc.get_leaf_names() if x.split("_")[0] in spTree.get_leaf_names()]
	additionalGenes = [x for x in allGenes if x not in proteins]
	mainInfo["addGenes"] = len(additionalGenes)
	if len(additionalGenes) != 0:
		#Check whether they surpass the taxonomic range
		found,scope = check_taxonomic_scope(additionalGenes,taxonomy,tax_thr)
		mainInfo["RangeTHROther"] = found
		mainInfo["taxScopeOther"] = scope[-1]
	else:
		mainInfo["RangeTHROther"] = "-"
		mainInfo["taxScopeOther"] = "-"
	found,scope = check_taxonomic_scope(fprots,taxonomy,tax_thr)
	mainInfo["RangeTHROwn"] = found
	mainInfo["taxScopeOwn"] = scope[-1]
	#Calculate the treeKO distance between the subtree and the species tree
	mainInfo["TK"] = "Unknown"
	#Check monophyly and number of losses inferred from the species tree
	species = [x.split("_")[0] for x in protCounter[prot]]
	species2 = [x for x in species if x in spTree.get_leaf_names()]
	if len(species2) > 1:
		monophyly,anc = check_monophyly(spTree,species2)
		numLosses = count_losses(spTree,species2)
		mainInfo["STmonophyly"] = monophyly
		mainInfo["Losses"] = numLosses
	elif len(species2) == 1:
		mainInfo["STmonophyly"] = True
		mainInfo["Losses"] = 0
	else:
		mainInfo["STmonophyly"] = "-"
		mainInfo["Losses"] = "-"
	return mainInfo

########################################################################
# Obtain tree distances
########################################################################

#Calculate treeKO distance between species tree and gene subtree
def calculate_tree_distance(tree,spTree):
	speciesR = spTree.get_species()
	speciesT = tree.get_species()
	common = speciesR.intersection(speciesT)
	if len(common) <= 2:
		print "Not enough species for comparison"
		distance = None
	else:
		for leaf in tree.iter_leaves():
			if leaf.species not in common:
				leaf.delete()
		tree = tree.collapse_lineage_specific_expansions()
		try:
			info = tree.compare(spTree,has_duplications=True,source_tree_attr='species')
		except:
			info = {}
		#print "TREE",tree.write()
		#~ print "SPTREE",spTree.write()
		#~ print "INFO",info
		if len(info) == 0:
			distance = "Unable to compute distance"
		else:
			if info["effective_tree_size"] != 0.0:
				distance = info["treeko_dist"]
			else:
				distance = None
	return distance

#Delete species specific duplications
def collapse_lineage_specific_expansions(t,seed):
	events = t.get_descendant_evol_events()
	taken = set([])
	for ev in events:
		if ev.etype == "D":
			if len(ev.node.get_species()) == 1:
				leaves = set(ev.node.get_leaf_names())
				common = leaves.intersection(taken)
				if len(common) == 0:
					for leaf in leaves:
						taken.add(leaf)
					if seed in leaves:
						chosen = seed
					else:
						chosen = list(leaves)[0]
					ev.node.add_feature("collapse",chosen)
	for node in t.traverse():
		if "collapse" in node.features:
			if node.is_root():
				t = None
			else:
				parent = node.up
				newName = node.collapse
				node.detach()
				parent.add_child(name=newName)
	return t


########################################################################
# Determine final calls
########################################################################

#Given two different factors and given a line of stats it decides the kind of evolution that has happened
#VE = Vertical evolution, CE = convergent evolution / independent formation

def call_tree_event(dades,fP,fL,tax_thr):
	monophylyGT = dades[4]
	monophylyST = dades[5]
	event_type = dades[6]
	try:
		numLosses = int(dades[7])
	except:
		numLosses = 0
	numAddGenes = int(dades[8])
	if len(tax_thr) == 0:
		taxThrOther = False
		taxThrOwn = False
	else:
		taxOther = dades[9]
		if taxOther in tax_thr:
			taxThrOther = True
		else:
			taxThrOther = False
		taxOwn = dades[10]
		if taxOwn in tax_thr:
			taxThrOwn = True
		else:
			taxThrOwn = False
	species = int(dades[3])
	proteins = int(dades[2])
	lossTHR = None
	conclusion = "FAILED"
	a = 0
	if monophylyGT == "True":
		if monophylyST == "True":
			conclusion = "VE"
		elif monophylyST == "False":
			lossTHR = int(species*fL)
			if numLosses > lossTHR:
				if taxThrOwn == False:
					conclusion = "HGT"
				else:
					conclusion = "VE"
			else:
				conclusion = "VE"
	else:
		addTHR = int(proteins*fP)
		if addTHR < numAddGenes:
			if taxThrOther == False:
				conclusion = "CE"
			else:
				conclusion = "VE"
		else:
			if monophylyST == "True":
				conclusion = "VE"
			elif monophylyST == "False":
				lossTHR = int(species*fL)
				if numLosses > lossTHR:
					if taxThrOwn == False:
						conclusion = "HGT"
					else:
						conclusion = "VE"
				else:
					conclusion = "VE"
	return conclusion

#Print tables and calculate final calling
def print_table(info,outfileName,distances):
	outfile = open(outfileName,"w")
	vertical_evol = 0
	print >>outfile,"#ClName\tVE\tCE\tHGT\tMGR_distance\tDecision"
	for cl in info:
		total = 0
		counter = {}
		for prot in info[cl]:
			for n in info[cl][prot]:
				t = info[cl][prot][n]
				if t not in counter:
					counter[t] = 0
				counter[t] += 1
				total += 1
		tags = ["VE","CE","HGT"]
		string = cl
		decision = None
		for t in tags:
			if t in counter:
				p = float(counter[t])/float(total) * 100.0
				if p > args.decision_threshold:
					decision = t
				string += "\t%.2f" % (p)
			else:
				string += "\t0.00"
		if cl in distances:
			string += "\t"+str(distances[cl])
			if distances[cl] <= args.mgr_threshold:
				if not decision:
					pass
				else:
					if decision == "CE":
						decision = None
		else:
			string += "\tNA"
		if not decision:
			string += "\tUndecided"
		else:
			if decision == "VE":
				vertical_evol += 1
			string += "\t"+decision
		print >>outfile,string
	outfile.close()
	return vertical_evol

def format_table(inFile,outFile):
	info = {}
	for line in open(inFile):
		line = line.strip()
		if "#" not in line:
			dades = line.split("\t")
			spe1,spe2 = dades[0].split("*")
			if spe1 not in info:
				info[spe1] = {}
			if spe2 not in info:
				info[spe2] = {}
			info[spe1][spe2] = dades[-1]
			info[spe2][spe1] = dades[-1]
	outfile = open(outFile,"w")
	species = info.keys()
	print >>outfile,"\t"+"\t".join(species)
	for s1 in species:
		string = s1
		for s2 in species:
			if s1 == s2:
				string +="\t-"
			elif s2 not in info[s1]:
				string += "\tUnrelated"
			else:
				string += "\t"+info[s1][s2]
		print >>outfile,string
	outfile.close()

#Print an image of the species tree
def draw_tree(species,spTree,outfileNameTree):
	anc = spTree.get_common_ancestor(species)
	ts = ete3.TreeStyle()
	ts.show_leaf_name = False
	ts.layout_fn = geneTree_layout
	ts.mode = "c"
	spTree.prune(species)
	return spTree

def draw_gene_trees(trees,allProteins,outFolder):
	for code in trees:
		#~ print code
		t = trees[code]
		common = set(t.get_leaf_names()).intersection(allProteins)
		if len(common) > 1:
			anc = t.get_common_ancestor(common)
			if not anc.is_root():
				anc = anc.up
			outfileName = outFolder+"/"+code+".pdf"
			ts = ete3.TreeStyle()
			ts.show_leaf_name = False
			ts.mode = "c"
			ts.layout_fn = geneTree_layout
			#anc.show(tree_style=ts)
			anc.render(outfileName,tree_style=ts,w=200)

def geneTree_layout(node):
	node.img_style["size"] = 0
	if node.is_leaf():
		if node.name in allProteins:
			nameFace = ete3.faces.AttrFace("name", fsize=14, fgcolor="blue")
		else:
			nameFace = ete3.faces.AttrFace("name", fsize=14, fgcolor="grey")
		ete3.faces.add_face_to_node(nameFace,node,column=0,position="aligned")

#From the first round of analysis, it groups together those clusters that are related through vertical evolution
def cluster_vertical_evol(fileName):
	outfile = open(args.dirName+"/t","w")
	for line in open(fileName):
		line = line.strip()
		dades = line.split("\t")
		if dades[-1] == "VE":
			string = dades[0].replace("*",",").replace("|",",")
			print >>outfile,string
		elif "|" in dades[0]:
			c1,c2 = dades[0].split("*")
			if "|" in c1:
				string = c1.replace("|",",")
				print >>outfile,string
			if "|" in c2:
				string = c2.replace("|",",")
				print >>outfile,string
	outfile.close()
	cmd = "python "+clusteringPath+" -i "+args.dirName+"/t -o "+args.dirName+"/cl -m combined -f 0.5 --prefix "+args.dirName+"/"
	run_command(cmd,False)
	groups = []
	for line in open(args.dirName+"/cl"):
		line = line.strip()
		dades = line.split()
		codes = dades[-1].split(",")
		groups.append(codes)
	return groups

def recount_events(fileName,spTree,outfileName):
	for node in spTree.traverse():
		node.add_feature("event",set([]))
	info = {}
	for line in open(fileName):
		line = line.strip()
		if "#" not in line:
			dades = line.split("\t")
			names = dades[0].split("*")
			for n in names:
				if "|" in n:
					codes = n.split("|")
					species = [x.split("_")[0] for x in codes]
					anc = spTree.get_common_ancestor(species)
					anc.event.add("VE")
					anc.event.add(dades[-1])
				else:
					if "_" in n:
						name = n.split("_")[0]
					else:
						name = n
					node = spTree.get_leaves_by_name(name)[0]
					node.event.add(dades[-1])
	outfile = open(outfileName,"w")
	for node in spTree.traverse():
		if len(node.event) == 0:
			pass
		else:
			name = "|".join(node.get_leaf_names())
			print >>outfile,name+"\t"+" ".join(list(node.event))
	outfile.close()

########################################################################
# Gene order distances as calculated by MGR 2.03
########################################################################

def calculate_distances(clusterFile,distances_file,frames):
	clusters = load_clusters_general(clusterFile)
	conversion = load_all_conversion(conversionFolder)
	outfileD = open(distances_file,"w")
	for code in clusters:
		species = set([x.split("_")[0] for x in clusters[code]])
		if len(species) == len(clusters[code]):
			families,size = get_families(clusters[code],conversion)
			if families:
				outfile = open("family.txt","w")
				speConv = {}
				a = 1
				for spe in families:
					speConv[spe] = str(a)
					a += 1
					print >>outfile,">"+speConv[spe]
					print >>outfile," ".join(families[spe])+" $\n"
				outfile.close()
				species = set([x.split("_")[0] for x in families])
				spTree = load_spTree(species_tree)
				spTree.prune(species)
				for leaf in spTree.iter_leaves():
					leaf.name = speConv[leaf.name]
				spTree.write(outfile="spT.txt",format=9)
				cmd = args.mgr_path+" -L -t spT.txt -f family.txt -o result.MGR.txt"
				run_command(cmd,False)
				for line in open("result.MGR.txt"):
					line = line.strip()
					if "Score" in line:
						val = line.split()[4]
				print >>outfileD,code+"\t"+val+"\t"+str(size)
			else:
				print >>outfileD,code+"\tNone"
	outfileD.close()
	
def calculate_distances_pairwise(cluster,outfileD,conversion,outfolder):
	families,size = get_families(cluster,conversion)
	val = "None"
	if families:
		outfile = open(outfolder+"/family.txt","w")
		speConv = {}
		a = 1
		for spe in families:
			speConv[spe] = str(a)
			a += 1
			print >>outfile,">"+speConv[spe]
			print >>outfile," ".join(families[spe])+" $\n"
		outfile.close()
		species = set([x.split("_")[0] for x in families])
		spTree = load_spTree(species_tree)
		spTree.prune(species)
		for leaf in spTree.iter_leaves():
			leaf.name = speConv[leaf.name]
		spTree.write(outfile=outfolder+"/spT.txt",format=9)
		cmd = args.mgr_path+" -L -t "+outfolder+"/spT.txt -f "+outfolder+"/family.txt -o "+outfolder+"/result.MGR.txt"
		run_command(cmd,False)
		for line in open(outfolder+"/result.MGR.txt"):
			line = line.strip()
			if "Score" in line:
				val = line.split()[4]
	return val,size
		
def get_families(cluster,conversion):
	proteins = {}
	info = {}
	cluster1 = {}
	size = 0
	for tag in cluster:
		spe = cluster[tag][0].split("_")[0]
		genes = cluster[tag]
		genes = sorted(genes,key=lambda x:int(x.split("_")[2]))
		codes = []
		for g in genes:
			if g in conversion[spe]:
				if frames[g] == "+":
					codes.append(conversion[spe][g])
				else:
					codes.append("-"+conversion[spe][g])
		info[tag] = codes
		for c in codes:
			noFrameC = c.replace("-","")
			if noFrameC not in proteins:
				proteins[noFrameC] = {}
			if tag not in proteins[noFrameC]:
				proteins[noFrameC][tag] = 0
			proteins[noFrameC][tag] += 1
	to_keep = []
	for p in proteins:
		if len(proteins[p]) == len(cluster):
			duplicates = 0
			for c in proteins[p]:
				if proteins[p][c] != 1:
					duplicates += 1
			if duplicates == 0:
				to_keep.append(p)
	if len(to_keep) <= 3:
		cluster1 = None
	else:
		for name in info:
			cl = []
			for c in info[name]:
				cNoFrame = c.replace("-","")
				if cNoFrame in to_keep:
					if "-" in c:
						cl.append("-"+str(to_keep.index(cNoFrame)+1))
					else:
						cl.append(str(to_keep.index(cNoFrame)+1))
			cluster1[name] = cl
			size = len(cl)
	if cluster1:
		cluster1 = turn_around_if_necessary(cluster1)
	return cluster1,size

def turn_around_if_necessary(clusters):
	species = clusters.keys()
	found = []
	cl1 = ";".join(clusters[species[0]])
	found.append(cl1)
	for s in species[1:]:
		cl = clusters[s]
		clA = ";".join(cl)
		cl2 = []
		for c in cl:
			if "-" not in c:
				n = "-"+c
				cl2.append(n)
			else:
				cl2.append(c.replace("-",""))
		cl2.reverse()
		clB = ";".join(cl2)
		if clA in found:
			pass
		elif clB in found:
			clusters[s] = clB.split(";")
		else:
			found.append(clA)
	return clusters
	
#Get all tree-based statistics from the gene and species trees (can be re-used)
def get_cluster_statistics(clusters,conversion,trees,spTree,outfileNameStats,method):
	#Open the results file
	outfile = open(outfileNameStats,"w")
	print >>outfile,"#ClusterName\tProtCode\tNumProt\tNumSpecies\tMonophylyGT\tMonophylyST\tAncNodeType\tLosses\tAdditionalGenes\tTaxScopeOwnGenes\tTaxScopeAddGenes"
	print "Calculate main stats table..."
	#For each cluster
	a = 0
	for clusterName in clusters:
		if a % 1000 == 0:
			p = str(int(float(a)/float(len(clusters))*100.0))
			print "Clusters processed: "+str(a)+" out of "+str(len(clusters))+" ("+p+"%)"
		a += 1
		#Obtain the proteins that belong to each protein family
		protCounter = count_protein_presence(clusters[clusterName],conversion)
		proteins = protCounter.keys()
		proteins = sorted(proteins,key=lambda x: len(protCounter[x]))
		proteins.reverse()
		#For each protein family present in the trees and with more than one member
		for prot in proteins:
			if prot in trees and len(protCounter[prot]) > 1:
				if method == "general":
					fprots = check_protein_presence(trees[prot],protCounter[prot])
					if len(fprots) > 1:
						mainInfo = analyse_protein_family(protCounter,fprots,trees[prot],fprots,prot)
						print >>outfile,clusterName+"\t"+prot+"\t"+str(mainInfo["prot"])+"\t"+str(mainInfo["species"])+"\t"+str(mainInfo["GTmonophyly"])+"\t"+str(mainInfo["STmonophyly"])+"\t"+mainInfo["EvolNode"]+"\t"+str(mainInfo["Losses"])+"\t"+str(mainInfo["addGenes"])+"\t"+mainInfo["taxScopeOwn"]+"\t"+mainInfo["taxScopeOther"]
				elif method == "detailed":
					fprots,foundA,foundB = check_protein_presence1(trees[prot],protCounter[prot],clusterName)
					if len(fprots) > 1:
						if foundA and foundB:
							mainInfo = analyse_protein_family(protCounter,fprots,trees[prot],allProteins,prot)
							print >>outfile,clusterName+"\t"+prot+"\t"+str(mainInfo["prot"])+"\t"+str(mainInfo["species"])+"\t"+str(mainInfo["GTmonophyly"])+"\t"+str(mainInfo["STmonophyly"])+"\t"+mainInfo["EvolNode"]+"\t"+str(mainInfo["Losses"])+"\t"+str(mainInfo["addGenes"])+"\t"+mainInfo["taxScopeOwn"]+"\t"+mainInfo["taxScopeOther"]				
	outfile.close()	
	p = str(int(float(a)/float(len(clusters))*100.0))
	print "Finished processing clusters: "+str(a)+" out of "+str(len(clusters))+" ("+p+"%)"

def load_distances(outfileName):
	distances = {}
	for line in open(outfileName):
		line = line.strip()
		dades = line.split("\t")
		if "None" not in line:
			distances[dades[0]] = int(dades[1])
	return distances


def obtain_all_calls(outfileNameStats,factorsP,factorsL,tax_thr):
	info = {}
	for line in open(outfileNameStats):
		line = line.strip()
		if "#" not in line:
			dades = line.split("\t")
			if dades[0] not in info:
				info[dades[0]] = {}
			if dades[1] not in info[dades[0]]:
				info[dades[0]][dades[1]] = {}
			for f1 in factorsP:
				for f2 in factorsL:
					name = str(f1)+"-"+str(f2)
					tag = call_tree_event(dades,f1,f2,tax_thr)
					info[dades[0]][dades[1]][name] = tag
	return info

def load_frames(inFile):
	frames = {}
	for line in open(inFile):
		line = line.strip()
		dades = line.split()
		frames[dades[0]] = dades[1]
	return frames

parser = argparse.ArgumentParser(description="Needed to calculate different statistics from the clusters")
parser.add_argument("-i","--cluster_file",dest="clFile",action="store",required=True,help="File where the clusters have been stored")
parser.add_argument("-f","--family",dest="famName",action="store",default=None,help="Name of the family to analyse")
parser.add_argument("-d","--dir_name",dest="dirName",action="store",required=True,help="Folder name where the results will be stored")
parser.add_argument("-c","--conversion",dest="convFile",action="store",required=True,help="Files that contain the conversion between protein codes and the protein family number")
parser.add_argument("-g","--gene_trees",dest="treeFile",action="store",required=True,help="List of newick trees for the protein families. Format needs to be protein_family_number<tab>newick")
parser.add_argument("-t","--taxonomy",dest="taxFile",action="store",required=True,help="Taxonomy table")
parser.add_argument("-s","--species_tree",dest="spTreeFile",action="store",required=True,help="Species tree file")
parser.add_argument("-m","--method",dest="method",action="store",choices=["general","detailed"],default="general",help="Method that will be used to calculate the evolutionary mechanisms of cluster families")
parser.add_argument("--frames_file",dest="framesFile",action="store",default="None",help="File containing the Frames of all genes")
parser.add_argument("--main_stats_file",dest="mainStatsFile",action="store",default=None,help="Pre-calculated tree statistics - Obtained from previous runs")
parser.add_argument("--rank_filter",dest="rank_filter",action="store",choices=["order","class","none"],default="none",help="Taxonomic filter setting")
parser.add_argument("--clustering_path",dest="clusteringPath",action="store",default="./cluster_arrays.py",help="Location of script to fuse clusters together - cluster_arrays.py")
parser.add_argument("--MGR_path",dest="mgr_path",action="store",default="/users/tg/mmarcet/HGT3_proteomes/analyses5/event_calling/MGR",help="Path where the MGR executable can be found")
parser.add_argument("--MGR_results",dest="mgr_results",action="store",default=None,help="File containing pre-computed results of MGR")
parser.add_argument("--MGR_threshold",dest="mgr_threshold",action="store",default=0,type=int,help="Internal gene order conservation below which CE is not trusted; ")
parser.add_argument("--P_threshold_lower",dest="lowPTHR",action="store",default=2.0,type=float,help="Lower multiplier for P threshold")
parser.add_argument("--P_threshold_upper",dest="upPTHR",action="store",default=5.0,type=float,help="Upper multiplier for P threshold")
parser.add_argument("--L_threshold_lower",dest="lowLTHR",action="store",default=1.0,type=float,help="Lower multiplier for L threshold")
parser.add_argument("--L_threshold_upper",dest="upLTHR",action="store",default=2.0,type=float,help="Upper multiplier for L threshold")
parser.add_argument("--max_combination_clusters",dest="maxComb",action="store",default=1000,type=int,help="Maximum number of pairwise cluster combinations that will be analysed by the detailed method")
parser.add_argument("--decision_threshold",dest="decision_threshold",action="store",default=50.0,type=float,help="Threshold needed to calculate when a mechanism is considered valid")
parser.add_argument("-r","--replace",dest="replace",action="store_true",help="Replace previous results")

args = parser.parse_args()

#Load cluster file
clusterFile = args.clFile

clusteringPath = args.clusteringPath

#Write down paths to stable information
#List of the trees calculated for protein families of more than three members. The format is seed<tab>number_protein_family<tab>newick_file
treesFile = args.treeFile
#Conversion files contain a relation between protein name and protein family code as obtained from the mcl clustering
conversionFolder = args.convFile
taxonomy_table = args.taxFile
species_tree = args.spTreeFile

#Create outfolder
create_folder(args.dirName)

#Load taxonomic information
taxonomy,tax_thr = load_taxonomy(taxonomy_table,args.rank_filter)

#Create the different multipliers
factorsP = numpy.arange(args.lowPTHR,args.upPTHR+0.25,0.25)
factorsL = numpy.arange(args.lowLTHR,args.upLTHR+0.25,0.25)

if args.method == "general":
	#Create outfiles
	if args.mainStatsFile:
		outfileNameStats = args.mainStatsFile
		if not os.path.exists(outfileNameStats):
			exit("Main stats file provided does not exists")
	else:
		outfileNameStats = args.dirName+"/main_stats.txt"
	
	if args.mgr_results:
		outfileNameMGR = args.mgr_results
	else:	
		outfileNameMGR = args.dirName+"/MGR.distances.txt"
	outfileName = args.dirName+"/evolclassifier.results.txt"


	#Get basic statistics for each gene in the cluster
	if not os.path.exists(outfileNameStats) or args.replace:
		#Load needed information for statistics
		print "Loading clusters..."
		clusters = load_clusters_general(clusterFile)
		print "Loading conversion..."
		conversion = load_all_conversion(conversionFolder)
		print "Loading trees..."
		trees = load_all_trees(treesFile)
		spTree = load_spTree(species_tree)
		get_cluster_statistics(clusters,conversion,trees,spTree,outfileNameStats,args.method)

	#Calculate number of rearrangements based on MGR software
	if not os.path.exists(outfileNameMGR):
		if not args.framesFile:
			exit("A file containing the frames of the different genes needs to be provided to calculate distances")
		frames = {}
		for line in open(args.framesFile):
			line = line.strip()
			dades = line.split()
			frames[dades[0]] = dades[1]
		calculate_distances(clusterFile,outfileNameMGR,frames)
	
	#Load MGR values into memory
	distances = load_distances(outfileNameMGR)
	
	print "Obtaining calls..."
	info = obtain_all_calls(outfileNameStats,factorsP,factorsL,tax_thr)
	vertical_evol = print_table(info,outfileName,distances)
	
	
elif args.method == "detailed":
	if not args.famName:
		exit("The detailed analysis is only performed for a specific family")
	#Create outfiles
	outfileNameStats1 = args.dirName+"/main_stats1.txt"
	outfileCluster = args.dirName+"/cluster_family.txt"
	outfileName1 = args.dirName+"/results1.txt"
	outfileNameMGR = args.dirName+"/MGR.distances1.txt"

	outfileNameFinalTable = args.dirName+"/evolclassifier.final_table.txt"
	outfileNameTree = args.dirName+"/spTree.pdf"
	outfileRecount = args.dirName+"/evolclassifier.recount.txt"
	outfolderTrees = args.dirName+"/geneTrees/"
	create_folder(outfolderTrees)

	#Get basic statistics for each gene in the cluster
	#Load a single cluster into memory and create all pairwise combinations - They will be loaded as if they were independent clusters
	print "Loading clusters..."
	clusters,proteins2spe = load_one_cluster(clusterFile,args.famName,outfileCluster)

	#If too many combinations are found, the cluster will not be analysed
	if len(clusters) > args.maxComb:
		outfile = open(args.dirName+"/note.txt","w")
		print >>outfile,"Unable to process so many clusters: "+str(len(clusters))
		outfile.close()
		exit()

	print "Loading conversion..."
	conversion,protFam = load_detailed_conversion(conversionFolder,proteins2spe)

	if not args.framesFile:
		exit("A file containing the frames of the different genes needs to be provided")
	
	#Calculate number of re-arrangements
	frames = load_frames(args.framesFile)
	
	outfileDistances = open(outfileNameMGR,"w")
	distances = {}
	for clName in clusters:
		distance,size = calculate_distances_pairwise(clusters[clName],outfileDistances,conversion,args.dirName)
		print >>outfileDistances,clName+"\t"+distance+"\t"+str(size)
		if distance != "None":
			distances[clName] = int(distance)
	outfileDistances.close()
	
	#Load gene and species trees into memory
	print "Loading trees..."
	trees = load_detailed_trees(treesFile,protFam)
	spTree = load_spTree(species_tree)
	
	#Remember all the proteins found in the cluster family
	allProteins = set([])
	for code in proteins2spe:
		for p in proteins2spe[code]:
			allProteins.add(p)
			
	#Call events
	get_cluster_statistics(clusters,conversion,trees,spTree,outfileNameStats1,args.method)
	print "Obtaining calls..."
	info = obtain_all_calls(outfileNameStats1,factorsP,factorsL,tax_thr)
	vertical_evol = print_table(info,outfileName1,distances)
	file2analyse = outfileName1

	a = 2
	#If there are 
	print "VE:",vertical_evol
	if vertical_evol == 0:
		pass
	else:
		rounds = 0
		while vertical_evol != 0 and rounds < 10:
			#Creates groups of clusters that evolved vertically
			groups = cluster_vertical_evol(file2analyse)
			#Loads the grouped clusters and makes combinations between groups
			clusters = combine_clusters(clusterFile,groups,args.famName)
			if len(clusters) != 0:
				#Create new outfiles so that the previous ones are not over-written
				outfileName2 = args.dirName+"/results"+str(a)+".txt"
				outfileNameStats2 = args.dirName+"/main_stats"+str(a)+".txt"
				outfileNameMGR2 = args.dirName+"/MGR.distances"+str(a)+".txt"
				#Keeps track of the file naming
				a += 1
				#Calculates distances
				outfileDistances = open(outfileNameMGR2,"w")
				distances = {}
				for clName in clusters:
					distance,size = calculate_distances_pairwise(clusters[clName],outfileDistances,conversion,args.dirName)
					print >>outfileDistances,clName+"\t"+distance+"\t"+str(size)
					if distance != "None":
						distances[clName] = int(distance)
				outfileDistances.close()
				#Calculates events between clusters
				print "Obtaining calls..."
				get_cluster_statistics(clusters,conversion,trees,spTree,outfileNameStats2,args.method)
				info = obtain_all_calls(outfileNameStats2,factorsP,factorsL,tax_thr)
				vertical_evol = print_table(info,outfileName2,distances)
				file2analyse = outfileName2
				print "VE:",vertical_evol
				rounds += 1
			else:
				vertical_evol = 0	

	format_table(file2analyse,outfileNameFinalTable)
	spTree = draw_tree(proteins2spe.keys(),spTree,outfileNameTree)
	recount_events(file2analyse,spTree,outfileRecount)
	draw_gene_trees(trees,allProteins,outfolderTrees)
