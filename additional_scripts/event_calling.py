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
import os
import subprocess as sp

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

#Load clusters into memory
def load_clusters(fileName):
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

#Correlate each protein to the cluster/s it belongs to
def load_prot2clusters(fileName):
	prot2clusters = {}
	for line in open(fileName):
		line = line.strip()
		if "##" in line:
			pass
		elif "#" in line:
			name = line.split()[1]
		else:
			codes = line.split()[1].split(";")
			for c in codes:
				if c not in prot2clusters:
					prot2clusters[c] = set([])
				prot2clusters[c].add(name)
	return prot2clusters

#Load trees into memory and correlate them to their family number
def speName(node):
	return node.split("_")[0]

def whole_name(node):
	return node

def load_trees(fileName,conversion):
	trees = {}
	for line in open(fileName):
		line = line.strip()
		dades = line.split()
		t = ete3.PhyloTree(dades[1],sp_naming_function=speName)
		code = dades[0]
		num = conversion[code.split("_")[0]][code]
		trees[num] = t
	return trees

#Loads the species tree
def load_spTree(fileName):
	spTree = ete3.PhyloTree(fileName,sp_naming_function=whole_name)
	for leaf in spTree.iter_leaves():
		if len(leaf.name) != 5:
			leaf.delete()
	return spTree

#Loads conversion into memory
def load_conversion(conversionFolder):
	conversion = {}
	for fileName in glob.glob(conversionFolder+"/*txt"):
		spe = fileName.split("/")[-1].split(".")[0]
		conversion[spe] = {}
		for line in open(fileName):
			line = line.strip()
			dades = line.split()
			conversion[spe][dades[0]] = dades[1]
	return conversion

#Load taxonomic information into memory
def load_taxonomy(fileName):
	taxonomy = {}
	classes = set([])
	for line in open(fileName):
		line = line.strip()
		dades = line.split("\t")
		taxonomy[dades[1]] = dades[3:]
		classes.add(dades[6])
	return taxonomy,classes

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
	for p in proteins:
		if p in t.get_leaf_names():
			proteins2.append(p)
	return proteins2

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
def check_taxonomic_scope(additionalGenes,taxonomy,classes):
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
		if b in classes:
			found = True
	return found,base_tax
	
#Obtain statistics of a protein family within a cluster		
def analyse_protein_family(protCounter,fprots,t,treeDist):
	mainInfo = {}
	mainInfo["prot"] = len(protCounter[prot])
	mainInfo["protInTree"] = len(fprots)
	#Check whether gene tree is monophyletic
	monophyly,anc = check_monophyly(t,fprots)
	mainInfo["GTmonophyly"] = monophyly
	#Check whether the common ancestor of the species involved is a speciation node or a duplication node
	species1 = set([x for x in anc.get_children()[0].get_species()])
	species2 = set([x for x in anc.get_children()[0].get_species()])
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
	additionalGenes = [x for x in allGenes if x not in fprots]
	mainInfo["addGenes"] = len(additionalGenes)
	if len(additionalGenes) != 0:
		#Check whether they surpass the taxonomic range
		found,scope = check_taxonomic_scope(additionalGenes,taxonomy,classes)
		mainInfo["RangeTHR"] = found
		mainInfo["taxScope"] = scope[-1]
	else:
		mainInfo["RangeTHR"] = "-"
		mainInfo["taxScope"] = "-"
	#Calculate the treeKO distance between the subtree and the species tree
	if treeDist:
		mainInfo["TK"] = calculate_tree_distance(anc,spTree)
	else:
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
#V = Vertical evolution, VL = Vertical evolution plus loss, IF = convergent evolution / independent formation

def call_tree_event(dades,f1,f2):
	monophylyGT = dades[4]
	monophylyST = dades[5]
	try:
		numLosses = int(dades[7])
	except:
		numLosses = 0
	numAddGenes = int(dades[8])
	taxThr = dades[9]
	species = int(dades[3])
	proteins = int(dades[2])
	conclusion = "FAILED"
	if monophylyGT == "True":
		if monophylyST == "True":
			conclusion = "V"
		elif monophylyST == "False":
			lossTHR = species*f1
			if numLosses > lossTHR:
				conclusion = "HGT"
			else:
				conclusion = "VL"
	else:
		addTHR = proteins*f2
		if addTHR < numAddGenes or taxThr == "False":
			conclusion = "IF"
		else:
			if monophylyST == "True":
				conclusion = "V"
			elif monophylyST == "False":
				lossTHR = species*f1
				if numLosses > lossTHR:
					conclusion = "HGT"
				else:
					conclusion = "VL"
	return conclusion

#Print tables and calculate final calling
def print_table(info,outfileName,outfileNameFinal,names):
	outfile = open(outfileName,"w")
	finalCall = {}
	print >>outfile,"#ClName\tProtFamily\t"+"\t".join(names)+"\tFinalCall"
	for cl in info:
		counterFam = {}
		totalFam = 0
		for prot in info[cl]:
			string = cl+"\t"+prot
			counter = {}
			for n in names:
				t = info[cl][prot][n]
				string += "\t"+t
				if t not in counter:
					counter[t] = 0
				counter[t] += 1
			decision = ["Undecided",0.00]
			for t in counter:
				p = float(counter[t]) / float(len(names)) * 100.00
				if p >= 50.0:
					decision = [t,p]
			if decision[0] not in counterFam:
				counterFam[decision[0]] = 0
			counterFam[decision[0]] += 1
			totalFam += 1
			print >>outfile,string+"\t"+decision[0]+"\t"+str(decision[1])
		finalCall[cl] = ["Undecided",0.00]
		for code in counterFam:
			p = float(counterFam[code]) / float(totalFam)*100.00
			if p >= 50.0:
				finalCall[cl] = [code,p]
	outfile.close()
	outfile = open(outfileNameFinal,"w")
	for code in finalCall:
		print >>outfile,code+"\t"+finalCall[code][0]+"\t"+str(finalCall[code][1])
	outfile.close()


					
parser = argparse.ArgumentParser(description="Needed to calculate different statistics from the clusters")
parser.add_argument("-i","--cluster_file",dest="clFile",action="store",required=True,help="File where the clusters have been stored")
parser.add_argument("-d","--dir_name",dest="dirName",action="store",required=True,help="Folder name where the results of --assess_tree_topology will be stored")
parser.add_argument("-c","--conversion",dest="convFile",action="store",required=True,help="Files that contain the conversion between protein codes and the protein family number")
parser.add_argument("-g","--gene_trees",dest="treeFile",action="store",required=True,help="List of newick trees for the protein families. Format needs to be seed<\t>newick")
parser.add_argument("-t","--taxonomy",dest="taxFile",action="store",required=True,help="Taxonomy table")
parser.add_argument("-s","--species_tree",dest="spTreeFile",action="store",required=True,help="Species tree file")
parser.add_argument("--tree_distances",dest="td",action="store_true",help="In addition obtains distance between gene trees and species trees")
parser.add_argument("-r","--replace",dest="replace",action="store_true",help="Replace previous results")

args = parser.parse_args()

#Load cluster file
clusterFile = args.clFile

#Write down paths to stable information
#List of the trees calculated for protein families of more than three members. The format is seed<tab>number_protein_family<tab>newick_file
treesFile = args.treeFile
#Conversion files contain a relation between protein name and protein family code as obtained from the mcl clustering
conversionFolder = args.convFile
taxonomy_table = args.taxFile
species_tree = args.spTreeFile

#Create outfiles
create_folder(args.dirName)
outfileNameStats = args.dirName+"/main_stats.txt"
outfileName = args.dirName+"/event_calling.stats.txt"
outfileNameFinal = args.dirName+"/event_calling.final.txt"

#Get basic statistics for each gene in the cluster
if not os.path.exists(outfileNameStats) or args.replace:
	#Load data into memory
	print "Loading clusters..."
	clusters = load_clusters(clusterFile)
	prot2clusters = load_prot2clusters(clusterFile)
	print "Loading conversion..."
	conversion = load_conversion(conversionFolder)
	taxonomy,classes = load_taxonomy(taxonomy_table)
	print "Loading trees..."
	trees = load_trees(treesFile,conversion)
	spTree = load_spTree(species_tree)
	#Open the results file
	outfile = open(outfileNameStats,"w")
	print >>outfile,"#ClusterName\tProtCode\tNumProt\tNumSpecies\tMonophylyGT\tMonophylyST\tAncNodeType\tLosses\tAdditionalGenes\tTaxTHR\tTaxScopeAddGenes\tTreeKo"
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
				fprots = check_protein_presence(trees[prot],protCounter[prot])
				if len(fprots) > 1:
					mainInfo = analyse_protein_family(protCounter,fprots,trees[prot],args.td)
					print >>outfile,clusterName+"\t"+prot+"\t"+str(mainInfo["prot"])+"\t"+str(mainInfo["species"])+"\t"+str(mainInfo["GTmonophyly"])+"\t"+str(mainInfo["STmonophyly"])+"\t"+mainInfo["EvolNode"]+"\t"+str(mainInfo["Losses"])+"\t"+str(mainInfo["addGenes"])+"\t"+str(mainInfo["RangeTHR"])+"\t"+mainInfo["taxScope"]+"\t"+str(mainInfo["TK"])
	outfile.close()
	if a % 1000 == 0:
		p = str(int(float(a)/float(len(clusters))*100.0))
		print "Clusters processed: "+str(a)+" out of "+str(len(clusters))+" ("+p+"%)"			

print "Obtaining calls..."
#We're going to run the analysis for different thresholds
factors = [0.5,0.75,1.0,1.25,1.5,1.75,2.0,2.25,2.5,2.75,3.0]
allNames = []
info = {}
for line in open(outfileNameStats):
	line = line.strip()
	if "#" not in line:
		dades = line.split("\t")
		if dades[0] not in info:
			info[dades[0]] = {}
		if dades[1] not in info[dades[0]]:
			info[dades[0]][dades[1]] = {}
		for f1 in factors:
			for f2 in factors:
				name = str(f1)+"-"+str(f2)
				if name not in allNames:
					allNames.append(name)
				tag = call_tree_event(dades,f1,f2)				
				info[dades[0]][dades[1]][name] = tag
print_table(info,outfileName,outfileNameFinal,allNames)
