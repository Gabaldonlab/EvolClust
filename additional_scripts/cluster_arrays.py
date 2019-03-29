#!/usr/bin/env python
import os
import sys
import time
import cPickle
import logging as log
from ete3 import Tree
from string import strip
from itertools import combinations
from argparse import ArgumentParser
from hcluster import linkage, to_tree

sys.setrecursionlimit(10000)

__DESCRIPTION__ = """
Performs a UPGMA clustering analysis on a selected list of vectors and produces
a reduced list in which overlaping vectors are removed.
Current version is based on jhuerta contribution.
Salvador Capella-Gutierrez (twitter/@sj_capella)
"""

def euc_dist(a, b):
  return  1 - (float(len(a & b)) / max(len(b), len(a)))
  #return (float(len(a.symmetric_difference(b)))) / max(len(a), len(b))

def dist_matrix(X, metric):
  matrix = []
  t1 = time.time()
  period = 10000.0
  expected = (len(X)**2)/2. - len(X)/2.

  print "array size", len(X), "matrix size", expected
  for i, (x,y) in enumerate(combinations(X, 2)):
    matrix.append(metric(x, y))

    i += 1
    if (i % period) == 0:
      t2 = time.time()
      rem_time = int((((expected - i) / period) * (t2 - t1))/60.)
      t1 = t2
      print "\r", i, "/", expected, "remaining time: %d min" % rem_time,

  print "\nmatrix length", len(matrix)
  return matrix

def cluster(items, cache_clustering_file = None, dist_fn = euc_dist, \
    prefix_output = None):

  if not cache_clustering_file:
    print "Generating distance matrix..."
    sys.stdout.flush()
    Y = dist_matrix(items, dist_fn)

    print "Linkage clustering..."
    sys.stdout.flush()
    Z = linkage(Y, "single") # average, complete = max, single = min ?

    print "Dumping clustering information into cache file";
    sys.stdout.flush()
    cPickle.dump([Y,Z], open(prefix_output + "clustering_dump.pkl", "w"))

  else:
    print "Loading clustering cache from '%s'" % cache_clustering_file.name
    Y, Z = cPickle.load(cache_clustering_file)

  print "Converting into ETE tree..."; sys.stdout.flush()
  T = to_tree(Z)

  root = Tree()
  root.dist = 0
  root.name = "root"
  item2node = {T: root}

  to_visit = [T]
  while to_visit:
    node = to_visit.pop()
    cl_dist = node.dist /2.0
    for ch_node in [node.left, node.right]:
      if ch_node:
        ch = Tree()
        #try:
        #  ch.add_features(content = str(items[ch_node.id]))
        #except IndexError:
        #  pass
        ch.dist = cl_dist
        ch.name = str(ch_node.id)
        item2node[node].add_child(ch)
        item2node[ch_node] = ch
        to_visit.append(ch_node)

  return root

def get_cluster_representatives(root, items, min_compact, merge_type):
  print "Merging clusters. Min compactness = ", min_compact,
  print "Merge type =", merge_type

  is_leaf = lambda x:  1 - sum([ch.dist for ch in x.children]) >= min_compact \
    or x.is_leaf()

  item_groups = [node for node in root.iter_leaves(is_leaf_fn = is_leaf)]
  if set([n.up for n in item_groups]) == set(root):
    item_groups = [root]

  representatives = []
  for n in item_groups:
    cluster_items = [items[int(item_idx)] for item_idx in n.get_leaf_names()]
    if merge_type == "largest":
      cluster_items.sort(lambda x, y: cmp(len(x), len(y)), reverse = True)
      representatives.append(cluster_items[0])
    elif merge_type == "smallest":
      cluster_items.sort(lambda x, y: cmp(len(x), len(y)), reverse = False)
      representatives.append(cluster_items[0])
    elif merge_type == "combined":
      representatives.append(set([i for memb in cluster_items for i in memb]))
  return representatives

log.basicConfig(level = log.INFO, format = "%(levelname)s - %(message)s")

if __name__ == "__main__":
  parser = ArgumentParser(description = __DESCRIPTION__)

  parser.add_argument("-i", "--input", dest = "input", type = open, required = \
    True, help = "Source file containing a list of element arrays")

  parser.add_argument("--min_size", dest = "min_size", type = int, default = 0,
    help = "Set the minimum cluster size for the dumped clusters")

  parser.add_argument("-c", "--column", dest = "column", type = int, default = \
    0, help = "Identify the column containing the array of each line")

  parser.add_argument("--delimiter", dest = "delimiter", default = "\t", type =\
    str, help = "How fields are delimited within lines")

  parser.add_argument("--array_delimiter", dest = "array_delimiter", type = str,
    default = ",", help = "How elements are separated within arrays")

  parser.add_argument("-o", "--outfile", dest = "outfile", type = str, required\
    = True, help = "Where to write the resulting clusters")

  parser.add_argument("-m", "--merging", dest = "merge_type", default = \
    "largest",choices = ["combined", "largest", "smallest"], help = "Clusters "
    + "are represented by the sum, largest or smallest element in the group")

  parser.add_argument("-f", "--compact_factor", dest = "compact_factor", type =\
    float, required = True, help = "Minimun compactness to collapse clusters")

  parser.add_argument("--cache_cluster", dest = "cache_cluster_file", type = \
    open, help = "A clustering file from previous executions")

  parser.add_argument("--prefix", dest = "output_prefix", type = str, default \
    = "", help = "A prefix for output files")

  args = parser.parse_args()

  in_items = []
  for line in args.input:
    if line.startswith("#"):
      continue
    f = map(strip, line.split(args.delimiter))
    in_items.append(set(map(str, f[args.column].split(args.array_delimiter))))
  print ("\rTotal: %d processed arrays") % (len(in_items))

  items = []
  slide = int(len(in_items)/100.) if len(in_items) > 100 else 1
  sorted_items = sorted([(len(array), array) for array in in_items])
  for pos in range(len(sorted_items)):
    size, array = sorted_items[pos]
    overlap = False
    for record in sorted_items[pos+1:]:
      if record[0] == size:
        continue
      if array - record[1] == set():
        overlap = True
        break
    if not overlap:
      items.append(array)
    if pos and (pos % slide) == 0:
      ratio = float(pos)/slide
      print ("\r%d processed arrays (%.2f%%)") % (pos, ratio),
      sys.stdout.flush()
  print ("\rAfter removing redundancy: %d arrays") % (len(items))
  sys.stdout.flush()

  ## Determine output file prefix, if given
  if args.output_prefix and not args.output_prefix[-1] in [".", "_", "-"]:
    args.output_prefix += "."

  if len(items) == 1:
    oFile = open(args.outfile, "w")
    cluster = ",".join(sorted(items[0]))
    print >> oFile, ("cluster_1\tsize_%d\t%s") % (len(items[0]), cluster)
    oFile.close()
    exit()

  tree = cluster(items, args.cache_cluster_file, prefix_output = \
    args.output_prefix)

  ofile = ("%sdump_tree.nw") % (args.output_prefix)
  print "Dumping clustering tree into '%s'" % (ofile)
  tree.write(outfile = ofile)

  repre = get_cluster_representatives(tree, items, args.compact_factor, \
    args.merge_type)
  print "Detected", len(repre), "clusters"
  sys.stdout.flush()

  ## New Code
  ## Collapse smaller clusters into bigger ones
  clusters = {}
  for rp in repre:
    clusters.setdefault(len(clusters), set(rp))

  order = sorted([(len(clusters[cl]), cl) for cl in clusters], reverse = True)

  n = 1
  discarded = []
  for current in order:
    members = clusters[current[1]]
    for oth_clust in [oth for oth in order[n:] if not oth[1] in discarded]:
      ## In case bigger cluster contains the smaller one, just discard
      ## the small one
      if clusters[oth_clust[1]] - members == set():
        discarded.append(oth_clust[1])
    n += 1

  cl = 0
  oFile = open(args.outfile, "w")
  for entry in [record[1] for record in order if not record[1] in discarded]:
    clusterSize = len(clusters[entry])
    if clusterSize >= args.min_size:
      cl += 1
      cluster = ",".join(sorted(map(str, clusters[entry])))
      print >> oFile, ("cluster_%s\tsize_%d\t%s") % (cl, clusterSize, cluster)
  oFile.close()

  print "Final dumped clusters", cl
