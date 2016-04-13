#!/usr/bin/python
#treeshrink_for_itol.py
'''
A small utility for summarizing protein presence/absence in summary phylogenetic trees.
Uses ETE for tree analysis and NCBI Taxonomy handling.
ETE will already do many of the things iTOL does...
but iTOL is friendlier to programming laymen.

Remaining changes:
Some downloaded NOGs don't appear to download correctly

'''
import glob, gzip, os, re, urllib2, sys
from ete2 import Tree, NCBITaxa

##Options

ncbi = NCBITaxa()	#Note that this stores in the home dir
#ncbi.update_taxonomy_database()

##Classes

##Functions

def parse_tree_file(file_name):
	#Returns parent taxids of nodes at a specific level
	#For now this is the Order level
	
	with open(file_name) as tree_file:
		treestring = tree_file.readline()
	
	leaf_taxids = []	#All taxids from leaves of the input tree
	parent_taxids = {}	#All highest-level parents of the leaf taxids
						#A dict with unique taxids as keys
						#And genome counts as values 
						
	fulltree = Tree(treestring)
	for leaf in fulltree:
		leaf_taxid = ((leaf.name).split("."))[0]
		leaf_taxids.append(leaf_taxid)
	#print(leaf_taxids)
	
	print("Retrieving parent nodes for each leaf...")
	for taxid in leaf_taxids:	
		#We're assuming taxids are all in the database here to save time
		#But the following can be used too
		#if taxid in all_nodes.keys():
		parent_taxid = find_parent(taxid, ncbi)
		
		if parent_taxid in parent_taxids.keys():
			parent_taxids[parent_taxid] = parent_taxids[parent_taxid] +1
		else:
			parent_taxids[parent_taxid] = 1
			
	total_count = len(leaf_taxids)
	
	output_taxidlist(parent_taxids)
	output_ann_file(parent_taxids, total_count)
	#print(parent_taxids)
			
def find_parent(taxid, db):
	lineage = ncbi.get_lineage(taxid)
	if len(lineage) > 1:	#Some lineages resolve to root only. Error?
		parent_taxid = lineage[4]	#Not the direct parent, more of the ancestor
	else:
		parent_taxid = lineage[0]
	#print(lineage)
	#names = ncbi.get_taxid_translator(lineage)
	#print [names[taxid] for taxid in lineage]
	return parent_taxid
	
def output_taxidlist(parent_taxids):
	out_filename = "output_taxidlist.txt"
	with open(out_filename, "w+b") as taxid_file:
		for taxid in parent_taxids.keys():
			taxid_file.write("%s\n" % taxid)
	print("Wrote new taxid list to %s." % out_filename)
	
def output_ann_file(parent_taxids, total):
	out_filename = "output_annotations.txt"
	
	with open(out_filename, "w+b") as ann_file:
		ann_file.write("DATASET_MULTIBAR\n" +
						"SEPARATOR TAB\n" +
						"DATASET_LABEL\tGenomes With OG Member\n" +
						"COLOR\t#006400\n" +
						"FIELD_COLORS\t#006400\t#00ff00\n" +
						"FIELD_LABELS\tAbsolute\tRelative\n" +
						"ALIGN_FIELDS\t1\n" +
						"DATA\n")
		for taxid in parent_taxids:
			count = parent_taxids[taxid] 
			relative_count = (count/float(total))*100	
				#Actually a percentage for display purposes
			ann_file.write("%s\t%s\t%s\n" % (taxid, count,
											relative_count))
	print("Wrote new annotation file to %s." % out_filename)
	
def getNOGTree(nog):
	base_url = "http://eggnogapi.embl.de/nog_data/text/tree/"
	nog = str(nog.rstrip())
	nog_tree_url = base_url + nog
	tree_file_name = nog + ".txt"
	
	if os.path.isfile(tree_file_name): 
		print("Found a matching tree file on disk: %s" % tree_file_name)
	else:
		print("Downloading tree for %s from eggNOG (%s)." % (nog, nog_tree_url))
		try:
			response = urllib2.urlopen(nog_tree_url)
			tree_file = open(os.path.basename(tree_file_name), "w") #Start local file
			while 1:
				for line in response:
					tree_file.write(line)
				print("\nDownload complete.")
				break
		except urllib2.HTTPError:
			print("Could not find that NOG on eggNOG.")
			return "NA" 

	return tree_file_name
	
##Main
def main():
	choice = raw_input("Please choose an option:\n" +
						"1\tSpecify a NOG\n" +
						"2\tProvide a tree file\n")
						
	if choice.rstrip() == str(1):
		nog = raw_input("Name of the NOG?\n")
		tree_file_name = getNOGTree(nog)
		if tree_file_name == "NA":
			sys.exit("Please try again.")
		
	elif choice.rstrip() == str(2):
		have_tree_file = False
		
		while not have_tree_file:
			tree_file_name = raw_input("Name of the input file?\n")
			if not os.path.isfile(tree_file_name):
				print("Couldn't find that file. Try again.")
			else:
				have_tree_file = True
				
	else:
		print("Did not recognize that choice. Please try again.")
		
	parse_tree_file(tree_file_name)
	
if __name__ == "__main__":
	sys.exit(main())
