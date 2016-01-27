#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
#  
#  Copyright 2012 Unknown <diogo@arch>
#  
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#  
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#  
#  Author: Diogo N. Silva
#  Version: 0.1
#  Last update: 11/02/14


import argparse
import sys
import subprocess
import pickle

####### #######
#
#	ARGUMENTS
#
###### #######

parser = argparse.ArgumentParser(description="Toolbox to parse and analyze the ouput file of OrthoMCL")

parser.add_argument("-in",dest="infile",nargs="+",help="Input orthoMCL cluster file(s)")
parser.add_argument("-out",dest="outfile",nargs=1,help="Output name file for the table")
parser.add_argument("-threshold",dest="threshold",nargs=2,default=[1,200],help="Provide the thresholds for the number of genes per species in each cluster AND for the minimum number of species per cluster. Default values are 1 and 200, respectively")
parser.add_argument("-f",dest="filt",action="store_const",const=True,help="Using this option creates new filtered group files instead of basic stats tables")
parser.add_argument("-s",dest="stat",choices=["1a","1b","2a","3a","3b"],help="Using this option creates basic stats tables from the information of the groups file. Option 1a can be used to obtain basic statistics comparing different orthoMCL group files (with different inflation values for example); Option 1b creates a csv table with the number of paralog containing groups per species")
parser.add_argument("-g2f",dest="groups2fasta",choices=["single_sp","multi_genes"],help="Use this option if you want to retrieve the protein sequences in the orthoMCL groups into fasta files. Use the 'single' option if you want to save all sequences in all sequences to a single fasta file (in which case you must also specify the name of the output file using the '-out' option), and the 'multi' option if you want to save sequences from different clusters into different files (in which case there is no need to specifiy the output file - the fasta files will have the name of the cluster)")
parser.add_argument("-db",dest="database",nargs=1,help="Provide the BLAST database from which the sequences are to be retrieved")
parser.add_argument("-p",dest="plot",choices=["1","1b","2","2b"],help="Option for plotting: 1 - Simple counting of sequences/species per cluster creating a single graphic with the data for all species; 1b - Same as 1 but creating a graph for each species")
parser.add_argument("-pout",dest="plot_outfile",nargs=1,help="Output name file for plot")
parser.add_argument("--with-outgroup",dest="outgroup",nargs="+",help="Use this option to fetch outgroup taxa from the complete dataset to the BasidioOnly. Please provide the names of the outgroup taxa that you wish to add")
parser.add_argument("--key",dest="key",nargs=1,help="Use this flag if you already have a compiled dictionary with the overlapping between two datasets. Provide this dictionary (in pickle format) using the -in option")
parser.add_argument("-taxa", dest="taxa_subset", nargs="+", help="Specifiy a specific taxa subset to be used in the filters")
parser.add_argument("-taxa-strict", dest="taxa_strict", const=True, action="store_const", help="Use this flag if you want to use the strict mode for the specified taxa only, i.e., filtered data sets will have gene copies below a given threshold only for the focal taxa subset")


arg = parser.parse_args()

####### #######
#
#	LISTS AND TRANSLATION
#
###### #######

Sp_translation = {'Naful': 'Nadsonia_fulvescens', 'Wococ': 'Wolfiporia cocos', 'Waseb': 'Wallemia sebi', 'Bjadu': 'Bjerkandera adusta', 'Crgat': 'Cryptococcus gattii', 'Fuoxy': 'Fusarium_oxysporum', 'Disqu': 'Dichomitus squalens', 'Putri': 'Puccinia triticina', 'Catro': 'Candida_tropicalis', 'Asfum': 'Aspergillus_fumigatus', 'Migyp': 'Microsporum_gypseum', 'Spthe': 'Sporotrichum_thermophile', 'Pisti': 'Pichia_stipitis', 'Gasp': 'Ganoderma sp', 'Ascla': 'Aspergillus_clavatus', 'Selac': 'Serpula lacrymans', 'Phbre': 'Phlebia brevispora', 'Setur': 'Setosphaeria_turcica', 'Chglo': 'Chaetomium_globosum', 'Unree': 'Uncinocarpus_reesii', 'Cosat': 'Cochliobolus_sativus', 'Sacer': 'Saccharomyces_cerevisiae', 'Cocin': 'Coprinopsis cinerea', 'Hecyl': 'Hebeloma cylindrosporum', 'Melar': 'Melampsora laricis-populina', 'Necra': 'Neurospora_crassa', 'Nefis': 'Neosartorya_fischeri', 'Bacom': 'Baudoinia_compniacensis', 'Cohet': 'Cochliobolus_heterostrophus', 'Dosep': 'Dothistroma_septosporum', 'Plost': 'Pleurotus ostreatus', 'Thter': 'Thielavia_terrestris', 'Bocin': 'Botrytis_cinerea', 'Labic': 'Laccaria bicolor', 'Phgig': 'Phlebiopsis gigantea', 'Cohig': 'Colletotrichum_higginsianum', 'Vealb': 'Verticillium_albo-atrum', 'Sacom': 'Saitoella_complicata', 'Trmes': 'Tremella mesenterica', 'Aster': 'Aspergillus_terreus', 'Calus': 'Candida_lusitaniae', 'Trton': 'Trichophyton_tonsurans', 'Netet': 'Neurospora_tetrasperma', 'Trvir': 'Trichoderma_virens', 'Mican': 'Microsporum_canis', 'Tratr': 'Trichoderma_atroviride', 'Audel': 'Auricularia delicate', 'Scjap': 'Schizosaccharomyces_japonicus', 'Fopin': 'Fomitopsis pinicola', 'Sel79': 'Serpula lacrymans S7.9', 'Sccry': 'Schizosaccharomyces_cryophilus', 'Pimem': 'Pichia_membranifaciens', 'Asory': 'Aspergillus_oryzae', 'Painv': 'Paxillus involutus', 'Fomed': 'Fomitiporia mediterranea', 'Yalip': 'Yarrowia_lipolytica', 'Gltra': 'Gloeophyllum trabeum', 'Blder': 'Blastomyces_dermatitidis', 'Lista': 'Lipomyces_starkeyi', 'Agbis': 'Agaricus bisporus', 'Asacu': 'Aspergillus_aculeatus', 'Crneo': 'Cryptococcus neoformans', 'Loelo': 'Lodderomyces_elongisporus', 'Phchr': 'Phanerochaete chrysosporium', 'Copos': 'Coccidioides_posadasii', 'Sccom': 'Schizophyllum commune', 'Poans': 'Podospora_anserina', 'Caalb': 'Candida_albicans', 'Sppas': 'Spathaspora_passalidarum', 'Cesub': 'Ceriporiopsis subvermispora', 'Asnid': 'Aspergillus_nidulans', 'Trver': 'Trametes versicolor', 'Asnig': 'Aspergillus_niger', 'Gedes': 'Geomyces_destructans', 'Ascar': 'Aspergillus_carbonarius', 'Sthir': 'Stereum hirsutum', 'Popla': 'Postia placenta', 'Hicap': 'Histoplasma_capsulatum', 'Fugra': 'Fusarium_graminearum', 'Heann': 'Heterobasidion annosum', 'Pugra': 'Puccinia graminis', 'Mapoa': 'Magnaporthe_poae', 'Acalc': 'Acremonium_alcalophilum', 'Phcar': 'Phanerochaete carnosa', 'Asfla': 'Aspergillus_flavus', 'Pytri': 'Pyrenophora_tritici-repentis', 'Hapol': 'Hansenula_polymorpha', 'Cagui': 'Candida_guilliermondii', 'Scscl': 'Sclerotinia_sclerotiorum', 'Pabra': 'Paracoccidioides_brasiliensis', 'Coput': 'Coniophora puteana', 'Wiano': 'Wickerhamomyces_anomalus', 'Magri': 'Magnaporthe_grisea', 'Gagra': 'Gaeumannomyces_graminis', 'Fuver': 'Fusarium_verticillioides', 'Pustr': 'Punctularia strigosozonata', 'Spros': 'Sporobolomyces roseus', 'Cogra': 'Colletotrichum_graminicola', 'Trrub': 'Trichophyton_rubrum', 'Crpar': 'Cryphonectria_parasitica', 'Vedah': 'Verticillium_dahliae', 'Coimm': 'Coccidioides_immitis', 'Nedis': 'Neurospora_discreta', 'Trequ': 'Trichophyton_equinum', 'Trree': 'Trichoderma_reesei', 'Dasp': 'Dacryopinax sp', 'Nehae': 'Nectria_haematococca', 'Mygra': 'Mycosphaerella_graminicola', 'Usmay': 'Ustilago maydis', 'Scoct': 'Schizosaccharomyces_octosporus', 'Cacas': 'Candida_caseinolytica', 'Phnod': 'Phaeosphaeria_nodorum', 'Patan': 'Pachysolen_tannophilus', 'Myfij': 'Mycosphaerella_fijiensis','Agcyl_EST':'Agrocybe cylindracea','Amare_EST':'Amylostereum areolatum','Ammus_EST':'Amanita muscaria','Artab_EST':'Armillaria tabescens','Ashum_EST':'Asterotremella humicola','Atrol_EST':'Athelia rolfsii','Auaur_EST':'Auricularia auricula-judae','Crlau_EST':'Cryptococcus laurentii','Crque_EST':'Cronartium quercuum f. sp. fusiforme','Crvis_EST':'Cryptococcus vishniacii','Flver_EST':'Flammulina velutipes','Fopal_EST':'Fomitopsis palustris','Galuc_EST':'Ganoderma lucidum','Hecyl_EST':'Hebeloma cylindrosporum','Hesp_EST':'Hericium sp','Hevas_EST':'Hemileia vastatrix','Laqui_EST':'Lactarius quietus','Leedo_EST':'Lentinula edodes','Lesco_EST':'Leucosporidium scottii','Memeddel_EST':'Melampsora medusae f. sp. deltoidis','Memedtre_EST':'Melampsora medusae f. sp. tremuloidis','Meocc_EST':'Melampsora occidentalis','Mivio_EST':'Microbotryum violaceum','Moper_EST':'Moniliophthora perniciosa','Moror_EST':'Moniliophthora roreri','Phchr_EST':'Phanerochaete chrysosporium','Phnam_EST':'Pholiota nameko','Phpac_EST':'Phakopsora pachyrhizi','Phsul_EST':'Phellinidium sulphurascens','Pialb_EST':'Pisolithus albus','Pimic_EST':'Pisolithus microcarpus','Pitin_EST':'Pisolithus tinctorius','Plsp_EST':'Pleurotus sp. "Florida"','Pucor_EST':'Puccinia coronata var. lolii','Pugratri_EST':'Puccinia graminis f. sp. tritici','Pustrtri_EST':'Puccinia striiformis f. sp. tritici','Rhsol_EST':'Rhizoctonia solani','Sulut':'Suillus luteus','Tacam_EST':'Taiwanofungus camphoratus','Thcuc_EST':'Thanatephorus cucumeris','Urapp_EST':'Uromyces appendiculatus','Urvifa_EST':'Uromyces viciae-fabae','Vovol_EST':'Volvariella volvacea'}

Basidiomycota_list = ['Trmes', 'Crneo', 'Phbre', 'Wococ', 'Popla', 'Phchr', 'Spros', 'Phcar', 'Audel', 'Phgig', 'Waseb', 'Sccom', 'Bjadu', 'Usmay', 'Painv', 'Pustr', 'Cesub', 'Crgat', 'Disqu', 'Fomed', 'Sel79', 'Sthir', 'Trver', 'Cocin', 'Hecyl', 'Coput', 'Gltra', 'Melar', 'Labic', 'Putri', 'Heann', 'Dasp', 'Plost', 'Gasp', 'Pugra', 'Fopin', 'Selac', 'Agbis']

Basidio_EST_list = ['64608','103385','41956','47431','5417','39291','29892','5418','136831_','89929','38945','186125','5315','76867','1007625','203904','111141','5353','5278','258770','374510','82102','5272','153609','221103','5306','61267','170000','175648','178870','178872','37468','188765','219187','56615','168172','456999','5384','196114','107832','5264','55588','36659']

Ascomycota_list = ['Naful', 'Fuoxy', 'Catro', 'Asfum', 'Migyp', 'Poans', 'Ascla', 'Setur', 'Chglo', 'Unree', 'Cosat', 'Sacer', 'Nefis', 'Bacom', 'Cohet', 'Dosep', 'Thter', 'Bocin', 'Necra', 'Vealb', 'Sacom', 'Aster', 'Calus', 'Trton', 'Netet', 'Trvir', 'Mican', 'Tratr', 'Scjap', 'Sccry', 'Pimem', 'Spthe', 'Asory', 'Yalip', 'Blder', 'Lista', 'Nedis', 'Loelo', 'Copos', 'Pisti', 'Caalb', 'Sppas', 'Asnid', 'Trequ', 'Asnig', 'Gedes', 'Ascar', 'Hicap', 'Fugra', 'Mapoa', 'Mygra', 'Cohig', 'Asfla', 'Pytri', 'Hapol', 'Cagui', 'Scscl', 'Pabra', 'Wiano', 'Magri', 'Gagra', 'Fuver', 'Cogra', 'Trrub', 'Crpar', 'Vedah', 'Coimm', 'Asacu', 'Trree', 'Nehae', 'Acalc', 'Scoct', 'Cacas', 'Phnod', 'Patan', 'Myfij']

####### #######
#
#	FUNCTIONS
#
###### #######

def loading (current_state,size,prefix,width,suffix):
	""" Function that prints the loading progress of the script """
	percentage = int(((current_state+1)/size)*100)
	complete = int(width*percentage*0.01)
	if percentage == 100:
		sys.stdout.write("\r%s [%s%s] %s%% -- Done!\n" % (prefix,"#"*complete,"."*(width-complete),percentage))
	else:
		sys.stdout.write("\r%s [%s%s] %s%% (%s)" % (prefix,"#"*complete,"."*(width-complete),percentage,suffix))
	sys.stdout.flush()

def load_dic (ortho_infile):
	""" Create dictionary with orthoMCL ID as key and protein clusters as values """
	ortho_dic,ortho_sort = {},[]
	ortho_read = open(ortho_infile)
	for line in ortho_read:
		ortho_dic[line.split(":")[0]] = line.split(":")[1]
		ortho_sort.append(line.split(":")[0])
	return ortho_dic, ortho_sort
	
def species_frequency (sp_list,groups_dictionary,current_group,threshold=1):
	" Function that creates a dictionary with species as keys, and their frequency in the current ortholog group as values "
	sp_freq_dict = dict((species,freq) for species,freq in zip(sp_list,map( lambda species: groups_dictionary[current_group].count(species),sp_list)) if freq >= threshold )
	return sp_freq_dict
	
def frequency_filter (species_frequency_dic,gene_threshold,sp_threshold, taxa_subset, taxa_strict):
	""" Function that filters clusters with gene copy numbers that exceed a given threshold. It also returns a list of the species with gene copy numbers that violate the threshold """
	gene_frequency_flag = 0 # If the value of this flag is changed, then the gene frequency of one or more species is above the threshold
	sp_frequency_flag = 0
	#Basidiomycota_only = [species for species in species_frequency_dic.keys() if species in Basidiomycota_list or Basidio_EST_list]

	if taxa_subset != None: # If the user has specified a taxa set, then the species threshold should apply to that subset only
		focal_species = [sp for sp in taxa_subset if sp in species_frequency_dic]
	else:
		focal_species = species_frequency_dic.keys()

	if len(focal_species) <= int(sp_threshold):
			sp_frequency_flag += 1
	for species,frequency in species_frequency_dic.items():
		if taxa_strict == None and species in taxa_subset:
			return 0, sp_frequency_flag
		else:
			if frequency > int(gene_threshold):
				gene_frequency_flag += 1
	return gene_frequency_flag,sp_frequency_flag
	
def cluster_manipulation (ortho_dic,gene_threshold,sp_threshold, taxa_subset=None,taxa_strict=None):
	""" Function that computes basic statistics on the orthoMCL groups file, such as the number of clusters, number of proteins number of clusters that contain only one gene copy per species """
	try:
		gene_threshold = int(gene_threshold)
		sp_threshold = int(sp_threshold)
	except:
		print ("TypeError: Arguments of the -threshold option must be intergers")
		raise SystemExit
	num_clusters = len(ortho_dic)
	num_sequences = 0 
	below_gene_threshold = 0 # Stores the number of clusters with gene copy numbers per species below a certain threshold
	gene_species_threshold = 0
	# Create a new dictionary to store the filtered groups, if the "f" option is applied. Otherwise, leave it alone.
	filtered_groups = {}
	filtered_order = []
	# Parameters for the "loading" function
	current_group = 0
	# Species list for the "species_frequency" function
	sp_list = list(Sp_translation.keys())
	for ID, cluster in ortho_dic.items():
		loading(current_group,num_clusters,"Parsing orthoMCL groups",50,"Processing group %s" % (ID))
		sequences = cluster.split()
		sp_freq_dic = species_frequency(sp_list,ortho_dic,ID)
		gene_flag, sp_flag = frequency_filter(sp_freq_dic,gene_threshold,sp_threshold, taxa_subset,taxa_strict)
		
		# Counting only the clusters with the number of genes per species below the given threshold
		if gene_flag == 0 and arg.filt == None:
			below_gene_threshold += 1
			
		# Counting (first "if") or selecting (first "elif") only the clusters with the number of genes per species below the given threshold AND the total number of species above a given species threshold
		if gene_flag == 0 and sp_flag == 0 and arg.filt == None:
			gene_species_threshold += 1
			
		if gene_flag == 0 and sp_flag == 0:
			filtered_groups[ID] = cluster
			filtered_order.append(ID)
		num_sequences += len(sequences)
		
		current_group += 1 # Parameter of the "loading" function
		
	return num_clusters, num_sequences,below_gene_threshold,gene_species_threshold,filtered_groups,filtered_order
	
	
def csv_basicStats (outfile,csv_template):
	""" This function grabs a number os statistics (stored in a dictionary variable) to be written on a csv file """
	outfile_handle = open("".join(outfile),"w")
	outfile_handle.write("file;# Clusters;# Sequences; # Clusters (single copy); # Clusters (single copy/ minimum Basidio sp)\n")
	for infile,stats in csv_template.items():
		outfile_handle.write("%s ; %s\n" % (infile,stats))
	else:
		outfile_handle.close()
		
def groups_dump (groups_dic,outfile,order=None):
	outfile_handle = open(outfile,"w")
	current_group = 0 # Parameter of the "loading" function
	if order == None:
		for ID, cluster in groups_dic.items():
			loading(current_group,len(groups_dic),"Dumping groups",50,"Processing group %s" % (ID))
			outfile_handle.write("%s : %s" % (ID,cluster))
			current_group += 1
	else:
		for ID in order:
			loading(current_group,len(order),"Dumping groups",50, "Processing group %s" % (ID))
			outfile_handle.write("%s : %s" % (ID,groups_dic[ID]))
			current_group += 1
	outfile_handle.close()
	
def groups2fasta (group_dic,database,multi_files,fasta_file=None):
	""" Function that converts the clusters contained in an OrthoMCL's group file ito fasta files. The fuction requires a parsed orthoMCL dictionary, and a BLAST database from which sequences are to be retrieved. Is has two option of return: it either creates a ew fasta file per cluster, or is creates a sigle flasta file for the entire groups file (e.g., to make a database). """
	current_group = 0
	database = "".join(database)
	if multi_files == "single_sp":
		for cluster in group_dic.values():
			loading(current_group,len(group_dic),"Retrieving clusters",50,"Processing cluster %s" % (current_group))
			proteins = cluster.split()
			[subprocess.Popen(["blastdbcmd -db %s -dbtype prot -entry '%s' >> %s" % (database,protein,fasta_file)],shell=True).wait() for protein in proteins]
			current_group += 1
	elif multi_files == "multi_genes":
		for ID, cluster in group_dic.items():
			loading(current_group,len(group_dic),"Retrieving clusters",50,"Processing group %s" % (ID))
			proteins = cluster.split()
			[subprocess.Popen(["blastdbcmd -db %s -dbtype prot -entry '%s' >> %s" % (database,protein,ID)],shell=True).wait() for protein in proteins]
			current_group += 1
			
def cross_ortho_dics (Basidio_dic, Complete_dic,key=None):
	""" This function compares both dictionaries and tries to find the Complete_dic group that contains each Basidio_dic group """
	groups_key = {} # This dictionary will be groups_key[group_in_complete] = group_in_basidio
	cur_group = 1
	for complete_group, complete_cluster in Complete_dic.items():
		loading(cur_group,len(Complete_dic),"Crossing dictionaries", 50, "Processing cluster %s" % (cur_group))
		complete_prots = complete_cluster.split()
		for basidio_group, basidio_cluster in Basidio_dic.items():
			basidio_prots = basidio_cluster.split()
			complement = [x for x in basidio_prots if x in complete_prots]
			if complement != []:
				groups_key[complete_group] = basidio_group
		cur_group += 1
	if key == None:
		pickle.dump(groups_key, open("Crossing_dic","wb"))
	print ("There were %s overlapping groups" % (len(groups_key)))
	return groups_key
	
def fetch_outgroups (Basidio_dic, Complete_dic, groups_key, outgroups):
	""" This function fetches the required outgroups from their corresponding groups in the Complete_dic to the Basidio_dic """
	groups_key_inv = dict((values,key) for key, values in groups_key.items())
	for taxa in outgroups:
		for group in Basidio_dic:
			if group in groups_key_inv.keys():
				complete_group = groups_key_inv[group]
				complete_prots = Complete_dic[complete_group].split()
				outgroup_seq = [seq for seq in complete_prots if taxa in seq]
				if outgroup_seq != []:
					outgroup_seq = "".join(outgroup_seq[0])
					basidio_prots = Basidio_dic[group].split()
					basidio_prots.append(outgroup_seq)
					basidio_prots.append("\n")
					Basidio_dic[group] = " ".join(basidio_prots)
	return Basidio_dic
	

####### #######
#
#	MAIN 
#
###### ######	

def main_multiOrtho (ortho_file_list):
	""" Main function for when multiple orthoMCL groups files are provided """
	csv_template = {} # Stores the basic stats for each ortho_file in the following format: ortho_file : BasicStat1;BasicStat2;...
	for ortho_file in ortho_file_list:
		ortho_dic,ortho_order = load_dic (ortho_file)
		num_clusters,num_sequences,gene_filtered,gene_sp_filtered,filtered_groups,filtered_order = cluster_manipulation (ortho_dic,arg.threshold[0],arg.threshold[1])
		csv_template[ortho_file] = "%s ; %s ; %s ; %s" % (num_clusters,num_sequences,gene_filtered,gene_sp_filtered)
	return csv_template

def main_filterGroups (ortho_file,gene_threshold,sp_threshold, taxa_subset=None, taxa_strict=None):
	""" Main function that parses a single orthomcl groups file and filters it accoding to the specified threshold """
	if len(ortho_file) != 1:
		print ("InputFileError: Only one orthoMCL groups file can be specified when using the '-f' option")
		raise SystemExit
	input_ortho_file = "".join(ortho_file)
	ortho_dic,ortho_order = load_dic(input_ortho_file)
	returned_values = cluster_manipulation (ortho_dic,gene_threshold,sp_threshold,taxa_subset, taxa_strict)
	filtered_groups,filtered_order = returned_values[-2],returned_values[-1]
	return filtered_groups,filtered_order
	
def main_groups2fasta (ortho_file=None,ortho_dic=None):
	""" Main function that parses a orthoMCL groups file and retrieves their sequences into fasta file(s) """
	# Checking if the BLAST database was provided
	if arg.database == None:
		print ("MissingFileError: Pleave provide the BLAST database file when using the '-g2f' option")
		raise SystemExit
	# Checking if if the single fasta output file was provided
	if arg.groups2fasta == "single_sp":
		if arg.outfile == None:
			print ("MissingArgumentError: Please provide the name of the single output fasta file")
			raise SystemExit
	# Checking if only one orthoMCL groups file is provided
	if len(ortho_file) != 1:
		print ("InputFileError: Only one orthoMCL groups file can be specified when using the '-f' option\n\tAborting...")
		raise SystemExit
	# If the argument ortho_dic is left None, the program assumes that the groups will be read from a orthoMCL groups FILE that will be give when running the script. If this option is selected, please ensure that the FILE provided is the final version (i.e., has the desired filters).
	if ortho_dic == None:
		input_ortho_file = "".join(ortho_file)
		ortho_dic,ortho_order = load_dic(input_ortho_file)
	if arg.outfile == None:
		groups2fasta(ortho_dic,arg.database,arg.groups2fasta)
	elif arg.outfile != None:
		groups2fasta(ortho_dic,arg.database,arg.groups2fasta,"".join(arg.outfile))
				
				
def paralogPerSp (filtered_dic,outfile):
	""" This function is a spin-off of the cluster_manipulation function. It parses an already filtered orthoMCL dic and produces a table with the number of clusters containing paralogs for each species """
	outfile_handle = open(outfile,"w")
	outfile_handle.write("Species;Number of paralog containing groups\n")
	# Parameters for the "loading" function
	num_clusters = len(filtered_dic)
	current_group = 0
	# Creating the initial dictionary that will store the number of paralogs for each species. The dictionary will be initiated with each species as key, and the intial values as 0
	paralog_frequency_storage = dict((species,initial_value) for species, initial_value in zip(Basidiomycota_list,[0]*len(Basidiomycota_list)))
	for ID, cluster in filtered_dic.items():
		loading(current_group,num_clusters,"Creating table",50, "Processing group %s" % (ID))
		sp_freq_dic = species_frequency(Basidiomycota_list,filtered_dic,ID)
		for species,value in sp_freq_dic.items():
			if value > 1:
				paralog_frequency_storage[species] += 1
		current_group += 1
	for species,value in paralog_frequency_storage.items():
		outfile_handle.write("%s;%s\n" % (Sp_translation[species],value))
		
def main_cluster_stats (input_file,output_file,max_threshold):
	input_handle = open(input_file)
	output_handle = open(output_file,"w")
	output_handle.write("Species threshold; Sequences\n")
	stats_storage = {}
	for i in range (5,max_threshold):
		filtered_groups, ortho_order = main_filterGroups(arg.infile,arg.threshold[0],i)
		stats_storage[str(i)+"_species"] = len(filtered_groups)
	[output_handle.write(str(threshold)+";"+str(stats_storage[threshold])+"\n") for threshold in stats_storage]
	
def main_outgroup (Complete_file, Basidio_file, output_file, outgroups, key=None):
	complete_dic,complete_order = load_dic (Complete_file)
	basidio_dic,basidio_order = load_dic(Basidio_file)
	if key == None:
		cross_dic = cross_ortho_dics (basidio_dic, complete_dic)
		modified_basidio_dic = fetch_outgroups (basidio_dic, complete_dic, cross_dic, outgroups)
		groups_dump (modified_basidio_dic,output_file)
	else:
		cross_dic = key
		modified_basidio_dic = fetch_outgroups (basidio_dic, complete_dic, cross_dic, outgroups)
		groups_dump (modified_basidio_dic,output_file)

def stats_per_sp (Sp_list, ortho_dic, outfile):
	output_handle = open(outfile,"w")
	output_handle.write("Species; # groups\n")
	for species in Sp_list:
		species_count = 0
		for cluster in ortho_dic.values():
			species_count += cluster.count(species)
		output_handle.write("%s; %s\n" % (Sp_translation[species], species_count))
	

####### #######
#
#	EXECUTION 
#
###### ######	


if arg.stat != None:
	if arg.stat == "1a":
		csv_template = main_multiOrtho(arg.infile)
		csv_basicStats (arg.outfile,csv_template)
	elif arg.stat == "1b":
		filtered_groups,ortho_order = main_filterGroups(arg.infile,arg.threshold[0],arg.threshold[1],arg.taxa_subset, arg.taxa_strict)
		paralogPerSp (filtered_groups,arg.outfile[0])
	elif arg.stat == "2a":
		main_cluster_stats ("".join(arg.infile),"".join(arg.outfile[0]),len(Basidiomycota_list))
	elif arg.stat == "3a":
		ortho_dic, ortho_order = load_dic("".join(arg.infile))
		stats_per_sp (Ascomycota_list, ortho_dic, "".join(arg.outfile))
	elif arg.stat == "3b":
		ortho_dic, ortho_order = load_dic("".join(arg.infile))
		stats_per_sp (Basidiomycota_list, ortho_dic, "".join(arg.outfile))
		
if arg.filt != None:
	filtered_groups,ortho_order = main_filterGroups(arg.infile,arg.threshold[0],arg.threshold[1],arg.taxa_subset,arg.taxa_strict)
	groups_dump (filtered_groups,"".join(arg.outfile),ortho_order)
	#orthoMCL_ploter.species_count(filtered_groups)
if arg.groups2fasta != None:
	main_groups2fasta (ortho_file=arg.infile)
if arg.plot != None:
	filtered_groups,ortho_order = main_filterGroups(arg.infile,arg.threshold[0],arg.threshold[1],arg.taxa_subset,arg.taxa_strict)
	# Checking if every thing is okay with the file name for the plot output
	if arg.plot_outfile == None:
		print("MissingArgumentError: Please provide the file name for the plot output (if creating a single plot) or the plot suffix (if creating a plot for each sepecies) using the '-pout' option")
		raise SystemExit
	if arg.plot_outfile != None and len(arg.plot_outfile) != 1:
		print("InputFileError: Only one file name for the plot output can be specified using the '-pout' option")
		raise SystemExit
	# Running the plot options
	if arg.plot == "1":
		orthoMCL_ploter.genes_per_species(filtered_groups,"".join(arg.plot_outfile),singlePlot=True)
	elif arg.plot == "1b":
		orthoMCL_ploter.genes_per_species(filtered_groups,"".join(arg.plot_outfile))
	elif arg.plot == "2":
		orthoMCL_ploter.genes_per_species(filtered_groups,"".join(arg.plot_outfile))
	elif arg.plot == "2b":
		orthoMCL_ploter.genes_per_species(filtered_groups,"".join(arg.plot_outfile),singlePlot=True)

if arg.outgroup != None:
	if len (arg.infile) != 2:
		print ("InputFileError: The -outgroup option required exactly two input files (The first is the complete dataset, and the second the basidioOnly dataset); %s arguments were given" % (len(arg.infile)))
		raise SystemExit
	if arg.outgroup == []:
		print ("MissingInputError: Please provide the at least one outgroup taxa to be added using the -outgroup option (e.g. '-outgroup Bjadu Audel')")
		raise SystemExit
	if arg.key != None:
		# In this case the input is a dictionary in pickle format
		key_dic = pickle.load(open("".join(arg.key),"rb"))
		main_outgroup (arg.infile[0], arg.infile[1],"".join(arg.outfile), arg.outgroup, key=key_dic)
	else:
		main_outgroup (arg.infile[0], arg.infile[1], "".join(arg.outfile), arg.outgroup)
