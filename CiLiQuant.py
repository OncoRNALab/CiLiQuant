import pandas as pd
import argparse
import os
import sys

def readjunctions(filename, fw_col, overlap_col, fw_filter):
	""" Read the (normal) junction file
	Start and stop position in some mappers reflects the start of the leftmost read and the stop of the rightmost read that contain the junction instead of actual splicing sites (see ifelse)
	"""
	junctions = pd.read_csv(filename,sep="\t",header=None,usecols=list(range(0,6,1)), names=["chromosome","start","stop","name","score","strand"],dtype={"chromosome":str,"start":float,"stop":float,"name":str,"score":str,"strand":str})
	junction_reads = pd.read_csv(filename,sep="\t",header=None,usecols=[int(fw_col)-1], names=["reads"],dtype={"reads":float}) #add read column (can be one of the first 6 columns, e.g. score column in TopHat)
	junctions = pd.concat([junctions, junction_reads], axis=1) #add read column to junction bed file
	if overlap_col!="":
		junctions_overlap = pd.read_csv(filename,sep="\t",header=None,usecols=[int(overlap_col)-1], names=["overlap"],dtype={"overlap":str})
		junctions_overlap[['overlapleft','overlapright']]=junctions_overlap['overlap'].str.split(',', expand=True)
		junctions['start'] = junctions['start']+junctions_overlap['overlapleft'].astype('int')-2 #add left overlap to get exact start site of junction, then subtract 2 to look for overlap within 2nt        
		junctions['stop'] = junctions['stop']-junctions_overlap['overlapright'].astype('int')+2 #subtract right overlap to get exact stop site of junction, then add 2 to look for overlap within 2nt   
	else: #if the exact start and stop site were already in df: just modify start and stop to look at overlap within 2nt of exact start/stop site (accounts for discrepancies due to annotation or coordinate system)
		junctions['start'] = junctions['start']-2
		junctions['stop'] = junctions['stop']+2
	junctions = junctions[(junctions['start']>0) & (junctions['stop']>0)]
	junctions['start']=junctions['start'].astype('int') #start of junction (or of leftmost read that contains the junction)
	junctions['stop']=junctions['stop'].astype('int') #stop of junction (or of rightmost read that contains the junction)
	junctions['name']=junctions['chromosome'].astype('str') + '_' + junctions['start'].astype('str') + '_' + junctions['stop'].astype('str') + "_" + junctions['strand'].astype('str')
	junctions['reads']=junctions['reads'].astype('int') #number of reads that contain the junction
	#filter out forward-splice junctions that have less than the requested reads
	junctions_mask=junctions['reads'] >= int(fw_filter)
	filtered_junctions = junctions[junctions_mask]
	return(filtered_junctions)

def readbacksplice(filename, bs_col, bs_filter):
	""" Read the backsplice junction file """
	backsplice_junctions = pd.read_csv(filename,sep="\t",header=None,usecols=list(range(0,6,1)),names=["chromosome","start","stop","name","score","strand"],dtype={"chromosome":str,"start":int,"stop":int,"name":str,"score":str,"strand":str})
	backsplice_reads = pd.read_csv(filename,sep="\t",header=None,usecols=[int(bs_col)-1], names=["reads"],dtype={"reads":int}) #add read column (can be one of the first 6 columns)
	backsplice_junctions = pd.concat([backsplice_junctions, backsplice_reads], axis=1) #add read column to backsplice junction bed file
	backsplice_junctions = backsplice_junctions.sort_values(["chromosome","start","stop","strand"])
	backsplice_junctions['circ_id'] = backsplice_junctions['chromosome'].astype('str') + '_' + backsplice_junctions['start'].astype('str') + '_' + backsplice_junctions['stop'].astype('str') + "_" + backsplice_junctions['strand'].astype('str')
	#filter out backsplice junctions that have less than the requested reads
	backsplice_mask = backsplice_junctions['reads'] >= int(bs_filter)
	filtered_bsjunctions = backsplice_junctions[backsplice_mask]
	return(filtered_bsjunctions)

def readexons(filename):
	""" Read the exon file """
	exons = pd.read_csv(filename,sep="\t",header=None,usecols=list(range(0,6,1)), names=['chromosome','start','stop','name','score','strand'],dtype={'chromosome':str,'start':int,'stop':int,'name':str,'score':str,'strand':str})
	exons['exon_id'] = exons['chromosome'].astype('str') + '_' + exons['start'].astype('str') + '_' + exons['stop'].astype('str') + '_' + exons['strand'].astype('str')
	return (exons)

def label_rflank(row):
	""" Substitute the 0 counts by NA if the only right flanking junction(s) is/are ambiguous """
	if (row['rfl_reads'] == 0) & (row['rfl_junctions_ambi'] > 0):
		return "NA"
	else:
		return row['rfl_reads']

def label_lflank(row):
	""" Substitute the 0 counts by NA if the only left flanking junction(s) is/are ambiguous """
	if (row['lfl_reads'] == 0) & (row['lfl_junctions_ambi'] > 0):
		return "NA"
	else:
		return row['lfl_reads']	

def genelevelfraction(df):
	"""
	Calculate for every gene the circRNA fraction based on average lin and circ junction reads in gene
	CI with adjusted Wald method (Agresti and Coull, 1998; Brown et al., 2001)
	"""
	
	# Divide sum of circular and linear reads by the respective nr of junctions (if no junctions present: set to 0)
	df['lin_reads_av'] = df['linear_reads'].divide(df['linear_junctions']).where( (df['linear_junctions']>0), 0)
	df['circ_reads_av'] = df['circ_reads'].divide(df['circ_junctions']).where( (df['circ_junctions']>0), 0)
	
	# Determine the abundance of circ (backsplice) junction reads (relative to the total nr of junction reads that are not ambiguous)
	df['circ_fraction'] = df['circ_reads_av']/(df['circ_reads_av'] + df['lin_reads_av'])
	
	## Wald method
	#df['ci_lower_Wald'] = df['circ_proportion'] - 1.96*( (df['circ_proportion'] * (1-df['circ_proportion'])) / (df['bs_reads']+df['linfl_reads_av']) )**(0.5)
	#df['ci_upper_Wald'] = df['circ_proportion'] + 1.96*( (df['circ_proportion'] * (1-df['circ_proportion'])) / (df['bs_reads']+df['linfl_reads_av']) )**(0.5)
	# Adjusted Wald method (Agresti and Coull, 1998; Brown et al., 2001)
	df['p_AC'] = (df['circ_reads_av'] + 1.96**2/2) / (df['circ_reads_av']+df['lin_reads_av'] + 1.96**2)
	df['ci_lower_AC'] = df['p_AC'] - 1.96*( (df['p_AC'] * (1-df['p_AC'])) / (df['circ_reads_av']+df['lin_reads_av']+1.96**2) )**(0.5)
	df['ci_upper_AC'] = df['p_AC'] + 1.96*( (df['p_AC'] * (1-df['p_AC'])) / (df['circ_reads_av']+df['lin_reads_av']+1.96**2) )**(0.5)
	return(df.round({'lin_reads_av':3,'circ_reads_av':3,'circ_fraction':5,'p_AC':5,'ci_lower_AC':5,'ci_upper_AC':5}))

def circlevelfraction(df):
	"""
	Calculate for every backsplice junction the circRNA fraction based on flanking linear only junction reads
	CI with adjusted Wald method (Agresti and Coull, 1998; Brown et al., 2001)
	"""
	
	# Take sum of lfl and rfl reads if they are both present, otherwise use nr of reads from lfl or rfl 
	df['linfl_reads'] = df[['rfl_reads', 'lfl_reads']].sum(axis=1).where((df['lfl_junctions']>0) & (df['rfl_junctions']>0), df['lfl_reads']) #take sum if both are present, otherwise, take lfl reads
	df['linfl_reads'] = df[['rfl_reads']].where((df['lfl_junctions']==0) & (df['rfl_junctions']>0), df['linfl_reads'], axis=0) #copy nr of reads from rfl (overwrites NaN from lfl if needed, otherwise it keeps the value that was already there)
	# Result: total column value is only NaN in case neither lfl or rfl reads are present
	df['linfl_reads_av'] = df['linfl_reads']/(df['lfl_junctions'] + df['rfl_junctions'])
	# If there were no flanking junctions (neither linear nor ambiguous, keep the average 0)
	df['linfl_reads_av'] = df[['lfl_reads']].where((df['lfl_reads']==0) & (df['rfl_reads']==0), df['linfl_reads_av'], axis=0)
	df['circ_fraction_fl'] = df['bs_reads']/(df['bs_reads'] + df['linfl_reads_av'])
	## Wald method
	#df['ci_lower_Wald'] = df['circ_proportion'] - 1.96*( (df['circ_proportion'] * (1-df['circ_proportion'])) / (df['bs_reads']+df['linfl_reads_av']) )**(0.5)
	#df['ci_upper_Wald'] = df['circ_proportion'] + 1.96*( (df['circ_proportion'] * (1-df['circ_proportion'])) / (df['bs_reads']+df['linfl_reads_av']) )**(0.5)
	# Adjusted Wald method (Agresti and Coull, 1998; Brown et al., 2001)
	df['p_AC_fl'] = (df['bs_reads'] + 1.96**2/2) / (df['bs_reads']+df['linfl_reads_av'] + 1.96**2)
	df['ci_lower_AC_fl'] = df['p_AC_fl'] - 1.96*( (df['p_AC_fl'] * (1-df['p_AC_fl'])) / (df['bs_reads']+df['linfl_reads_av']+1.96**2) )**(0.5)
	df['ci_upper_AC_fl'] = df['p_AC_fl'] + 1.96*( (df['p_AC_fl'] * (1-df['p_AC_fl'])) / (df['bs_reads']+df['linfl_reads_av']+1.96**2) )**(0.5)
	return(df.round({'linfl_reads_av':3,'circ_fraction_fl':5,'p_AC_fl':5,'ci_lower_AC_fl':5,'ci_upper_AC_fl':5}))

def nonflankcirclevelfraction(df_circ, df_gene):
	"""
	Calculate for every backsplice junction the circRNA fraction based on linear only junction reads in gene
	CI with adjusted Wald method (Agresti and Coull, 1998; Brown et al., 2001)
	"""
	#merge the average of linear only splice junction reads in gene (from gene table) with the individual circ counts (circ table)
	df = df_circ.merge(df_gene[['gene_id','lin_reads_av']], how='left',on='gene_id')
	df.rename(columns={'lin_reads_av':'lin_reads_gene_av'}, inplace=True) #modify name for clarity    
	# If there were no flanking junctions (neither linear nor ambiguous, keep the average 0)
	df['circ_fraction_all'] = df['bs_reads']/(df['bs_reads'] + df['lin_reads_gene_av'])
	## Wald method
	#df['ci_lower_Wald'] = df['circ_proportion'] - 1.96*( (df['circ_proportion'] * (1-df['circ_proportion'])) / (df['bs_reads']+df['linfl_reads_av']) )**(0.5)
	#df['ci_upper_Wald'] = df['circ_proportion'] + 1.96*( (df['circ_proportion'] * (1-df['circ_proportion'])) / (df['bs_reads']+df['linfl_reads_av']) )**(0.5)
	# Adjusted Wald method (Agresti and Coull, 1998; Brown et al., 2001)
	df['p_AC_all'] = (df['bs_reads'] + 1.96**2/2) / (df['bs_reads']+df['lin_reads_gene_av'] + 1.96**2)
	df['ci_lower_AC_all'] = df['p_AC_all'] - 1.96*( (df['p_AC_all'] * (1-df['p_AC_all'])) / (df['bs_reads']+df['lin_reads_gene_av']+1.96**2) )**(0.5)
	df['ci_upper_AC_all'] = df['p_AC_all'] + 1.96*( (df['p_AC_all'] * (1-df['p_AC_all'])) / (df['bs_reads']+df['lin_reads_gene_av']+1.96**2) )**(0.5)
	return(df.round({'lin_reads_gene_av':3,'circ_fraction_all':5,'p_AC_all':5,'ci_lower_AC_all':5,'ci_upper_AC_all':5}))

def splitread(junctions, backsplice_junctions, exons, strandedness):
	# only keep the exons whose chromosomes also occur in normal junction or backsplice junction file (others are not relevant here)
	all_chr = list(set().union(junctions.chromosome.unique(), backsplice_junctions.chromosome.unique()))
	exons_filt = exons[(exons['chromosome'].isin(all_chr))].copy()
	del exons #remove larger exon df
	
	unique_genes = exons_filt["name"].unique()
	output_df = pd.DataFrame(columns = ['gene_id','circ_id','bs_reads', 'lfl_reads','lfl_junctions','rfl_reads','rfl_junctions','lfl_junctions_ambi','rfl_junctions_ambi'])
	#output_df = pd.DataFrame(columns = ['gene_id','circ_id','bs_reads', 'lfl_reads','lfl_junctions','rfl_reads','rfl_junctions','lfl_junctions_ambi','lfl_reads_ambi','rfl_junctions_ambi','rfl_reads_ambi'])
	output_df2 = pd.DataFrame(columns = ['gene_id','linear_reads','linear_junctions','ambiguous_reads','ambiguous_junctions','circ_reads','circ_junctions'])
	output_exons = pd.DataFrame(columns = ['chromosome','start','stop','name','score','strand','exon_id'])
	output_fwjunctions = pd.DataFrame(columns = ['chromosome','start','stop','name','score','strand','reads','gene','pos_to_circ'])
	
	for gene in unique_genes:
		DF_per_gene = exons_filt[exons_filt["name"]== gene] #all exons of the gene		
		# Retrieve backsplice junction reads (same chr, start of BS between min start and max stop of all exons of this gene) 
		# Retrieve junction reads (same chr, same strand, start of junction between min start and max stop (inclusive) of all exons of this gene)
		if strandedness == "yes": #run in strand-specific mode
			backsplice_junctions_gene = backsplice_junctions[(backsplice_junctions['chromosome'].isin(DF_per_gene['chromosome'])) &
				(backsplice_junctions['start'].between(DF_per_gene['start'].min(), DF_per_gene['stop'].max())) &
				(backsplice_junctions['stop'].between(DF_per_gene['start'].min(), DF_per_gene['stop'].max())) &
				(backsplice_junctions['strand'].isin(DF_per_gene['strand']))].copy() 
			junctions_gene = junctions[(junctions['chromosome'].isin(DF_per_gene['chromosome'])) &
				(junctions['start'].between(DF_per_gene['start'].min(), DF_per_gene['stop'].max())) &
				(junctions['stop'].between(DF_per_gene['start'].min(), DF_per_gene['stop'].max())) &
				(junctions['strand'].isin(DF_per_gene['strand']))].copy()
		else: #run in unstranded mode
			backsplice_junctions_gene = backsplice_junctions[(backsplice_junctions['chromosome'].isin(DF_per_gene['chromosome'])) & 
				(backsplice_junctions['start'].between(DF_per_gene['start'].min(), DF_per_gene['stop'].max())) &
				(backsplice_junctions['stop'].between(DF_per_gene['start'].min(), DF_per_gene['stop'].max()))].copy()
			junctions_gene = junctions[(junctions['chromosome'].isin(DF_per_gene['chromosome'])) &
				(junctions['start'].between(DF_per_gene['start'].min(), DF_per_gene['stop'].max())) &
				(junctions['stop'].between(DF_per_gene['start'].min(), DF_per_gene['stop'].max()))].copy()
		junctions_gene['gene'] = gene ## add gene name
		junctions_gene['pos_to_circ'] = 'outside' ## fill in position relative to circ (change in case it is flanking or inside any circ)
		# Make the structure of the output_df: one line per circRNA (BS junction)
		if len(backsplice_junctions_gene)>0: #empty dataframe will have length 0
			#print(gene)
			for i in list(backsplice_junctions_gene['circ_id']):
				output_df = output_df.append({'gene_id':gene,'circ_id':i,'bs_reads':int(0),'lfl_junctions':int(0),'lfl_reads':int(0),'rfl_junctions':int(0),'rfl_reads':int(0),'lfl_junctions_ambi':int(0),'rfl_junctions_ambi':int(0)}, ignore_index=True)
				#output_df = output_df.append({'gene_id':gene,'circ_id':i,'bs_reads':int(0),'lfl_junctions':int(0),'lfl_reads':int(0),'rfl_junctions':int(0),'rfl_reads':int(0),'lfl_junctions_ambi':int(0),'lfl_reads_ambi':int(0),'rfl_junctions_ambi':int(0),'rfl_reads_ambi':int(0)}, ignore_index=True)

		# For every gene, start with counters for linear and ambiguous junctions at 0
		ambiguous_junc, ambiguous_reads = 0,0 
		linear_junc, linear_reads = 0,0
		circ_junc = backsplice_junctions_gene['circ_id'].count() #count total nr of BS junctions in that gene
		circ_reads = backsplice_junctions_gene['reads'].sum() #sum up all the BS junction reads in that gene

		# If there are both linear and backsplice junctions:
		if (len(junctions_gene) > 0) & (len(backsplice_junctions_gene) > 0): 
			# Go over every normal junction & check whether it is ambiguous/linear_only/flanking circles
			for fw in junctions_gene['name'].unique():
				fwjunc = junctions_gene[junctions_gene['name'] == fw].sort_values('reads', ascending=False).drop_duplicates(['name']).copy() #in case there are duplicates, keep only the one with the max nr of reads
				inside,outside,unknown,bs_counter = int(0),int(0),int(0),int(0)
				# Go over all backsplice junctions in the gene and check whether the normal junction falls in between BS sites or is completely outside BS sites
				for bs in backsplice_junctions_gene['circ_id'].unique():
					bsjunc = backsplice_junctions_gene[backsplice_junctions_gene['circ_id'] == bs].copy()
					# if entirely inside circ:
					if (int(bsjunc['start']) <= int(fwjunc['start']) < int(fwjunc['stop']) <= int(bsjunc['stop'])): #junction lies inside a circRNA
						inside =1
						junctions_gene.loc[junctions_gene['name'] == fw,'pos_to_circ'] = 'inside'
						break #if junction lies in between BS sites: do not continue, fw junction can only be ambiguous (not clear if derived from linear or circRNA) -> continue outside the "for bs loop"                    
					# else the junction is flanking the circ or entirely outside the circ
					else:
						if (int(fwjunc['stop']) <= int(bsjunc['start'])) | (int(fwjunc['start']) > int(bsjunc['stop'])) | (int(fwjunc['start']) < int(bsjunc['start']) < int(bsjunc['stop']) < int(fwjunc['stop'])): #junction lies entirely before OR entirely after OR partially before+partially after BS -> outside!
							outside += 1
						elif (int(fwjunc['start']) < int(bsjunc['start']) < int(fwjunc['stop']) <= int(bsjunc['stop'])) | (int(bsjunc['start']) <= int(fwjunc['start']) < int(bsjunc['stop']) < int(fwjunc['stop'])): #junction flanks circle at one side (left or right)
							outside += 1
							junctions_gene.loc[junctions_gene['name'] == fw,'pos_to_circ'] = 'flanking'
						else:
							unknown +=1
					bs_counter +=1 #add 1 to BS junction counter

				# Count total number of ambigous junctions and sum number of reads 
				if inside==1: #as soon as the junction lies inside one circ -> ambiguous
					ambiguous_junc += 1
					ambiguous_reads += int(fwjunc['reads'])
					
				# Count total number of "linear only" junctions and sum number of reads
				elif outside == bs_counter: #if junction lies outside every BS junction, the outside counter will be equal to the BS junction counter
					linear_junc += 1
					linear_reads += int(fwjunc['reads'])
				else:
					print('problem to assign to ambiguous vs linear')
					print(junctions_gene)
					print(backsplice_junctions_gene)
								
				# Find whether the junction flanks certain circles
				# Note1: one junction can be flanking a circle on the left and another one on the right
				# Note2: can be an ambiguous flanking junction (keep counter of this as well to discriminate abscence of flanking junctions from flanking junctions part of another circRNA!)
				lfl_present,rfl_present=0,0			
				for bs in backsplice_junctions_gene['circ_id'].unique(): #go over all backsplice junctions in the gene and check if normal junction flanks BS (and on which side)
					bsjunc = backsplice_junctions_gene[backsplice_junctions_gene['circ_id'] == bs]
					pos = output_df[(output_df['gene_id'] == gene) & (output_df['circ_id'] == bs)].index.values[0] 	#get index of line that contains the specific circle in output_df
					output_df.loc[pos,'bs_reads'] = int(bsjunc['reads']) 	#fill in the backsplice read counts
					
					lfl_present = int(fwjunc['start']) < int(bsjunc['start']) < int(fwjunc['stop']) <= int(bsjunc['stop'])
					rfl_present = int(bsjunc['start']) <= int(fwjunc['start']) < int(bsjunc['stop']) < int(fwjunc['stop'])
					if lfl_present:
						if inside:
							output_df.loc[pos,'lfl_junctions_ambi'] +=1
							#output_df.loc[pos,'lfl_reads_ambi'] += int(fwjunc['reads'])
						else:
							output_df.loc[pos,'lfl_junctions'] +=1	
							output_df.loc[pos,'lfl_reads'] += int(fwjunc['reads']) 	#fill in leftflanking read counts
					elif rfl_present:
						if inside:
							output_df.loc[pos,'rfl_junctions_ambi'] += 1
							#output_df.loc[pos,'rfl_reads_ambi'] += int(fwjunc['reads'])	
						else:
							output_df.loc[pos,'rfl_junctions'] +=1
							output_df.loc[pos,'rfl_reads'] += int(fwjunc['reads'])	#fill in right flanking read counts
					elif inside | (outside>0): 
						pass 
					else: #Just as a check: if junction is not marked as inside, outside or flanking -> problem!
						print("NOT OK! inside:",inside, "outside:",outside, "flanking:",lfl_present,rfl_present)
						print(fwjunc)
						print(bsjunc)
						
			
			
		# if there are only backsplice junctions
		elif (len(backsplice_junctions_gene) > 0): #len(junctions_gene)==0 (so no normal junctions)
			for bs in backsplice_junctions_gene['circ_id'].unique():
				bsjunc = backsplice_junctions_gene[backsplice_junctions_gene['circ_id'] == bs]
				pos = output_df[(output_df['gene_id'] == gene) & (output_df['circ_id'] == bs)].index.values[0]  #get index of line that contains the specific circle in output_df
				output_df.loc[pos,'bs_reads'] = int(bsjunc['reads'])    #fill in the backsplice read counts
				#linear_junc and linear_reads can stay 0
				
		# if there are only linear junctions
		elif (len(junctions_gene) > 0): #len(backsplice_junctions_gene)==0 (no backsplice junctions)
			#assign all linear junction counts to linear_junctions in output_df2
			linear_junc = junctions_gene['name'].count() ##
			linear_reads = junctions_gene['reads'].sum() ##
		
		output_fwjunctions = output_fwjunctions.append(junctions_gene, ignore_index = True) #add the position and gene name of every forward junction to new df
		#print("INTERMEDIATE RES: ambiguous:",ambiguous_junc,"linear:",linear_junc)
		output_df2 = output_df2.append({'gene_id':gene,'linear_junctions':int(linear_junc),'linear_reads':int(linear_reads),'ambiguous_junctions':int(ambiguous_junc),\
										'ambiguous_reads':int(ambiguous_reads), 'circ_junctions':int(circ_junc), 'circ_reads':int(circ_reads)},ignore_index=True)
		
	#if (len(output_df2) > 0): #if the output_df contains at least one gene (otherwise error)
	output_df2 = output_df2.infer_objects() #infer type of columns
	output_df2_ci = genelevelfraction(output_df2)
	
	#Replace 0 by NA in flanking reads columns when there are no flanking junctions except ambiguous ones
	if (len(output_df) > 0): #if the output_df contains at least one circle (otherwise this gives an error)
		output_df = output_df.infer_objects() #infer type of columns
		output_df_ci = circlevelfraction(output_df) #calculate fraction and ci based on flanking linear junction reads
		output_df_ci_flnonfl = nonflankcirclevelfraction(output_df_ci, output_df2_ci) #calculate fraction and ci based on all linear only junction reads
		output_df_ci_flnonfl["lfl_reads"] = output_df_ci_flnonfl.apply(lambda row: label_lflank(row), axis=1)
		output_df_ci_flnonfl["rfl_reads"] = output_df_ci_flnonfl.apply(lambda row: label_rflank(row), axis=1)
	else: #just return the empty df at bs level
                output_df_ci_flnonfl = output_df
	return (output_df_ci_flnonfl, output_df2_ci, output_fwjunctions)

if __name__ == '__main__':
	ap = argparse.ArgumentParser()
	ap.add_argument("-j","--junctions", required=True, help="Normal (forward) junction file (BED format)")
	ap.add_argument("-b","--backsplice", required=True, help="Backsplice junction file (BED format)")
	ap.add_argument("-e","--exons", required=True, help="Exon or gene file (BED format)")
	ap.add_argument("-bc","--bsreads_column",required=True, help="Column in backsplice junction file that contains the number of junction reads")
	ap.add_argument("-fc","--fsreads_column", required=True, help="Column in forward junction file that contains the number of junction reads")
	ap.add_argument("-v","--overlap_column", default="", help="In case the start and stop columns contain the max spanning read positions instead of exact sites, indicate which column contains overlap at left and right side of junction (e.g. column 11 in TopHat junction file, not needed in STAR junction file)") 
	ap.add_argument("-o","--output", default=".", help="Optional output directory")
	ap.add_argument("-n","--name_prefix", default ="", help="Optional pefix (e.g. sample name) for output files")
	ap.add_argument("-ff","--fsfilter", default="1", help="Filter out forward junctions with less reads than this filter (default=1)")
	ap.add_argument("-bf","--bsfilter", default="1", help="Filter out backsplice junctions with less reads than this filter (default=1)")
	ap.add_argument("-s","--strand", choices=['yes','no'], default="yes", help="Run in strand-specific mode or not? (default=yes)")
	args = vars(ap.parse_args())

	outputdir = os.path.normpath(args["output"])
	os.makedirs(outputdir, exist_ok=True) # if outputdir does not exist yet, make one
	nameprefix = args["name_prefix"].strip()
	if nameprefix != "": #if the user specified a certain prefix, add an underscore
		nameprefix += "_"

	#read backsplice junctions
	backsplice_junctions = readbacksplice(args["backsplice"], args["bsreads_column"], args["bsfilter"])
	#read all junctions
	junctions = readjunctions(args["junctions"], args["fsreads_column"], args["overlap_column"], args["fsfilter"])
	#read exons
	exons = readexons(args["exons"])
	
	# execute LinCircSplit and save output files in output directory
	(output_circ,output_gene, output_fwjunctions) =splitread(junctions,backsplice_junctions,exons, args["strand"])
	
	output_circ.to_csv(outputdir +"/"+ nameprefix + "CiLiQuant_circ.txt",sep="\t",index=False)
	output_gene.to_csv(outputdir +"/"+ nameprefix + "CiLiQuant_gene.txt",sep="\t",index=False)
	#output_fwjunctions.to_csv(outputdir +"/"+ nameprefix + "LinCircSplit_FWjunctions.bed",sep="\t",index=False)
