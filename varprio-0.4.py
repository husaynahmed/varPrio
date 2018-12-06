##############################################
##############################################
##############################################
############# v a r P r i o 0.4 ##############
##############################################
##############################################
##############################################

#####################################################################
### varPrio is a tool for the prioritization of genetic variants  ###
### from WES/WGS data. Variants which are relevant and associated ###
### to the disease phenotype are prioritized based on in silico   ###
### predictions of damaging mutations and based on occurrence or  ###
### frequency across pedigrees and in the population. varPrio is  ###
### developed as part of the The Accelerator program for Discovery###
### in Brain disorders using Stem cells (ADBS) at NCBS. 	      ###
### This program is tailor-made for large-scale analysis of 	###
### pedigrees recruited in ADBS. The input formats recognized by 	###
### this tool is based on the files generated in ADBS. This tool	###
### is not generalized to read any type of annotated VCFs.  	###
### Please read the README file before using this program. 		###
#####################################################################

# packages import
import pandas as pd
import numpy as np
import os
import glob
import argparse


###
##### Argument for varPrio
###

parser = argparse.ArgumentParser(description='varPrio version 0.4', epilog='Please give absolute(full) path to all the files. Refer the README at https://github.com/husaynahmed/varPrio for more information.')
parser.add_argument('-T', '--typeofvariant',
                     help='Type of variant to prioritize {snp,indel}',
                     required='True',
                     choices=["snp","indel"])

parser.add_argument('-I', '--inputfileinfo',
                     help='Path to the text file containing 3 rows. 1st row - Sample identifier of the affected individuals; 2nd row - Family identifier; 3rd row - Path to the annotated file (ANNOVAR tab delimmited TXT files).',
                     required='True')

parser.add_argument('-PC', '--populationcontrol',
                     help='Path to population control variant data file.',
                     required='True')

parser.add_argument('-AFC', '--allfamilycontrol',
                     help='Path to all familial control variant data file of multiple families.',
			   required='True')

parser.add_argument('-O', '--outdir',
                     help='Path to the output directory where the varprio results will be written.',
                     required='True')


args = parser.parse_args()
typevar = args.typeofvariant
inputFI = args.inputfileinfo
pcfile = args.populationcontrol
afcfile = args.allfamilycontrol
outdir = args.outdir
#####

###
##### Removing all VPR files before analysis in the Output directory
###
os.chdir(outdir)

print "||| Deleting all vpr files in the output directory |||"
for f in glob.glob("*.vpr"):
    os.remove(f)
#####

###
##### Reading INPUT File information 
###

input = pd.read_table(inputFI, header=None, prefix="C",comment='#', usecols=(0,1,2))
print "||| Reading cases input files information |||"
sample_num = input["C0"].tolist()
family_num = input["C1"].tolist()
file_names = input["C2"].tolist()

len1 = len(sample_num)
len2 = len(family_num)
len3 = len(file_names)
#####


###
##### Reading population control variant data 
###
print "||| Reading population control variant data |||"
df_pc = pd.read_table(pcfile,header=None,dtype=str,comment='#',low_memory=False)
#####

###
##### Reading familial control variant data 
###
print "||| Reading familial	 control variant data |||"
df_afc = pd.read_table(afcfile,header=None,dtype=str,comment='#',low_memory=False)
#####

###
##### VARIANT TYPE :  SNP
###

# Creating the names of the dataframes 
def frNameCr(nd):
	"Create the names of dataframes"
	frName = []
	for i in range(0,nd):
		frName.append("df_"+str(i))
	return frName

###
##### SNP PRIORITIZATION
###
	

if(typevar == "snp"):
	print "||| Performing SNP Prioritization |||"
	
	# Reading the files into the dataframes
	j = 0
	frP = frNameCr(len1)
	print "||| Reading all annotated files - annovar txt files tab demilimited|||"
	for i in file_names:
		df = pd.read_table(i, header=None,comment='#',low_memory=False,skiprows=1)
		frP[j] = df	
		#print df
		j=j+1
	print "||| Creating dataframes |||"
	
	###
	##### Step - 01 #####
	### Any variant shared by >= n-1 affected individuals

	print "||| Executing n-1 prioritization approach |||" 
	
	print "||| Writing temp1_allVar.vpr |||"	
	j=0
	for i in sample_num:
		dfx = frP[j]
		dfx["C_ind",] = i
		dfx["C_fam",] = family_num[j] 
		j = j+1	
		dfx.to_csv('temp1_allVar.vpr', mode='a', header=False,sep='\t',index=False)
	
	dfx_1 = pd.read_table("temp1_allVar.vpr",header=None,dtype=str, comment='#',low_memory=False)
	dfx_1.sort_values(by=[0,1],inplace=True)

	
	print "||| STEP 1: Any variant shared by >= n-1 affected individuals |||"
	
	n = len1 - 1
	step_1 = dfx_1.groupby([0,1]).filter(lambda x: len(x) >= n)
	step_1.sort_values(by=[0,1],inplace=True)

	step_1.to_csv('step1_upto1GTmissingallowed.vpr',mode='a',header=False, sep='\t',index=False)
	print "||| Writing step1_upto1GTmissingallowed.vpr |||"
	
	#######################################################################################

	###
	##### STEP - 02 #####
	### Filter : 1-5 predictors 
	print "||| STEP 2: Applying filter - Any of the 5 predictors |||"
	
	dfy = step_1.loc[(step_1[25] == "D") | ((step_1[27] == "D") | (step_1[27] == "P")) | (step_1[31] == "D") | ((step_1[33] == "A") | (step_1[33] == "D")) | ((step_1[35] == "H") | (step_1[35] == "M"))]
	
	dfy.to_csv('step2_1to5P_damaging.vpr', mode='a', header=False,sep='\t',index=False)
	print "||| Writing step2_1to5P_damaging.vpr |||"
	
	#######################################################################################	
	### Count the number of predictors for each variant (np)
	dfy = dfy.reset_index(drop=True)
	len_dfy = len(dfy)
	for i in range(0,len_dfy):
		count = 0
		if (dfy.loc[i,25] == "D"):
			count = count + 1
		if ( (dfy.loc[i,27] == "D") | (dfy.loc[i,27] == "P") ):
			count = count + 1
		if (dfy.loc[i,31] == "D"):
			count = count + 1
		if ( (dfy.loc[i,33] == "A") | (dfy.loc[i,33] == "D") ):
			count = count + 1
		if ( (dfy.loc[i,35] == "H") | (dfy.loc[i,35] == "M") ):
			count = count + 1
		dfy.loc[i,"np"] = str(count)

	dfy.sort_values(by=["np",0,1],inplace=True,ascending=[False,True,True])
	dfy.to_csv('step2_1to5P_sorted.vpr', mode='a', header=False,sep='\t',index=False)
	print "||| Writing step2_1to5P_sorted.vpr |||"

	#######################################################################################

	###
	##### STEP - 03 #####
	### Addition of PC and AFC info
	print "||| Writing LIST1B_step3_upto1GTmissingallowed_withPCAFC.vpr |||"
	step_1_pc = step_1.merge(df_pc,how='left',on=[0,1])
	step_1_pc_afc = step_1_pc.merge(df_afc,how='left',on=[0,1])
	step_1_pc_afc.to_csv('LIST1B_step3_upto1GTmissingallowed_withPCAFC.vpr', mode='a', header=False,sep='\t',index=False,na_rep="0") 
	print "||| Writing LIST1A_step3_upto1GTmissingallowed_withPCAFC.vpr |||"
	step_1_pc_afc.to_csv('LIST1A_step3_upto1GTmissingallowed_withPCAFC.vpr', mode='a', columns=(0,1,3,4,5,6,8,9,13,15,22,23,25,27,31,33,35,71,72,"2_y",2), header=False,sep='\t',index=False,na_rep="0") 
	

	print "||| Writing LIST2B_step3_1to5P_withPCAFC.vpr |||"
	dfy_pc = dfy.merge(df_pc,how='left',on=[0,1])
	dfy_pc_afc = dfy_pc.merge(df_afc,how='left',on=[0,1])
	dfy_pc_afc.to_csv('LIST2B_step3_1to5P_withPCAFC.vpr', mode='a', header=False,sep='\t',index=False,na_rep="0") 
	print "||| Writing LIST2A_step3_1to5P_withPCAFC.vpr |||"
	dfy_pc_afc.to_csv('LIST2A_step3_1to5P_withPCAFC.vpr', mode='a', columns=(0,1,3,4,5,6,8,9,13,15,22,23,25,27,31,33,35,71,72,"np","2_y",2), header=False,sep='\t',index=False,na_rep="0")
	

	print "||| COMPLETED SNP Prioritization !!! |||"

	#######################################################################################
	#######################################################################################
	#######################################################################################


	
###
##### INDEL PRIORITIZATION
###

if(typevar == "indel"):
	print "||| Performing INDEL Prioritization |||"
	# Reading the files into the dataframes
	j = 0
	frP = frNameCr(len1)
	print "||| Reading all annotated files - annovar txt files tab demilimited|||"
	for i in file_names:
		df = pd.read_table(i, header=None,comment='#',low_memory=False,skiprows=1)
		frP[j] = df	
		j=j+1
	print "||| Creating dataframes |||"
	
	###
	##### STEP - 01 #####
	### Filter : Frameshift / stopgain / stoploss 
	print "||| STEP 1: Applying Frameshift / stopgain / stoploss  filter |||"
	j=0
	for i in sample_num:
		dfx = frP[j]
		dfy = dfx.loc[(dfx[8] == "frameshift insertion") | (dfx[8] == "frameshift deletion") | (dfx[8] == "stopgain") | (dfx[8] == "stoploss")]
		dfy["C_ind",] = i
		dfy["C_fam",] = family_num[j] 
		j = j+1	
		dfy.to_csv('step1_filtered_INDEL.vpr', mode='a', header=False,sep='\t',index=False)
	print "||| Writing step1_filtered_INDEL.vpr |||"
	# Sorting
	list5 = pd.read_table("step1_filtered_INDEL.vpr",header=None,dtype=str, comment='#',low_memory=False)
	list5.sort_values(by=[0,1],inplace=True)
	
	###
	##### STEP - 02 #####
	### Adding population control details
	print "||| STEP 2: Presence/Absence in Population control |||"
	df_7 = pd.merge(list5, df_pc, how='inner',on=[0,1,2])
	keys = [0,1,2]
	i5 = list5.set_index(keys).index
	i6 = df_7.set_index(keys).index
	df_8 = list5[~i5.isin(i6)]
	df_7["PC",] = "Present in PC"
	df_8["PC",] = "Absent in PC"
	df_9 = pd.concat([df_8, df_7])
	df_9.to_csv('step2_prioritized_INDEL_LIST3.vpr',mode='a',header=False, sep='\t',index=False)
	print "||| Writing step2_prioritized_INDEL_LIST3.vpr |||"
	print "||| SUCCESSFULLY performed INDEL prioritization |||"


	#######################################################################################
	#######################################################################################
	#######################################################################################


#######################################################################################
###  Developed: Husayn Ahmed P 
###  https://github.com/husaynahmed
###  At ADBS, National Centre for Biological Sciences, Bengaluru
###
###  Please acknowledge the usage of varPrio by citing the following :
###  Suhas Ganesh,  Husayn Ahmed P,  Ravi K Nadella, Ravi P More, Manasa Sheshadri, Biju Viswanath, Mahendra Rao, Sanjeev Jain, The ADBS consortium, Odity Mukherjee. 2018. Exome sequencing in families with severe mental illness identifies novel and rare variants in genes implicated in Mendelian neuropsychiatric syndromes. Psychiatry and Clinical Neurosciences. doi:10.1111/pcn.12788
###
#######################################################################################


