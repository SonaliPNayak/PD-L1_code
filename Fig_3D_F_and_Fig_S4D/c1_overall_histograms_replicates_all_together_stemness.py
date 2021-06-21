import os
import itertools
import numpy as np
import collections
import pandas as pd
import seaborn as sns
import matplotlib as mpl
import scipy.stats as ss
from textwrap import wrap
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans
from gap_statistic import OptimalK
from scipy.stats import gaussian_kde
import scipy.cluster.hierarchy as shc
from sklearn.decomposition import PCA
from scipy.signal import argrelextrema
from scipy.cluster.hierarchy import fcluster
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import AgglomerativeClustering
from matplotlib.colors import LinearSegmentedColormap
from scipy.cluster.hierarchy import dendrogram, linkage

from sklearn.mixture import GaussianMixture

#cmap = mpl.colors.ListedColormap(['#7F7F7F','#FF7F7F','#FFD17F'])
#shc.set_link_color_palette([mpl.colors.rgb2hex(rgb[:3]) for rgb in cmap])

def reading_the_cfg_file(path_to_folder):
	for filename in os.listdir(path_to_folder):
		#print(filename[-4:])
		if filename[-4:] == ".cfg":
			name = filename[:-4]
	flag = -1
	d_genes = {}
	with open(path_to_folder+name+".cfg") as f:
		for line in f:
			a = line[:-1].split("\t")
			if a[0] == "NumberOfRACIPEModels":
				num_models = float(a[1])
			if a[0] == "NumberOfGenes":
				nodes = int(a[1])
				flag = 0
				continue
			if flag >= 0 and flag < nodes:
				flag += 1
				d_genes[a[0]] = a[1]

	return name,num_models,nodes,d_genes

def collating_all_runs_for_z_score_calculation(path_to_dat_files,network_name,genes,num_stability_to_consider):

    state_dataframe = pd.DataFrame(columns = np.sort(genes))
    sdf_1 = pd.DataFrame(columns = np.sort(genes))
    sdf_2 = pd.DataFrame(columns = np.sort(genes))
    sdf_3 = pd.DataFrame(columns = np.sort(genes))

    # iterating over all the files with mono, bi, tri and tetra stable cases
    for i in range(1,num_stability_to_consider+1):
        # reading each file separately, getting the sub data structures and appending it to the state dataframe
        data = pd.read_csv(path_to_dat_files+network_name+"_solution_"+str(i)+"_z_norm.dat",delimiter="\t",header=None)
        data_1 = data[data[0] < 10000]
        data_2 = data[(data[0] >= 10000) & (data[0] < 20000)]
        data_3 = data[(data[0] >= 20000) & (data[0] < 30000)]
        for j in range(0,i):
            sub_dataframe = data[data.columns[len(genes)*j+2:len(genes)*j+(len(genes)+2)]]
            sub_dataframe.columns = np.sort(genes)
            state_dataframe = state_dataframe.append(sub_dataframe,ignore_index = True)
        for j in range(0,i):
            sub_dataframe = data_1[data_1.columns[len(genes)*j+2:len(genes)*j+(len(genes)+2)]]
            sub_dataframe.columns = np.sort(genes)
            sdf_1 = sdf_1.append(sub_dataframe,ignore_index = True)
        for j in range(0,i):
            sub_dataframe = data_2[data_2.columns[len(genes)*j+2:len(genes)*j+(len(genes)+2)]]
            sub_dataframe.columns = np.sort(genes)
            sdf_2 = sdf_2.append(sub_dataframe,ignore_index = True)
        for j in range(0,i):
            sub_dataframe = data_3[data_3.columns[len(genes)*j+2:len(genes)*j+(len(genes)+2)]]
            sub_dataframe.columns = np.sort(genes)
            sdf_3 = sdf_3.append(sub_dataframe,ignore_index = True)
    return state_dataframe,sdf_1,sdf_2,sdf_3

def plot_EMT_score_hist(state_dataframe,genes,network_name,path_to_plots):
	EMT_score = (state_dataframe['ZEB1'] + state_dataframe['SLUG'] - state_dataframe['CDH1'] - state_dataframe['miR200'])/4
	ax = sns.distplot(EMT_score,hist=True,kde_kws={"color": "black"})
	plt.title(network_name+"_EMT_score")
	#plt.xlim(-3.5,3)
	plt.xlabel("EMT Score")
	plt.ylabel("Frequency")
	plt.savefig(path_to_plots+network_name+"_EMT_score_histogram.png")
	plt.close()
	stem_score = (state_dataframe['OCT4'] + state_dataframe['LIN28'] - state_dataframe['miR145'] - state_dataframe['let7'])/4
	ax = sns.distplot(stem_score,hist=True,kde_kws={"color": "black"})
	plt.title(network_name+"_stem_score")
	#plt.xlim(-3.5,3)
	plt.xlabel("Stem Score")
	plt.ylabel("Frequency")
	plt.savefig(path_to_plots+network_name+"_stem_score_histogram.png")
	plt.close()

def plot_scatters(state_dataframe,genes,network_name,path_to_plots):
	state_dataframe["EMT_score"] = (state_dataframe['ZEB1'] + state_dataframe['SLUG'] - state_dataframe['CDH1'] - state_dataframe['miR200'])/4
	state_dataframe["stem_score"] = (state_dataframe['OCT4'] + state_dataframe['LIN28'] - state_dataframe['miR145'] - state_dataframe['let7'])/4
	sns.regplot(x="stem_score", y="PDL1", data=state_dataframe, scatter_kws={'s':2})
	plt.axhline(y=0,c="black")
	#plt.axhline(y=-0.5,c="black")
	plt.axvline(x=-0.25,c="black")
	plt.axvline(x=0.5,c="black")
	#plt.axhline(y=+0.5)
	plt.show()

def plot_clustermap(df):
    #g = sns.clustermap(df, cmap="vlag", standard_scale=1)
    #plt.show()
    df1 = pd.DataFrame(columns = ["EMT_score","PDL1"])
    #g = sns.clustermap(df1, cmap="vlag")
    df["EMT_score"] = (df['ZEB1'] + df['SLUG'] - df['CDH1'] - df['miR200'])/4
    sns.distplot(df["EMT_score"],kde = True,color="black")
    plt.axvline(x=-0.25,c="red")
    plt.axvline(x=0.5,c="red")
    plt.savefig("stemness_EMT_score_hist.png",dpi=800)
    plt.close()
    
    df["Stemness_score"] = (df['OCT4'] + df['LIN28'] - df['let7'] - df['miR145'])/4
    sns.distplot(df["Stemness_score"],kde = True,color="black")
    plt.axvline(x=-0.7,c="red")
    plt.axvline(x=0.6,c="red")
    plt.savefig("stemness_stemness_score_hist.png",dpi=800)
    plt.close()
    
    sns.distplot(df["PDL1"],kde = True,color="black")
    plt.axvline(x=0,c="red")
    plt.savefig("stemness_PDL1_hist.png",dpi=800)
    plt.close()
    
    print("EMT_pDL1",ss.spearmanr(df["EMT_score"],df["PDL1"])[0],ss.spearmanr(df["EMT_score"],df["PDL1"])[1])
    plt.scatter(df["EMT_score"],df["PDL1"],c="black",s=1)
    plt.axhline(y=0,c="red")
    plt.axvline(x=-0.25,c="red")
    plt.axvline(x=0.5,c="red")
    plt.xlabel("EMT_score")
    plt.ylabel("Normalised Frequency")
    plt.savefig("stemness_EMT_score_scatter.png",dpi=800)
    plt.close()
    
    print(ss.spearmanr(df["EMT_score"],df["Stemness_score"])[0],ss.spearmanr(df["EMT_score"],df["Stemness_score"])[1])
    plt.scatter(df["EMT_score"],df["Stemness_score"],c="black",s=1)
    plt.axhline(y=-0.7,c="red")
    plt.axhline(y=0.6,c="red")
    plt.axvline(x=-0.25,c="red")
    plt.axvline(x=0.5,c="red")
    plt.xlabel("EMT_score")
    plt.ylabel("Stemness")
    plt.savefig("stemness_EMT_score_stemness_scatter.png",dpi=800)
    plt.close()
    
    print("stemness_pdl1",ss.spearmanr(df["Stemness_score"],df["PDL1"])[0],ss.spearmanr(df["Stemness_score"],df["PDL1"])[1])
    plt.scatter(df["Stemness_score"],df["PDL1"],c="black",s=1)
    plt.axhline(y=0,c="red")
    plt.axvline(x=-0.7,c="red")
    plt.axvline(x=0.6,c="red")
    plt.xlabel("stemness_score")
    plt.ylabel("Normalised Frequency")
    plt.savefig("stemness_stemness_pdl1_score_scatter.png",dpi=800)
    plt.close()
    
    print(ss.spearmanr(df["CDH1"],df["PDL1"])[0],ss.spearmanr(df["CDH1"],df["PDL1"])[1])
    plt.scatter(df["CDH1"],df["PDL1"],c="black",s=1)
    plt.xlabel("CDH1")
    plt.ylabel("PDL1")
    plt.savefig("stemness_CDH1_PDL1.png",dpi=800)
    plt.close()
    print(ss.spearmanr(df["SLUG"],df["PDL1"])[0],ss.spearmanr(df["SLUG"],df["PDL1"])[1])
    plt.scatter(df["SLUG"],df["PDL1"],c="black",s=1)
    plt.xlabel("SLUG")
    plt.ylabel("PDL1")
    plt.savefig("stemness_SLUG_PDL1.png",dpi=800)
    plt.close()
    print(ss.spearmanr(df["ZEB1"],df["PDL1"])[0],ss.spearmanr(df["ZEB1"],df["PDL1"])[1])
    plt.scatter(df["ZEB1"],df["PDL1"],c="black",s=1)
    plt.xlabel("ZEB1")
    plt.ylabel("PDL1")
    plt.savefig("stemness_ZEB1_PDL1.png",dpi=800)
    plt.close()

def stats_EMT_stemness(df):
    E_pos=E_neg=H_pos=H_neg=M_pos=M_neg=0
    sco = list((df['ZEB1'] + df['SLUG'] - df['CDH1'] - df['miR200'])/4)
    z = list((df['OCT4'] + df['LIN28'] - df['let7'] - df['miR145'])/4)
    for i,j in enumerate(sco):
        score = j
        stem = z[i]
        if (score <= -0.25) & ((stem < -0.7) | (stem > 0.6)):
            E_neg += 1
        if (score <= -0.25) & ((stem >= -0.7) & (stem <= 0.6)):
            E_pos += 1
        if (score > 0.5) & ((stem < -0.7) | (stem > 0.6)):
            M_neg += 1
        if (score > 0.5) & ((stem >= -0.7) & (stem <= 0.6)):
            M_pos += 1
        if (score > -0.25) & (score <= 0.5) & ((stem < -0.7) | (stem > 0.6)):
            H_neg += 1
        if (score > -0.25) & (score <= 0.5) & ((stem >= -0.7) & (stem <= 0.6)):
            H_pos += 1
    print(E_pos,E_neg,H_pos,H_neg,M_pos,M_neg)

def stats(df):
    E_pos=E_neg=H_pos=H_neg=M_pos=M_neg=0
    sco = list((df['ZEB1'] + df['SLUG'] - df['CDH1'] - df['miR200'])/4)
    z = list(df["PDL1"])
    for i,j in enumerate(sco):
        score = j
        pdl1 = z[i]
        if (score <= -0.25) & (pdl1 < 0):
            E_neg += 1
        if (score <= -0.25) & (pdl1 > 0):
            E_pos += 1
        if (score > 0.5) & (pdl1 < 0):
            M_neg += 1
        if (score > 0.5) & (pdl1 > 0):
            M_pos += 1
        if (score > -0.25) & (score <= 0.5) & (pdl1 < 0):
            H_neg += 1
        if (score > -0.25) & (score <= 0.5) & (pdl1 > 0):
            H_pos += 1
    print(E_pos,E_neg,H_pos,H_neg,M_pos,M_neg)
    
def stats_PDL1_stemness(df):
    S_pos=S_neg=NS_pos=NS_neg=0
    sco = list((df['OCT4'] + df['LIN28'] - df['let7'] - df['miR145'])/4)
    z = list(df["PDL1"])
    for i,j in enumerate(sco):
        score = j
        pdl1 = z[i]
        if ((score >= -0.7) & (score <= 0.6)) & (pdl1 < 0):
            S_neg += 1
        if ((score >= -0.7) & (score <= 0.6)) & (pdl1 > 0):
            S_pos += 1
        if ((score < -0.7) | (score > 0.6)) & (pdl1 < 0):
            NS_neg += 1
        if ((score < -0.7) | (score > 0.6)) & (pdl1 > 0):
            NS_pos += 1
    print(S_pos,S_neg,NS_pos,NS_neg)    

core_path = "./../../"
num_stability_to_consider = 5
network_name = "core"
path_to_dat_files = core_path+"stemness/"
path_to_output_z_norm = path_to_dat_files+"Z_normed/"
path_to_plots = path_to_dat_files+"plots/histograms/"
if not os.path.exists(path_to_plots):
	os.makedirs(path_to_plots)
name,num_models,nodes,d_genes = reading_the_cfg_file(path_to_dat_files)
genes = list(d_genes.values())
state_dataframe,sdf_1,sdf_2,sdf_3 = collating_all_runs_for_z_score_calculation(path_to_output_z_norm,network_name,genes,num_stability_to_consider)
plot_clustermap(sdf_1)
stats(sdf_1)
stats(sdf_2)
stats(sdf_3)
stats_EMT_stemness(sdf_1)
stats_EMT_stemness(sdf_2)
stats_EMT_stemness(sdf_3)
stats_PDL1_stemness(sdf_1)
stats_PDL1_stemness(sdf_2)
stats_PDL1_stemness(sdf_3)
plot_scatters(state_dataframe,genes,network_name,path_to_plots)
plot_EMT_score_hist(state_dataframe,genes,network_name,path_to_plots)
