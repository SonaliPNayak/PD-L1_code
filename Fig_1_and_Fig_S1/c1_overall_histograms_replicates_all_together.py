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

def fancy_dendrogram(*args, **kwargs):
    max_d = kwargs.pop('max_d', None)
    if max_d and 'color_threshold' not in kwargs:
        kwargs['color_threshold'] = max_d
    annotate_above = kwargs.pop('annotate_above', 0)

    ddata = dendrogram(*args, **kwargs)

    if not kwargs.get('no_plot', False):
        plt.title('Hierarchical Clustering Dendrogram (truncated)')
        plt.xlabel('sample index or (cluster size)')
        plt.ylabel('distance')
        for i, d, c in zip(ddata['icoord'], ddata['dcoord'], ddata['color_list']):
            x = 0.5 * sum(i[1:3])
            y = d[1]
            if y > annotate_above:
                plt.plot(x, y, 'o', c=c)
                #plt.annotate("%.3g" % y, (x, y), xytext=(0, -5),textcoords='offset points',va='top', ha='center')
        if max_d:
            plt.axhline(y=max_d, c='k')
    return ddata

def plotting_the_dendrogram(state_dataframe,network_name,path_to_plots,genes):
	max_d = 85
	Z = shc.linkage(state_dataframe, method='ward')
	shc.set_link_color_palette(['#7F7F7F','#FFD17F','#FF7F7F'])
	dend = fancy_dendrogram(Z,show_leaf_counts=True,leaf_rotation=45.,leaf_font_size=12.,max_d=max_d,truncate_mode='lastp',p=25,show_contracted=False,above_threshold_color="black")
	clusters = fcluster(Z, max_d, criterion='distance')
	#plotting_overall_heatmaps(state_dataframe,network_name,path_to_plots,pd.Series(clusters))
	state_dataframe["cluster_label"] = pd.Series(clusters)
	#for cluster_number in ["1","2","3","4","5","6"]:
	#	cluster_chacteristics(cluster_number,state_dataframe,genes)
	plt.savefig("dendrogram.png")
	plt.close()

def plot_clustermap(df):
    #g = sns.clustermap(df, cmap="vlag", standard_scale=1)
    #plt.show()
    df1 = pd.DataFrame(columns = ["EMT_score","PDL1"])
    #g = sns.clustermap(df1, cmap="vlag")
    df["EMT_score"] = (df['ZEB1'] + df['SLUG'] - df['CDH1'] - df['miR200'])/4
    sns.distplot(df["EMT_score"],kde = True,color="black")
    plt.axvline(x=-0.25,c="red")
    plt.axvline(x=0.5,c="red")
    plt.savefig("EMT_score_hist.png",dpi=800)
    plt.close()
    sns.distplot(df["PDL1"],kde = True,color="black")
    plt.axvline(x=0,c="red")
    plt.savefig("PDL1_hist.png",dpi=800)
    plt.close()
    print(ss.spearmanr(df["EMT_score"],df["PDL1"])[0],ss.spearmanr(df["EMT_score"],df["PDL1"])[1])
    plt.scatter(df["EMT_score"],df["PDL1"],c="black",s=1)
    plt.axhline(y=0,c="red")
    plt.axvline(x=-0.25,c="red")
    plt.axvline(x=0.5,c="red")
    plt.xlabel("EMT_score")
    plt.ylabel("Normalised Frequency")
    plt.savefig("EMT_score_scatter.png",dpi=800)
    plt.close()
    print(ss.spearmanr(df["CDH1"],df["PDL1"])[0],ss.spearmanr(df["CDH1"],df["PDL1"])[1])
    plt.scatter(df["CDH1"],df["PDL1"],c="black",s=1)
    plt.xlabel("CDH1")
    plt.ylabel("PDL1")
    plt.savefig("CDH1_PDL1.png",dpi=800)
    plt.close()
    print(ss.spearmanr(df["SLUG"],df["PDL1"])[0],ss.spearmanr(df["SLUG"],df["PDL1"])[1])
    plt.scatter(df["SLUG"],df["PDL1"],c="black",s=1)
    plt.xlabel("SLUG")
    plt.ylabel("PDL1")
    plt.savefig("SLUG_PDL1.png",dpi=800)
    plt.close()
    print(ss.spearmanr(df["ZEB1"],df["PDL1"])[0],ss.spearmanr(df["ZEB1"],df["PDL1"])[1])
    plt.scatter(df["ZEB1"],df["PDL1"],c="black",s=1)
    plt.xlabel("ZEB1")
    plt.ylabel("PDL1")
    plt.savefig("ZEB1_PDL1.png",dpi=800)
    plt.close()

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

def cluster_wise_expression(df):
    E_ZEB1=[]
    E_SLUG=[]
    E_CDH1=[]
    E_PDL1=[]
    H_ZEB1=[]
    H_SLUG=[]
    H_CDH1=[]
    H_PDL1=[]
    M_ZEB1=[]
    M_SLUG=[]
    M_CDH1=[]
    M_PDL1=[]
    sco = list((df['ZEB1'] + df['SLUG'] - df['CDH1'] - df['miR200'])/4)
    z = list(df["ZEB1"])
    s = list(df["SLUG"])
    c = list(df["CDH1"])
    p = list(df["PDL1"])
    for i,j in enumerate(sco):
        score = j
        if (score <= -0.25):
            E_ZEB1.append(z[i])
            E_SLUG.append(s[i])
            E_CDH1.append(c[i])
            E_PDL1.append(p[i])
        if (score > 0.5):
            M_ZEB1.append(z[i])
            M_SLUG.append(s[i])
            M_CDH1.append(c[i])
            M_PDL1.append(p[i])
        if (score > -0.25) & (score <= 0.5):
            H_ZEB1.append(z[i])
            H_SLUG.append(s[i])
            H_CDH1.append(c[i])
            H_PDL1.append(p[i])

    e_means, e_std = (np.mean(E_CDH1), np.mean(E_SLUG), np.mean(E_ZEB1), np.mean(E_PDL1)), (np.std(E_CDH1), np.std(E_SLUG), np.std(E_ZEB1), np.std(E_PDL1))
    h_means, h_std = (np.mean(H_CDH1), np.mean(H_SLUG), np.mean(H_ZEB1), np.mean(H_PDL1)), (np.std(H_CDH1), np.std(H_SLUG), np.std(H_ZEB1), np.std(H_PDL1))
    m_means, m_std = (np.mean(M_CDH1), np.mean(M_SLUG), np.mean(M_ZEB1), np.mean(M_PDL1)), (np.std(M_CDH1), np.std(M_SLUG), np.std(M_ZEB1), np.std(M_PDL1))

    ind = np.arange(len(e_means))  # the x locations for the groups
    width = 0.25  # the width of the bars

    fig, ax = plt.subplots()
    rects1 = ax.bar(ind - width, e_means, width, yerr=e_std, alpha=0.5, color = 'black', capsize=2)
    rects2 = ax.bar(ind, h_means, width, yerr=h_std, alpha=0.5, color = 'orange', capsize=2)
    rects2 = ax.bar(ind + width, m_means, width, yerr=m_std, alpha=0.5, color = 'red', capsize=2)

    # Add some text for labels, title and custom x-axis tick labels, etc.
    ax.set_ylabel('Scores')
    ax.set_title('Individual Gene Expression')
    ax.set_xticks(ind)
    ax.set_xticklabels(('CDH1', 'SLUG', 'ZEB1', 'PDL1'))
    #ax.legend()
    ax.set_ylim([-1.6,1.6])
    plt.savefig("phenotype_wise_expression.png",dpi=1000)
    plt.close()

def PCA_analysis(state_dataframe_PCA,genes,network_name,path_to_plots):
    num_clusters = [2,3,4,5,6,7,8]
    #scaling the normalized z dataframe to reduce mean to 0 and std to 1
    scaled_state_dataframe = StandardScaler().fit_transform(state_dataframe_PCA)

    #performing pca
    pca = PCA(n_components = 2)
    principalComponents = pca.fit_transform(scaled_state_dataframe)
    principalDf = pd.DataFrame(data = principalComponents
             , columns = ['principal component 1', 'principal component 2'])

    #kmeans_silhouette(principalDf,range(8),path_to_plots)

    cluster = AgglomerativeClustering(n_clusters=3, affinity='euclidean', linkage='ward')
    cluster.fit_predict(state_dataframe_PCA)

    #gene_levels(state_dataframe_PCA,genes)

    elements_count = collections.Counter(cluster.labels_)
    # printing the element and the frequency
    for key, value in elements_count.items():
       print(f"{key}   {value}")

    PrincipalComponents = principalDf

    #to find the variance of each prinicpal component and the contribuition of each gene in the PCA
    explained_variance=pca.explained_variance_
    explained_variance_ratio=pca.explained_variance_ratio_

    pca_components = pd.DataFrame(pca.components_,columns=state_dataframe_PCA.columns,index=['PC-1','PC-2'])
    print(pca_components)

    fig, ax = plt.subplots()
    cmap = mpl.colors.ListedColormap(['#7F7F7F','#FF7F7F','#FFD17F'])
    fig = principalDf.plot(kind='scatter',x='principal component 1', y='principal component 2', c=cluster.labels_,cmap=cmap,s=0.8,ax=ax)
    ax.set_xlabel('PC-1('+str(round(explained_variance_ratio[0],4)*100)+'% variance)')
    ax.set_ylabel('PC-2('+str(round(explained_variance_ratio[1],4)*100)+'% variance)')
    #plt.show()
    plt.savefig("PCA_scatter_colored_via_clustering.png",dpi=1000)
    plt.close()

    print(principalDf)
    fig, ax = plt.subplots()
    cmap = mpl.colors.ListedColormap(['#7F7F7F','#FF7F7F'])
    x = state_dataframe_PCA['PDL1'].apply(np.sign)
    fig = principalDf.plot(kind='scatter',x='principal component 1', y='principal component 2', c=x,cmap=cmap,s=0.8,ax=ax)
    ax.set_xlabel('PC-1('+str(round(explained_variance_ratio[0],4)*100)+'% variance)')
    ax.set_ylabel('PC-2('+str(round(explained_variance_ratio[1],4)*100)+'% variance)')
    #plt.show()
    plt.savefig("PCA_scatter_colored_via_PDL1.png",dpi=1000)
    plt.close()

core_path = "./../../"
num_stability_to_consider = 5
network_name = "core"
path_to_dat_files = core_path+"core/"
path_to_output_z_norm = path_to_dat_files+"Z_normed/"
path_to_plots = path_to_dat_files+"plots/histograms/"
if not os.path.exists(path_to_plots):
	os.makedirs(path_to_plots)
name,num_models,nodes,d_genes = reading_the_cfg_file(path_to_dat_files)
genes = list(d_genes.values())
state_dataframe,sdf_1,sdf_2,sdf_3 = collating_all_runs_for_z_score_calculation(path_to_output_z_norm,network_name,genes,num_stability_to_consider)
stats(sdf_1)
stats(sdf_2)
stats(sdf_3)
cluster_wise_expression(sdf_1)
plot_EMT_score_hist(state_dataframe,genes,network_name,path_to_plots)
plotting_the_dendrogram(sdf_1,network_name,path_to_plots,genes)
PCA_analysis(sdf_1,genes,network_name,path_to_plots)
plot_clustermap(sdf_1)



