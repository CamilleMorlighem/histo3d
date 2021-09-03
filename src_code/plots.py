import pandas as pd 
import geopandas as gpd 
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.patches as mpatches

"""
This script is used to generate plots for all the quality metrics of the building plots after they were (i) extracted from the historical maps, (ii) vectorised and (iii) generalised 
"""

df = pd.read_csv(r"C:\Users\Camille\Documents\mem_docs\3_ANALYSIS\final_workflow\metrics_bxls.csv", sep =";")

metrics = [("avg_area", "Average building plot area"), ("med_area", "Median building plot area"), ("med_perim", "Median building plot perimeter"), 
("med_n_vert", "Median number of vertices per building plot"),  ("avg_convex", "Average building plot convexity"),
("med_rect", "Median building plot rectangularity"), ("nber_plots", "Number of building plots"), ("nber_concave_poly", "Number of concave building plots")]
df=df.sort_values(by=['year'])
df = df.set_index('year')

"""
fig, axes = plt.subplots(4,2)
lst=[]

for i in range(0,4): 
    for j in range(0,2):
        lst.append([i,j])

idx=0
"""

for metric, metric_name in metrics: 
    subset0 = df.loc[df["file"] == metric]
    subset = subset0.drop(columns=["map"])
    ymax=subset.drop(columns=["file"]).to_numpy().max()
   
    """             
    ax = subset.plot.bar(rot=0,  color={"raw_plots": "xkcd:light salmon", "alpha_plots": "xkcd:pastel red", "generalized_plots":"xkcd:dark red"})
    ax = subset.plot.bar(rot=0,  color={"raw_plots": "xkcd:baby blue", "alpha_plots": "xkcd:sky blue", "generalized_plots":"xkcd:medium blue"})
    ax = subset.plot.bar(rot=0,  color={"raw_plots": "xkcd:baby blue", "alpha_plots": "xkcd:light salmon", "generalized_plots":"xkcd:light grey green"})
    ax = subset.plot.bar(rot=0,  color={"raw_plots": '#96ceb4', "alpha_plots": '#ffcc5c', "generalized_plots":'#ff6f69'})
    ax = subset.plot.bar(rot=0,  color={"raw_plots": '#96ceb4', "alpha_plots": '#F8C471', "generalized_plots":"xkcd:medium blue"})
    ax = subset.plot.bar(rot=0,  color={"raw_plots": "#2EB76C", "alpha_plots": "#338FCF", "generalized_plots":"#F7AA3B"}, title = "Average building plot area")
    """

    ax = subset.plot.bar(rot=0,  figsize =(5.8,3), color={"raw_plots": "#1A5276", "alpha_plots": "#6F903F", "generalized_plots":"#C0392B"}, 
    title = metric_name,use_index=False, edgecolor='white', linewidth=1, width=0.5, ylabel= "Metric value", xlabel="Historical map", legend=False)#, xticks={"Brussels 1700", "Brussels 1890", "Brussels 1924"})

    nan_idx = np.where(subset['generalized_plots'].isnull())[0]
    plt.axvspan(nan_idx+0.09, nan_idx+0.25, facecolor='white', edgecolor="grey", hatch='//')
    ax.set_xticklabels(subset0.map)
    ax.tick_params(axis="y", labelsize=9)
    ax.tick_params(axis="x", labelsize=9)

    plt.subplots_adjust(right=0.68,bottom=0.2)
    nan_legend = mpatches.Patch(facecolor='white', edgecolor='gray', hatch='//', label='Not available')
    vect_leg= mpatches.Patch(color='#1A5276',label='Vectorisation')
    gene1_leg= mpatches.Patch(color="#6F903F", label="$1^{st}$ generalisation")
    gene2_leg= mpatches.Patch(color="#C0392B", label="$2^{nd}$ generalisation")
    plt.xlabel("Historical map", labelpad=8, fontsize=9.5)
    plt.ylabel("Metric value", fontsize=9.5)
    plt.title(metric_name, fontsize=11,fontweight="bold")
    axes = plt.gca()
    plt.legend(handles=[vect_leg, gene1_leg, gene2_leg, nan_legend], frameon=False, prop={'size': 10}, bbox_to_anchor=(1, 0.5), loc=6)
    #plt.legend(handles=[vect_leg, gene1_leg, gene2_leg],frameon=False, prop={'size': 10}, bbox_to_anchor=(1, 0.5), loc=6)
    #plt.savefig("C:\\Users\\Camille\\Documents\\mem_docs\\3_ANALYSIS\\final_workflow\\graphs\\bxls" + metric + ".png")
    #idx+=1

#plt.tight_layout() 
#plt.legend(handles=[vect_leg, gene1_leg, gene2_leg],frameon=False, prop={'size': 9}, loc="best")# bbox_to_anchor=(1, 0.5), loc=6)
#plt.subplots_adjust(wspace=0.3, hspace=0.6)
plt.show()
    
