import matplotlib.patches as mpatches
import pandas as pd
import datetime
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from shutil import copyfile
from matplotlib.patches import Rectangle
import matplotlib.patches as mpatches
import copy
from colors_and_countries import *
from clusters import clusters as cluster_definitions
from helpers import *

figure_path = "./figures/"
grey_color = "#cccccc"
fmt = "pdf"

if __name__=="__main__":

    clusters = {}
    for clus in cluster_definitions.keys():
        clusters[clus] = {}
        print(f"Running cluster {clus}")
        cluster_data, total_data, summary = load_cluster(clus)
        clusters[clus].update(cluster_definitions[clus])
        clusters[clus]['cluster_data'] = cluster_data
        clusters[clus]['country_info'] = summary

    print("\n\nYou can check on the country_info with: clusters['S477']['country_info'] ")

    for clus in clusters.keys():
        c_i = clusters[clus]['country_info']
        c_i[c_i['num_seqs'] > 10]
        print(f"Countries with >10 seqs in cluster {clus}:")
        print("\t", ", ".join(c_i[c_i['num_seqs'] > 10].index))
        print("\n")

    #fix cluster order in a list so it's reliable
    clus_keys = list(clusters.keys())

    my_df = [clusters[x]["country_info"] for x in clus_keys]
    all_num_seqs = pd.concat([x.loc[:,"num_seqs"] for x in my_df],axis=1)
    all_num_seqs.columns = clus_keys

    has10 = []
    has10_countries = []
    for index, row in all_num_seqs.iterrows():
        if any(row > 10) and index not in uk_countries:
            has10.append(True)
            has10_countries.append(index)
        else:
            has10.append(False)

    all_num_seqs["has_10"] = has10

    print("Countries who have more than 10 in any cluster:", has10_countries, "\n")
    print(all_num_seqs)

    ############## Plot

    #Use the S477 data to decide what to plot.
    countries_to_plot = ["France", "United Kingdom", "Netherlands",
        "Switzerland", "Belgium", "Spain", "Norway", "Ireland", "Denmark", "Czech Republic"]
    #Remember to adjust the number of axes if needed below....


    country_week = {clus: {} for clus in clusters}

    #fig, ax1 = plt.subplots(nrows=1,figsize=(10,7))
    fs = 14
    #fig, (ax1, ax2, ax3, ax4, ax5, ax6, ax7) = plt.subplots(nrows=7, sharex=True,figsize=(9,9),
    #                                    gridspec_kw={'height_ratios':[1,1,1,1,1,1,1]})
    fig, axs = plt.subplots(nrows=len(countries_to_plot)+1, sharex=True,figsize=(9,11))
    smoothing = get_smoothing(width=1)

    week_as_dates = {}

    #for coun in [x for x in countries_to_plot]:
    for coun, ax in zip(countries_to_plot, axs[1:]):
        i=0
        first_clus_count = []
        ptchs = []
        #for cluster_data in [cluster_data_S477, cluster_data_S222]:
        for clus in clusters.keys():
            cluster_data = clusters[clus]['cluster_data']
            week_as_date, cluster_count, total_count = non_zero_counts(cluster_data, total_data, coun, smoothing=smoothing)

            week_as_dates[coun] = week_as_date
            country_week[clus][coun] = cluster_count/total_count

            linesty = '-'
            lab = clusters[clus]["display_name"] #f"{clus}"
            if i == 0:
                first_clus_count = [0] * len(cluster_count)
            #if i == 1:
            #    linesty = ':'
            cluster_count = first_clus_count + cluster_count #unindented
            ax.fill_between(week_as_date, first_clus_count/total_count, cluster_count/total_count, facecolor=clusters[clus]['col'])
            patch = mpatches.Patch(color=clusters[clus]['col'], label=lab)

            #exclude 501 for now (not present)
            if clus is not "S501":
                ptchs.append(patch)
            if i == len(clusters)-1 :
                ax.fill_between(week_as_date, cluster_count/total_count, 1, facecolor=grey_color)
                patch = mpatches.Patch(color=grey_color, label=f"other")
                ptchs.append(patch)
            #if i == 0:
            first_clus_count = cluster_count # unindented
            i+=1
        ax.text(datetime.datetime(2020,6,1), 0.7, coun, fontsize=fs)
        ax.tick_params(labelsize=fs*0.8)
        #ax.set_ylabel('frequency')
        #ax.legend(ncol=1, fontsize=fs*0.8, loc=2)

    axs[0].legend(handles=ptchs, loc=3, fontsize=fs*0.7, ncol=4)
    axs[0].axis('off')
    fig.autofmt_xdate(rotation=30)
    plt.show()
    plt.tight_layout()

    plt.savefig(figure_path+f"EUClusters_compare.{fmt}")
    trends_path = figure_path+f"EUClusters_compare.{fmt}"
    copypath = trends_path.replace("compare", "compare-{}".format(datetime.date.today().strftime("%Y-%m-%d")))
    copyfile(trends_path, copypath)


#for clus in clusters.keys():
#    for coun in countries_to_plot:
#        print(clus, coun, len(country_week[clus][coun]))