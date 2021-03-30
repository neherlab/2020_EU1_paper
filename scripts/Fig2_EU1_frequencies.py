
import pandas as pd
import datetime
import json
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from shutil import copyfile
from collections import defaultdict
from matplotlib.patches import Rectangle
from colors_and_countries import country_list, uk_countries, all_countries
from travel_data import *
from helpers import load_cluster, non_zero_counts, trim_last_data_point
from paths import *
from clusters import clusters as cluster_definitions

def marker_size(n):
    if n>200:
        return 150
    elif n>100:
        return 75
    elif n>50:
        return 35
    elif n>10:
        return 15
    else:
        return 5


if __name__=="__main__":
    clus = 'S222'
    figure_path = 'figures/'
    fmt = 'pdf'
    print_files = True
    # Convert into dataframe
    cluster_data, total_data, country_info = load_cluster(clus)

    width = 1
    smoothing = np.exp(-np.arange(-10,11)**2/2/width**2)
    smoothing /= smoothing.sum()

    cutoff_num_seqs = 200
    hasMinNumber = []
    hasMinNumber_countries = []
    for index, row in country_info.iterrows():
        if row.num_seqs > cutoff_num_seqs and index not in uk_countries:
            hasMinNumber.append(True)
            hasMinNumber_countries.append(index)
        else:
            hasMinNumber.append(False)

    country_info["has_min"] = hasMinNumber

    print(f"{len(hasMinNumber_countries)} countries have more than {cutoff_num_seqs} in the cluster.")
    print(f"Countries who have more than {cutoff_num_seqs} in the cluster:", hasMinNumber_countries, "\n")
    print(country_info)
    json_output = {}
    #plot those with >=200 seqs, plus a few extra to try:
    countries_to_plot = [
        'France',
        'Norway',
        'United Kingdom',
        'Netherlands',
        'Spain',
        'Belgium',
        'Switzerland',
        'Ireland',
        'Italy',
        'Denmark'
    ]

    clus_display  = cluster_definitions[clus]["display_name"]
    json_output[clus_display] = {}

    # Make a plot
    #fig = plt.figure(figsize=(10,5))
    #fig, axs=plt.subplots(1,1, figsize=(10,5))
    fs = 14
    #fig, (ax1, ax2, ax3) = plt.subplots(nrows=3, sharex=True,figsize=(10,7),
    #                                    gridspec_kw={'height_ratios':[1,1,3]})
    # Change to just show Travel to spain only. see above for old 3 panel version
    fig, (ax1, ax3) = plt.subplots(nrows=2, sharex=True,figsize=(10,6),
                                        gridspec_kw={'height_ratios':[1, 3]})
    ax1.grid(True, axis='x', which='major')
    ax3.grid(True, axis='x', which='major')
    ax1.set_axisbelow(True)
    ax3.set_axisbelow(True)
    i=0
    #for coun in [x for x in countries_to_plot]:
    for coun in travel_order:
        if coun in q_free_to_spain:
            q_times = q_free_to_spain[coun]
            strt = datetime.datetime.strptime(q_times["start"], "%Y-%m-%d")
            end = datetime.datetime.strptime(q_times["end"], "%Y-%m-%d")
            y_start = i*0.022
            height = 0.02
            ax1.add_patch(Rectangle((strt,y_start), end-strt, height,
                        ec=country_styles[coun]['c'], fc=country_styles[coun]['c']))
            #ax1.text(strt, y_start+0.002, q_times["msg"], fontsize=fs*0.8)
            ax1.text(strt, y_start+0.003, q_times["msg"], fontsize=fs*0.8)
            if coun == "Denmark":
                strt = datetime.datetime.strptime(q_free_to_spain["Denmark2"]["start"], "%Y-%m-%d")
                end = datetime.datetime.strptime(q_free_to_spain["Denmark2"]["end"], "%Y-%m-%d")
                ax1.add_patch(Rectangle((strt,y_start), end-strt, height,
                        ec=country_styles[coun]['c'], fc="none", hatch="/"))
                ax1.text(strt, y_start+0.003, q_free_to_spain["Denmark2"]["msg"], fontsize=fs*0.8)
        i=i+1
    ax1.set_ylim([0,y_start+height])
    ax1.text(datetime.datetime.strptime("2020-05-03", "%Y-%m-%d"), y_start,
            "Quarantine-free", fontsize=fs)
    ax1.text(datetime.datetime.strptime("2020-05-03", "%Y-%m-%d"), y_start-height-0.005,
            "Travel to/from Spain", fontsize=fs)
    ax1.text(datetime.datetime.strptime("2020-05-03", "%Y-%m-%d"), y_start-height-height-0.01,
            "(on return)", fontsize=fs)
    ax1.get_yaxis().set_visible(False)

    #for a simpler plot of most interesting countries use this:
    for coun in [x for x in countries_to_plot]:
        week_as_date, cluster_count, total_count = non_zero_counts(cluster_data, total_data, coun, smoothing=smoothing)
        # remove last data point if that point as less than frac sequences compared to the previous count
        week_as_date, cluster_count, total_count = trim_last_data_point(week_as_date, cluster_count, total_count, frac=0.1, keep_count=10)

        ax3.plot(week_as_date, cluster_count/total_count,
                color=country_styles[coun]['c'],
                linestyle=country_styles[coun]['ls'], label=coun)
        ax3.scatter(week_as_date, cluster_count/total_count, s=[marker_size(n) for n in total_count],
                color=country_styles[coun]['c'],
                linestyle=country_styles[coun]['ls'])

        json_output[clus_display][coun] = {}
        json_output[clus_display][coun]["week"] = [datetime.datetime.strftime(x, "%Y-%m-%d") for x in week_as_date]
        json_output[clus_display][coun]["smoothed_total_sequences"] = [int(x) for x in total_count]
        json_output[clus_display][coun]["smoothed_cluster_sequences"] = [int(x) for x in cluster_count]

    for ni,n in enumerate([0,10,50,100,200]):
        ax3.scatter([week_as_date[0]], [0.08+ni*0.07], s=marker_size(n+0.1), edgecolor='k', facecolor='w')
        ax3.text(week_as_date[1], 0.06+ni*0.07, f"n>{n}" if n else "n<=10")
        #          color=country_styles[coun]['c'], linestyle=country_styles[coun]['ls'], label=coun)


    plt.legend(ncol=2, fontsize=fs*0.8, loc=2)
    fig.autofmt_xdate(rotation=30)
    ax3.tick_params(labelsize=fs*0.8)
    ax3.set_ylim(0,1)
    ax3.set_ylabel('frequency', fontsize=fs)
    ax3.set_xlim(datetime.datetime(2020,5,1), datetime.datetime(2020,12,1))
    plt.show()
    plt.tight_layout()

    #spain opens borders
    ax3.text(datetime.datetime.strptime("2020-06-21", "%Y-%m-%d"), 0.05,
            "Spain opens borders", rotation='vertical', fontsize=fs*0.8)

    if print_files:
        trends_path = figure_path+f"Fig2_{clus}_overall_trends.{fmt}"
        plt.savefig(trends_path)
        copypath = trends_path.replace("trends", "trends-{}".format(datetime.date.today().strftime("%Y-%m-%d")))
        copyfile(trends_path, copypath)

    if print_files:
        with open(figure_path+f'{clus_display}_data.json', 'w') as fh:
            json.dump(json_output[clus_display], fh)
