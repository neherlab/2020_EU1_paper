# THIS SCRIPT IS TO BE RUN FROM THE MAIN REPO LEVEL AS
# python scripts/generate_cluster_counts.py
###
# It converts gisaid metadata and sequence data into tables that
# aggregate sequence counts in different clusters by country and week
import pandas as pd
import datetime
import numpy as np
from collections import defaultdict
from colors_and_countries import country_list, all_countries, uk_countries
from helpers import CW_to_date, date_to_CW
from clusters import clusters as cluster_definitions
#read in bad sequences to be excluded
from bad_sequences import bad_seqs

# Give a cutoff date to allow back-fill of metadata but no sequences past end Nov
max_cutoff_date = "2020-11-30"
cluster_path = "../clusters"

if __name__=="__main__":
    # Get diagnostics file - used to get list of SNPs of all sequences, to pick out seqs that have right SNPS
    diag_file = "gisaid_data/sequence-diagnostics_2021-01-20.tsv"
    diag = pd.read_csv(diag_file, sep='\t', index_col=False)
    # Read metadata file
    input_meta = "gisaid_data/metadata_2021-01-20.tsv"
    meta = pd.read_csv(input_meta, sep='\t', index_col=False)
    meta = meta.fillna('')

    has_valid_date = set()
    for ri, row in meta.iterrows():
        try:
            date = datetime.datetime.strptime(row.date, '%Y-%m-%d')
            if date <= datetime.datetime.strptime(max_cutoff_date, "%Y-%m-%d"):
                has_valid_date.add(row.strain)
        except:
            pass

    clus_to_run = cluster_definitions.keys()

    for clus in clus_to_run:
        print(f"\nRunning cluster {clus}\n")

        snps = cluster_definitions[clus]['snps']

        # get the sequences that we want - which are 'part of the cluster:
        wanted_seqs = []
        for index, row in diag.iterrows():
            strain = row['strain']
            if strain not in has_valid_date:
                continue
            snplist = row['all_snps']
            if not pd.isna(snplist):
                intsnp = [int(x) for x in snplist.split(',')]
                if all(x in intsnp for x in snps):
                    wanted_seqs.append(row['strain'])

        # Exclude sequences that seem to have bad dates or are overdiverged compared to date
        for key, value in bad_seqs.items():
            bad_seq = meta[meta['strain'].isin([key])]
            if not bad_seq.empty and bad_seq.date.values[0] == value and key in wanted_seqs:
                wanted_seqs.remove(key)

        print("Sequences found: ")
        print(len(wanted_seqs)) # how many are there?
        print("\n")

        # get metadata for these sequences
        cluster_meta = meta[meta['strain'].isin(wanted_seqs)]
        observed_countries = [x for x in cluster_meta['country'].unique()]

        # What countries do sequences in the cluster come from?
        print(f"The cluster is found in: {observed_countries}\n")
        if clus != "S222":
            print("Remember, countries are not set for clusters other than S222")

        if len(observed_countries) > len(country_list) and clus=="S222":
            print("\nWARNING!! Appears a new country has come into the cluster!")
            print([x for x in observed_countries if x not in country_list])

        # Let's get some summary stats on number of sequences, first, and last, for each country.
        country_info = pd.DataFrame(index=all_countries, columns=['first_seq', 'num_seqs', 'last_seq', "sept_to_nov_freq"])
        country_dates = {}
        cutoffDate = datetime.datetime.strptime("2020-09-01", '%Y-%m-%d')
        for coun in all_countries:
            if coun in uk_countries:
                temp_meta = cluster_meta[cluster_meta['division'].isin([coun])]
            else:
                temp_meta = cluster_meta[cluster_meta['country'].isin([coun])]
            country_info.loc[coun].first_seq = temp_meta['date'].min()
            country_info.loc[coun].last_seq = temp_meta['date'].max()
            country_info.loc[coun].num_seqs = len(temp_meta)

            country_dates[coun] = [datetime.datetime.strptime(dat, '%Y-%m-%d') for dat in temp_meta['date']]

            herbst_dates = [x for x in country_dates[coun] if x >= cutoffDate]
            if coun in uk_countries:
                temp_meta = meta[meta['division'].isin([coun])]
            else:
                temp_meta = meta[meta['country'].isin([coun])]
            all_dates = [datetime.datetime.strptime(x, '%Y-%m-%d') for x in temp_meta["date"] if len(x) is 10 and "-XX" not in x and datetime.datetime.strptime(x, '%Y-%m-%d') >= cutoffDate]
            if len(all_dates) == 0:
                country_info.loc[coun].sept_to_nov_freq = 0
            else:
                country_info.loc[coun].sept_to_nov_freq = round(len(herbst_dates)/len(all_dates),2)

        country_info.sort_index(inplace=True)
        country_info.to_csv(f'cluster_tables/{clus}_summary.tsv', sep='\t')
        print(country_info)
        print("\n")

        country_info_df = pd.DataFrame(data=country_info)
        if clus == "S222":
            print("\nOrdered list:")
            print(country_info_df.sort_values(by="first_seq"))
            print("\n")

        # We want to look at % of samples from a country that are in this cluster
        # To avoid the up-and-down of dates, bin samples into weeks
        countries_to_plot = country_list
        acknowledgement_table = []
        # Get counts per week for sequences in the cluster
        clus_week_counts = {}
        for coun in all_countries:
            counts_by_week = defaultdict(int)
            for dat in country_dates[coun]:
                counts_by_week[date_to_CW(dat)]+=1  #returns ISO calendar week
            clus_week_counts[coun] = counts_by_week

        cluster_data = pd.DataFrame(data=clus_week_counts).sort_index()
        cluster_data.to_csv(f'cluster_tables/{clus}_counts.tsv', sep='\t')

    # Get counts per week for sequences regardless of whether in the cluster or not - from week 20 only.
    total_week_counts = {}
    for coun in all_countries:
        counts_by_week = defaultdict(int)
        if coun in uk_countries:
            temp_meta = meta[meta['division'].isin([coun])]
        else:
            temp_meta = meta[meta['country'].isin([coun])]
        #week 20
        for ri, row in temp_meta.iterrows():
            strain = row['strain']
            if strain not in has_valid_date:
                continue
            dat = row.date
            if len(dat) is 10 and "-XX" not in dat: # only take those that have real dates
                dt = datetime.datetime.strptime(dat, '%Y-%m-%d')
                #exclude sequences with identical dates & underdiverged
                if coun == "Ireland" and dat == "2020-09-22":
                    continue

                wk = date_to_CW(dt) #returns ISO calendar week
                if wk >= 20:
                    counts_by_week[wk]+=1
                    acknowledgement_table.append([row.strain, row.gisaid_epi_isl, row.originating_lab, row.submitting_lab, row.authors])
        total_week_counts[coun] = counts_by_week

    with open(f'acknowledgment_tables/{clus}_acknowledgement_table.tsv', 'w') as fh:
        fh.write('#strain\tEPI_ISOLATE_ID\tOriginating lab\tsubmitting lab\tauthors\n')
        for d in acknowledgement_table:
            fh.write('\t'.join(d)+'\n')


    # Convert into dataframe and write to file
    total_data = pd.DataFrame(data=total_week_counts).sort_index()
    total_data.to_csv(f'cluster_tables/total_counts.tsv', sep='\t')
