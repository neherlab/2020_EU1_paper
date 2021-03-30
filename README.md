# [Spread of a SARS-CoV-2 variant through Europe in the summer of 2020](https://www.medrxiv.org/content/10.1101/2020.10.25.20219063v3)

**by Emma B. Hodcroft, Moira Zuber, Sarah Nadeau, Timothy G. Vaughan, Katharine H. D. Crawford, Christian L. Althaus, Martina L. Reichmuth, John E. Bowen, Alexandra C. Walls, Davide Corti, Jesse D. Bloom, David Veesler, David Mateo, Alberto Hernando, Iñaki Comas, Fernando González Candelas, SeqCOVID-SPAIN consortium, Tanja Stadler, Richard A. Neher**

## Abstract

_Following its emergence in late 2019, severe acute respiratory syndrome coronavirus 2 (SARS-CoV-2) has caused a global pandemic resulting in unprecedented efforts to reduce transmission and develop therapies and vaccines (WHO Emergency Committee, 2020; Zhu et al., 2020). Rapidly generated viral genome sequences have allowed the spread of the virus to be tracked via phylogenetic analysis (Worobey et al., 2020; Hadfield et al., 2018; Pybus et al., 2020). While the virus spread globally in early 2020 before borders closed, intercontinental travel has since been greatly reduced, allowing continent-specific variants to emerge. However, within Europe travel resumed in the summer of 2020, and the impact of this travel on the epidemic is not well understood. Here we report on a novel SARS-CoV-2 variant, 20E (EU1), that emerged in Spain in early summer, and subsequently spread to multiple locations in Europe. We find no evidence of increased transmissibility of this variant, but instead demonstrate how rising incidence in Spain, resumption of travel across Europe, and lack of effective screening and containment may explain the variant’s success. Despite travel restrictions and quarantine requirements, we estimate 20E (EU1) was introduced hundreds of times to countries across Europe by summertime travellers, likely undermining local efforts to keep SARS-CoV-2 cases low. Our results demonstrate how a variant can rapidly become dominant even in absence of a substantial transmission advantage in favorable epidemiological settings. Genomic surveillance is critical to understanding how travel can impact SARS-CoV-2 transmission, and thus for informing future containment strategies as travel resumes._

## Scripts and resources

This repository contains the relevant code and data to reproduce the analysis presented in the manuscript.
The analysis can be performed in the standard nextstrain [conda environment](https://github.com/nextstrain/conda) with the additional package `geopandas`.

### Main text

### Generating summary data

the script `scripts/generate_cluster_counts.py` ingests metadata and a file listing mutations for each strain and outputs the number of times a particular sequence cluster is observed by country and calendar week.
The raw data can't be shared due to GISAIDs data access agreement, but all sequences that go into this analysis are listed in the acknowledgement table.
The output of this scripts are the tables in the directory `cluster_tables` which are committed to this repository and on which the subsequent analysis are based.

#### Figure 1

Figure 1 is an annotated screen shot using the [nextstrain/auspice](https://github.com/nextstrain/auspice).
It is based on the auspice json file [EU1_cluster.json](XXX).

#### Figure 2

Figure 2 summarizes the frequency dynamics of the EU1 cluster in different countries in Europe.
These include data available on GISAID until 2021-01-20.
The full list of accession numbers is provided in [acknowledgement_tables/acknowledgement_tables.tsv](acknowledgement_tables/acknowledgement_tables.tsv).

The data aggregated by calendar week and country are available in XXX and serve as alternative entry point for the script [scripts/Fig2_EU1_frequencies.py](scripts/Fig2_EU1_frequencies.py).

### Extended data
