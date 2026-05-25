# The Impact of Wildfires on Rental Prices: Evidence from Canada

**Robin Masson** — Université Paris 1 Panthéon-Sorbonne, M2 Économie  
Master's Research Thesis, 2025–2026

## Research Question

Does climate-induced displacement raise rental prices in receiving cities? This paper examines whether the Williams Lake wildfire (British Columbia, July 2017) increased rents in Canadian cities with pre-existing migration links to the affected area.

## Main Result

Cities with stronger pre-fire migration ties to Williams Lake experienced higher rent growth after 2017. The baseline estimate (city and year fixed effects, 131 cities, 2008–2024) is β = 0.818 (p < 0.01), implying a 19 percent rent increase in the most exposed destination (Kamloops, BC). Results should be interpreted as reduced-form evidence consistent with a wildfire-induced housing demand shock. Pre-trend concerns documented in the paper preclude a causal interpretation.

## Identification Strategy

Bartik-style network exposure index motivated by Goldsmith-Pinkham, Sorkin, and Swift (AER 2020):

Z_WL = share_wl × Post2017

where share_wl is the pre-fire (2016/2017) migration share from Williams Lake to each destination city (in decimal units, 0 to 1), and Post2017 is an indicator equal to 1 from 2017 onward.

## Data

- **Internal migration**: Statistics Canada, Table 17-10-0141-01
- **Rental prices**: CMHC via Statistics Canada, Table 34-10-0133-01
- **Coordinates**: CensusMapper API (Census 2021), supplemented manually

Raw data files are not included in this repository due to Statistics Canada licensing restrictions. See Part I of `replication_code.R` for data download instructions.

## Repository Contents

```
replication_code.R          # Full replication script (17 sections)
master_thesis.pdf           # Final submitted PDF
city_migration_exposure.csv # City-level exposure shares for spatial map
README.md
```

## Replication

Set the working directory to this folder, ensure all packages listed in Part I of the script are installed, and run `replication_code.R` sequentially. Raw data must be downloaded from Statistics Canada before running the script.

**R packages required**: dplyr, stringr, tidyr, cancensus, sf, geosphere, fixest, ggplot2, ggrepel, scales, modelsummary, purrr, broom, rnaturalearth, rnaturalearthdata
