# ==============================================================================
# Title:     Retrospective analysis of clinical and environmental genotyping results revealing persistence of Pseudomonas aeruginosa in the water system of a large tertiary children’s hospital in England
# Author:    Esha Sheth, Yu Wan
# Affiliation: Department of Infection Prevention and Control/ Alder Hey Children's NHS Foundation Trust/ Department of Infection Prevention and Control, Alder Hey Children's Hospital Trust/ Department of Pathology, Alder Hey Children's Hospital Trust/ Department of Estates and Facilities, Alder Hey Children's Hospital Trust/ David Price Evans Global Health and Infectious diseases Research Group, University of Liverpool
# Date:      March 2025
# Version:   1.0
# 
# Description:
#   This script processes extracted VNTR data from UKHSA PDF reports.
#   Steps include: data import, source classification from filenames,
#   date harmonisation, anonymisation, deduplication, and creation of a
#   unique profile summary table with cluster labels.
# 
# Input:     summary_water_180825.xlsx (or .csv)
# Output:    vntr_deduplicated.csv, unique_vntr_summary_clusters.csv
# 
# Dependencies: tidyverse, lubridate
# 
# License:   GNU General Public License v3.0
# Contact:   [e.d.m.sheth@liverpool.ac.uk]
# ==============================================================================

# Install and load packages

install.packages("tidyverse")
install.packages(lubridate")
install.packages("dplyr")
install.packages("tidyr")
install.packages("readxl")
library(tidyverse)
library(lubridate)
library(dplyr)
library(tidyr)
library(readxl)


# Import UKHSA summary results 

setwd("C:/path/to/your/folder")

water_df <- read_excel("filename.xlsx", sheet = "SheetName")


# Clean imported data frame
water_df <- water_df %>%
  drop_na(VNTR_Profile)


# Initial cleaning & date parsing

df <- water_df %>%
  mutate(
    Date_Received = dmy(Date_Received),
    Date_Collected = dmy(Date_Collected),
    
    Harmonised_Date = coalesce(Date_Collected, Date_Received),
    
    # Assign Source and clean any 'Unknown' in Excel later
    Source = case_when(
      str_detect(Patient_Name, regex("ENVIRON|TAP|FEEDWATER|THEATRE|WARD|ENVIRONMENTAL|WATER|OUTLET|"), ignore_case = TRUE) ~ "Environmental",
      str_detect(Patient_Name, regex("^[A-Za-z]+\\s*[A-Za-z+", ignore_case = TRUE) ~ "Clinical"),
      TRUE ~ "Unknown"
      ),
    
    # Assign Isolation_Class, clean further on Excel and investigate gaps from IPC team
    Isolation_Class = case_when(
      str_detect(Isolation_Site, regex("Blood|csf|BAL|alveolar|Bronchial", ignore.case = TRUE)) ~"Invasive",
      str_detect(Isolation_Site, regex("sputum|urine|wound|swab|faeces|aspirate", ignore_case = TRUE)) ~ "TRUE",
      Source == "Environmental" ~ "-",
      TRUE ~ "Unknown"
    )
  )

# Check Source and Isolation_site assignment on excel


write_xlsx(df, "filename_cleaned.xlsx")


# Data anonymisation #####

df_anon <- df %>%
  mutate(
    #Remove all identifiable columns
    across(c(Sender_Ref, PHE_Ref, Hospital_ID, Patient_Name, Date_of_birth), ~ NA_character_),
    AHPa_ID = paste0("AHPa", sprintf("%03d", row_number())),
    across(c(Harmonised_Date, VNTR_Profile, Organism, Isolation_Site, Source, Isolation_Class, AHPa_ID))
  )

write_xlsx(df_anon, "filename_cleaned_anonymised.xlsx")

# Deduplication #####

df_dedup <- water_df %>%
  distinct(Harmonised_collection_dates, VNTR_Profile, Patient_Name, .keep_all = TRUE)

nrow(df_dedup)

write_xlsx(df_dedup, "filename_deduplicated.xlsx")

# Create unique summary table (with cluster labels) #####
# Make sure all 'Unknowns' are clarified manually and import consolidated clean file for next steps


unique_summary <- df_dedup %>%
  group_by(VNTR_Profile) %>%
  summarise(
    Count             = n(),
    AHPA_IDs          = paste0(AHPA_ID, collapse = ", "),
    Num_Clinical      = sum(Source == "Clinical", na.rm = TRUE),
    Num_Environmental = sum(Source == "Environmental", na.rm = TRUE),
    Num_Invasive      = sum(Isolation_Class == "Invasive", na.rm = TRUE),
    Num_Non_invasive  = sum(Isolation_Class == "Non-invasive", na.rm = TRUE),
    
    Isolate_Date      = paste(sort(unique(Harmonised_Date)), collapse = ", "),
    First_Date        = min(Harmonised_Date, na.rm = TRUE),
    Last_Date         = max(Harmonised_Date, na.rm = TRUE),
    
    .groups = "drop"
  ) %>%
  filter(Count >= 2) %>%
  arrange(desc(Count)) %>%
  mutate(
    Cluster_Label = paste0("Cluster ", row_number(),
                           " (n=", Count, ")")
  ) %>%
  select(Cluster_Label, VNTR_Profile, Count, everything())

write_xlsx(unique_summary, "unique_vntr_summary.xlsx")

unique_summary <- unique_summary %>%
  mutate(Cluster_ID = paste0("Cluster ", row_number(),
                             " (n=", Count, ")"))

# Create datasets for iToL ####
# Clean and sort imported dataframe

vntr_df <- water_df %>%
  separate(VNTR_Profile, into = paste0("Gene", 1:9), sep = ",", convert = TRUE) %>%
  select(Gene1:Gene9) %>%
  mutate(across(Gene1:Gene9, ~ifelse(. == "-", 0, .)))

# Generate a distance matrix; n is the number of isolates

Isolate_ID <- paste0("AHPa", sprintf("%03d", 1 : n))

D <- matrix(0, nrow = n, ncol = n, dimnames = list(Isolate_ID, Isolate_ID))

for (i in 1 : n-1) {
  for (j in 2: n) {
    d <- sum(vntr_df[i, ] != vntr_df[j, ])
    D[i, j] <- d
    D[j, i] <- d
  }
}

write.csv(D, file = "VNTR_distance_matrix.csv")

# Tree reconstruction

install.packages("ape")
install.packages("phangorn")

library(ape)
library(phangorn)

tr_nj <- nj(as.dist(D))
write.tree(tr_nj, file = "VNTR_NJ.newick")


# Comparison of distances
D_tr <- cophenetic.phylo(tr_nj_mid)
D_tr_original <- cophenetic.phylo(tr_nj)

r <- cor(as.vector(D), as.vector(D_tr))
r_tr <- cor(as.vector(D_tr), as.vector(D_tr_original))

# Dendogram
dnd <- hclust(d = as.dist(D), method ="complete")
tr_dnd <- as.phylo(dnd)
D_dnd <- cophenetic.phylo(tr_dnd)
write.tree(tr_dnd, file = "VNTR_hclust_complete.newick")
r_dnd_D <- cor(as.vector(D), as.vector(D_dnd))
r_dnd_tr <- cor(as.vector(D_tr), as.vector(D_dnd))


water_df <- water_df %>%
  select(-Isolate_ID) %>%                     # remove old column
  mutate(Isolate_ID = paste0("AHPa", 
                             sprintf("%03d", row_number())))
itol_source <- water_df %>%
  select(Isolate_ID, Isolation_Class)
itol_source <- itol_source %>%
  mutate(color = case_when(
    Isolation_Class == "Non-invasive"    ~ "#FF9999",   # light red
    Isolation_Class == "Invasive" ~ "#990000",  # dark red
    Isolation_Class == "-"     ~ "#00AA00",  # light green
    TRUE          ~ "#CCCCCC" 
  )) # grey for anything else)

# Sort columns to make it copy-paste ready for the iToL colorstrip template, which can found on iToL help page

itol_source <- itol_source[, c(1, 3, 2)]


write.table(itol_source, file = "vntr_classification_colorstrip.txt",
            sep = "\t",
            row.names = FALSE,
            col.names = FALSE,
            quote = FALSE)

# Itol year colorstrip

itol_date <- water_df %>%
  select(Isolate_ID, Harmonised_collection_dates)

itol_date <- itol_date %>%
  mutate(Year = format(as.Date(Harmonised_collection_dates), "%Y"))

years <- 2016:2024

# Soft rainbow using HCL (low chroma = muted)
year_colors <- setNames(
  hcl(
    h = seq(0, 360, length.out = length(years) + 1)[1:length(years)],
    c = 40, # chroma (colour intensity) - lower = more muted
    l = 70 # luminance (brightness) - mid-range = soft
  ),
  years
)

install.packages("colorspace")
library(colorspace)

colors_rainbow <- qualitative_hcl(9, palette = "Set 3")

year_colors <- setNames(colors_rainbow, years)

itol_date <- itol_date %>%
  mutate(color = year_colors[Year])

write.table(itol_date, file = "vntr_year_colorstrip_.txt",
            sep = "\t",
            row.names = FALSE,
            col.names = FALSE,
            quote = FALSE)



