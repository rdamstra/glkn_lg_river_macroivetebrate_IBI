##################### NPS GLKN LARGE RIVER IBI CALCULATIONS ######################
# Rick Damstra, Aquatic Ecologist (used DOI Version of ChatGPT for help on some sections of code)
# v0.1 (Beta) September 2024

#Built to work with GLKN Large River Macroinvertebrate data from NPS Data Store:

#Damstra RA and Hester CM. 2023. Saint Croix National Scenic Riverway Invertebrate Data Package by the Great Lakes Inventory and Monitoring Network. 
#National Park Service. Fort Collins CO https://doi.org/10.57830/2301852

#Link: https://irma.nps.gov/DataStore/Reference/Profile/2301852

#Calculations are based on:

#Brian M. Weigel and Jeffrey J. Dimick "Development, validation, and application of a macroinvertebrate-based Index of Biotic Integrity for 
#nonwadeable rivers of Wisconsin," Journal of the North American Benthological Society 30(3), 665-679, (31 May 2011). https://doi.org/10.1899/10-161.1

#Link: https://www.jstor.org/stable/10.1899/10-161.1?item_view=read_online

library(tidyverse)
#set working directory and import local copy of data downloaded from NPS Data Store. May need to change filename of .csv file, in this case it is "glkn_inverts1.csv"
setwd("C:/Users/rdamstra/OneDrive - DOI/Desktop/GLKN RV WQ Invert R Code")
i <- read.csv("C:/Users/rdamstra/OneDrive - DOI/Desktop/GLKN RV WQ Invert R Code/glkn_inverts1.csv")

i_old <- read.csv("C:/Users/rdamstra/OneDrive - DOI/Desktop/GLKN RV WQ Invert R Code/glkn_inverts.csv")
# Create the Duplicate_Samples data frame by filtering rows with "Duplicate" in STATION_NAME
Duplicate_Samples <- i %>%
  filter(grepl("Duplicate", STATION_NAME, ignore.case = TRUE))

# remove duplicates from the original data frame 'i'
i <- i %>%
  filter(!grepl("Duplicate", STATION_NAME, ignore.case = TRUE))

#create column to indicate if the taxon is unique. "FALSE" in DUPLICATE_T_4_SITEDATE means that the result is unique. 
#This will be used in calculations to avoid double counts of the same taxa when some animals can't be identified past genus and there are others IDed to species.
#This is pertinent to below metrics "INSECT_T" and "EPT_T", which are calculated by counting numbers of taxa present in 
#class insecta and orders "Ephemeroptera", "Plecoptera" and "Trichoptera". Otherwise, scores of INSECT_T and EPT_T will be inflated.   
# Function to determine unique taxon


#chiros <- i %>% 
#  group_by(STATION_ID) %>% 
#  filter(grepl("Chironomidae", FAMILY, ignore.case = TRUE)) %>% 
#  unique(i$TAXON_ID_EXPERT)





# Clean the TAXON_ID_EXPERT column and determine uniqueness based on exact matches
i <- i %>%
  mutate(TSN_LOOKUP = trimws(as.character(TSN_LOOKUP))) %>%  # Remove whitespace
  group_by(STATION_ID, YEAR) %>%
  mutate(
    # Create a new column to identify unique taxa 
    UNIQUE_TAXON = ifelse(n_distinct(na.omit(TSN_LOOKUP)) == 1, "Y" , "N")
  ) %>% 
  ungroup()


################## IBI Metrics ###################
ibi_calc_metrics <- i %>% 
  group_by(STATION_ID, YEAR) %>% 
  summarise(
    INSECT_T = sum(CLASS == "Insecta"),                              #Count Insecta Taxa for INSECT_T
    Insect_I = sum(TOTAL[CLASS == "Insecta"], na.rm = TRUE), #Sum total for all Insecta
    Total = sum(TOTAL, na.rm = TRUE),  
    INSECT_PERCENT_I = (Insect_I / Total) * 100, # Calculate percentage CHANGE TO INDIVIDUALS using "TOTAL" column (Insect_Percent_I)
    #Total of all results in TOTAL
    EPT_T = sum(ORDER %in% c("Ephemeroptera", "Plecoptera", "Trichoptera")),  # Count EPT orders
    DOM3_P_I = {
      #calculate the % of the top 3 abundant taxa based on TOTAL (individual organisms)
      top_taxa <- cur_data() %>% 
        group_by(TSN_LOOKUP) %>% 
        summarise(total_abundance = sum(TOTAL, na.rm = TRUE), .groups = 'drop') %>% 
        arrange(desc(total_abundance)) %>% 
        head(3) #get the top 3 taxa
      
      #Calculate the total count of the top 3 taxa
      total_top_taxa <- sum(top_taxa$total_abundance)
      (total_top_taxa / Total) * 100 #Calculate percentage of top 3 taxa
    },
    MPTV = sum(TVAL, na.rm = TRUE) / sum(!is.na(TVAL)),  # Calculate MPTV 
    INTOL_EPT2_P_I = {
      #calculate the sum of TOTAL for interolerant EPT individuals with TVAL<=2
      intolerant_EPT_count <- sum(TOTAL[ORDER %in% c("Ephemeroptera", "Plecoptera", "Trichoptera") & TVAL <= 2], na.rm = TRUE)
      if (Total > 0) {
        (intolerant_EPT_count / Total) * 100 #Calculate percentage of intolerant EPT individuals
      } else {
          0 #avoids division by zero
      }
      },
    TOLERANT_CHIRO8_P_I = {
      #calculate the sum of TOTAL for tolerant chiro individuals with TVAL>=8
      tolerant_chiro_count <- sum(TOTAL[FAMILY %in% c("Chironomidae") & TVAL >= 8], na.rm = TRUE)
      if (Total > 0) {
        (tolerant_chiro_count / Total) * 100 #Calculate percentage of tolerant chiro individuals
      } else {
        0 #avoids division by zero
      }
    },
    ECO_FTN_UNIQUE = n_distinct(EcoFTN),  # Count of unique values in EcoFTN
    GATH_P_I = {
      #calculate the sum of TOTAL for in individuals in the Gatherer guild (TFUNCT contains "gatherer")
      gatherer_count <- sum(TOTAL[grepl("gatherer", TFUNCT)], na.rm = TRUE)
      if (Total > 0) {
        (gatherer_count / Total) * 100 #Calculate percentage of gatherers among all individuals
      } else {
        0 #avoids division by zero
      }
    },
    SCRAPE_P_I = {
      #calculate the sum of TOTAL for in individuals in the scraper guild (TFUNCT contains "scraper")
      scraper_count <- sum(TOTAL[grepl("scraper", TFUNCT)], na.rm = TRUE)
      if (Total > 0) {
        (scraper_count / Total) * 100 #Calculate percentage of gatherers among all individuals
      } else {
        0 #avoids division by zero
      }
    },
    .groups = 'drop' # drop grouping by station ID and year after summarizing
  )
############# IBI Metrics for Duplicate samples ################
dup_ibi_calc_metrics <- Duplicate_Samples %>% 
  group_by(STATION_ID, YEAR) %>% 
  summarise(
    INSECT_T = sum(CLASS == "Insecta"),                              #Count Insecta Taxa for INSECT_T
    Insect_I = sum(TOTAL[CLASS == "Insecta"], na.rm = TRUE), #Sum total for all Insecta
    Total = sum(TOTAL, na.rm = TRUE),  
    INSECT_PERCENT_I = (Insect_I / Total) * 100, # Calculate percentage CHANGE TO INDIVIDUALS using "TOTAL" column (Insect_Percent_I)
    #Total of all results in TOTAL
    EPT_T = sum(ORDER %in% c("Ephemeroptera", "Plecoptera", "Trichoptera")),  # Count EPT orders
    DOM3_P_I = {
      #calculate the % of the top 3 abundant taxa based on TOTAL
      top_taxa <- cur_data() %>% 
        group_by(TSN_LOOKUP) %>% 
        summarise(total_abundance = sum(TOTAL, na.rm = TRUE), .groups = 'drop') %>% 
        arrange(desc(total_abundance)) %>% 
        head(3) #get the top 3 taxa
      
      #Calculate the total count of the top 3 taxa
      total_top_taxa <- sum(top_taxa$total_abundance)
      (total_top_taxa / Total) * 100 #Calculate percentage of top 3 taxa
    },
    MPTV = sum(TVAL, na.rm = TRUE) / sum(!is.na(TVAL)),  # Calculate MPTV 
    INTOL_EPT2_P_I = {
      #calculate the sum of TOTAL for interolerant EPT individuals with TVAL<=2
      intolerant_EPT_count <- sum(TOTAL[ORDER %in% c("Ephemeroptera", "Plecoptera", "Trichoptera") & TVAL <= 2], na.rm = TRUE)
      if (Total > 0) {
        (intolerant_EPT_count / Total) * 100 #Calculate percentage of intolerant EPT individuals
      } else {
        0 #avoids division by zero
      }
    },
    TOLERANT_CHIRO8_P_I = {
      #calculate the sum of TOTAL for tolerant chiro individuals with TVAL>=8
      tolerant_chiro_count <- sum(TOTAL[FAMILY %in% c("Chironomidae") & TVAL >= 8], na.rm = TRUE)
      if (Total > 0) {
        (tolerant_chiro_count / Total) * 100 #Calculate percentage of tolerant chiro individuals
      } else {
        0 #avoids division by zero
      }
    },
    ECO_FTN_UNIQUE = n_distinct(EcoFTN),  # Count of unique values in EcoFTN
    GATH_P_I = {
      #calculate the sum of TOTAL for in individuals in the Gatherer guild (TFUNCT contains "gatherer")
      gatherer_count <- sum(TOTAL[grepl("gatherer", TFUNCT)], na.rm = TRUE)
      if (Total > 0) {
        (gatherer_count / Total) * 100 #Calculate percentage of gatherers among all individuals
      } else {
        0 #avoids division by zero
      }
    },
    SCRAPE_P_I = {
      #calculate the sum of TOTAL for in individuals in the scraper guild (TFUNCT contains "scraper")
      scraper_count <- sum(TOTAL[grepl("scraper", TFUNCT)], na.rm = TRUE)
      if (Total > 0) {
        (scraper_count / Total) * 100 #Calculate percentage of gatherers among all individuals
      } else {
        0 #avoids division by zero
      }
    },
    .groups = 'drop' # drop grouping by station ID and year after summarizing
  )

#change year in ibi_calc_metrics and ibi_dup_metrics to char.
ibi_calc_metrics <- ibi_calc_metrics %>%
  mutate(YEAR = as.character(YEAR))

dup_ibi_calc_metrics <- dup_ibi_calc_metrics %>% 
  mutate(YEAR = as.character(YEAR))

# Join IBI.INSECT_T with ibi_calc_metrics from Cyrus' code.
ibi_calc_metrics <- ibi_calc_metrics %>%
  left_join(IBI.INSECT_T, by = c("STATION_ID", "YEAR")) %>%
  mutate(INSECT_T = coalesce(INSECT_T.y, INSECT_T.x)) %>%
  select(-INSECT_T.x, -INSECT_T.y)

# Join IBI.EPT_T with ibi_calc_metrics from Cyrus' code.
dup_ibi_calc_metrics <- dup_ibi_calc_metrics %>%
  left_join(IBI.EPT_T, by = c("STATION_ID", "YEAR")) %>%
  mutate(EPT_T = coalesce(EPT_T.y, EPT_T.x)) %>%
  select(-EPT_T.x, -EPT_T.y)

dup_ibi_calc_metrics <- dup_ibi_calc_metrics %>%
  left_join(IBI.INSECT_T, by = c("STATION_ID", "YEAR")) %>%
  mutate(INSECT_T = coalesce(INSECT_T.y, INSECT_T.x)) %>%
  select(-INSECT_T.x, -INSECT_T.y)


ibi_calc_metrics <- ibi_calc_metrics %>% select(STATION_ID, YEAR, INSECT_T, INSECT_PERCENT_I, EPT_T, DOM3_P_I, MPTV, INTOL_EPT2_P_I, TOLERANT_CHIRO8_P_I, ECO_FTN_UNIQUE, GATH_P_I, SCRAPE_P_I)

dup_ibi_calc_metrics <- dup_ibi_calc_metrics %>% select(STATION_ID, YEAR, INSECT_T, INSECT_PERCENT_I, EPT_T, DOM3_P_I, MPTV, INTOL_EPT2_P_I, TOLERANT_CHIRO8_P_I, ECO_FTN_UNIQUE, GATH_P_I, SCRAPE_P_I)

# Now, ibi_calc_metrics will have the updated INSECT_T and EPT_T columns

######## Scoring for individual metrics (in ibi_calc_metrics) I TWEAKED MANY OF THESE RANGES FROM THE PAPER... there were gaps in some of the ranges #########
# Poor (0)
# Fair (5)
# Good (10)
###### INSECT_T #####
# < 22 = 0
# 22-31 = 5
# > 31 = 10
###### INSECT_PERCENT_I #########
# < 90 = 0
# 90-95 = 5
# > 95 = 10
######## EPT_T ########
# < 6 = 0
# 6-15 = 5
# > 15 = 10
########  DOM3_P_I #######
# > 66 = 0
# 41-66 = 5
# < 40 = 10
####### MPTV ########
# > 6.440 = 0
# 5.876 - 6.440 = 5
# < 5.876 = 5
######## INTOL_EPT2_P_I ##########
# < 0.1 = 0
# 0.1 - 3 = 5
# >3 = 10
######## TOLERANT_CHIRO8_P_I #########
# > 16 = 0
# 2.4 - 16 = 5
# < 2.4 = 10
######## ECO_FTN_UNIQUE #########
# < 8 = 0
# 9 - 12 = 5
# > 12 = 12
######## GATH_P_I ############
# > 54 = 0
# 16-54 = 5
# < 16 = 10
######### SCRAPE_P_I ##########
# < 0.1 = 0
# 0.1 - 7.4 = 5
# > 7.4 = 10

ibi_score_parameters <- function(data) {
  data %>% 
        mutate(
          INSECT_T_Score = case_when(
          INSECT_T < 22 ~ 0,
          INSECT_T >= 22 & INSECT_T <= 31 ~5,
          INSECT_T >31 ~ 10
        ),
        INSECT_PERCENT_I_Score = case_when(
          INSECT_PERCENT_I < 90 ~ 0,
          INSECT_PERCENT_I >= 90 & INSECT_PERCENT_I <= 95 ~ 5,
          INSECT_PERCENT_I > 95 ~ 10
        ),
        EPT_T_Score = case_when(
          EPT_T < 6 ~ 0,
          EPT_T >= 6 & EPT_T <= 15 ~ 5,
          EPT_T > 15 ~ 10
        ),
        DOM3_P_I_Score = case_when(
          DOM3_P_I > 66 ~ 0,
          DOM3_P_I >= 41 & DOM3_P_I <= 66 ~ 5,
          DOM3_P_I < 41 ~ 10
        ),
        MPTV_Score = case_when(
          MPTV > 6.440 ~ 0,
          MPTV >= 5.876 & MPTV <= 6.440 ~ 5,
          MPTV < 5.876 ~ 10
        ),
        INTOL_EPT2_P_I_Score = case_when(
          INTOL_EPT2_P_I < 0.1 ~ 0,
          INTOL_EPT2_P_I >= 0.1 & INTOL_EPT2_P_I <= 3 ~ 5,
          INTOL_EPT2_P_I > 3 ~ 10
        ),
        TOLERANT_CHIRO8_P_I_Score = case_when(
          TOLERANT_CHIRO8_P_I > 16 ~ 0,
          TOLERANT_CHIRO8_P_I >= 2.4 & TOLERANT_CHIRO8_P_I <= 16 ~ 5,
          TOLERANT_CHIRO8_P_I < 2.4 ~ 10
        ),
        ECO_FTN_UNIQUE_Score = case_when(
          ECO_FTN_UNIQUE < 8 ~ 0,
          ECO_FTN_UNIQUE >= 8 & ECO_FTN_UNIQUE <= 12 ~ 5,
          ECO_FTN_UNIQUE > 12 ~ 10
        ),
        GATH_P_I_Score = case_when(
          GATH_P_I > 54 ~ 0,
          GATH_P_I >= 16 & GATH_P_I <= 54 ~ 5,
          GATH_P_I < 16 ~ 10
        ),
        SCRAPE_P_I_Score = case_when(
          SCRAPE_P_I < 0.1 ~ 0,
          SCRAPE_P_I >= 0.1 & SCRAPE_P_I < 7.4 ~ 5,
          SCRAPE_P_I >= 7.4 ~ 10
            )
        )
}

# Apply the scoring function to the ibi_calc_metrics dataframe
scored_metrics <- ibi_score_parameters(ibi_calc_metrics)

print(scored_metrics)
########### CALCULATE FINAL NUMERICAL IBI SCORES BY SUMMING SCORES OF ALL METRICS ################## 
scored_metrics <- scored_metrics %>% 
  mutate(FINAL_IBI_SCORE_NUM = (
    INSECT_T_Score + INSECT_PERCENT_I_Score + EPT_T_Score + DOM3_P_I_Score + MPTV_Score + INTOL_EPT2_P_I_Score + TOLERANT_CHIRO8_P_I_Score + ECO_FTN_UNIQUE_Score +
      GATH_P_I_Score +  SCRAPE_P_I_Score))

############ FINAL IBI Score Criteria ##############
# Overall score ranges for final ranking for each site (in the pub, listed as <=19; might have to revist these since some numbers may be left out based on code behavior, 
#since there are gaps in the ranges... possibly deal with by rounding.
# 0=Worst, 100=Best
# Very Poor <=19 
# Poor = >19-39
# Fair = >39-59
# Good = >59-80
# Excellent >=80
scored_metrics <- scored_metrics %>% 
  mutate(IBI_FINAL_SCORE_CHR = case_when(FINAL_IBI_SCORE_NUM <= 19 ~ "VERY POOR",
                                         FINAL_IBI_SCORE_NUM  > 19 & FINAL_IBI_SCORE_NUM <= 39 ~ "POOR",
                                         FINAL_IBI_SCORE_NUM > 39 & FINAL_IBI_SCORE_NUM <= 59 ~ "FAIR",
                                         FINAL_IBI_SCORE_NUM > 59 & FINAL_IBI_SCORE_NUM < 80 ~ "GOOD",
                                         FINAL_IBI_SCORE_NUM >= 80 ~ "EXCELLENT"))




# Calculate mean, min, and max for specified metrics
mean_metrics <- scored_metrics %>%
  group_by(STATION_ID) %>%
  summarise(
    INSECT_T_Score = mean(INSECT_T_Score, na.rm = TRUE),
    INSECT_PERCENT_I_Score = mean(INSECT_PERCENT_I_Score, na.rm = TRUE),
    EPT_T_Score = mean(EPT_T_Score, na.rm = TRUE),
    DOM3_P_I_Score = mean(DOM3_P_I_Score, na.rm = TRUE),
    MPTV_Score = mean(MPTV_Score, na.rm = TRUE),
    INTOL_EPT2_P_I_Score = mean(INTOL_EPT2_P_I_Score, na.rm = TRUE),
    TOLERANT_CHIRO8_P_I_Score = mean(TOLERANT_CHIRO8_P_I_Score, na.rm = TRUE),
    ECO_FTN_UNIQUE_Score = mean(ECO_FTN_UNIQUE_Score, na.rm = TRUE),
    GATH_P_I_Score = mean(GATH_P_I_Score, na.rm = TRUE),
    SCRAPE_P_I_Score = mean(SCRAPE_P_I_Score, na.rm = TRUE),
    .groups = 'drop'
  )

# Calculate min and max for error bars with unique names
range_metrics <- scored_metrics %>%
  group_by(STATION_ID) %>%
  summarise(
    INSECT_T_min = min(INSECT_T_Score, na.rm = TRUE),
    INSECT_T_max = max(INSECT_T_Score, na.rm = TRUE),
    INSECT_PERCENT_I_min = min(INSECT_PERCENT_I_Score, na.rm = TRUE),
    INSECT_PERCENT_I_max = max(INSECT_PERCENT_I_Score, na.rm = TRUE),
    EPT_T_min = min(EPT_T_Score, na.rm = TRUE),
    EPT_T_max = max(EPT_T_Score, na.rm = TRUE),
    DOM3_P_I_min = min(DOM3_P_I_Score, na.rm = TRUE),
    DOM3_P_I_max = max(DOM3_P_I_Score, na.rm = TRUE),
    MPTV_min = min(MPTV_Score, na.rm = TRUE),
    MPTV_max = max(MPTV_Score, na.rm = TRUE),
    INTOL_EPT2_P_I_min = min(INTOL_EPT2_P_I_Score, na.rm = TRUE),
    INTOL_EPT2_P_I_max = max(INTOL_EPT2_P_I_Score, na.rm = TRUE),
    TOLERANT_CHIRO8_P_I_min = min(TOLERANT_CHIRO8_P_I_Score, na.rm = TRUE),
    TOLERANT_CHIRO8_P_I_max = max(TOLERANT_CHIRO8_P_I_Score, na.rm = TRUE),
    ECO_FTN_UNIQUE_min = min(ECO_FTN_UNIQUE_Score, na.rm = TRUE),
    ECO_FTN_UNIQUE_max = max(ECO_FTN_UNIQUE_Score, na.rm = TRUE),
    GATH_P_I_min = min(GATH_P_I_Score, na.rm = TRUE),
    GATH_P_I_max = max(GATH_P_I_Score, na.rm = TRUE),
    SCRAPE_P_I_min = min(SCRAPE_P_I_Score, na.rm = TRUE),
    SCRAPE_P_I_max = max(SCRAPE_P_I_Score, na.rm = TRUE),
    .groups = 'drop'
  )

# Reshape the range metrics for plotting
range_metrics_long <- range_metrics %>%
  pivot_longer(cols = -STATION_ID, 
               names_to = c("Metric", ".value"), 
               names_pattern = "(.+)_(min|max)")

# Reshape mean metrics for all metrics
mean_metrics_long <- mean_metrics %>%
  pivot_longer(cols = -STATION_ID, names_to = "Metric", values_to = "Mean_Score")

mean_metrics_long %>%
  group_by(STATION_ID, Metric) %>%
  summarise(count = n()) %>%
  filter(count > 1)

range_metrics_long %>%
  group_by(STATION_ID, Metric) %>%
  summarise(count = n()) %>%
  filter(count > 1)

str(mean_metrics_long)
str(range_metrics_long)

unique(mean_metrics_long$STATION_ID)
unique(range_metrics_long$STATION_ID)

unique(mean_metrics_long$Metric)
unique(range_metrics_long$Metric)

combined_metrics_debug <- mean_metrics_long %>%
  full_join(range_metrics_long, by = c("STATION_ID", "Metric"))

sum(is.na(mean_metrics_long$STATION_ID))
sum(is.na(mean_metrics_long$Metric))
sum(is.na(range_metrics_long$STATION_ID))
sum(is.na(range_metrics_long$Metric))

mean_metrics_long$STATION_ID <- trimws(mean_metrics_long$STATION_ID)
range_metrics_long$STATION_ID <- trimws(range_metrics_long$STATION_ID)

mean_metrics_long$Metric <- trimws(mean_metrics_long$Metric)
range_metrics_long$Metric <- trimws(range_metrics_long$Metric)

# Remove "_score" from the Metric column in mean_metrics_long
mean_metrics_long <- mean_metrics_long %>%
  mutate(Metric = sub("_Score$", "", Metric))

unique(mean_metrics_long$Metric)
unique(range_metrics_long$Metric)


# Combine mean and range data
combined_metrics <- mean_metrics_long %>%
  left_join(range_metrics_long, by = c("STATION_ID", "Metric"))

# Create the horizontal bar plot for all metrics
ggplot(combined_metrics, aes(x = Mean_Score, y = Metric)) +
  geom_bar(stat = "identity", position = "dodge", aes(fill = STATION_ID), 
           width = 0.7) +
  geom_errorbar(aes(ymin = min, ymax = max), 
                width = 0.2, position = position_dodge(0.7)) +
  labs(title = "Mean Scores by Metric and Station",
       x = "Mean Score",
       y = "Metric") +
  theme_minimal() +
  facet_wrap(~ STATION_ID, scales = "free_y")  # Separate plots for each STATION_ID




# Define the output directory
output_dir <- "C:/Users/rdamstra/OneDrive - DOI/Desktop/GLKN RV WQ Invert R Code/plots"

# Loop through each unique STATION_ID
for (station_id in unique(combined_metrics$STATION_ID)) {
  
  # Filter data for the current STATION_ID
  station_data <- combined_metrics %>%
    filter(STATION_ID == station_id)
  
  # Create the plot
  p <- ggplot(station_data, aes(x = Mean_Score, y = Metric)) +
    geom_bar(stat = "identity", position = "dodge", aes(fill = STATION_ID), 
             width = 0.7) +
    geom_errorbar(aes(xmin = min, xmax = max), 
                  width = 0.2, position = position_dodge(0.7), color = "black") +  # Set error bars to black
    labs(title = paste("Mean Scores by Metric for Station:", station_id),
         x = "Mean Score",
         y = "Metric") +
    theme(
      panel.background = element_rect(fill = "white"),  # White background for the panel
      plot.background = element_rect(fill = "white"),   # White background for the entire plot
      panel.grid.major = element_blank(),  # Remove major gridlines
      panel.grid.minor = element_blank(),  # Remove minor gridlines
      axis.line = element_line(color = "black"),  # Black lines for x and y axes
      legend.position = "none",  # Position the legend
      text = element_text(size = 15)  # Set text size
    )
  
  # Define the file path
  file_path <- file.path(output_dir, paste0(station_id, ".png"))
  
  # Save the plot to a PNG file with a white page background
  ggsave(filename = file_path, plot = p, width = 10, height = 6, dpi = 300, bg = "white")
}


# Calculate mean, min, and max for FINAL_IBI_SCORE_NUM by STATION_ID
summary_metrics <- scored_metrics %>%
  group_by(STATION_ID) %>%
  summarise(
    Mean_Score = mean(FINAL_IBI_SCORE_NUM, na.rm = TRUE),
    Min_Score = min(FINAL_IBI_SCORE_NUM, na.rm = TRUE),
    Max_Score = max(FINAL_IBI_SCORE_NUM, na.rm = TRUE)
  )

# Define the desired order for STATION_ID
station_order <- c("SACN_NAKA_84.6", "SACN_NAKA_41.3", "SACN_NAKA_4.8", 
                   "SACN_STCR_138.9", "SACN_STCR_104.0", "SACN_STCR_89.7", 
                   "SACN_STCR_63.8", "SACN_STCR_51.9", "SACN_STCR_43.7", 
                   "SACN_STCR_26.0_A", "SACN_STCR_26.0_B")

# Convert STATION_ID to a factor with the specified order
summary_metrics <- summary_metrics %>%
  mutate(STATION_ID = factor(STATION_ID, levels = station_order))

# Assign colors based on Mean_Score criteria
summary_metrics <- summary_metrics %>%
  mutate(Color = case_when(
    Mean_Score <= 19 ~ "Very Poor",
    Mean_Score > 19 & Mean_Score <= 39 ~ "Poor",
    Mean_Score > 39 & Mean_Score <= 59 ~ "Fair",
    Mean_Score > 59 & Mean_Score <= 80 ~ "Good",
    Mean_Score > 80 ~ "Excellent"
  ))

# Create a legend data frame
legend_data <- data.frame(
  Category = c("Very Poor", "Poor", "Fair", "Good", "Excellent"),
  Color = c("red", "yellow", "green", "blue", "violet")
)

# Create the vertical bar plot with error bars
p <- ggplot(summary_metrics, aes(x = STATION_ID, y = Mean_Score, fill = Color)) +
  geom_bar(stat = "identity", width = 0.7) +
  geom_errorbar(aes(ymin = Min_Score, ymax = Max_Score), 
                width = 0.2, color = "black") +
  labs(title = "Mean FINAL IBI Score by Station",
       x = "Station ID",
       y = "Mean FINAL IBI Score") +
  scale_fill_manual(values = setNames(legend_data$Color, legend_data$Category)) +  # Use the legend data for fill colors
  theme(
    panel.background = element_rect(fill = "white"),  # White background for the panel
    plot.background = element_rect(fill = "white"),   # White background for the entire plot
    panel.grid.major = element_blank(),  # Remove major gridlines
    panel.grid.minor = element_blank(),  # Remove minor gridlines
    axis.line = element_line(color = "black"),  # Black lines for x and y axes
    legend.position = c(0.9, 0.9),  # Position the legend on the top right
    legend.title = element_blank(),  # Remove legend title
    axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels for better readability
    text = element_text(size = 15)  # Set text size
  ) +
  guides(fill = guide_legend(title = "Categories"))  # Add title to the legend

# Print the plot
print(p)

# Optionally, save the plot
ggsave("mean_final_ibi_score_by_station.png", plot = p, width = 10, height = 6, dpi = 300, bg = "white")

# Optionally, save the plot
ggsave("mean_final_ibi_score_by_station.png", plot = p, width = 10, height = 6, dpi = 300, bg = "white")

# Optionally, save the plot
ggsave("C:/Users/rdamstra/OneDrive - DOI/Desktop/GLKN RV WQ Invert R Code/plots/mean_final_ibi_score_by_station.png", plot = p, width = 10, height = 6, dpi = 300, bg = "white")

# chiros only
chiros <- i %>% 
  filter(grepl("Chironomidae", FAMILY, ignore.case = TRUE)) %>% 
  select(STATION_ID, TAXON_ID_EXPERT, TAXON_ID_NOTES, TSN_LOOKUP, TSN_SOURCE, TSN, TVAL) %>% 
  group_by(STATION_ID, TAXON_ID_EXPERT) %>%
  distinct() %>% 
  group_by(STATION_ID, TAXON_ID_EXPERT) %>% 
  summarise(
    tsn_id = first(TSN_LOOKUP), 
    taxonomist_notes = first(TAXON_ID_NOTES),
    tsn_source = first(TSN_SOURCE),
    tsn_number = first(TSN),
    tolerance_value = first(TVAL),
    ,.groups = 'drop'
  )

unique_count <- n_distinct(chiros$TAXON_ID_EXPERT, chiros$taxonomist_notes)

print(unique_count)

write.csv(chiros, "C:/Users/rdamstra/OneDrive - DOI/Desktop/GLKN RV WQ Invert R Code/unique_chiros_2016-2023.csv")





