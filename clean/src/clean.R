#!/usr/bin/env Rscript --vanilla
# set expandtab ft=R ts=4 sw=4 ai fileencoding=utf-7
#
# Author: JR
# Maintainer(s): JR
# License:  2020, GPL v3 or later
#
# -----------------------------------------------------------
# epimodel_practice/clean/src/clean.R

pacman::p_load("tidyverse", "lubridate", "readr", 
               "here", "assertr", "janitor", "forcats")

files <- list(
  cen_data = here("epimodel_practice/clean/input/QuickFacts Jun-29-2020.csv"),
  ga_countycases = here("epimodel_practice/clean/input/countycases.csv"), 
  ga_deaths = here("epimodel_practice/clean/input/deaths.csv"), 
  ga_demo = here("epimodel_practice/clean/input/demographics.csv"), 
  
  cen_clean = here("epimodel_practice/model/input/cendata_clean.csv"), 
  ga_cc_clean = here("epimodel_practice/model/input/countycases_clean.csv"), 
  ga_d_clean = here("epimodel_practice/model/input/deaths_clean.csv"), 
  ga_demo_clean = here("epimodel_practice/model/input/demo_clean.csv")
  
)

stopifnot(length(files) == 8)

# cen_data come from here https://www.census.gov/quickfacts/fact/table/GA/PST045219

cen_df <- as.data.frame(read_delim(files$cen_data, 
                                    delim = ",")) %>%
  clean_names() %>%
  select(-c("value_note_for_georgia", "fact_note")) %>%
  drop_na(georgia) %>%
  mutate(ga_num = as.numeric(parse_number(georgia, trim_ws = TRUE)))

# unit tests

cen_df <- cen_df %>%
  verify(ncol(cen_df) == 3 & (nrow(cen_df) == 62)) %>%
  verify(is.na(ga_num) == FALSE & ga_num[1] == 10617423.00) %>%
  write_delim(files$cen_clean, delim = "|")

#done  