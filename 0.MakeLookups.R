#devtools::load_all()

### 1. Install and load required packages
# The following tidyverse packages need to be installed (but do not need to be loaded): dplyr, haven, purrr, stringr, labelled, simPop
# The magrittr package should also be installed. Once installed, it should be loaded:
library(magrittr)

### 2. Load all functions
read_fns <- function(fns_dir_chr = "data-raw/fns/"){

  fns_path_chr_vec <- list.files(fns_dir_chr,

                                 pattern="*.R$",

                                 full.names=TRUE,

                                 ignore.case=TRUE) %>%

    purrr::walk(~source(.x))

  return(fns_path_chr_vec)

}

fns_path_chr_vec <- read_fns()

## 3. Set up blank data import lookup tibble
db_var_name_matches_lup <- readr::read_csv("Data/lookup_variable_names.csv")

look_up_column_for_bin <- readr::read_csv("Data/look_up_column_for_bin.csv")

LSAC_lookup_parent_variables <- readr::read_csv("Data/LSAC_lookup_parent_variables.csv")

empty_sp_data_import_tb <- data_import_make_sp_lookup_tb()

## 4. Add data to data import lookup tibble
sp_data_import_read_tb <- read.csv("data-raw/sp_import_data.csv")  %>%
  tibble::as.tibble()
sp_data_import_read_tb <- sp_data_import_read_tb %>%
  cbind(list(inc_files_to_rename = rep(NA_character_, nrow(sp_data_import_read_tb)))) %>%
  cbind(list(new_names_for_inc_files = rep(NA_character_, nrow(sp_data_import_read_tb)))) %>%
  dplyr::mutate(inc_files_to_rename = ifelse(name=="aus_sa2_nat_cor_2011_2016",
                                             "CG_SA2_2011_SA2_2016.xls",
                                             inc_files_to_rename),
                new_names_for_inc_files = ifelse(name=="aus_sa2_nat_cor_2011_2016",
                                                 "sa2correspondence.xls",
                                                 new_names_for_inc_files))
sp_data_import_tb <- data_import_add_tb_rows_to_sp_tb(sp_data_import_tb,
                                                      sp_data_import_read_tb ,
                                                      1:nrow(sp_data_import_read_tb))

## 5. Add Parental Variables Lookup Table
LSAC_lookup_parent_variables <- readr::read_csv("~/Readyforwhatsnext/DATASETS/LSAC_lookup_parent_variables.csv")

# save data lookup table as internal file
sp_data_import_tb <- dplyr::arrange(sp_data_import_tb,
                                    name)
## 6. Import census data for sex by state/territory and create data table for calibration
auspop_sex_state_2009 <- readxl::read_excel("Data/auspop_sex_state_2009.xls",
                                            sheet = "Table_8", skip=4)

#filter out age-specific data and retain only totals by sex
auspop_sex_age_state <-  auspop_sex_state_2009 %>%
  dplyr::filter(`Age (years)` %in% c("0–4",
                                     "5–9",
                                     "10–14",
                                     "15–19",
                                     "20–24",
                                     "25–29",
                                     "30–34",
                                     "35–39",
                                     "40–44",
                                     "45–49",
                                     "50–54",
                                     "55–59",
                                     "60–64",
                                     "65–69",
                                     "70–74",
                                     "75–79",
                                     "80–84",
                                     "85–89",
                                     "90–94",
                                     "95–99",
                                     "100 and over"))

# Rename and convert to factor
auspop_sex_age_state <-  auspop_sex_age_state %>%
  dplyr::rename(rfwn_age_cat = `Age (years)`) %>%
  dplyr::mutate(rfwn_age_cat = factor(rfwn_age_cat))

# Drop total persons rows ; create and bind sex column
auspop_sex_age_state_male <- auspop_sex_age_state[1:21, 1:9]
auspop_sex_age_state_female <- auspop_sex_age_state[22:42, 1:9]
auspop_sex_age_state_male <- cbind(data.frame(rfwn_demographics_sex="male"), auspop_sex_age_state_male)
auspop_sex_age_state_female <- cbind(data.frame(rfwn_demographics_sex="female"), auspop_sex_age_state_female)
auspop_sex_age_state <- rbind(auspop_sex_age_state_female, auspop_sex_age_state_male)

# Reshape into long data, gathering by State/Tty
auspop_sex_age_state_2009 <- auspop_sex_age_state %>%
  tidyr::gather(key = "rfwn_spatial_home_state", value="freq", -rfwn_demographics_sex, -rfwn_age_cat)

# Recode state/tty into factored numeric code matching HILDA (NB. same code for LSAC data)
auspop_sex_age_state_2009$rfwn_spatial_home_state <-  factor(dplyr::recode(auspop_sex_age_state_2009$rfwn_spatial_home_state,
                                                                      "New South Wales" = 1,
                                                                      "Victoria" = 2,
                                                                      "Queensland"= 3,
                                                                      "South Australia"= 4,
                                                                      "Western Australia"= 5,
                                                                      "Tasmania"= 6,
                                                                      "Northern Territory"= 7,
                                                                      "Australian Capital Territory"= 8))

#recode sex into numeric factors as per HILDA baseline microdata
auspop_sex_age_state_2009$rfwn_demographics_sex <-  factor(dplyr::recode(auspop_sex_age_state_2009$rfwn_demographics_sex,
                                                                    "male" = 1,
                                                                    "female" = 2))
## 7. Import ABS ERP data and set up ERP by LGA table as per simPop requirements
ERP_LGA_2009 <- readr::read_csv("Data/ABS_ERP_LGA2016_03092019155543058.csv")

#Select columns of interest and rename; mutate State_Tty col
ERP_LGA_2009 <- ERP_LGA_2009 %>%
  dplyr:: select(LGA_2016, Value) %>%
  dplyr::rename(LGA_code_2016 = LGA_2016, Freq = Value) %>%
  dplyr::mutate(State_Tty = stringr::str_sub(LGA_code_2016, 1,1)) %>%
  #filter out special LGA codes with 0 residents and state/tty totals
  dplyr::filter(Freq>0) %>%
  dplyr::filter(LGA_code_2016>10) %>%
  dplyr:: select(State_Tty, LGA_code_2016, Freq) %>%
  # reclass vars as factors
  dplyr::mutate(State_Tty = factor(State_Tty),
                LGA_code_2016 = factor(LGA_code_2016))


