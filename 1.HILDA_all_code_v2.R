## For this script to execute correctly:
## The db_var_name_matches_lup object created in 0.MakeLookups must be loaded 
## The functions from https://github.com/ready4-dev/fakefolk must be loaded.
##
## Create included variables list
# Create a vector that specifies the generic (non-wave specific) variable names
# Note: vector should include the cross wave identifier "xwaveid"
included_var_list <- c("xwaveid",
                       "hhstate",
                       "hhsos",
                       "hhlga",
                       "hhslga",
                       "hhssa2",
                       "hhsad10",
                       "hhmvehk",
                       "mvncar",
                       "hgage",
                       "anatsi",
                       "anlote",
                       "anengf",
                       "anbcob",
                       "hgsex",
                       "lssexor",
                       "hhrhid",
                       "hhtype",
                       "hhfty",
                       "hhrih",
                       "hhmxid",
                       "hhfxid",
                       "hhpxid",
                       "mrcurr",
                       "tcr",
                       "tcnr",
                       "hstenr",
                       "hsllord",
                       "hsfrea",
                       "edceei",
                       "caeft",
                       "caept",
                       "edhists",
                       "edhigh1",
                       "lshremp",
                       "wscei",
                       "hifeftn",
                       "hifeftp",
                       "esbrd",
                       "esdtl",
                       "fmfo6s",
                       "fmmo6s",
                       "jbmo6s",
                       "pjoto6s",
                       "ujljo6s",
                       "bncoth1",
                       "bnfhave",
                       "bncalli",
                       "losat",
                       "ghsf6an",
                       "ghsf6ap",
                       "hglth",
                       "helth",
                       "heany",
                       "lejls",
                       "lejlf",
                       "pdk10s",
                       "henec",
                       "hemirh",
                       "hedep",
                       "heomi",
                       "lshedep",
                       "pdk10rc",
                       "lsdrex",
                       "lspact",
                       "apcat",
                       "fmpdiv",
                       "ledsc",
                       "psyodf",
                       "psyodm",
                       "lesep",
                       "lslaha",
                       "ffcake",
                       "ffconf",
                       "ffprocm",
                       "ffsnack",
                       "ffspud"
)

## Create merged dataset (all waves)
  # Note: specify required file extension type
HILDA_merged_dataset <- read_merge_waves_HILDA(waves = c(1:17),
                                   path_to_source = "~/data_projects/confidential_data/HILDA/STATA 170u_Combined Data Files/", # UPDATE THIS TO PATH ON YOUR MACHINE
                                   included_vars_generic = included_var_list,
                                   included_vars_xwave = "xwaveid",
                                   data_type = "Combined_",
                                   file_extension  = ".dta"
)

# ## Create vector of ids to use if balancing dataset
# ## Note: specify required file extension type
# balanced_ids <- balanced_dataset_ids(path_to_source = "~/data_projects/confidential_data/HILDA/STATA 170u_Other Data Files/",
#                                      file_extension = ".dta")
#
# ## Create balanced dataset
# HILDA_merged_dataset_balanced <- HILDA_merged_dataset %>%
#   dplyr::filter(xwaveid %in% balanced_ids)

## Import, select and append HILDA death variables from wave 17 Master file
Master_q170u <- haven::read_dta("~/data_projects/confidential_data/HILDA/STATA 170u_Other Data Files/Master_q170u.dta") # UPDATE TO PATH ON YOUR MACHINE

master_w17_death <- Master_q170u[c("xwaveid", "yodeath", "aadeath" )]

HILDA_merged_dataset <- dplyr::left_join(HILDA_merged_dataset, master_w17_death, by="xwaveid")

## Create single wave raw baseline dataset
waveI <- HILDA_merged_dataset %>%
  dplyr::select(xwaveid, dplyr::starts_with("i"),
                #Set of ACE-related variables from all preceding waves
                amrcurr,
                bmrcurr,
                cmrcurr,
                dmrcurr,
                emrcurr,
                fmrcurr,
                gmrcurr,
                hmrcurr,
                glsdrex,
                blejlf,
                clejlf,
                dlejlf,
                elejlf,
                flejlf,
                glejlf,
                hlejlf,
                blejls,
                clejls,
                dlejls,
                elejls,
                flejls,
                glejls,
                hlejls,
                alslaha,
                blslaha,
                clslaha,
                dlslaha,
                flslaha,
                hlslaha,
                hpsyodf,
                hpsyodm)

# Drop observations without data (i.e. survey dropouts)
waveI <- waveI %>%
  dplyr::filter(!is.na(ihhrhid))

##  Recode missing values for numeric columns
# Note: The following function only recodes numeric values. Factor and boolean variables imported from STATA files appear as numeric variables once imported into R.
HILDA_wave_i_clean <- rc_miss_num(waveI)

# Recode missing values for character columns (i.e. mother, father, partner ID columns)
HILDA_wave_i_clean <-  HILDA_wave_i_clean %>%
  dplyr::mutate(ihhmxid = ifelse(ihhmxid == "",
                          NA_character_,
                          ihhmxid)) %>%
  dplyr::mutate(ihhfxid = ifelse(ihhfxid == "",
                          NA_character_,
                          ihhfxid)) %>%
  dplyr::mutate(ihhpxid = ifelse(ihhpxid == "",
                          NA_character_,
                          ihhpxid))

## Generating derived variables not requiring targeting of parental variables:
HILDA_wave_i_bl  <-  HILDA_wave_i_clean %>%
  dplyr::select(xwaveid, ihhrhid, ihhfxid, ihhmxid, ihhpxid, dplyr::everything()) %>%
  dplyr::mutate(demographics_CALD_status = ifelse((ianlote==1 | ianbcob==3),
                                                  TRUE,
                                                  FALSE)) %>%
  dplyr::mutate(demographics_ATSI_status = ifelse(is.na(ianatsi),
                                                  NA,
                                                  ifelse(ianatsi %in% c(2:4),
                                                         TRUE,
                                                         FALSE))) %>%
  dplyr::mutate(demographics_ATSI_status = ifelse(ianbcob %in% c(2:3),
                                                  FALSE,
                                                  demographics_ATSI_status)) %>%

  #dplyr::mutate(demographics_same_sex_attracted = ifelse((llssexor==2 | llssexor==3),TRUE,FALSE)) %>%
  dplyr::mutate(socio_econ_accom_status_bl= ifelse(((ihstenr %in% c(1,3)) |
                                                      (ihsllord %in% c(1,3)) |
                                                      ihsfrea %in% 8),
                                                      "secure",
                                                      "at-risk")) %>%
  dplyr::mutate(socio_econ_education_status_bl= ifelse((icaept==0 & icaeft==0),
                                                              "No",
                                                       ifelse(icaept==1,
                                                              "PT",
                                                              "FT"))) %>%
  dplyr::mutate(socio_econ_welfare_bl = ifelse((ibncoth1==1 | ibnfhave==1),
                                               TRUE,
                                               FALSE))

# Calculate HH gross income and its standardised z-score (for household SES calculation)
HILDA_wave_i_bl  <-  HILDA_wave_i_bl %>%
  # Calculate HH gross income and recode ~45 observations with negative income to 0
  dplyr::mutate(socio_econ_household_gross_income_bl = ifelse(ihifeftp - ihifeftn >=0,
                                                              ihifeftp - ihifeftn,
                                                              0)) %>%
  # Cube root-transform raw score to more normal distribution and convert to standardised Z-score
  dplyr::mutate(socio_econ_household_gross_income_bl_tf = (socio_econ_household_gross_income_bl)^(1/3)) %>%
  dplyr::mutate(household_income_z = (socio_econ_household_gross_income_bl_tf - mean(socio_econ_household_gross_income_bl_tf, na.rm=TRUE))/sd(socio_econ_household_gross_income_bl_tf, na.rm=TRUE))

# Calculate highest occupational prestige at HH level and standardised z-score (for household SES calculation)
  # Note: evaluates current job or most recent job if unemployed or NILF
  # Note: this may be inaccurate in group or multifamily households which are expected to be very small minority
HILDA_wave_i_bl  <-  HILDA_wave_i_bl %>%
  dplyr::mutate(occ_status = pmax(ijbmo6s, ipjoto6s, iujljo6s, na.rm=TRUE)) %>%
  dplyr::group_by(ihhrhid) %>%
  dplyr::mutate(household_occ_status = ifelse(all(is.na(occ_status)), # otherwise -Inf generated if all NULL
                                              NA,
                                              max(occ_status, na.rm = TRUE))) %>%
  dplyr::ungroup() %>%
  # Convert raw occupational prestige score to standardised Z-score
  dplyr::mutate(household_occ_z = (household_occ_status - mean(household_occ_status, na.rm=TRUE))/sd(household_occ_status, na.rm=TRUE))

# Calculate average z-score for HH-income and occupational prestige and convert to factor based on z-score quintiles
HILDA_wave_i_bl <-  HILDA_wave_i_bl %>%
  dplyr::rowwise() %>%
  dplyr::mutate(socio_econ_household_ses  = ifelse(mean(c(household_occ_z, household_income_z), na.rm=TRUE) <= -0.842,
                                                   "Q1",
                                                   ifelse(mean(c(household_occ_z, household_income_z), na.rm=TRUE) <= -0.253,
                                                          "Q2",
                                                          ifelse(mean(c(household_occ_z, household_income_z), na.rm=TRUE) <= 0.253,
                                                                 "Q3",
                                                                 ifelse(mean(c(household_occ_z, household_income_z), na.rm=TRUE) <= 0.842,
                                                                        "Q4",
                                                                        "Q5"))))) %>%
  dplyr::ungroup()

HILDA_wave_i_bl <-  HILDA_wave_i_bl %>%
  dplyr::mutate(socio_econ_incarcerated_hh_mbr_bl  = ifelse((ilejlf==2  | ilejls==2),
                                                            TRUE,
                                                            FALSE)) %>%
  dplyr::mutate(socio_econ_ever_incarcerated_hh_mbr = ifelse((ilejlf==2 |
                                                                ilejls==2 |
                                                                blejlf==2 |
                                                                clejlf==2 |
                                                                dlejlf==2 |
                                                                elejlf==2 |
                                                                flejlf==2 |
                                                                glejlf==2 |
                                                                hlejlf==2 |
                                                                blejls==2 |
                                                                clejls==2 |
                                                                dlejls==2 |
                                                                elejls==2 |
                                                                flejls==2 |
                                                                glejls==2 |
                                                                hlejls==2),
                                                             TRUE,
                                                             FALSE)) %>%
  dplyr::mutate(mental_health_curr_prev_adult_md_bl = ifelse(ihgage < 18, NA,
                                                             ifelse((ihenec %in% c(0,NA) &
                                                                       ihemirh %in% c(0,NA) &
                                                                       ihedep %in% c(0,NA) &
                                                                       iheomi %in% c(0,NA)),FALSE, TRUE))) %>%
  dplyr::mutate(mental_health_curr_prev_adult_md_bl = ifelse(is.na(ilshedep),
                                                             mental_health_curr_prev_adult_md_bl,
                                                             ifelse(mental_health_curr_prev_adult_md_bl==FALSE & ilshedep==1,
                                                                    TRUE,
                                                                    mental_health_curr_prev_adult_md_bl))) %>%
  dplyr::mutate(mental_health_substance_abuse_bl =ifelse(ilsdrex >5,
                                                         TRUE,
                                                         FALSE)) %>%
  dplyr::mutate(mental_health_ever_substance_abuse =ifelse((ilsdrex >5 |
                                                              glsdrex >5),
                                                           TRUE,
                                                           FALSE)) %>%
  dplyr::mutate(mental_health_mental_disorder_bl = ifelse(ipdk10rc == 1,
                                                          "No/Low",
                                                          ifelse(ipdk10rc == 2,
                                                                 "Subthreshold",
                                                                 "Disorder"))) %>%
  dplyr::mutate(health_utility_sf_6d_bl  = ighsf6ap - ighsf6an) %>%
  dplyr::mutate(ace_low_ses_bl = ifelse(ihgage >=18,
                                        NA,
                                        ifelse(socio_econ_household_ses == "Q1",
                                               TRUE,
                                               FALSE))) %>%
  dplyr::mutate(family_ever_widowed = ifelse((amrcurr==5 |
                                                bmrcurr==5 |
                                                cmrcurr==5 |
                                                dmrcurr==5 |
                                                emrcurr==5 |
                                                fmrcurr==5 |
                                                gmrcurr==5 |
                                                hmrcurr==5 |
                                                imrcurr==5),
                                             TRUE,
                                             FALSE)) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(ace_exposure_community_violence_bl = ifelse(ihgage >=18,
                                                            NA,
                                                            ifelse(all(is.na(alslaha),
                                                                       is.na(blslaha),
                                                                       is.na(clslaha),
                                                                       is.na(dlslaha),
                                                                       is.na(flslaha),
                                                                       is.na(hlslaha)),
                                                                   NA,
                                                                   any(alslaha >3,
                                                                       blslaha >3,
                                                                       clslaha >3,
                                                                       dlslaha >3,
                                                                       flslaha >3,
                                                                       hlslaha >3,
                                                                       na.rm=TRUE)))) %>%
  # Estimate history of MI in total 15+ respondent population (i.e. not just adults; does not include history of SU)
  dplyr::mutate(rpf_psychosocial_history_poor_mh_bl = ifelse(all(is.na(ihemirh),
                                                                 is.na(ihedep),
                                                                 is.na(iheomi),
                                                                 is.na(ilshedep)),
                                                             NA,
                                                             any(ihemirh==1,
                                                                 ihedep==1,
                                                                 iheomi==1,
                                                                 ilshedep==1, na.rm=TRUE)))  %>%
  dplyr::ungroup() %>%
  dplyr::mutate(rpf_physiological_disability_illness_bl = ifelse(ihgage <15,
                                                                 ifelse(ihglth==1,
                                                                        TRUE,
                                                                        FALSE),
                                                                 ifelse((ihelth==1 | iheany==1),
                                                                        TRUE,
                                                                        FALSE)))%>%
  dplyr::mutate(rpf_lifestyle_low_physical_act_bl= ifelse((ilspact<4),
                                                          TRUE,
                                                          FALSE)) %>%
  #Calculate raw score for poor diet (unhealthy food frequency)
  dplyr::mutate(diet_score = iffcake + iffconf + iffprocm + iffsnack + iffspud) %>%
  # Convert raw score to standardised Z-score
  dplyr::mutate(diet_score_z = (diet_score - mean(diet_score, na.rm=TRUE))/sd(diet_score, na.rm=TRUE)) %>%
  #Calculate those in lowest quintile as at risk of poor diet
  dplyr::mutate(rpf_lifestyle_poor_diet_bl = ifelse(diet_score_z <= -0.842,
                                                    TRUE,
                                                    FALSE))

##Scripts to generate variables targeting parent data
  # NB.order important for logic sequence
HILDA_wave_i_bl  <-  HILDA_wave_i_bl %>%
  dplyr::mutate(demographics_country_of_birth = ifelse(ihgage>=15,
                                                       ianbcob,
                                                       gen_cob_pop(db = .,
                                                                   target_variable = "ianbcob",
                                                                   agt_age_vect = ihgage,
                                                                   moth_id_vect = ihhmxid,
                                                                   fath_id_vect = ihhfxid,
                                                                   agt_id_col = "xwaveid")))
HILDA_wave_i_bl  <-  HILDA_wave_i_bl %>%
  dplyr::mutate(demographics_CALD_status = ifelse(ihgage>=15,
                                                  demographics_CALD_status,
                                                  gen_parental_CALD_pop(db = .,
                                                                        target_variable = "demographics_CALD_status",
                                                                        agt_age_vect = ihgage,
                                                                        moth_id_vect = ihhmxid,
                                                                        fath_id_vect = ihhfxid,
                                                                        agt_id_col = "xwaveid")))

HILDA_wave_i_bl  <-  HILDA_wave_i_bl %>%
  dplyr::mutate(demographics_ATSI_status = ifelse(ihgage>=15,
                                                  demographics_ATSI_status,
                                                  gen_parental_ATSI_pop(db = .,
                                                                        target_variable = "demographics_ATSI_status",
                                                                        agt_age_vect = ihgage,
                                                                        moth_id_vect = ihhmxid,
                                                                        fath_id_vect = ihhfxid,
                                                                        agt_id_col = "xwaveid")))

HILDA_wave_i_bl  <-  HILDA_wave_i_bl %>%
  dplyr::mutate(demographics_generation = ifelse(demographics_country_of_birth %in% c(2,3),
                                                 "FIRST",
                                                 ifelse(is.na(ihhmxid) & is.na(ihhfxid),
                                                        NA,
                                                        gen_ausgen_pop(db = .,
                                                                       target_variable = "demographics_country_of_birth",
                                                                       agt_id_vect = xwaveid,
                                                                       moth_id_vect = ihhmxid,
                                                                       fath_id_vect = ihhfxid,
                                                                       agt_id_col = "xwaveid"))))

HILDA_wave_i_bl  <-  HILDA_wave_i_bl %>%
  dplyr::mutate(ace_parental_SU_bl = gen_parental_SU_pop(db = .,
                                                         target_variable = "mental_health_ever_substance_abuse",
                                                         agt_age_vect = ihgage,
                                                         moth_id_vect = ihhmxid,
                                                         fath_id_vect = ihhfxid,
                                                         agt_id_col = "xwaveid"))


HILDA_wave_i_bl  <-  HILDA_wave_i_bl %>%
  dplyr::mutate(ace_parental_MI_bl = gen_parental_MI_pop(db = .,
                                                         target_variable = "mental_health_curr_prev_adult_md_bl",
                                                         agt_age_vect = ihgage,
                                                         moth_id_vect = ihhmxid,
                                                         fath_id_vect = ihhfxid,
                                                         agt_id_col = "xwaveid"))

HILDA_wave_i_bl  <-  HILDA_wave_i_bl %>%
  dplyr::mutate(ace_parental_death_bl = ifelse(ihgage>=18,
                                               NA,
                                               ifelse(ihgage<16,
                                                      gen_parental_death_pop(db = .,
                                                                             target_variable = "family_ever_widowed",
                                                                             agt_age_vect = ihgage,
                                                                             moth_id_vect = ihhmxid,
                                                                             fath_id_vect = ihhfxid,
                                                                             agt_id_col = "xwaveid"),
                                                      ifelse((HILDA_wave_i_bl$ihgage %in% c(16:17) & (!is.na(HILDA_wave_i_bl$hpsyodf) | !is.na(HILDA_wave_i_bl$hpsyodm))),
                                                             TRUE,
                                                             FALSE))))

HILDA_wave_i_bl  <-  HILDA_wave_i_bl %>%
  dplyr::mutate(ace_incarcerated_household_mbr_bl = gen_incarcerated_hh_mbr_pop(db = .,
                                                                                target_variable = "socio_econ_ever_incarcerated_hh_mbr",
                                                                                agt_age_vect = ihgage,
                                                                                moth_id_vect = ihhmxid,
                                                                                fath_id_vect = ihhfxid,
                                                                                agt_id_col = "xwaveid"))

##  Generate overall 'ACE present' variable
HILDA_wave_i_bl <-  HILDA_wave_i_bl %>%
  dplyr::rowwise() %>%
  dplyr::mutate(ace_present_bl = ifelse(all(is.na(ace_parental_MI_bl),
                                            is.na(ace_parental_SU_bl),
                                            is.na(ace_incarcerated_household_mbr_bl),
                                            is.na(ace_parental_death_bl),
                                            is.na(ace_low_ses_bl),
                                            is.na(ace_exposure_community_violence_bl)),
                                        NA,
                                        any(ace_parental_MI_bl,
                                            ace_parental_SU_bl,
                                            ace_incarcerated_household_mbr_bl,
                                            ace_parental_death_bl,
                                            ace_low_ses_bl,
                                            ace_exposure_community_violence_bl,
                                            na.rm=TRUE))) %>%
  dplyr::ungroup()

## Script for generating sibling IDs and child IDs (list data type)
  # NB. First need to strip STATA labels and attributes which interfere with 'character(0)' recoding
HILDA_wave_i_bl <- labelled::remove_labels(HILDA_wave_i_bl)
HILDA_wave_i_bl <- labelled::remove_attributes(HILDA_wave_i_bl, "format.stata")

HILDA_wave_i_bl  <-  HILDA_wave_i_bl %>%
  dplyr::mutate(family_sibling_id = get_sibling_id_pop(db = .,
                                                       target_variable = "xwaveid",
                                                       agt_id_vec = xwaveid,
                                                       hh_id_col = "ihhrhid",
                                                       agt_rel_type_col = "ihhrih",
                                                       agt_id_col = "xwaveid"))

HILDA_wave_i_bl  <-  HILDA_wave_i_bl %>%
  dplyr::mutate(family_child_id = get_child_id_pop(db = .,
                                                   target_variable = "xwaveid",
                                                   agt_id_vec = xwaveid,
                                                   hh_id_col = "ihhrhid",
                                                   agt_rel_type_col = "ihhrih",
                                                   agt_id_col = "xwaveid"))

## Script for determining parent and sibling of child with ACE
HILDA_wave_i_bl  <-  HILDA_wave_i_bl %>%
  dplyr::mutate(family_sibling_of_child_with_ace  = purrr::map_lgl(family_sibling_id,
                                                                   ~ sibling_ace_ind(data_tb = HILDA_wave_i_bl ,
                                                                                     sibling_ids = .x[1])))

HILDA_wave_i_bl  <-  HILDA_wave_i_bl %>%
  dplyr::mutate(family_parent_of_child_with_ace  = purrr::map_lgl(family_child_id,
                                                                  ~ parent_ace_ind(data_tb = HILDA_wave_i_bl ,
                                                                                   child_ids = .x)))


## 8. Rename original HILDA names to rfwn variable names:
HILDA_wave_i_bl_final <- HILDA_wave_i_bl %>%
  rename_to_rfwn_vars(lookup_tb = db_var_name_matches_lup,
                      old_var_name_col = "HILDA_var_name",
                      old_var_stub_col = "HILDA_var_stub",
                      wave_prefix = "i")

## Create weights for HILDA households/individuals and join to baseline dataset
# Import 'eperson' file for wave I and select weight vars
enumerated_i <- haven::read_dta("~/data_projects/confidential_data/HILDA/STATA 170u_Other Data Files/Eperson_i170u.dta")

i_weights_ind <- enumerated_i %>%
  dplyr::select(xwaveid, ihhrhid, ihhwte, ihhwtes)

HILDA_wave_i_bl_final <- dplyr::left_join(HILDA_wave_i_bl_final, i_weights_ind, by  = c("xwaveid", "rfwn_family_household_id" = "ihhrhid"))

# rename weight vars to rfwn style
HILDA_wave_i_bl_final <-  HILDA_wave_i_bl_final %>%
  dplyr::rename(rfwn_hh_pop_weight = ihhwte, rfwn_hh_samp_weight = ihhwtes)

## Drop all raw HILDA variables no longer relevant to rfwn dataset
HILDA_wave_i_bl_final <-  HILDA_wave_i_bl_final %>% dplyr::select(dplyr::starts_with("rfwn"), xwaveid)

rm(HILDA_wave_i_bl, HILDA_wave_i_clean, enumerated_i, i_weights_ind, waveI)

