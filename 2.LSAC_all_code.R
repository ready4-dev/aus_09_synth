## For this script to execute correctly:
## The db_var_name_matches_lup object created in 0.MakeLookups must be loaded 
## The functions from https://github.com/ready4-dev/fakefolk must be loaded.
## Rename default data files to merge-friendly file names:
rename_default_LSAC(old_path = "~/data_projects/confidential_data/LSAC/Wave_7/Original_files/",
                    new_path ="~/data_projects/confidential_data/LSAC/Wave_7/Renamed_files_for_merge/",
                    age_from = 0,
                    age_to = 18,
                    age_interval = 2,
                    letters = letters[seq(from = 1, to = 10)],
                    files_countdown = files_countdown)

## Create included variables list (a vector that specifies the generic (non-wave specific) variable names)
## Note: vector should include the cross-wave "hicid", "cohort" and "wave" identifiers
## Note: include static variables  with z-prefix, e.g.: "zf12_","zf09_" and "zf02_" (indigenous status, country of birth and sex, respectively)
## Note: Include variables for Parent 1 and 2 with wave-specific suffix
included_var_list <- c("hicid",
                       "cohort",
                       "wave",
                       "state",
                       "sla2011",
                       "sa22011",
                       "sa22016",
                       "pcoden",
                       "sos",
                       "cnfsad",
                       "cnfsad2",
                       "f03m1",
                       "f03m2",
                       "f03m3",
                       "f03bp1",
                       "f03bp2",
                       "f03cp1",
                       "f03cp2",
                       "f03dp1",
                       "f03dp2",
                       "f03ep1",
                       "f03ep2",
                       "f03fp1",
                       "f03fp2",
                       "f03gp1",
                       "f03gp2",
                       "f03hp1",
                       "f03hp2",
                       "f03ip1",
                       "f03ip2",
                       "zf12m1",
                       "zf12m2",
                       "zf12m3",
                       "zf12bp1",
                       "zf12bp2",
                       "zf12cp1",
                       "zf12cp2",
                       "zf12dp1",
                       "zf12dp2",
                       "zf12ep1",
                       "zf12ep2",
                       "zf12fp1",
                       "zf12fp2",
                       "zf12gp1",
                       "zf12gp2",
                       "zf12hp1",
                       "zf12hp2",
                       "zf12ip1",
                       "zf12ip2",
                       "f11m1",
                       "f11m2",
                       "f11m3",
                       "f11bp1",
                       "f11bp2",
                       "f11cp1",
                       "f11cp2",
                       "f11dp1",
                       "f11dp2",
                       "f11ep1",
                       "f11ep2",
                       "f11fp1",
                       "f11fp2",
                       "f11gp1",
                       "f11gp2",
                       "f11hp1",
                       "f11hp2",
                       "f11ip1",
                       "f11ip2",
                       "zf09m1",
                       "zf09m2",
                       "zf09m3",
                       "zf09bp1",
                       "zf09bp2",
                       "zf09cp1",
                       "zf09cp2",
                       "zf09dp1",
                       "zf09dp2",
                       "zf09ep1",
                       "zf09ep2",
                       "zf09fp1",
                       "zf09fp2",
                       "zf09gp1",
                       "zf09gp2",
                       "zf09hp1",
                       "zf09hp2",
                       "zf09ip1",
                       "zf09ip2",
                       "fd14a",
                       "re23a",
                       "zf02m1",
                       "zf02m2",
                       "zf02m3",
                       "zf02bp1",
                       "zf02bp2",
                       "zf02cp1",
                       "zf02cp2",
                       "zf02dp1",
                       "zf02dp2",
                       "zf02ep1",
                       "zf02ep2",
                       "zf02fp1",
                       "zf02fp2",
                       "zf02gp1",
                       "zf02gp2",
                       "zf02hp1",
                       "zf02hp2",
                       "zf02ip1",
                       "zf02ip2",
                       "prel",
                       "mmn",
                       "p1mn",
                       "fmn",
                       "p2mn",
                       "nyngsib",
                       "noldsib",
                       "parpart",
                       "relstat",
                       "ho04a5",
                       "ho11a2",
                       "ho11a1a",
                       "fd09a2",
                       "fd09b2",
                       "fd09a1",
                       "fd09b1",
                       "fd08a1",
                       "fd08b1",
                       "fd08a3a",
                       "fd08b3a",
                       "pw09a",
                       "pw09b",
                       "fn13a",
                       "fn13b",
                       "fn05",
                       "hinci",
                       "aemp",
                       "bemp",
                       "awork",
                       "bwork",
                       "sep",
                       "sep2",
                       "pe23d",
                       "pe07a",
                       "pe07b",
                       "pe07b1",
                       "fn04",
                       "fn02a5",
                       "fn02b5",
                       "smfq",
                       "aemot",
                       "cemot",
                       "asdqtb",
                       "asdqta",
                       "ak6s",
                       "bk6s",
                       "hs25a1",
                       "hs25a2",
                       "hs25b1",
                       "f18cm2",
                       "f18cm3",
                       "f18dm2",
                       "f18dm3",
                       "hs48a32",
                       "hs48a33",
                       "hs48a34",
                       "hs48a35",
                       "hs48a36",
                       "hs48a7",
                       "hs48b32",
                       "hs48b32b",
                       "hs48b32c",
                       "hs48b33",
                       "hs48b33b",
                       "hs48b33c",
                       "hs48b34",
                       "hs48b34b",
                       "hs48b34c",
                       "hs48b35",
                       "hs48b35b",
                       "hs48b35c",
                       "hs48b36",
                       "hs48b36b",
                       "hs48b36c",
                       "hs48b7",
                       "hs48b7b",
                       "hs48b7c",
                       "smfqc",
                       "ak6g",
                       "bk6g",
                       "hb16c10",
                       "calcharm",
                       "aalcp",
                       "balcp",
                       "hs48a30",
                       "hs48a31",
                       "hs48b30",
                       "hs48b30b",
                       "hs48b30c",
                       "hs48b31",
                       "hs48b31b",
                       "hs48b31c",
                       "se21b1",
                       "pa06b1",
                       "pa06b2",
                       "pa06b3",
                       "awarm",
                       "bwarm",
                       "aang",
                       "bang",
                       "ahostc",
                       "bhostc",
                       "hs54e",
                       "hs54b",
                       "hs54d",
                       "chu9d",
                       "re08c",
                       "re26c2",
                       "pa21a3",
                       "pa21b3",
                       "pa06c",
                       "pa23",
                       "pa23a1",
                       "parsit",
                       "f15m2",
                       "f15m3",
                       "f15bp1",
                       "f15bp2",
                       "f15dp1",
                       "f15dp2",
                       "f15cp1",
                       "f15cp2",
                       "f15ep1",
                       "f15ep2",
                       "f15gp1",
                       "f15gp2",
                       "f15fp1",
                       "f15fp2",
                       "f15hp1",
                       "f15hp2",
                       "f15ip1",
                       "f15ip2",
                       "ho11a3h",
                       "re15a5",
                       "re15b5",
                       "aarga",
                       "barga",
                       "pc46",
                       "se03c5d",
                       "se03t5d",
                       "se03a5d",
                       "pc33c",
                       "marshpr",
                       "apeer",
                       "tpeer",
                       "cpeer",
                       "sc26c7",
                       "sc26c2",
                       "sc26c5",
                       "sc26c1",
                       "sc26c8",
                       "sc26c9",
                       "marshgs",
                       "se21a1",
                       "se21a2",
                       "se21a3",
                       "se21a4",
                       "se21a5",
                       "academic",
                       "tlit",
                       "tlitb",
                       "tlitc",
                       "tmath",
                       "tmathb",
                       "tmathc",
                       "f18cm1",
                       "f18dm1",
                       "bulimia",
                       "anorexia",
                       "eatingdis",
                       "f13m1",
                       "f13am1",
                       "f13bm1",
                       "shcn",
                       "mresponse",
                       "mautonomy",
                       "mdemand",
                       "fresponse",
                       "fautonomy",
                       "fdemand",
                       "aparmon",
                       "aparmonb",
                       "bparmon",
                       "bparmonb",
                       "scheng",
                       "csrsla",
                       "csrsavc",
                       "pssm",
                       "hb14c3",
                       "hb14c4",
                       "hb14c5",
                       "he06b2",
                       "he06b1",
                       "he06c2",
                       "he06c1",
                       "re23c",
                       "hb28c3",
                       "hb26c2",
                       "hb16c12",
                       "hb27c2",
                       "ahfat",
                       "chfatb",
                       "ahsdrnk",
                       "chsdrnkb",
                       "re09a",
                       "re09b",
                       "re09c",
                       "re09d",
                       "re09e",
                       "re09f",
                       "re09g",
                       "re09h",
                       "hb28c1",
                       "hb27c1",
                       "hb29a1",
                       "hb29b1",
                       "wnsco",
                       "wseoi"
)

## Create merged (all waves) datasets
LSAC_merged_B_cohort <- read_merge_waves_LSAC(waves = c(1:7),
                                         path_to_source = "~/data_projects/confidential_data/LSAC/Wave_7/Renamed_files_for_merge/",
                                         included_vars_generic = included_var_list,
                                         included_vars_xwave = c("hicid", "cohort",
                                                                 included_var_list[startsWith(included_var_list,"z")],
                                                                 "bf15m2",	"bf15m3",
                                                                 "df15m2",	"df15m3",
                                                                 "cf15bp1",	"cf15bp2",
                                                                 "ef15dp1",	"ef15dp2",
                                                                 "df15cp1",	"df15cp2",
                                                                 "ff15ep1",	"ff15ep2"),
                                         data_type = "lsac",
                                         cohort= "b_",
                                         file_extension  = ".dta"
)
LSAC_merged_K_cohort <- read_merge_waves_LSAC(waves = c(3:9),
                                         path_to_source = "~/data_projects/confidential_data/LSAC/Wave_7/Renamed_files_for_merge/",
                                         included_vars_generic = included_var_list,
                                         included_vars_xwave = c("hicid", "cohort",
                                                                 included_var_list[startsWith(included_var_list,"z")],
                                                                 "bf15m2",	"bf15m3",
                                                                 "df15m2",	"df15m3",
                                                                 "cf15bp1",	"cf15bp2",
                                                                 "ef15dp1",	"ef15dp2",
                                                                 "df15cp1",	"df15cp2",
                                                                 "ff15ep1",	"ff15ep2"),
                                         data_type = "lsac",
                                         cohort= "k_",
                                         file_extension  = ".dta"
)

## Create household IDs:
LSAC_merged_B_cohort <- dplyr::mutate(LSAC_merged_B_cohort, hhid= paste0("hh_",hicid)) %>%
  dplyr::select("hhid", dplyr::everything())

LSAC_merged_K_cohort <- dplyr::mutate(LSAC_merged_K_cohort, hhid= paste0("hh_",hicid)) %>%
  dplyr::select("hhid", dplyr::everything())

## Remove any duplicated hicids due to changes in "static/historical" variables:
##NB 4 duplicates remove from CohortB and 8 from Cohort K. Duplicates typically reflected changes in country of birth code or ATSI status
LSAC_merged_B_cohort <- LSAC_merged_B_cohort[!duplicated(LSAC_merged_B_cohort$hicid), ]
LSAC_merged_K_cohort <- LSAC_merged_K_cohort[!duplicated(LSAC_merged_K_cohort$hicid), ]

##Create B cohort "cross-sectional" database using random sample for each wave
set.seed(2302)
index_vec_B <- sample(1:nrow(LSAC_merged_B_cohort))

list_index_vec_B <- split(index_vec_B, ceiling(seq_along(index_vec_B)/(nrow(LSAC_merged_B_cohort)/7)
))

LSAC_merged_B_cohort_filtered <- purrr::map(list_index_vec_B,
                                            ~ LSAC_merged_B_cohort %>%
                                              dplyr::slice(.x))

LSAC_B_cross_sec <- purrr::map2(LSAC_merged_B_cohort_filtered,
                                letters[1:7],
                                ~ .x %>%
                                  dplyr::select(dplyr::starts_with(.y),
                                                hicid, hhid, cohort, dplyr::starts_with("z")))

LSAC_B_cross_sec <- purrr::map(LSAC_B_cross_sec,
                               ~ .x %>%
                                 dplyr::rename_at(dplyr::vars(names(.x)[!names(.x)%in% c("hicid","hhid", "cohort")]),
                                                  ~(stringr::str_sub(.,start=2))))

LSAC_B_cross_sec <- purrr::map(LSAC_B_cross_sec,
                               ~ .x %>%
                                 dplyr::mutate(pcoden = as.character(pcoden)))

LSAC_cohortB_baseline <- purrr::reduce(LSAC_B_cross_sec,
                                       ~ dplyr::full_join(.x,.y))

##Mutate Wave ID and wave letter based on child age:
LSAC_cohortB_baseline <- LSAC_cohortB_baseline %>% dplyr::mutate(wave=ifelse(is.na(f03m1),NA,
                                                                             ifelse(f03m1 %in% c(0,1),1,
                                                                                    ifelse(f03m1 %in% c(2,3),2,
                                                                                           ifelse(f03m1 %in% c(4,5),3,
                                                                                                  ifelse(f03m1 %in% c(6,7),4,
                                                                                                         ifelse(f03m1 %in% c(8,9),5,
                                                                                                                ifelse(f03m1 %in% c(10,11),6,7)))))))) %>%
  dplyr::mutate(wave_letter= dplyr::recode(wave, '1' = "a", '2'="b", '3'="c", '4'="d", '5' = "e", '6'="f", '7'="g"))

##Create K cohort "cross-sectional" database using random sample for each wave


index_vec_K <- sample(1:nrow(LSAC_merged_K_cohort))

list_index_vec_K <- split(index_vec_K, ceiling(seq_along(index_vec_K)/(nrow(LSAC_merged_K_cohort)/7)
))

LSAC_merged_K_cohort_filtered <- purrr::map(list_index_vec_K,
                                            ~ LSAC_merged_K_cohort %>%
                                              dplyr::slice(.x))

LSAC_K_cross_sec <- purrr::map2(LSAC_merged_K_cohort_filtered,
                                letters[3:9],
                                ~ .x %>%
                                  dplyr::select(dplyr::starts_with(.y),
                                                hicid, hhid, cohort, dplyr::starts_with("z")))


LSAC_K_cross_sec <- purrr::map(LSAC_K_cross_sec,
                               ~ .x %>%
                                 dplyr::rename_at(dplyr::vars(names(.x)[!names(.x)%in% c("hicid","hhid", "cohort" )]),
                                                  ~(stringr::str_sub(.,start=2))))

LSAC_K_cross_sec <- purrr::map(LSAC_K_cross_sec,
                               ~ .x %>%
                                 dplyr::mutate(pcoden = as.character(pcoden)))

LSAC_cohortK_baseline <- purrr::reduce(LSAC_K_cross_sec,
                                       ~ dplyr::full_join(.x,.y))

##Mutate Wave ID and wave letter based on child age:
LSAC_cohortK_baseline <- LSAC_cohortK_baseline %>% dplyr::mutate(wave=ifelse(is.na(f03m1),NA,
                                                                             ifelse(f03m1 %in% c(4,5),1,
                                                                                    ifelse(f03m1 %in% c(6,7),2,
                                                                                           ifelse(f03m1 %in% c(8,9),3,
                                                                                                  ifelse(f03m1 %in% c(10,11),4,
                                                                                                         ifelse(f03m1 %in% c(12,13),5,
                                                                                                                ifelse(f03m1 %in% c(14,15),6,7)))))))) %>%
  dplyr::mutate(wave_letter=dplyr::recode(wave, '1' = "c", '2'="d", '3'="e", '4'="f", '5' = "g", '6'="h", '7'="i"))

##Create merged B and K baseline dataset and drop rows with no wave-specific data (i.e. drop outs with hicid only)
LSAC_mergedcohorts_baseline <- dplyr::full_join(LSAC_cohortB_baseline, LSAC_cohortK_baseline) %>%
  dplyr::select(hicid, hhid, cohort, wave, wave_letter, mmn, fmn, p1mn, p2mn, dplyr::everything()) %>%
  dplyr::filter(!is.na(wave))

# Recode missing values
LSAC_mergedcohorts_baseline[LSAC_mergedcohorts_baseline== -9 |
                              LSAC_mergedcohorts_baseline== -4 |
                              LSAC_mergedcohorts_baseline== -3 |
                              LSAC_mergedcohorts_baseline== -2 |
                              LSAC_mergedcohorts_baseline== -1] <- NA

# Remove observations with missing data from key variables age and sex:
LSAC_mergedcohorts_baseline <- LSAC_mergedcohorts_baseline %>% dplyr::filter(!is.na(f03m1)) %>% dplyr::filter(!is.na(f02m1))

rm(list_index_vec_B, list_index_vec_K, LSAC_B_cross_sec, LSAC_K_cross_sec, LSAC_cohortB_baseline, LSAC_cohortK_baseline,
   LSAC_merged_B_cohort, LSAC_merged_B_cohort_filtered, LSAC_merged_K_cohort, LSAC_merged_K_cohort_filtered)

## Generate derived variables applicable to PARENT 1:
LSAC_baseline <- LSAC_mergedcohorts_baseline %>%
  # Using stricter MI definition (i.e. without hs25a1)
  dplyr::mutate(mental_health_curr_prev_adult_md_bl_P1 = ifelse((hs25a2==1 |
                                                                   f18cm2==1 |
                                                                   f18dm2==1 |
                                                                   hs48a32 %in% c(3:4) |
                                                                   hs48a33 %in% c(3:4) |
                                                                   hs48a34 %in% c(3:4) |
                                                                   hs48a35 %in% c(3:4) |
                                                                   hs48a36 %in% c(3:4) |
                                                                   hs48a7 %in% c(3:4) |
                                                                   ak6g==1),
                                                                TRUE,
                                                                FALSE)) %>%
  dplyr::mutate(mental_health_curr_prev_adult_SU_bl_P1 = ifelse((aalcp ==1 |
                                                                 hs48a30 %in% c(3:4) |
                                                                 hs48a31 %in% c(3:4)),
                                                                    TRUE,
                                                                    FALSE)) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(behaviours_parenting_style_issues_bl_P1 = ifelse(all(is.na(awarm),
                                                                     is.na(aang),
                                                                     is.na(ahostc)),
                                                                 NA,
                                                                 any(awarm <3,
                                                                     aang >3,
                                                                     ahostc >5.5,
                                                                     na.rm=TRUE))) %>%
  dplyr::ungroup()
rm(LSAC_mergedcohorts_baseline)

## Generate derived variables applicable to PARENT 2:
LSAC_baseline <- LSAC_baseline %>%
  #Using stricter MI definition (i.e. without hs25b1)
  dplyr::mutate(mental_health_curr_prev_adult_md_bl_P2 = ifelse((f18cm3==1 |
                                                                   f18dm3==1 |
                                                                   hs48b32 %in% c(3:4) |
                                                                   hs48b32b %in% c(3:4) |
                                                                   hs48b32c %in% c(3:4) |
                                                                   hs48b33 %in% c(3:4) |
                                                                   hs48b33b %in% c(3:4) |
                                                                   hs48b33c %in% c(3:4) |
                                                                   hs48b34 %in% c(3:4) |
                                                                   hs48b34b %in% c(3:4) |
                                                                   hs48b34c %in% c(3:4) |
                                                                   hs48b35 %in% c(3:4) |
                                                                   hs48b35b %in% c(3:4) |
                                                                   hs48b35c %in% c(3:4) |
                                                                   hs48b36 %in% c(3:4) |
                                                                   hs48b36b %in% c(3:4) |
                                                                   hs48b36c %in% c(3:4) |
                                                                   hs48b7 %in% c(3:4) |
                                                                   hs48b7b %in% c(3:4) |
                                                                   hs48b7c %in% c(3:4) |
                                                                   bk6g==1),
                                                                TRUE,
                                                                FALSE)) %>%
  dplyr::mutate(mental_health_curr_prev_adult_SU_bl_P2 = ifelse((balcp ==1 |
                                                                   hs48b30  %in% c(3:4) |
                                                                   hs48b30b  %in% c(3:4) |
                                                                   hs48b30c  %in% c(3:4) |
                                                                   hs48b31  %in% c(3:4) |
                                                                   hs48b31b  %in% c(3:4) |
                                                                   hs48b31c  %in% c(3:4)),
                                                                TRUE,
                                                                FALSE)) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(behaviours_parenting_style_issues_bl_P2 = ifelse(all(is.na(bwarm),
                                                                     is.na(bang),
                                                                     is.na(bhostc)),
                                                                 NA,
                                                                 any(bwarm <3,
                                                                     bang >3,
                                                                     bhostc >5.5,
                                                                     na.rm=TRUE))) %>%
  dplyr::ungroup() %>%
  # For P2 (only), additional data exists on ongoing drug and alcohol abuse in waves 5 and 7:
dplyr::rowwise() %>%
  dplyr::mutate(mental_health_substance_abuse_bl_P2 = ifelse(wave %in% c(1:4, 6),
                                                             ifelse(balcp ==1,
                                                                    TRUE,
                                                                    FALSE),
                                                             ifelse(all(is.na(balcp),
                                                                        is.na(hs48b30b),
                                                                        is.na(hs48b30c),
                                                                        is.na(hs48b31b),
                                                                        is.na(hs48b31c)),
                                                                    NA,
                                                                    any(balcp==1,
                                                                        hs48b30b==4,
                                                                        hs48b30c==4,
                                                                        hs48b31b ==4,
                                                                        hs48b31c ==4, na.rm=TRUE)))) %>%
  dplyr::ungroup()

## Generate derived variables applicable to STUDY CHILD or HOUSEHOLD level:
  # NB. many logic checks across multiple vars require rowwise() or equivalent group_by(hicid)
LSAC_baseline <- LSAC_baseline %>%
  dplyr::rowwise() %>%
  dplyr::mutate(demographics_CALD_status_SC = ifelse(all(is.na(f11m1), is.na(fd14a)),
                                                     NA,
                                                     any(f11m1 !=1201, fd14a==1, na.rm=T))) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(demographics_ATSI_status_SC = ifelse(f12m1 %in% c(2:4),
                                                     TRUE,
                                                     FALSE)) %>%
  dplyr::mutate(demographics_generation_SC = ifelse(f09m1 !=1101,
                                                    "FIRST",
                                                    ifelse((f09m1==1101 & (f09m2!=1101 | f09m3 !=1101)),
                                                           "SECOND",
                                                           "THIRD_PLUS"))) %>%
  # recode Country of Birth to match brief HILDA categories
  dplyr::mutate(demographics_country_of_birth_SC = ifelse(f09m1 == 1101,
                                                          1,
                                                          ifelse(f09m1 %in% c(2100,
                                                                              1201,
                                                                              8102,
                                                                              8104,
                                                                              2201,
                                                                              9225),
                                                                 2,
                                                                 3))) %>%
  dplyr::mutate(demographics_same_sex_attracted = ifelse(((f02m1==1 & re23a==2) |
                                                            (f02m1==2 & re23a==1) |
                                                            re23a==3),
                                                         TRUE,
                                                         FALSE)) %>%
  dplyr::mutate(spatial_home_area_seifa = ifelse(wave %in% c(1:4),
                                                 cnfsad,
                                                 cnfsad2)) %>%
  dplyr::mutate(socio_econ_accom_status_bl = ifelse((ho04a5 %in% c(6,7) | ho11a2==1),
                                                    "at-risk",
                                                    "secure"))

#Recode fn05 (yearly household gross income) factor variable to numeric. [NB. labelled class data]
library(labelled)
LSAC_baseline$fn05 <- dplyr::recode(LSAC_baseline$fn05, `0` = 0,
                                    `1`	=	124800,
                                    `2`	=	119600,
                                    `3`	=	109200,
                                    `4`	=	91000,
                                    `5`	=	65000,
                                    `6`	=	46800,
                                    `7`	=	39000,
                                    `8`	=	33800,
                                    `9`	=	28600,
                                    `10`	=	23400,
                                    `11`	=	18200,
                                    `12`	=	13000,
                                    `13`	=	7800,
                                    `14`	=	3900,
                                    `15`	=	1300,
                                    `-99` = 0)

LSAC_baseline <- LSAC_baseline %>%
  dplyr::mutate(socio_econ_household_gross_income_bl = ifelse(wave==1,
                                                              fn05,
                                                              hinci/7*365.25)) %>%
  dplyr::mutate(socio_econ_household_ses = calculate_sep(wave = wave,
                                                         sep = sep,
                                                         sep2 = sep2)) %>%
# Note incarceration data only collected for "parent living elsewhere" subsample (small pop)
  dplyr::mutate(socio_econ_incarcerated_hh_mbr_bl = ifelse(wave %in% 1:5,
                                                           ifelse(pe07b==3,
                                                                  TRUE,
                                                                  FALSE),
                                                           ifelse(pe07b1==3,
                                                                  TRUE,
                                                                  FALSE))) %>%
  dplyr::group_by(hicid) %>%
  dplyr::mutate(socio_econ_welfare_bl = any(fn02a5==1, fn02b5==1, na.rm = T)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(behaviours_recent_suicide_attempt_bl = ifelse(hs54e>0,
                                                              TRUE,
                                                              FALSE)) %>%
  dplyr::mutate(mental_health_SDQ_emot_bl = ifelse(f03m1<12,
                                                   aemot,
                                                   cemot)) %>%
    dplyr::mutate(mental_health_SDQ_bl = ifelse(wave==1,
                                                asdqta,
                                                asdqtb)) %>%
  dplyr::mutate(mental_health_substance_abuse_bl_SC = ifelse(is.na(hb16c10),
                                                             NA,
                                                             ifelse(f03m1 <16,
                                                                    ifelse(hb16c10>0,
                                                                           TRUE,
                                                                           FALSE),
                                                                    ifelse((f02m1==1 & hb16c10>14) | (f02m1==2 & hb16c10>7),
                                                                           TRUE,
                                                                           FALSE)))) %>%
  # LSAC-specific (assumes <18 years): evaluates LSAC Social-Emotional Index in 0-3 year olds, parent-reported SDQ in 4-11 years olds,
  # & self-reported SDQ and Branch Eating Disorder Test in 12-17 year olds [note different cut-offs for parent vs child report SDQ]
  dplyr::rowwise() %>%
  dplyr::mutate(mental_health_mental_disorder_bl_SC = ifelse(f03m1<4,
                                                             ifelse(wnsco == 0,
                                                                    "No/Low",
                                                                    "Disorder"),
                                                             ifelse(f03m1<12,
                                                                    ifelse(aemot <4,
                                                                           "No/Low",
                                                                           ifelse(aemot == 4,
                                                                                  "Subthreshold",
                                                                                  "Disorder")),
                                                                    ifelse(any(cemot >=7 , eatingdis==1, na.rm=T),
                                                                           "Disorder",
                                                                           ifelse(cemot == 6,
                                                                                  "Subthreshold",
                                                                                  "No/Low"))))) %>%
  # use parent-reported SDQ for older children 12+ with missing self-report (i.e. self-report takes precedence)
  dplyr::mutate(mental_health_mental_disorder_bl_SC = ifelse(f03m1>=12 & is.na(mental_health_mental_disorder_bl_SC),
                                                             ifelse(aemot >= 5,
                                                                    "Disorder",
                                                                    ifelse(aemot == 4,
                                                                           "Subthreshold",
                                                                           "No/Low")),
                                                             mental_health_mental_disorder_bl_SC)) %>%
  dplyr::ungroup() %>%
   dplyr::mutate(mental_health_mental_disorder_bl_P1 = ifelse(is.na(ak6g),
                                                              NA,
                                                              ifelse(ak6g==1,
                                                                     "Disorder",
                                                                     "No/Low"))) %>%
  dplyr::mutate(mental_health_mental_disorder_bl_P2 = ifelse(is.na(bk6g),
                                                              NA,
                                                              ifelse(bk6g==1,
                                                                     "Disorder",
                                                                     "No/Low"))) %>%

  dplyr::rowwise() %>%
  dplyr::mutate(ace_emotional_abuse_bl  = ifelse((is.na(pa06b1)& is.na(pa06b2)& is.na(pa06b3)),
                                                 NA,
                                                 any(pa06b1 %in% c(2,4), pa06b2==3, pa06b3==3, na.rm = TRUE))) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(ace_sexual_abuse_bl  = ifelse(re23c<14,
                                              TRUE,
                                              FALSE)) %>%
  dplyr::mutate(ace_neglect_bl = ifelse(is.na(pa21b3),
                                          ifelse((pa21a3==1 |
                                                  pa06c==1),
                                                   TRUE,
                                                   FALSE),
                                          ifelse(((pa21b3==1 & pa21a3 ==1) |
                                                  pa06c==1),
                                                    TRUE,
                                                    FALSE))) %>%
    dplyr::rowwise() %>%
  dplyr::mutate(ace_parental_death_bl = ifelse(all(is.na(parsit),
                                                   is.na(f15m2),
                                                   is.na(f15m3),
                                                   is.na(f15bp1),
                                                   is.na(f15bp2),
                                                   is.na(f15cp1),
                                                   is.na(f15cp2),
                                                   is.na(f15dp1),
                                                   is.na(f15dp2),
                                                   is.na(f15ep1),
                                                   is.na(f15ep2),
                                                   is.na(f15gp1),
                                                   is.na(f15gp2),
                                                   is.na(f15fp1),
                                                   is.na(f15fp2),
                                                   is.na(f15hp1),
                                                   is.na(f15hp2),
                                                   is.na(f15ip1),
                                                   is.na(f15ip2)),
                                               NA,
                                               any(parsit %in% c(6,8,10,11), #includes parental death BEFORE birth
                                                   f15m2==10,
                                                   f15m3==10,
                                                   f15bp1==10,
                                                   f15bp2==10,
                                                   f15cp1==10,
                                                   f15cp2==10,
                                                   f15dp1==10,
                                                   f15dp2==10,
                                                   f15ep1==10,
                                                   f15ep2==10,
                                                   f15gp1==10,
                                                   f15gp2==10,
                                                   f15fp1==10,
                                                   f15fp2==10,
                                                   f15hp1==10,
                                                   f15hp2==10,
                                                   f15ip1==10,
                                                   f15ip2==10,
                                                   na.rm=TRUE))) %>%
    dplyr::ungroup() %>%
  dplyr::mutate(ace_parental_abandonment_bl = ifelse(((parsit %in% c(18:20)) |
                                                        pe07a %in% c(5, 10, 12) |
                                                        pe07b %in% c(2, 7, 9) |
                                                        pe07b1 %in% c(2, 7, 9)),
                                                     TRUE,
                                                     FALSE)) %>%
  dplyr::mutate(ace_domestic_violence_bl = ifelse(wave %in% 5:7,
                                                  ifelse((ho11a3h==1 |
                                                            re15a5 >3 |
                                                            re15b5 >3),
                                                         TRUE,
                                                         FALSE),
                                                  ifelse((re15a5 >3 |
                                                            re15b5 >3),
                                                         TRUE,
                                                         FALSE))) %>%
  dplyr::mutate(ace_parental_conflict_bl = ifelse(is.na(barga),
                                                  ifelse(aarga >3,
                                                         TRUE,
                                                         FALSE),
                                                  ifelse((aarga >3 |
                                                            barga>3),
                                                         TRUE,
                                                         FALSE))) %>%
  dplyr::group_by(hicid) %>%
  dplyr::mutate(ace_parental_MI_bl = any(mental_health_curr_prev_adult_md_bl_P1,
                                                             mental_health_curr_prev_adult_md_bl_P2, na.rm=TRUE)) %>%
  dplyr::mutate(ace_parental_SU_bl = any(mental_health_curr_prev_adult_SU_bl_P1,
                                         mental_health_curr_prev_adult_SU_bl_P2, na.rm=TRUE)) %>%
  dplyr::ungroup() %>%
  dplyr::rowwise() %>%
  dplyr::mutate(ace_incarcerated_household_mbr_bl= ifelse((is.na(pe07a) &
                                                             is.na(socio_econ_incarcerated_hh_mbr_bl)),
                                                          NA,
                                                          any(pe07a==6,
                                                              socio_econ_incarcerated_hh_mbr_bl,
                                                              na.rm = TRUE))) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(ace_low_ses_bl = ifelse(socio_econ_household_ses == "Q1",
                                        TRUE,
                                        FALSE)) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(ace_bullying_bl = ifelse(all(is.na(pc46),
                                             is.na(se03c5d),
                                             is.na(se03t5d),
                                             is.na(se03a5d)),
                                         NA,
                                         any(pc46==1,
                                             se03c5d==3,
                                             se03t5d==3,
                                             se03a5d==3,
                                             na.rm=TRUE))) %>%

  dplyr::mutate(ace_peer_isolation_bl = ifelse(all(is.na(marshpr),
                                                    is.na(tpeer),
                                                    is.na(cpeer),
                                                    is.na(apeer)),
                                                NA,
                                                any(marshpr<3,
                                                    tpeer >5,
                                                    cpeer >5,
                                                    apeer >5,
                                                    na.rm=TRUE))) %>%
  dplyr::ungroup() %>%
  # Recode discrimination variables which are numeric (1=Yes and 2=No) to logicals
  dplyr::mutate_at(dplyr::vars(c(sc26c5, sc26c1, sc26c7, sc26c2, sc26c8, sc26c9)), ~(dplyr::recode(., `1`=TRUE, `2`=FALSE))) %>%
  dplyr::group_by(hicid) %>%
  dplyr::mutate(ace_ethnic_discrimination_bl = any(sc26c5, sc26c1, na.rm=TRUE)) %>% #OPTIONAL TBC
  dplyr::mutate(ace_present_bl = any(ace_emotional_abuse_bl,
                                     ace_sexual_abuse_bl,
                                     ace_neglect_bl,
                                     ace_parental_death_bl,
                                     ace_parental_abandonment_bl,
                                     ace_domestic_violence_bl,
                                     ace_parental_conflict_bl,
                                     ace_parental_MI_bl,
                                     ace_parental_SU_bl,
                                     ace_incarcerated_household_mbr_bl,
                                     ace_low_ses_bl,
                                     ace_bullying_bl,
                                     ace_peer_isolation_bl,
                                     sc26c7,
                                     sc26c2,
                                     ace_ethnic_discrimination_bl,
                                     sc26c8,
                                     sc26c9,
                                     na.rm=TRUE)) %>%
  dplyr::ungroup()

LSAC_baseline[, c('se21a1','se21a2','se21a3','se21a4', 'se21a5')] <- lapply(LSAC_baseline[, c('se21a1','se21a2','se21a3','se21a4', 'se21a5')], as.numeric)
LSAC_baseline <- LSAC_baseline %>%
  dplyr::mutate(sum_self_esteem_low = ifelse((se21a1 + se21a2 + se21a3 + se21a4 + se21a5<15),
                                             TRUE,
                                             FALSE))  %>%
  dplyr::rowwise() %>%
  dplyr::mutate(rpf_psychosocial_low_self_esteem_bl = ifelse((is.na(marshgs) & is.na(sum_self_esteem_low)),
                                                             NA,
                                                             any(marshgs<3, sum_self_esteem_low,
                                                                 na.rm = TRUE))) %>%
  dplyr::ungroup() %>%
  dplyr::rowwise() %>%
  dplyr::mutate(rpf_psychosocial_poor_grades_bl = ifelse((is.na(academic) &
                                                            is.na(tlit) &
                                                            is.na(tlitb) &
                                                            is.na(tlitc) &
                                                            is.na(tmath) &
                                                            is.na(tmathb) &
                                                            is.na(tmathc)),
                                                         NA,
                                                         any(academic >2,
                                                             tlit <3,
                                                             tlitb <3,
                                                             tlitc <3,
                                                             tmath <3 ,
                                                             tmathb <3,
                                                             tmathc <3,
                                                             na.rm = TRUE))) %>%
  dplyr::ungroup() %>%
  dplyr::rowwise() %>%
  dplyr::mutate(rpf_psychosocial_history_poor_mh_bl = any(f18cm1==1,
                                                          f18dm1==1,
                                                          bulimia==1,
                                                          anorexia==1,
                                                          eatingdis==1,
                                                          smfqc==1,
                                                          na.rm = T)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(rpf_physiological_disability_illness_bl = ifelse(f13m1==2,
                                                                 TRUE,
                                                                 FALSE)) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(rpf_relational_parenting_style_bl = ifelse((is.na(mresponse) &
                                                              is.na(mautonomy) &
                                                              is.na(mdemand) &
                                                              is.na(fresponse) &
                                                              is.na(fautonomy) &
                                                              is.na(fdemand) &
                                                              is.na(behaviours_parenting_style_issues_bl_P1) &
                                                              is.na(behaviours_parenting_style_issues_bl_P2)),
                                                           NA,
                                                           any(mresponse >3,
                                                               mautonomy >3,
                                                               mdemand <3, #HIGH 'TRUE' RATE -  reconsider if appropriate measure
                                                               fresponse >3,
                                                               fautonomy >3,
                                                               fdemand <3, #HIGH 'TRUE' RATE -  reconsider if appropriate measure
                                                               behaviours_parenting_style_issues_bl_P1,
                                                               behaviours_parenting_style_issues_bl_P2,
                                                               na.rm = TRUE))) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(rpf_relational_low_parental_trust_bl = ifelse(((re09a +
                                                                  re09b +
                                                                  re09c +
                                                                  re09d +
                                                                  re09e +
                                                                  re09f +
                                                                  re09g +
                                                                  re09h) <20),
                                                              TRUE,
                                                              FALSE)) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(rpf_relational_low_school_bonding_bl = ifelse((is.na(scheng) &
                                                                 is.na(csrsla) &
                                                                 is.na(csrsavc) &
                                                                 is.na(pssm)),
                                                              NA,
                                                              any(scheng <3,
                                                                  csrsla >2,
                                                                  csrsavc <2,
                                                                  pssm <36,
                                                                  na.rm = TRUE))) %>%
  dplyr::ungroup() %>%
  dplyr::rowwise() %>%
  dplyr::mutate(rpf_lifestyle_low_physical_act_bl = ifelse((is.na(hb14c3) &
                                                              is.na(hb14c4) &
                                                              is.na(hb14c5)),
                                                           NA,
                                                           any(hb14c3 <3,
                                                               hb14c4==1,
                                                               hb14c5 >2,
                                                               na.rm = TRUE))) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(rpf_lifestyle_risky_behaviour_bl = ifelse((calcharm>10 |
                                                             hb28c1==1 |
                                                             hb27c1==1 |
                                                             hb29a1==1 |
                                                             hb29b1==1),TRUE,FALSE)) %>%
  # recode diet vars to logicals
  dplyr::mutate(poor_diet_parent_report = ifelse((ahfat>6 & ahsdrnk>3),
                                                 TRUE,
                                                 FALSE)) %>%
  dplyr::mutate(poor_diet_child_report = ifelse((chfatb >6 & chsdrnkb>3),
                                                TRUE,
                                                FALSE)) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(rpf_lifestyle_poor_diet_bl = ifelse((is.na(poor_diet_parent_report) &
                                                       is.na(poor_diet_child_report)),
                                                    NA,
                                                    any(poor_diet_parent_report,
                                                        poor_diet_child_report,
                                                        na.rm = TRUE))) %>%
  dplyr::ungroup() %>%

  dplyr::mutate(family_parent_of_child_with_ace = ifelse(ace_present_bl == TRUE,
                                                         TRUE,
                                                         FALSE)) %>%
  dplyr::mutate(family_sibling_of_child_with_ace = ifelse(ace_present_bl == TRUE,
                                                          TRUE,
                                                          FALSE))

#Scripts to generate appropriate wave-specific parent variables

LSAC_baseline <- LSAC_baseline %>%
  add_wave_chort_dep_col(target_in_lookup = "Age_P1",
                         new_col_name = "demographics_age_P1")%>%
  add_wave_chort_dep_col(target_in_lookup = "Age_P2",
                         new_col_name = "demographics_age_P2")%>%

  add_wave_chort_dep_col(target_in_lookup = "Sex_P1",
                         new_col_name = "demographics_sex_P1")%>%
  add_wave_chort_dep_col(target_in_lookup = "Sex_P2",
                         new_col_name = "demographics_sex_P2")

LSAC_baseline <- LSAC_baseline %>%
  add_wave_chort_dep_col(target_in_lookup = "ATSI_status_P1",
                         new_col_name = "demographics_ATSI_status_P1")%>%
  add_wave_chort_dep_col(target_in_lookup = "ATSI_status_P2",
                         new_col_name = "demographics_ATSI_status_P2") %>%
  dplyr::mutate_at(dplyr::vars(demographics_ATSI_status_P1, demographics_ATSI_status_P2), ~(ifelse(is.na(.), NA,
                                                                                                             ifelse(. %in% c(2:4), TRUE,FALSE))))
LSAC_baseline <- LSAC_baseline %>%
  add_wave_chort_dep_col(target_in_lookup = "LOTE_P1",
                         new_col_name = "demographics_CALD_status_P1")%>%
  add_wave_chort_dep_col(target_in_lookup = "LOTE_P2",
                         new_col_name = "demographics_CALD_status_P2") %>%
  dplyr::mutate_at(dplyr::vars(demographics_CALD_status_P1, demographics_CALD_status_P2), ~(ifelse(is.na(.), NA,
                                                                                                             ifelse(. !=1201, TRUE, FALSE))))
LSAC_baseline <- LSAC_baseline %>%
  add_wave_chort_dep_col(target_in_lookup = "Country_of_birth_P1",
                         new_col_name = "demographics_country_of_birth_P1")%>%
  add_wave_chort_dep_col(target_in_lookup = "Country_of_birth_P2",
                         new_col_name = "demographics_country_of_birth_P2") %>%
  dplyr::mutate(demographics_generation_P1 = ifelse(is.na(demographics_country_of_birth_P1), NA,
                                                    ifelse(demographics_country_of_birth_P1 !=1101,
                                                           "FIRST",
                                                           "SECOND_PLUS"))) %>%
  dplyr::mutate(demographics_generation_P2 = ifelse(is.na(demographics_country_of_birth_P2), NA,
                                                    ifelse(demographics_country_of_birth_P2 !=1101,
                                                           "FIRST",
                                                           "SECOND_PLUS"))) %>%
# recode Country of Birth to match brief HILDA categories
dplyr::mutate(demographics_country_of_birth_P1 = ifelse(demographics_country_of_birth_P1 == 1101,
                                                        1,
                                                        ifelse(demographics_country_of_birth_P1 %in% c(2100,
                                                                                                       1201,
                                                                                                       8102,
                                                                                                       8104,
                                                                                                       2201,
                                                                                                       9225),
                                                               2,
                                                               3))) %>%

  dplyr::mutate(demographics_country_of_birth_P2 = ifelse(demographics_country_of_birth_P2 == 1101,
                                                          1,
                                                          ifelse(demographics_country_of_birth_P2 %in% c(2100,
                                                                                                         1201,
                                                                                                         8102,
                                                                                                         8104,
                                                                                                         2201,
                                                                                                         9225),
                                                                 2,
                                                                 3)))

## Rename original LSAC variable names to rfwn variable names:

LSAC_baseline_final <- LSAC_baseline %>% rename_to_rfwn_vars(lookup_tb = db_var_name_matches_lup,
                                                             old_var_name_col = "LSAC_var_name",
                                                             old_var_stub_col = "LSAC_var_stub",
                                                             wave_prefix = "")

##Drop all raw LSAC variables no longer relevant to rfwn dataset
LSAC_baseline_final <- LSAC_baseline_final %>% dplyr::select(dplyr::starts_with("rfwn"), hicid)

rm(LSAC_baseline)

################################################################################################
## OPTIONAL: Create agent-specific observations in baseline dataset

# 1. Create agent ids
LSAC_baseline_ind_agents <- LSAC_baseline_final %>%
  dplyr::mutate(sc_id=paste0(rfwn_family_household_id,"_","sc")) %>%
  dplyr::mutate(p1_id=paste0(rfwn_family_household_id,"_","p1")) %>%
  dplyr::mutate(p2_id=paste0(rfwn_family_household_id,"_","p2")) %>%
  dplyr::select(rfwn_family_household_id, sc_id, p1_id, p2_id, dplyr::everything())

# 2. Create long dataset by gathering agent type
LSAC_baseline_ind_agents <- LSAC_baseline_ind_agents %>%
  tidyr::gather(key="agent_type", value="agent_id", sc_id, p1_id, p2_id) %>%
  dplyr::select(rfwn_family_household_id, hicid, agent_type, agent_id, dplyr::everything())

# 2.1 Remove "phantom" parents from single parent households
LSAC_baseline_ind_agents <- LSAC_baseline_ind_agents[!(LSAC_baseline_ind_agents$agent_type=="p2_id" &
                                                         (is.na(LSAC_baseline_ind_agents$rfwn_family_mother_id) |
                                                            is.na(LSAC_baseline_ind_agents$rfwn_family_father_id))),]

## Create sibling agents based on study child's number of siblings
LSAC_generate_sibs <- LSAC_baseline_ind_agents[rep(row.names(LSAC_baseline_ind_agents),
                                              times = (ifelse(is.na(LSAC_baseline_ind_agents$rfwn_family_nbr_younger_siblings),
                                                              0,
                                                              LSAC_baseline_ind_agents$rfwn_family_nbr_younger_siblings) +
                                                      ifelse(is.na(LSAC_baseline_ind_agents$rfwn_family_nbr_older_siblings),
                                                              0,
                                                              LSAC_baseline_ind_agents$rfwn_family_nbr_older_siblings))),
                                         ] %>%
  dplyr::group_by(agent_id) %>%
  dplyr::mutate(child_nbr = 1:dplyr::n()) %>%
  dplyr::ungroup() %>%
  dplyr::filter(agent_type=="sc_id") %>%
  dplyr::mutate(agent_type = "sib_id") %>%
  dplyr::mutate(agent_id = paste0(agent_id %>% stringr::str_sub(end=-3),
                                  "sib_",
                                  child_nbr)) %>%
  dplyr::select(-child_nbr)

## Merge sibling table with main agent table
LSAC_baseline_ind_agents <- dplyr::bind_rows(LSAC_baseline_ind_agents, LSAC_generate_sibs)
rm(LSAC_generate_sibs)
