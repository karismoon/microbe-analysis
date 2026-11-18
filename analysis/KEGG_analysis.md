KEGG_analysis
================
Karis Moon
2026-11-18

- [Import data](#import-data)
- [Clean metabolite deep dive data](#clean-metabolite-deep-dive-data)
- [Filter for unique modules in
  KEGGModuleStepHit](#filter-for-unique-modules-in-keggmodulestephit)
- [Find KEGG modules using
  metabolites](#find-kegg-modules-using-metabolites)
  - [Pull compound IDs](#pull-compound-ids)
  - [Pull relevant KEGG modules](#pull-relevant-kegg-modules)

# Import data

``` r
# Libraries
library(readxl)
library(tidyverse)
library(janitor)
library(dplyr)

metabolic_result_full_community_path <- "../data/METABOLIC_result_full_community.xlsx"
metabolite_deep_dive_path <- "../data/Metabolite deep dive.xlsx"
metabolites_of_interest_path <- "../data/metabolites_of_interest.csv"

df_metabolites_of_interest <- read_csv(metabolites_of_interest_path)
df_metabolite_deep_dive <- read_excel(metabolite_deep_dive_path)
df_KEGGModuleStepHit <- read_excel(metabolic_result_full_community_path, sheet = "KEGGModuleStepHit")
df_dbCAN2Hit <- read_excel(metabolic_result_full_community_path, sheet = "dbCAN2Hit")
df_MEROPSHit <- read_excel(metabolic_result_full_community_path, sheet = "MEROPSHit")
```

Define microbe names for later use throughout

``` r
microbes <- c(
  "chitinophaga",
  "koribacteraceae",
  "parvibaculum",
  "propionicimonas",
  "rhodomicrobium",
  "rhodopseudomonas",
  "sporomusa",
  "terracidiphilus"
)
```

# Clean metabolite deep dive data

``` r
df_metabolite_deep_dive_longer <-
  df_metabolite_deep_dive %>%
  # Split comma-separated organisms into separate rows (creates all Produced x Consumed pairs)
  separate_longer_delim(`Produced by`, delim = ",") %>%
  separate_longer_delim(`Consumed by`, delim = ",") %>%
  
  # Split multiline pathways into separate rows
  separate_longer_delim(`Relevant KEGG Pathways`, delim = "\n") %>%
  
  # Trim leading/trailing whitespace across all text columns
  mutate(across(where(is.character), str_trim)) %>%
  
  # Remove empty pathway rows caused by blank lines (keep rows where pathway was originally NA)
  filter(is.na(`Relevant KEGG Pathways`) | `Relevant KEGG Pathways` != "") %>%
  
  # Standardize organism capitalization (e.g., "rhodomicrobium" vs "Rhodomicrobium")
  mutate(
    `Produced by` = str_to_sentence(`Produced by`),
    `Consumed by` = str_to_sentence(`Consumed by`)
  ) %>%
  
  # Extract Pathway ID (e.g., "map01110") and Pathway Name into separate columns
  separate_wider_regex(
    `Relevant KEGG Pathways`,
    patterns = c(
      `Pathway ID` = "^map\\d+",
      "\\s+",
      `Pathway Name` = ".*"
    ),
    too_few = "align_start" # Safely keeps NAs for metabolites with no pathways
  ) %>%
  
  # Remove the trailing " (1)" count from the Pathway Name
  mutate(`Pathway Name` = str_remove(`Pathway Name`, "\\s*\\(\\d+\\)$"))

df_metabolite_deep_dive_longer
```

    ## # A tibble: 364 × 6
    ##    Metabolite `Produced by`    `Consumed by` `KEGG Compound ID` `Pathway ID`
    ##    <chr>      <chr>            <chr>         <chr>              <chr>       
    ##  1 sulfate    Rhodopseudomonas Chitinophaga  C00059             map01110    
    ##  2 sulfate    Rhodopseudomonas Chitinophaga  C00059             map00261    
    ##  3 sulfate    Rhodopseudomonas Chitinophaga  C00059             map01320    
    ##  4 sulfate    Rhodopseudomonas Chitinophaga  C00059             map02010    
    ##  5 sulfate    Rhodopseudomonas Chitinophaga  C00059             map00230    
    ##  6 sulfate    Rhodopseudomonas Chitinophaga  C00059             map02024    
    ##  7 sulfate    Rhodopseudomonas Chitinophaga  C00059             map01120    
    ##  8 sulfate    Rhodopseudomonas Chitinophaga  C00059             map01100    
    ##  9 sulfate    Rhodopseudomonas Chitinophaga  C00059             map00920    
    ## 10 sulfate    Rhodopseudomonas Chitinophaga  C00059             map00270    
    ## # ℹ 354 more rows
    ## # ℹ 1 more variable: `Pathway Name` <chr>

# Filter for unique modules in KEGGModuleStepHit

For each microbe, find what modules it has that no other microbe has.

``` r
get_unique_kegg <- function(df, microbe) {

  target_col <- paste0(microbe, ".Module.step.presence")

  df %>% 
    filter(
      .data[[target_col]] == "Present",
      if_all(
        ends_with(".Module.step.presence") & 
          !all_of(target_col),
        ~ .x == "Absent"
      )
    )
}

unique_kegg <- 
  set_names(microbes) %>% 
  map(~ get_unique_kegg(df_KEGGModuleStepHit, .x))

unique_kegg_all <- 
  imap_dfr(
    unique_kegg,
    ~ .x %>% 
      select(Module.step, Module, KO.id, Module.Category) %>% 
      mutate(microbe = .y)
  ) %>% 
  select(microbe, Module.step, Module, KO.id, Module.Category) %>% 
  arrange(microbe, Module.step)

unique_kegg_all %>% 
  select(microbe, Module, everything()) %>% 
  write.csv("../data/unique_kegg_pathways.csv", row.names = FALSE)
```

# Find KEGG modules using metabolites

## Pull compound IDs

``` r
# 1. Download the full KEGG compound synonym table
kegg_cpds <- read_tsv(
  "https://rest.kegg.jp/list/compound",
  col_names = c("kegg_id", "synonyms"),
  show_col_types = FALSE
) %>%
  mutate(kegg_id = str_remove(kegg_id, "^cpd:")) %>%
  separate_rows(synonyms, sep = "; ") %>%
  mutate(syn_lower = str_to_lower(str_trim(synonyms)))
```

``` r
# 2. Match your metabolite list (case-insensitive exact synonym match)
df_metabolites_of_interest <-
  df_metabolites_of_interest %>% 
  mutate(syn_lower = str_to_lower(str_trim(metabolite)))

df_compound_matched <- df_metabolites_of_interest %>%
  left_join(kegg_cpds, by = "syn_lower") %>% 
  select(metabolite, producer, consumer, analyses_matched, kegg_id)

df_compound_matched
```

    ## # A tibble: 21 × 5
    ##    metabolite                         producer consumer analyses_matched kegg_id
    ##    <chr>                              <chr>    <chr>    <chr>            <chr>  
    ##  1 NH3                                parviba… koribac… unique_consumer  C00014 
    ##  2 NH3                                rhodomi… koribac… unique_consumer  C00014 
    ##  3 NH3                                rhodops… koribac… unique_consumer  C00014 
    ##  4 NH3                                sporomu… koribac… unique_consumer  C00014 
    ##  5 [SoxY protein]-S-(2-sulfodisulfan… rhodomi… rhodops… unique_consumer  C21899 
    ##  6 [SoxY protein]-S-sulfanyl-L-cyste… rhodops… rhodomi… unique_producer  C21898 
    ##  7 [SoxY protein]-S-sulfosulfanyl-L-… rhodomi… rhodops… unique_consumer  C21815 
    ##  8 adenylyl sulfate                   chitino… parviba… unique_consumer  C00224 
    ##  9 adenylyl sulfate                   rhodomi… parviba… unique_consumer  C00224 
    ## 10 adenylyl sulfate                   rhodops… parviba… unique_consumer  C00224 
    ## # ℹ 11 more rows

## Pull relevant KEGG modules

``` r
# 1. Pull the KEGG Compound <-> Module lookup table from the REST API
kegg_cpd_mod <- read_tsv(
  "https://rest.kegg.jp/link/module/compound",
  col_names = c("kegg_id", "module_id"),
  show_col_types = FALSE
) %>%
  mutate(
    kegg_id   = str_remove(kegg_id, "^cpd:"),
    module_id = str_remove(module_id, "^md:")
  )
```

``` r
# 2. Summarize module completeness per microbe from KEGGModuleStepHit
microbe_modules <- df_KEGGModuleStepHit %>%
  # Extract clean Module ID (e.g., "M00001" from "M00001+01")
  mutate(module_id = str_extract(Module.step, "^M\\d+")) %>%
  # Pivot all microbe presence columns into long format
  pivot_longer(
    cols = ends_with(".Module.step.presence"),
    names_to = "microbe",
    values_to = "presence"
  ) %>%
  # Clean microbe name (e.g., "chitinophaga.Module.step.presence" -> "chitinophaga")
  mutate(microbe = str_remove(microbe, "\\.Module\\.step\\.presence$")) %>%
  # Calculate how many steps of each module each microbe has
  group_by(microbe, module_id, Module, Module.Category) %>%
  summarise(
    steps_present = sum(presence == "Present", na.rm = TRUE),
    total_steps   = n(),
    completeness  = round(steps_present / total_steps, 2),
    .groups = "drop"
  ) %>%
  # Keep modules where the microbe has at least 1 step (or filter completeness >= 0.75 for near-complete)
  filter(steps_present > 0)

# 3. Join your metabolites -> KEGG Modules -> Microbe Module Presence
df_metabolite_modules <- df_compound_matched %>%
  filter(!is.na(kegg_id)) %>%
  inner_join(kegg_cpd_mod, by = "kegg_id", relationship = "many-to-many") %>%
  inner_join(microbe_modules, by = "module_id", relationship = "many-to-many")

# 4. Optional: Filter specifically to the producer or consumer microbes in each row
# (Uses case-insensitive matching in case capitalization differs between tables)
df_producer_consumer_modules <- df_metabolite_modules %>%
  mutate(
    role_in_link = case_when(
      str_to_lower(microbe) == str_to_lower(producer) &
        str_to_lower(microbe) == str_to_lower(consumer) ~ "Both",
      str_to_lower(microbe) == str_to_lower(producer) ~ "Producer",
      str_to_lower(microbe) == str_to_lower(consumer) ~ "Consumer",
      TRUE ~ "Other community member"
    )
  ) %>%
  filter(role_in_link != "Other community member") %>%
  arrange(metabolite, module_id, role_in_link)

df_producer_consumer_modules
```

    ## # A tibble: 191 × 13
    ##    metabolite producer       consumer analyses_matched kegg_id module_id microbe
    ##    <chr>      <chr>          <chr>    <chr>            <chr>   <chr>     <chr>  
    ##  1 NH3        parvibaculum   koribac… unique_consumer  C00014  M00029    koriba…
    ##  2 NH3        rhodomicrobium koribac… unique_consumer  C00014  M00029    koriba…
    ##  3 NH3        rhodopseudomo… koribac… unique_consumer  C00014  M00029    koriba…
    ##  4 NH3        sporomusa      koribac… unique_consumer  C00014  M00029    koriba…
    ##  5 NH3        parvibaculum   koribac… unique_consumer  C00014  M00029    parvib…
    ##  6 NH3        rhodomicrobium koribac… unique_consumer  C00014  M00029    rhodom…
    ##  7 NH3        rhodopseudomo… koribac… unique_consumer  C00014  M00029    rhodop…
    ##  8 NH3        sporomusa      koribac… unique_consumer  C00014  M00029    sporom…
    ##  9 NH3        parvibaculum   koribac… unique_consumer  C00014  M00175    koriba…
    ## 10 NH3        rhodomicrobium koribac… unique_consumer  C00014  M00175    koriba…
    ## # ℹ 181 more rows
    ## # ℹ 6 more variables: Module <chr>, Module.Category <chr>, steps_present <int>,
    ## #   total_steps <int>, completeness <dbl>, role_in_link <chr>
