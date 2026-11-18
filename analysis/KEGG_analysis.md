KEGG_analysis
================
Karis Moon
2026-11-18

- [Import data](#import-data)
- [Filter for unique modules in
  KEGGModuleStepHit](#filter-for-unique-modules-in-keggmodulestephit)

# Import data

``` r
# Libraries
library(readxl)
library(tidyverse)
library(janitor)
library(dplyr)

metabolic_result_full_community_path <- "../data/METABOLIC_result_full_community.xlsx"

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
  select(microbe, Module, everything())
```

    ## # A tibble: 105 × 5
    ##    microbe         Module                      Module.step KO.id Module.Category
    ##    <chr>           <chr>                       <chr>       <chr> <chr>          
    ##  1 chitinophaga    D-Glucuronate degradation,… M00061+03   (K01… Other carbohyd…
    ##  2 chitinophaga    Heparan sulfate degradation M00078+06   K012… Glycosaminogly…
    ##  3 chitinophaga    Pectin degradation          M00081+01   K010… Other carbohyd…
    ##  4 chitinophaga    Coenzyme M biosynthesis     M00358+01   K080… Methane metabo…
    ##  5 chitinophaga    Isoleucine biosynthesis, p… M00535+01   K090… Branched-chain…
    ##  6 chitinophaga    Menaquinone biosynthesis, … M00930+04   K117… Cofactor and v…
    ##  7 chitinophaga    Proline metabolism          M00972+03   K003… Arginine and p…
    ##  8 koribacteraceae N-glycan precursor biosynt… M00055+10   K007… Glycan biosynt…
    ##  9 koribacteraceae Thiamine biosynthesis, pro… M00127+01   K031… Cofactor and v…
    ## 10 koribacteraceae Formaldehyde assimilation,… M00346+02   K008… Methane metabo…
    ## # ℹ 95 more rows
