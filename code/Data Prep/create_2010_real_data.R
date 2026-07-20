# ACS 2009, PEP 2009 (vintage precensus), Census 2010
# Builds full-pop and AIAN datasets matching 2020 wDensity schema (no WP, no lag years).

library(tidycensus)
library(tidyverse)
library(tigris)
library(sf)
library(stringr)
library(rstudioapi)

# Working directory: this script's folder (RStudio or Rscript)
tryCatch({
  current_path <- getActiveDocumentContext()$path
  if (!is.null(current_path) && nzchar(current_path)) {
    setwd(dirname(current_path))
  }
}, error = function(e) invisible(NULL))
if (!file.exists("../../data/pep_2000_2010_vintage_precensus.csv")) {
  # Fallback when run via Rscript from repo root
  if (file.exists("code/data prep/create_2010_real_data.R")) {
    setwd("code/data prep")
  } else if (file.exists("create_2010_real_data.R")) {
    # already in script dir
  } else {
    stop("Cannot find data/ relative to working directory: ", getwd())
  }
}

data_dir <- "../../data"
options(tigris_use_cache = TRUE)

pad_geoid <- function(x) {
  sprintf("%05d", as.integer(as.character(x)))
}

build_adjacency_for_df <- function(df, adj_csv = file.path(data_dir, "countyadj.csv"),
                                   abb_csv = file.path(data_dir, "state_abbreviations.csv")) {
  D <- read.csv(adj_csv)
  abb <- read.csv(abb_csv)
  for (i in seq_len(nrow(D))) {
    nn <- D$X[i]
    st <- substr(nn, nchar(nn) - 1, nchar(nn))
    ind <- which(abb$Code == st)
    if (length(ind) == 1) {
      D$X[i] <- gsub(substr(nn, nchar(nn) - 2, nchar(nn)), paste0(", ", abb$State[ind]), nn)
    }
  }

  cat("Names in df but not adjacency:", length(setdiff(df$NAME, D$X)), "\n")
  cat("Names in adjacency but not df:", length(setdiff(D$X, df$NAME)), "\n")
  if (length(setdiff(df$NAME, D$X)) > 0) {
    print(setdiff(df$NAME, D$X))
  }

  df <- df[df$NAME %in% D$X, , drop = FALSE]
  D2 <- D
  rownames(D2) <- D2$X
  D2 <- D2[, -1, drop = FALSE]
  ind <- match(df$NAME, rownames(D2))
  D2 <- D2[ind, ind, drop = FALSE]
  # Align dimnames to df$NAME (CSV colnames are abbrev/check.names forms)
  colnames(D2) <- df$NAME
  rownames(D2) <- df$NAME
  stopifnot(identical(rownames(D2), df$NAME))
  stopifnot(identical(colnames(D2), df$NAME))
  list(df = df, adjacency = D2)
}

print_geoid_overlap <- function(geoid_2009, geoid_2019, label = "full") {
  g09 <- unique(as.character(geoid_2009))
  g19 <- unique(as.character(geoid_2019))
  both <- intersect(g09, g19)
  only_2009 <- setdiff(g09, g19)
  only_2019 <- setdiff(g19, g09)
  cat("\n==== GEOID overlap (", label, ") ====\n", sep = "")
  cat("n_2009:", length(g09), "\n")
  cat("n_2019:", length(g19), "\n")
  cat("n_both:", length(both), "\n")
  cat("n_only_2009:", length(only_2009), "\n")
  if (length(only_2009) > 0) cat("only_2009:", paste(only_2009, collapse = ", "), "\n")
  cat("n_only_2019:", length(only_2019), "\n")
  if (length(only_2019) > 0) cat("only_2019:", paste(only_2019, collapse = ", "), "\n")
  invisible(list(both = both, only_2009 = only_2009, only_2019 = only_2019))
}

compare_decades <- function(df_2010, df_2020, label = "full",
                            cols = c("acs", "pep", "census")) {
  m <- merge(
    df_2010 %>% select(GEOID, all_of(cols)),
    df_2020 %>% select(GEOID, all_of(cols)),
    by = "GEOID",
    suffixes = c("_2010", "_2020")
  )
  cat("\n==== Cross-decade comparison (", label, ") ====\n", sep = "")
  cat("n matched GEOIDs:", nrow(m), "\n")

  # Rows = 2009/10 series; columns = 2019/20 series
  mat <- matrix(NA_real_, nrow = length(cols), ncol = length(cols),
                dimnames = list(
                  paste0(cols, "_2009_10"),
                  paste0(cols, "_2019_20")
                ))
  for (i in seq_along(cols)) {
    for (j in seq_along(cols)) {
      x <- m[[paste0(cols[i], "_2010")]]
      y <- m[[paste0(cols[j], "_2020")]]
      ok <- is.finite(x) & is.finite(y)
      mat[i, j] <- cor(x[ok], y[ok])
    }
  }
  cat("\nCross-decade correlation table (rows=2009/10, cols=2019/20):\n")
  print(round(mat, 4))

  cat("\nOverall proportional sum change ((sum_2020 - sum_2010) / sum_2010):\n")
  for (col in cols) {
    s0 <- sum(m[[paste0(col, "_2010")]], na.rm = TRUE)
    s1 <- sum(m[[paste0(col, "_2020")]], na.rm = TRUE)
    pct <- 100 * (s1 - s0) / s0
    cat(sprintf("  %-8s sum_2010=%.0f  sum_2020=%.0f  change=%+.2f%%\n", col, s0, s1, pct))
  }
  invisible(list(merged = m, cor_table = mat))
}


#### Census 2010 ####
census <- get_decennial(
  geography = "county",
  variables = "P001001",
  year = 2010
) %>%
  select(GEOID, NAME, census = value) %>%
  mutate(GEOID = pad_geoid(GEOID))

#### ACS 2009 ####
ACS <- get_acs(
  geography = "county",
  variables = "B01003_001",
  year = 2009
) %>%
  select(GEOID, NAME, acs = estimate) %>%
  mutate(GEOID = pad_geoid(GEOID))

#### PEP 2009 (vintage precensus) ####
pep_vintage_raw <- read_csv(
  file.path(data_dir, "pep_2000_2010_vintage_precensus.csv"),
  show_col_types = FALSE
) %>%
  filter(SUMLEV == "050")

PEP <- pep_vintage_raw %>%
  mutate(
    GEOID = paste0(STATE, COUNTY),
    NAME = paste0(CTYNAME, ", ", STNAME),
    pep = POPESTIMATE072009
  ) %>%
  select(GEOID, NAME, pep)

if (FALSE) {
  # Sanity check: compare precensus vs postcensus 2009 estimates
  pep_postcensus <- read_csv(
    file.path(data_dir, "pep_2000_2010_postcensus.csv"),
    show_col_types = FALSE
  ) %>%
    filter(SUMLEV == 50)
  cat("Precensus 2009 total:", sum(pep_vintage_raw$POPESTIMATE072009), "\n")
  cat("Postcensus 2009 total:", sum(pep_postcensus$POPESTIMATE2009), "\n")
  cat("Precensus 2010 total:", sum(pep_vintage_raw$POPESTIMATE072010), "\n")
  cat("Postcensus 2010 total:", sum(pep_postcensus$POPESTIMATE2010), "\n")
}

#### Bring together (full pop) — no WP ####
# Merge on GEOID only — NAME spellings can differ (e.g. Doña Ana)
df <- census %>%
  select(GEOID, NAME, census) %>%
  merge(ACS %>% select(GEOID, acs), by = "GEOID", all = TRUE) %>%
  merge(PEP %>% select(GEOID, pep), by = "GEOID", all = TRUE)

#### Filter to mainland US ####
tg <- tigris::fips_codes %>%
  mutate(fips = paste0(state_code, county_code)) %>%
  filter(!(state_code %in% c("02", "15", "60", "66", "69", "72", "74", "78")))

df <- df %>% filter(GEOID %in% tg$fips)
df <- tg %>%
  select(state = state_name, GEOID = fips) %>%
  distinct() %>%
  merge(df, by = "GEOID")

NA_rows <- apply(df[, c("acs", "census", "pep")], 1, function(xx) any(is.na(xx)))
cat("Full-pop rows with NA in acs/census/pep:", sum(NA_rows), "\n")
if (sum(NA_rows) > 0) print(df[NA_rows, c("GEOID", "NAME", "acs", "census", "pep")])

# Drop incomplete rows before overlap / adjacency
df <- df[!NA_rows, ]

cat("Vintage PEP county sum check:\n")
cat("  filtered pep sum:", sum(df$pep), "\n")
cat("  raw POPESTIMATE072009 mainland-ish (all counties in CSV):",
    sum(pep_vintage_raw$POPESTIMATE072009), "\n")


#### GEOID overlap with 2019/20 full-pop; keep intersection ####
df_2010_pre <- df
load(file.path(data_dir, "census_ACS_PEP_WP_wDensity_11152024.RData"))
df_2020_full <- df
df <- df_2010_pre

overlap_full <- print_geoid_overlap(df$GEOID, df_2020_full$GEOID, label = "full-pop")
df <- df %>% filter(GEOID %in% overlap_full$both)
# Use 2019/20 NAME spellings so adjacency matching aligns (e.g. Doña Ana)
df <- df %>%
  select(-NAME) %>%
  merge(df_2020_full %>% select(GEOID, NAME), by = "GEOID")
cat("Full-pop rows after intersection filter:", nrow(df), "\n")


#### Densities from 2009/10 pops over 2010 land area ####
counties_2010 <- counties(cb = FALSE, year = 2010)
land_area <- counties_2010 %>%
  st_set_geometry(NULL) %>%
  transmute(
    GEOID = pad_geoid(GEOID10),
    ALAND = ALAND10,
    ALAND_sqmi = ALAND10 * 0.0000003861
  )

df <- merge(df, land_area, by = "GEOID")
df$acs_density <- df$acs / df$ALAND_sqmi
df$pep_density <- df$pep / df$ALAND_sqmi

# Column order: 2020 wDensity minus wp/wp_un
df <- df %>%
  select(GEOID, state, acs, census, pep, NAME, ALAND, ALAND_sqmi, acs_density, pep_density)

#### Adjacency matched to final 2009-aligned set ####
adj_out <- build_adjacency_for_df(df)
df <- adj_out$df
adjacency <- adj_out$adjacency

save(df, adjacency, file = file.path(data_dir, "census_ACS_PEP_2010era_wDensity.RData"))
cat("Saved full-pop:", file.path(data_dir, "census_ACS_PEP_2010era_wDensity.RData"),
    " nrow=", nrow(df), "\n")

df_2010_full <- df
adjacency_full <- adjacency


#### American Indian data ####
ACS_AIAN <- get_acs(
  geography = "county",
  variables = "B02001_004",
  year = 2009
) %>%
  select(GEOID, ACS = estimate) %>%
  mutate(GEOID = pad_geoid(GEOID))

census_AIAN <- get_decennial(
  geography = "county",
  variables = "P003004",  # AIAN alone, 2010 SF1
  year = 2010
) %>%
  select(GEOID, census = value) %>%
  mutate(GEOID = pad_geoid(GEOID))

# Vintage PEP race/ethnicity demographics for AIAN checks + production pep
pep_demo <- read_csv(
  file.path(data_dir, "pep_2000_2010_vintage_precensus_demographics.csv"),
  show_col_types = FALSE
)

pep_demo_cty_2009 <- pep_demo %>%
  filter(
    as.character(SUMLEV) %in% c("50", "050"),
    AGEGRP == 0,
    YEAR == 12
  ) %>%
  mutate(GEOID = sprintf("%02d%03d", as.integer(STATE), as.integer(COUNTY)))

# Race-alone groups should sum to TOT_POP (modified race; no SOR)
pep_demo_cty_2009 <- pep_demo_cty_2009 %>%
  mutate(
    pep_IA_alone = IA_MALE + IA_FEMALE,
    pep_IAC = IAC_MALE + IAC_FEMALE,
    pep_NHIA_alone = NHIA_MALE + NHIA_FEMALE,
    pep_HIA_alone = HIA_MALE + HIA_FEMALE,
    race_alone_sum = (WA_MALE + WA_FEMALE + BA_MALE + BA_FEMALE +
                        IA_MALE + IA_FEMALE + AA_MALE + AA_FEMALE +
                        NA_MALE + NA_FEMALE + TOM_MALE + TOM_FEMALE)
  )

cat("\n==== AIAN / PEP vintage demographic totals (YEAR=12, AGEGRP=0, all counties in CSV) ====\n")
cat(sprintf("  TOT_POP:              %s\n", format(sum(pep_demo_cty_2009$TOT_POP), big.mark = ",")))
cat(sprintf("  race_alone_sum:       %s\n", format(sum(pep_demo_cty_2009$race_alone_sum), big.mark = ",")))
cat(sprintf("  TOT_POP - race_alone: %s\n",
            format(sum(pep_demo_cty_2009$TOT_POP) - sum(pep_demo_cty_2009$race_alone_sum), big.mark = ",")))
cat(sprintf("  IA alone (IA_M+F):    %s\n", format(sum(pep_demo_cty_2009$pep_IA_alone), big.mark = ",")))
cat(sprintf("  IAC (alone or combo): %s\n", format(sum(pep_demo_cty_2009$pep_IAC), big.mark = ",")))
cat(sprintf("  NHIA alone:           %s\n", format(sum(pep_demo_cty_2009$pep_NHIA_alone), big.mark = ",")))
cat(sprintf("  HIA alone:            %s\n", format(sum(pep_demo_cty_2009$pep_HIA_alone), big.mark = ",")))
cat(sprintf("  NHIA+HIA (check=IA):  %s\n",
            format(sum(pep_demo_cty_2009$pep_NHIA_alone + pep_demo_cty_2009$pep_HIA_alone), big.mark = ",")))
cat(sprintf("  max |TOT_POP - race_alone_sum| over counties: %s\n",
            format(max(abs(pep_demo_cty_2009$TOT_POP - pep_demo_cty_2009$race_alone_sum)), big.mark = ",")))

PEP_AIAN <- pep_demo_cty_2009 %>%
  select(GEOID, PEP = pep_IA_alone)

cat("AIAN vintage PEP YEAR==12 county rows:", nrow(PEP_AIAN),
    " national sum (IA alone):", sum(PEP_AIAN$PEP), "\n")

df_AIAN <- merge(
  df_2010_full %>% select(GEOID, state, NAME),
  census_AIAN,
  by = "GEOID"
) %>%
  merge(PEP_AIAN, by = "GEOID") %>%
  merge(ACS_AIAN, by = "GEOID") %>%
  select(GEOID, state, acs = ACS, census, pep = PEP, NAME)

# Same GEOID universe as full-pop intersection; use matching NAME spellings
overlap_aian <- print_geoid_overlap(df_AIAN$GEOID, df_2020_full$GEOID, label = "AIAN vs 2019 full GEOIDs")
df_AIAN <- df_AIAN %>%
  filter(GEOID %in% df_2010_full$GEOID) %>%
  select(-NAME) %>%
  merge(df_2010_full %>% select(GEOID, NAME), by = "GEOID")
cat("AIAN rows after restricting to full-pop 2010era GEOIDs:", nrow(df_AIAN), "\n")
cat("AIAN GEOIDs dropped vs full-pop set:",
    length(setdiff(df_2010_full$GEOID, df_AIAN$GEOID)), "\n")

# Densities from 2009/10 AIAN pops / 2010 land; proportions vs full 2009 pops
df_AIAN <- merge(df_AIAN, land_area, by = "GEOID")
df_AIAN$acs_density <- df_AIAN$acs / df_AIAN$ALAND_sqmi
df_AIAN$pep_density <- df_AIAN$pep / df_AIAN$ALAND_sqmi

df2 <- merge(
  df_2010_full %>%
    select(GEOID, acs_full = acs, pep_full = pep,
           acs_density_full = acs_density, pep_density_full = pep_density),
  df_AIAN %>% select(GEOID, acs, pep),
  by = "GEOID"
) %>%
  mutate(
    acs_AIAN_proportion = acs / acs_full,
    pep_AIAN_proportion = pep / pep_full
  )

df_AIAN <- merge(
  df_AIAN,
  df2 %>% select(GEOID, acs_AIAN_proportion, pep_AIAN_proportion,
                 acs_density_full, pep_density_full),
  by = "GEOID"
)

df <- df_AIAN %>%
  select(
    GEOID, state, acs, census, pep, NAME,
    ALAND, ALAND_sqmi, acs_density, pep_density,
    acs_AIAN_proportion, pep_AIAN_proportion,
    acs_density_full, pep_density_full
  )

# Restrict to the same SW subset counties used in the 2020 AIAN analysis (111 counties).
# Reuse the exact GEOID universe from the saved 2020 subset so the panel matches across decades.
subset_env <- new.env()
load(file.path(data_dir, "census_ACS_PEP_WP_AIANsubset_wDensity_and2018_05212025.RData"),
     envir = subset_env)
subset_geoids <- pad_geoid(subset_env$df$GEOID)
cat("2020 AIAN subset GEOIDs:", length(subset_geoids), "\n")

df <- df %>% filter(GEOID %in% subset_geoids)
cat("2010 AIAN rows kept after SW subset filter:", nrow(df), "\n")
missing_geoids <- setdiff(subset_geoids, df$GEOID)
if (length(missing_geoids) > 0) {
  cat("WARNING: subset GEOIDs not present in 2010 AIAN data:",
      paste(missing_geoids, collapse = ", "), "\n")
}

adj_aian <- build_adjacency_for_df(df)
df <- adj_aian$df
adjacency <- adj_aian$adjacency

save(df, adjacency, file = file.path(data_dir, "census_ACS_PEP_AIANsubset_2010era_wDensity.RData"))
cat("Saved AIAN subset:", file.path(data_dir, "census_ACS_PEP_AIANsubset_2010era_wDensity.RData"),
    " nrow=", nrow(df), "\n")

df_2010_aian <- df
adjacency_aian <- adjacency


#### Schema validation vs 2020 wDensity (excluding wp) ####
cat("\n==== Schema validation (wp intentionally omitted from 2010era) ====\n")
load(file.path(data_dir, "census_ACS_PEP_WP_wDensity_11152024.RData"))
df_2020_full <- df
expected_full <- setdiff(names(df_2020_full), c("wp", "wp_un"))
cat("Full colnames match 2020 minus wp/wp_un:",
    identical(names(df_2010_full), expected_full), "\n")
cat("  2010era:", paste(names(df_2010_full), collapse = ", "), "\n")
cat("  expected:", paste(expected_full, collapse = ", "), "\n")

load(file.path(data_dir, "census_ACS_PEP_WP_AIAN_wDensity_11152024.RData"))
df_2020_aian <- df
expected_aian <- setdiff(names(df_2020_aian), "wp")
cat("AIAN colnames match 2020 minus wp:",
    identical(names(df_2010_aian), expected_aian), "\n")
cat("  2010era:", paste(names(df_2010_aian), collapse = ", "), "\n")
cat("  expected:", paste(expected_aian, collapse = ", "), "\n")

stopifnot(identical(rownames(adjacency_full), df_2010_full$NAME))
stopifnot(identical(rownames(adjacency_aian), df_2010_aian$NAME))


#### Within-era correlations + cross-decade ACS/PEP/census table ####
cat("\n==== Within-era correlation: 2009/10 full-pop ====\n")
print(round(cor(df_2010_full[, c("acs", "pep", "census")], use = "complete.obs"), 4))
cat("\n==== Within-era correlation: 2019/20 full-pop ====\n")
print(round(cor(df_2020_full[, c("acs", "pep", "census")], use = "complete.obs"), 4))

cat("\n==== Within-era correlation: 2009/10 AIAN ====\n")
print(round(cor(df_2010_aian[, c("acs", "pep", "census")], use = "complete.obs"), 4))
cat("\n==== Within-era correlation: 2019/20 AIAN ====\n")
print(round(cor(df_2020_aian[, c("acs", "pep", "census")], use = "complete.obs"), 4))

cat("\n==== 2009/10 AIAN totals on saved modeling set (n=", nrow(df_2010_aian), ") ====\n", sep = "")
cat(sprintf("  census AIAN alone (P003004): %s\n", format(sum(df_2010_aian$census), big.mark = ",")))
cat(sprintf("  ACS AIAN alone (B02001_004): %s\n", format(sum(df_2010_aian$acs), big.mark = ",")))
cat(sprintf("  PEP AIAN alone (IA_M+F):     %s\n", format(sum(df_2010_aian$pep), big.mark = ",")))
pep_chk <- pep_demo_cty_2009 %>%
  filter(GEOID %in% df_2010_aian$GEOID)
cat(sprintf("  PEP IAC (alone+combo) same GEOIDs: %s\n", format(sum(pep_chk$pep_IAC), big.mark = ",")))
cat(sprintf("  PEP NHIA alone same GEOIDs:        %s\n", format(sum(pep_chk$pep_NHIA_alone), big.mark = ",")))
cat(sprintf("  PEP as %% of full pep:              %.3f%%\n",
            100 * sum(df_2010_aian$pep) / sum(df_2010_full$pep)))
cat(sprintf("  ACS as %% of full acs:              %.3f%%\n",
            100 * sum(df_2010_aian$acs) / sum(df_2010_full$acs)))

compare_decades(df_2010_full, df_2020_full, label = "full-pop")
compare_decades(df_2010_aian, df_2020_aian, label = "AIAN")

cat("\nDone.\n")
