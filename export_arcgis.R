library(dplyr)
library(stringr)
library(tidyr)

setwd("~/Desktop/MASTER THESIS")

# ── Chargement données brutes ──────────────────────────────────────────────
data1 <- read.csv("raw/17100141-eng/17100141.csv") %>%
  select(REF_DATE, GEO, Geography.of.destination, VALUE) %>%
  filter(!grepl("Ontario part|Quebec part|Alberta part|Saskatchewan part|New Brunswick part", GEO))

data2_raw <- read.csv("raw/34100133-eng/34100133.csv")

data2_two <- data2_raw %>%
  filter(`Type.of.structure` == "Apartment structures of three units and over",
         `Type.of.unit`      == "Two bedroom units")

# ── Crosswalk ──────────────────────────────────────────────────────────────
crosswalk <- data.frame(
  GEO_cmhc = c(
    "Kitchener-Cambridge-Waterloo, Ontario",
    "St. Catharines-Niagara, Ontario",
    "Ottawa-Gatineau, Ontario/Quebec",
    "Ottawa-Gatineau, Ontario part, Ontario/Quebec",
    "Ottawa-Gatineau, Quebec part, Ontario/Quebec",
    "Hawkesbury, Ontario part, Ontario/Quebec",
    "Hawkesbury, Quebec part, Ontario/Quebec",
    "Abbotsford-Mission, British Columbia",
    "Campbellton, New Brunswick part, New Brunswick/Quebec",
    "Campbellton, Quebec part, New-Brunswick/Quebec",
    "Portage La Prairie, Manitoba",
    "Lloydminster, Saskatchewan part, Saskatchewan/Alberta",
    "Lloydminster, Alberta part, Saskachewan/Alberta"
  ),
  GEO_statcan = c(
    "Kitchener - Cambridge - Waterloo (CMA), Ontario",
    "St. Catharines - Niagara (CMA), Ontario",
    "Ottawa - Gatineau (CMA), Ontario/Quebec",
    "Ottawa - Gatineau (CMA), Ontario part, Ontario",
    "Ottawa - Gatineau (CMA), Quebec part, Quebec",
    "Hawkesbury (CA), Ontario part, Ontario",
    "Hawkesbury (CA), Quebec part, Quebec",
    "Abbotsford - Mission (CMA), British Columbia",
    "Campbellton (CA), New Brunswick part, New Brunswick",
    "Campbellton (CA), Quebec part, Quebec",
    "Portage la Prairie (CA), Manitoba",
    "Lloydminster (CA), Saskatchewan part, Saskatchewan",
    "Lloydminster (CA), Alberta part, Alberta"
  )
)

# ── Liste des 145 villes communes ─────────────────────────────────────────
data_geos <- data1 %>%
  distinct(GEO) %>%
  filter(!grepl("Area outside", GEO)) %>%
  mutate(GEO_clean = GEO %>% str_remove("\\s*\\(C[MA]+\\)") %>% str_trim())

cities_geos <- data.frame(GEO_cmhc = unique(data2_two$GEO)) %>%
  left_join(crosswalk, by = "GEO_cmhc") %>%
  mutate(GEO_clean = ifelse(
    is.na(GEO_statcan),
    GEO_cmhc %>% str_remove(",.*") %>% str_trim() %>%
      paste0(", ", str_extract(GEO_cmhc, "[^,]+$") %>% str_trim()),
    GEO_statcan %>% str_remove("\\s*\\(C[MA]+\\)") %>% str_trim()
  ))

in_both_clean <- cities_geos %>%
  filter(GEO_clean %in% data_geos$GEO_clean) %>%
  select(GEO_cmhc, GEO_clean) %>%
  mutate(name_cancensus = GEO_clean %>% str_remove(",.*") %>% str_trim())

# ── Coordonnées ────────────────────────────────────────────────────────────
coords <- readRDS("processed/coords.rds")

wl_coord <- coords %>% filter(name_cancensus == "Williams Lake")
wl_lon   <- wl_coord$longitude
wl_lat   <- wl_coord$latitude

# ── Parts migratoires WL 2016/2017 ────────────────────────────────────────
shares_wl <- data1 %>%
  filter(GEO == "Williams Lake (CA), British Columbia",
         REF_DATE == "2016/2017",
         !grepl("Area outside", Geography.of.destination)) %>%
  mutate(
    total_out = sum(VALUE, na.rm = TRUE),
    share_wl  = VALUE / total_out,
    name_cancensus = Geography.of.destination %>%
      str_remove("\\s*\\(C[MA]+\\)") %>% str_remove(",.*") %>% str_trim()
  ) %>%
  select(name_cancensus, share_wl) %>%
  distinct(name_cancensus, .keep_all = TRUE)

# ── Panel instrument WL ────────────────────────────────────────────────────
panel_instrument_wl <- expand.grid(
  name_cancensus = unique(in_both_clean$name_cancensus),
  REF_DATE       = unique(data2_two$REF_DATE)
) %>%
  as_tibble() %>%
  left_join(shares_wl, by = "name_cancensus") %>%
  mutate(
    Post2017  = as.integer(REF_DATE >= 2017),
    share_wl  = ifelse(name_cancensus == "Williams Lake", 0, share_wl),
    share_wl  = ifelse(name_cancensus %in% c("Campbellton", "Hawkesbury"), 0, share_wl),
    Z_hat_wl  = replace_na(share_wl, 0) * Post2017
  )

# ── Cylindrage (pour obtenir exactement les 131 villes) ───────────────────
name_bridge <- in_both_clean %>% select(name_cancensus, GEO_cmhc)

rent_biprov <- data2_two %>%
  filter(grepl("Campbellton|Hawkesbury", GEO)) %>%
  mutate(name_cancensus = ifelse(grepl("Campbellton", GEO), "Campbellton", "Hawkesbury")) %>%
  group_by(name_cancensus, REF_DATE) %>%
  summarise(rent = mean(VALUE, na.rm = TRUE), .groups = "drop") %>%
  mutate(rent = ifelse(is.nan(rent), NA, rent))

panel_merged <- panel_instrument_wl %>%
  left_join(name_bridge, by = "name_cancensus") %>%
  left_join(data2_two %>% select(REF_DATE, GEO, VALUE) %>% rename(GEO_cmhc = GEO, rent = VALUE),
            by = c("GEO_cmhc", "REF_DATE")) %>%
  left_join(rent_biprov, by = c("name_cancensus", "REF_DATE")) %>%
  mutate(rent = ifelse(is.na(rent.x), rent.y, rent.x)) %>%
  select(-rent.x, -rent.y) %>%
  filter(!name_cancensus %in% c("Wasaga Beach", "Whitehorse")) %>%
  filter(REF_DATE >= 2008, REF_DATE <= 2024) %>%
  mutate(rent = ifelse(rent == 0, NA, rent), log_rent = log(rent))

villes_na <- panel_merged %>%
  group_by(name_cancensus) %>%
  summarise(n_NA = sum(is.na(log_rent)), .groups = "drop") %>%
  filter(n_NA > 0)

panel_balanced_wl <- panel_merged %>%
  filter(!name_cancensus %in% villes_na$name_cancensus[
    villes_na$name_cancensus != "Leamington"])

cat("Villes dans le panel cylindré :", length(unique(panel_balanced_wl$name_cancensus)), "\n")

# ── Export CSV pour ArcGIS ─────────────────────────────────────────────────
fig_d4_data <- panel_balanced_wl %>%
  distinct(name_cancensus, share_wl) %>%
  left_join(coords %>% select(name_cancensus, longitude, latitude),
            by = "name_cancensus") %>%
  filter(!is.na(longitude), !is.na(latitude)) %>%
  mutate(
    share_wl       = as.numeric(replace_na(share_wl, 0)),
    share_wl_pct   = round(as.numeric(replace_na(share_wl, 0)) * 100, 2),
    exposed        = as.integer(share_wl > 0),
    exposure_label = ifelse(share_wl > 0, "Exposed", "Non-exposed")
  )

# Williams Lake (origine — point spécial dans ArcGIS)
wl_row <- tibble(
  name_cancensus = "Williams Lake",
  share_wl       = NA_real_,
  longitude      = wl_lon,
  latitude       = wl_lat,
  exposed        = NA_integer_,
  exposure_label = "Origin (Williams Lake)"
)

arcgis_export <- bind_rows(fig_d4_data, wl_row)

write.csv(arcgis_export, "cities_exposure_arcgis.csv", row.names = FALSE)

cat("Fichier exporté : cities_exposure_arcgis.csv\n")
cat("Lignes :", nrow(arcgis_export), "(", sum(arcgis_export$exposed == 1, na.rm=TRUE),
    "exposées,", sum(arcgis_export$exposed == 0, na.rm=TRUE), "non-exposées, 1 origine)\n")
