# =============================================================================
# THE IMPACT OF CLIMATE-INDUCED DISPLACEMENT ON RENTAL PRICES
# Evidence from Canada
# Robin Masson -- Université Paris 1 Panthéon-Sorbonne
# Mémoire M2 Économie
# =============================================================================
#
# QUESTION DE RECHERCHE :
# Le feu de forêt de Williams Lake (Colombie-Britannique, juillet 2017) a-t-il
# provoqué une hausse des loyers dans les villes canadiennes d'accueil des
# déplacés ?
#
# STRATÉGIE EMPIRIQUE :
# Bartik-style network exposure index, motivated by shift-share logic
# (Goldsmith-Pinkham, Sorkin & Swift, AER 2020).
#   Z_WL = share_wl × Post2017
# où share_wl = part migratoire pré-feu de chaque ville de destination
# (flux Williams Lake → ville d, période 2016/2017) et Post2017 = 1(année ≥ 2017).
#
# Le choc agrégé (l'incendie) est plausiblement exogène au calendrier des
# marchés locatifs de destination. L'interprétation de Z_WL comme mesure
# d'exposition reduced-form exige toutefois que les parts pré-feu ne soient pas
# corrélées aux trajectoires contrefactuelles de loyer.
#
# RÉSULTAT PRINCIPAL :
# β = 0.818*** (FE ville + année, 131 villes, 2008-2024).
# Les placebos révèlent une tendance pré-existante en BC → interprétation
# reduced-form uniquement (cf. DIAGNOSTIC FINAL, PARTIE IX).
#
# NOTE D'INTERPRÉTATION CAUSALE :
# "The coefficient does not estimate the effect of one additional migrant on
#  rents. Instead, it captures the reduced-form effect of exposure to the
#  wildfire shock through pre-existing migration networks."
#
# =============================================================================
# CONDITIONS D'INTERPRÉTATION DU DESIGN D'EXPOSITION
# =============================================================================
#
# 1. ACTIVATION DU CANAL MIGRATOIRE :
#    Les parts pré-feu doivent être économiquement pertinentes comme mesure de
#    l'exposition potentielle au choc. Intuition : elles captent les liens
#    migratoires historiques vers les villes de destination.
#
# 2. EXOGÉNÉITÉ CONDITIONNELLE DES PARTS PRÉ-FEU :
#    Conditional on city and year fixed effects, share_wl ne doit pas capter
#    des trajectoires de loyer contrefactuelles différentes pour des raisons
#    indépendantes du feu. LIMITE : les villes exposées sont concentrées en BC,
#    qui connaissait une hausse des loyers supérieure à la moyenne canadienne
#    avant 2017 → confounding provincial possible.
#
# 3. PARALLEL TRENDS (pour l'interprétation causale) :
#    En l'absence du feu, les villes plus et moins exposées auraient dû suivre
#    des trajectoires de loyer similaires. Les placebos (PARTIE VII) montrent
#    que cette hypothèse est violée dans la spécification baseline →
#    reduced-form uniquement.
#
# 4. EXPOSITIONS WL vs FM :
#    Z_WL utilise des parts OBSERVÉES, strictement pré-feu.
#    Z_FM utilise des parts PRÉDITES via un modèle de gravité (Wood Buffalo 2016)
#    pour éviter la contamination temporelle (feu mai 2016 ∈ période 2015/2016).
#
# =============================================================================
# STRUCTURE DU SCRIPT
# =============================================================================
#   PARTIE I     -- PACKAGES ET CONFIGURATION
#   PARTIE II    -- CHARGEMENT DES DONNÉES
#   PARTIE III   -- NETTOYAGE ET CROSSWALK
#   PARTIE IV    -- ÉCHANTILLON ANALYTIQUE + GÉOGRAPHIES
#   PARTIE V     -- CONSTRUCTION DES INDICES D'EXPOSITION (WL + FM)
#   PARTIE VI    -- CONSTRUCTION DES PANELS
#   PARTIE VII   -- VÉRIFICATIONS D'IDENTIFICATION (event study, placebos)
#   PARTIE VIII  -- RÉSULTATS PRINCIPAUX
#   PARTIE IX    -- VÉRIFICATIONS DE ROBUSTESSE
#     IX.1    Expositions alternatives (FM)
#     IX.2    Tendances spécifiques et échantillon BC (+ timing t>=2018)
#     IX.3    Hard Checks Vancouver + trimming
#     IX.4    Code préparé (données externes requises)
#     IX.4bis Échantillon tronqué 2019 (robustesse COVID)
#     IX.4ter Tendances estimées sur période pré-traitement uniquement
#     IX.5    Placebo Quesnel
#     IX.5bis Placebo distant — Greater Sudbury
#     IX.6    Gravity shares pour WL
#     IX.7    Tests complémentaires (LOO, clustering, hétérogénéité, permutation,
#             poids Rotemberg, exports LaTeX)
#   PARTIE X     -- ANALYSE D'HÉTÉROGÉNÉITÉ (4×4)
#   PARTIE XI    -- CALIBRATION ÉCONOMIQUE
#   PARTIE XI.bis-- TABLES 4 ET 5 (calibration + market tightness)
#   PARTIE XII   -- STATISTIQUES DESCRIPTIVES
#   PARTIE XII.bis -- FIGURES DESCRIPTIVES (D1–D4)
#   PARTIE XII.5 -- CARTE D'EXPOSITION SPATIALE (fig_spatial_exposure_qgis)
#   PARTIE XIII  -- EXPORT TABLEAUX ET FIGURES
#   PARTIE XIV   -- INFÉRENCE HONNÊTE RAMBACHAN-ROTH [mis en commentaire]
#   PARTIE XV    -- ABSORPTION RÉGIONALE + TIMING t>=2018 (Section 7.4–7.5)
#   PARTIE XVI   -- VALIDATION DU CANAL MIGRATOIRE (Section 5)
#   PARTIE XVII  -- 2SLS EXPLORATOIRE SUR FENÊTRE LIMITÉE
# =============================================================================


# =============================================================================
# PARTIE I -- PACKAGES ET CONFIGURATION
# =============================================================================
# Chargement de tous les packages nécessaires à l'analyse. install.packages()
# n'est à exécuter qu'une seule fois ; commenter ensuite pour accélérer le
# rechargement de la session.

# install.packages(c("dplyr", "stringr", "tidyr", "cancensus",
#                    "sf", "geosphere", "fixest",
#                    "ggplot2", "ggrepel", "scales",
#                    "modelsummary", "purrr", "broom",
#                    "rnaturalearth", "rnaturalearthdata"))

library(dplyr)
library(stringr)
library(tidyr)
library(cancensus)
library(sf)
library(geosphere)
library(fixest)
library(ggplot2)
library(ggrepel)
library(scales)
library(modelsummary)
library(purrr)
library(broom)
library(rnaturalearth)
library(rnaturalearthdata)

# Répertoire de travail -- ajuster si nécessaire
setwd("~/Desktop/MASTER THESIS")

# Clé API CensusMapper (inscription gratuite : https://censusmapper.ca)
# → décommenter uniquement pour régénérer les coordonnées géographiques
# set_cancensus_api_key(Sys.getenv("CANCENSUS_API_KEY"),
#                       install   = TRUE,
#                       overwrite = TRUE)


# =============================================================================
# PARTIE II -- CHARGEMENT DES DONNÉES
# =============================================================================
# data1 : Migrations inter-cités, StatCan 17-10-0141-01
#         173 géographies (CMA/CA), périodes juillet-juin 2016/2017→2020/2021
# data2 : Loyers SCHL, StatCan 34-10-0133-01
#         244 villes (CMA, CA, CSD 10 000+), 1987-2025, annuel
# Ces deux sources sont les seules disponibles pour construire un panel
# ville × année reliant migration et loyer au niveau sub-national canadien.

data1 <- read.csv("raw/17100141-eng/17100141.csv")
data2 <- read.csv("raw/34100133-eng/34100133.csv")


# =============================================================================
# PARTIE III -- NETTOYAGE ET CROSSWALK
# =============================================================================

# -----------------------------------------------------------------------------
# III.1 Nettoyage de data1 (migrations)
# -----------------------------------------------------------------------------
# On supprime les variables internes StatCan (constantes ou identifiants)
# et on exclut les sous-géographies biprovinciates pour éviter la
# double-comptabilisation avec les versions agrégées.

data1 <- data1 %>%
  select(REF_DATE, GEO, `Geography.of.destination`, VALUE)

# Exclusion des sous-géographies biprovinciates
data1 <- data1 %>%
  filter(!grepl("Ontario part|Quebec part|Alberta part|Saskatchewan part|New Brunswick part", GEO))

# Vérification : Lloydminster ne doit conserver que la version agrégée
data1 %>% distinct(GEO) %>% filter(grepl("Lloydminster", GEO))
# → "Lloydminster (CA), Alberta/Saskatchewan" uniquement ✓


# -----------------------------------------------------------------------------
# III.2 Nettoyage de data2 (loyers SCHL)
# -----------------------------------------------------------------------------
# CHOIX DU TYPE DE STRUCTURE (baseline) :
# "Apartment structures of three units and over"
# Justification : Saiz (2006, IZA DP 2189, section Data p.11) utilise les
# HUD Fair Market Rents définis comme le loyer d'une unité 2 chambres vacante
# au 45ème percentile -- mesure du marché locatif primaire standard.
# Parmi les 4 structures SCHL, "Apartment structures of three units and over"
# est le segment le plus comparable : immeubles d'appartements de taille
# standard, meilleure couverture géographique, excluant les maisons en rangée.
#
# CHOIX DU TYPE D'UNITÉ (baseline) :
# 2 chambres -- même référence que Saiz (2006) qui utilise explicitement
# le loyer d'une unité 2 chambres (FMR, 45ème percentile).
#
# Les 3 autres structures et les 4 types d'unité servent aux tests
# d'hétérogénéité (PARTIE X) selon la grille 4 structures × 4 types.
#
# Les 4 structures disponibles dans SCHL :
#   S1 : "Apartment structures of three units and over"  ← BASELINE
#   S2 : "Row and apartment structures of three units and over"
#   S3 : "Row structures of three units and over"
#   S4 : "Apartment structures of six units and over"

data2 <- data2 %>%
  select(REF_DATE, GEO, `Type.of.structure`, `Type.of.unit`, VALUE)

# Fonction de création d'un sous-dataset par structure × type d'unité
# Garantit une observation unique par ville-année (filtre double)
make_data2 <- function(structure, unit) {
  data2 %>%
    filter(`Type.of.structure` == structure,
           `Type.of.unit`      == unit)
}

# --- BASELINE (S1 × 2 chambres) ---
data2_two <- make_data2("Apartment structures of three units and over",
                        "Two bedroom units")
stopifnot(nrow(data2_two) == nrow(distinct(data2_two, GEO, REF_DATE)))  # ✓

# --- Structure S1 : Apartment structures of three units and over ---
data2_s1_bachelor <- make_data2("Apartment structures of three units and over", "Bachelor units")
data2_s1_one      <- make_data2("Apartment structures of three units and over", "One bedroom units")
data2_s1_two      <- data2_two   # déjà créé
data2_s1_three    <- make_data2("Apartment structures of three units and over", "Three bedroom units")

# --- Structure S2 : Row and apartment structures of three units and over ---
data2_s2_bachelor <- make_data2("Row and apartment structures of three units and over", "Bachelor units")
data2_s2_one      <- make_data2("Row and apartment structures of three units and over", "One bedroom units")
data2_s2_two      <- make_data2("Row and apartment structures of three units and over", "Two bedroom units")
data2_s2_three    <- make_data2("Row and apartment structures of three units and over", "Three bedroom units")

# --- Structure S3 : Row structures of three units and over ---
data2_s3_bachelor <- make_data2("Row structures of three units and over", "Bachelor units")
data2_s3_one      <- make_data2("Row structures of three units and over", "One bedroom units")
data2_s3_two      <- make_data2("Row structures of three units and over", "Two bedroom units")
data2_s3_three    <- make_data2("Row structures of three units and over", "Three bedroom units")

# --- Structure S4 : Apartment structures of six units and over ---
data2_s4_bachelor <- make_data2("Apartment structures of six units and over", "Bachelor units")
data2_s4_one      <- make_data2("Apartment structures of six units and over", "One bedroom units")
data2_s4_two      <- make_data2("Apartment structures of six units and over", "Two bedroom units")
data2_s4_three    <- make_data2("Apartment structures of six units and over", "Three bedroom units")


# -----------------------------------------------------------------------------
# III.3 Crosswalk data1 ↔ data2
# -----------------------------------------------------------------------------
# Les deux sources utilisent des conventions de nommage différentes pour les
# mêmes géographies. Ce crosswalk manuel assure une jointure propre.

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
    "Lloydminster, Alberta part, Saskachewan/Alberta"  # sic -- faute de frappe dans source SCHL
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


# =============================================================================
# PARTIE IV -- ÉCHANTILLON ANALYTIQUE + GÉOGRAPHIES
# =============================================================================

# -----------------------------------------------------------------------------
# IV.1 Intersection data1 ∩ data2 (145 villes)
# -----------------------------------------------------------------------------
# Seules les villes présentes dans les deux sources sont utilisables :
# il faut à la fois la variable de traitement (flux de migration depuis WL)
# et la variable dépendante (loyer SCHL).

data_geos <- data1 %>%
  distinct(GEO) %>%
  filter(!grepl("Area outside", GEO)) %>%
  mutate(GEO_clean = GEO %>%
           str_remove("\\s*\\(C[MA]+\\)") %>%
           str_trim())

cities_geos <- data.frame(GEO_cmhc = unique(data2$GEO)) %>%
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

cat("Échantillon analytique final :", nrow(in_both_clean), "villes\n")
# → 145 villes à l'intersection de data1 (migrations) et data2 (loyers SCHL) ✓


# -----------------------------------------------------------------------------
# IV.2 Coordonnées géographiques (CensusMapper / Recensement 2021)
# -----------------------------------------------------------------------------
# Les coordonnées sont nécessaires pour le modèle de gravité (exposition FM gravitaire).
# Le bloc API est commenté : les coordonnées ont été générées une seule fois
# et sauvegardées dans processed/coords.rds.
# =========================================================================
# BLOC API CANCENSUS -- COMMENTÉ POUR NE PAS GRILLER LES TOKENS
# coords a été généré via get_census() (CA21, niveau CMA, geo_format = "sf")
# puis enrichi manuellement pour 3 villes hors CA21 (Bay Roberts, Leamington,
# Cold Lake) et 2 villes biprovinciates (Campbellton, Hawkesbury).
# =========================================================================

# regions_cma <- get_census(
#   dataset    = "CA21",
#   regions    = list(C = "01"),
#   level      = "CMA",
#   geo_format = "sf"
# ) %>%
#   mutate(
#     centroid   = st_centroid(geometry),
#     longitude  = st_coordinates(centroid)[, 1],
#     latitude   = st_coordinates(centroid)[, 2],
#     name_clean = name %>% str_remove("\\s*\\([A-Z]\\)$") %>% str_trim()
#   )
#
# regions_cma_matched <- regions_cma %>%
#   mutate(name_cancensus = case_when(
#     name_clean == "Belleville - Quinte West"
#     ~ "Belleville",
#     name_clean == "Greater Sudbury / Grand Sudbury"
#     ~ "Greater Sudbury",
#     name_clean == "Campbellton (New Brunswick part / partie du Nouveau-Brunswick)"
#     ~ "Campbellton, New Brunswick part",
#     name_clean == "Campbellton (partie du Québec / Quebec part)"
#     ~ "Campbellton, Quebec part",
#     name_clean == "Hawkesbury (Ontario part / partie de l'Ontario)"
#     ~ "Hawkesbury, Ontario part",
#     name_clean == "Hawkesbury (partie du Québec / Quebec part)"
#     ~ "Hawkesbury, Quebec part",
#     name_clean == "Lloydminster (Saskatchewan part / partie de la Saskatchewan)"
#     ~ "Lloydminster, Saskatchewan part",
#     TRUE ~ name_clean
#   )) %>%
#   filter(name_cancensus %in% in_both_clean$name_cancensus) %>%
#   select(name_cancensus, longitude, latitude, Population)
#
# coords_manual : 3 villes non couvertes au niveau CMA dans CA21
# (subdivisions de recensement isolées). Coordonnées et populations relevées
# manuellement depuis les profils du Recensement 2021 de Statistique Canada
# (catalogue 98-316-X2021001, profils communautaires).
# Source -- Bay Roberts: Statistics Canada, 2021 Census Community Profile.
# Source -- Leamington: Statistics Canada, 2021 Census Community Profile.
# Source -- Cold Lake: Statistics Canada, 2021 Census Community Profile.
# Vérifié le 2026-05-15.
# coords_manual <- data.frame(
#   name_cancensus = c("Bay Roberts", "Leamington", "Cold Lake"),
#   longitude      = c(-53.2648,     -82.5997,     -110.1825),
#   latitude       = c(47.5953,       42.0534,       54.4642),
#   Population     = c(12993,         28403,         17990)
# )
#
# coords_biprov : 2 agglomérations biprovinciales. Le polygone CA21 existe
# au niveau de chaque partie provinciale, pas pour l'agglomération fusionnée
# utilisée dans la table 17-10-0141. Coordonnées des centroïdes
# d'agglomération calculées manuellement à partir des composantes provinciales
# de la géographie du Recensement 2021.
# Source -- Campbellton: centroid reconstructed from the provincial CA21
# geography components to match the aggregated migration-table unit.
# Source -- Hawkesbury: centroid reconstructed from the provincial CA21
# geography components to match the aggregated migration-table unit.
# coords_biprov <- data.frame(
#   name_cancensus = c("Campbellton", "Hawkesbury"),
#   longitude      = c(-66.75121,    -74.60900),
#   latitude       = c(47.83301,      45.61337),
#   Population     = c(13330,         12010)
# )
#
# coords <- bind_rows(regions_cma_matched, coords_manual, coords_biprov) %>%
#   st_drop_geometry() %>%
#   select(name_cancensus, longitude, latitude, Population)
# dir.create("processed", showWarnings = FALSE)
# saveRDS(coords, "processed/coords.rds")

# → Charger les coordonnées depuis le fichier sauvegardé (usage normal)
coords <- readRDS("processed/coords.rds")

# Coordonnées de Wood Buffalo (origine Fort McMurray -- exposition FM gravitaire)
wb_coords <- coords %>% filter(name_cancensus == "Wood Buffalo")
wb_lon <- wb_coords$longitude
wb_lat <- wb_coords$latitude


# =============================================================================
# PARTIE V -- CONSTRUCTION DES INDICES D'EXPOSITION (WL + FM)
# =============================================================================

# -----------------------------------------------------------------------------
# V.1 Parts migratoires pré-feu BC 2017 (Williams Lake) -- share_wl
# -----------------------------------------------------------------------------
# Période pré-feu : 2016/2017 (juillet 2016 – juin 2017)
# Strictement antérieure au feu de juillet 2017 → s_d^WL = s_d^pre ✓
# (Goldsmith-Pinkham, Sorkin & Swift 2020)
# Exclusion du résidu rural "Area outside" + renormalisation
# (Borusyak, Hull & Jaravel, ReStud 2022)
#
# Nommage : share_wl = part pré-feu fixe ; Z_wl alias Z_hat_wl = share_wl × Post2017

dataset_wl_clean <- data1 %>%
  filter(GEO == "Williams Lake (CA), British Columbia",
         REF_DATE == "2016/2017",
         !grepl("Area outside", `Geography.of.destination`)) %>%
  mutate(
    total_out = sum(VALUE, na.rm = TRUE),
    share_wl  = VALUE / total_out
  )

sum(dataset_wl_clean$share_wl, na.rm = TRUE)  # doit être 1 ✓

shares_wl <- dataset_wl_clean %>%
  mutate(name_cancensus = `Geography.of.destination` %>%
           str_remove("\\s*\\(C[MA]+\\)") %>% str_remove(",.*") %>% str_trim()) %>%
  select(name_cancensus, share_wl) %>%
  distinct(name_cancensus, .keep_all = TRUE)


# -----------------------------------------------------------------------------
# V.1b Validation gravitaire directe -- Williams Lake 2016/2017
# -----------------------------------------------------------------------------
# Modèle : ln(Flux_WL->d) = β0 + β1·ln(D_WL,d) + β2·Border + β3·ln(Pop_d)
# OLS classique + PPML (Santos Silva & Tenreyro 2006)
# Objectif : confirmer que distance + population prédisent les parts s_WL,d

wl_coord <- coords %>% filter(name_cancensus == "Williams Lake")
wl_lon   <- wl_coord$longitude
wl_lat   <- wl_coord$latitude

# Données WL avec coordonnées destination + distance
gravity_wl_data <- dataset_wl_clean %>%
  mutate(name_cancensus = `Geography.of.destination` %>%
           str_remove("\\s*\\(C[MA]+\\)") %>% str_remove(",.*") %>% str_trim()) %>%
  left_join(coords, by = "name_cancensus") %>%
  filter(!is.na(longitude), !is.na(latitude), name_cancensus != "Williams Lake") %>%
  rowwise() %>%
  mutate(dist_km = distHaversine(c(wl_lon, wl_lat), c(longitude, latitude)) / 1000) %>%
  ungroup() %>%
  mutate(
    province_dest = `Geography.of.destination` %>% str_extract("(?<=, ).*$") %>% str_trim(),
    border = as.integer(province_dest != "British Columbia")
  )

# --- OLS log-linéaire ---
grav_wl_ols <- lm(log(VALUE) ~ log(dist_km) + border + log(Population),
                  data = gravity_wl_data %>% filter(dist_km > 0, VALUE > 0))
summary(grav_wl_ols)

# --- PPML (gère les zéros) ---
grav_wl_ppml <- fepois(VALUE ~ log(dist_km) + border + log(Population),
                       data = gravity_wl_data %>% filter(dist_km > 0, !is.na(Population)))
summary(grav_wl_ppml)

# --- Validation : prédit vs observé (OLS) ---
grav_wl_valid <- gravity_wl_data %>%
  filter(!is.na(dist_km), dist_km > 0, !is.na(Population)) %>%
  mutate(
    fitted_ols    = exp(predict(grav_wl_ols, newdata = .)),
    total_fitted  = sum(fitted_ols, na.rm = TRUE),
    share_wl_grav = fitted_ols / total_fitted
  ) %>%
  filter(!is.na(share_wl))

r2_wl <- cor(grav_wl_valid$share_wl, grav_wl_valid$share_wl_grav,
             use = "complete.obs")^2
n_wl  <- nrow(grav_wl_valid)
cat(sprintf("WL Gravity OLS R² = %.3f  (N = %d)\n", r2_wl, n_wl))

# --- Figure scatter WL ---
p_val_wl <- ggplot(grav_wl_valid,
                   aes(x = share_wl_grav, y = share_wl)) +
  geom_point(alpha = 0.7, color = "#0072B2", size = 1.8) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
  annotate("text", x = Inf, y = -Inf,
           label = paste0("R² = ", round(r2_wl, 2), "  N = ", n_wl),
           hjust = 1.1, vjust = -0.5, size = 3.2) +
  ggrepel::geom_text_repel(
    data = grav_wl_valid %>% filter(share_wl > 0.05),
    aes(label = name_cancensus), size = 2.5, max.overlaps = 10) +
  labs(x = "Gravity-predicted migration share (Williams Lake)",
       y = "Observed migration share (WL, 2016/2017)") +
  theme_classic(base_size = 11)

ggsave("figures/validation_williams_lake.pdf", p_val_wl, width = 6.5, height = 4.5)
cat("  validation_williams_lake.pdf exported.\n")


# -----------------------------------------------------------------------------
# V.2 Modèle de gravité pour Fort McMurray 2016 -- exposition FM gravitaire
# -----------------------------------------------------------------------------
# Parts PRÉDITES car :
#   (1) Contamination temporelle : feu mai 2016 ∈ période "2015/2016"
#   (2) Endogénéité : Wood Buffalo dominée par sables bitumineux
# Modèle : ln(Flux_od) = β0 + β1·ln(D_od) + β2·Border_od + β3·ln(Pop_d)

gravity_data <- data1 %>%
  filter(!grepl("Area outside", GEO),
         !grepl("Area outside", `Geography.of.destination`),
         GEO != "Wood Buffalo (CA), Alberta",
         !is.na(VALUE), VALUE > 0) %>%
  mutate(
    GEO_clean_orig = GEO %>%
      str_remove("\\s*\\(C[MA]+\\)") %>% str_remove(",.*") %>% str_trim(),
    GEO_clean_dest = `Geography.of.destination` %>%
      str_remove("\\s*\\(C[MA]+\\)") %>% str_remove(",.*") %>% str_trim()
  ) %>%
  left_join(coords %>% rename(GEO_clean_orig = name_cancensus,
                              lon_orig = longitude, lat_orig = latitude,
                              pop_orig = Population),
            by = "GEO_clean_orig") %>%
  left_join(coords %>% rename(GEO_clean_dest = name_cancensus,
                              lon_dest = longitude, lat_dest = latitude,
                              pop_dest = Population),
            by = "GEO_clean_dest") %>%
  filter(!GEO_clean_orig %in% c("Arnprior", "Carleton Place", "Winkler", "Nelson"),
         !GEO_clean_dest %in% c("Arnprior", "Carleton Place", "Winkler", "Nelson")) %>%
  rowwise() %>%
  mutate(dist_km = distHaversine(c(lon_orig, lat_orig),
                                 c(lon_dest, lat_dest)) / 1000) %>%
  ungroup() %>%
  mutate(
    province_orig = GEO %>% str_extract("(?<=, ).*$") %>% str_trim(),
    province_dest = `Geography.of.destination` %>%
      str_extract("(?<=, ).*$") %>% str_trim(),
    border = as.integer(province_orig != province_dest)
  )

gravity_model <- lm(log(VALUE) ~ log(dist_km) + border + log(pop_dest),
                    data = gravity_data %>% filter(dist_km > 0))
summary(gravity_model)
# Résultats OLS (N = 49 042) :
# β1 log(dist_km)   = -0.506 (p<0.001) ✓ friction géographique standard
#    → une hausse de 1% de la distance réduit les flux de 0.5%
# β2 border         = -0.339 (p<0.001) ✓ effet frontière interprovinciale
#    → franchir une frontière provinciale réduit les flux de ~29% (exp(-0.339)-1)
#    → cohérent avec Helliwell (1997) sur les barrières internes au Canada
# β3 log(pop_dest)  = +0.450 (p<0.001) ✓ attraction gravitationnelle
#    → une ville 1% plus grande attire 0.45% de migrants supplémentaires
#    → élasticité sous-unitaire : pas de rendements croissants à l'attraction
# R² = 0.32 -- satisfaisant pour un modèle de gravité en coupe transversale

gravity_ppml <- fepois(VALUE ~ log(dist_km) + border + log(pop_dest),
                       data = gravity_data %>% filter(dist_km > 0))
summary(gravity_ppml)
# Résultats PPML (N = 49 046) :
# β1 log(dist_km)   = -0.734 (p<0.001) -- plus fort qu'en OLS
# β2 border         = -0.198 (p<0.001) -- plus faible qu'en OLS (~-18%)
# β3 log(pop_dest)  = +0.603 (p<0.001) -- plus fort qu'en OLS
# Pseudo R² = 0.42 > R² OLS = 0.32 ✓ meilleur ajustement global
# PPML provides a robustness check to the log-linear gravity specification.
# However, because zero flows are excluded from the estimation sample
# (!is.na(VALUE), VALUE > 0), this PPML specification should not be interpreted
# as solving the zero-flow problem.
# → OLS retenu pour la prédiction hors-échantillon (interprétabilité des β)


# -----------------------------------------------------------------------------
# V.3 Parts prédites Fort McMurray 2016 -- share_fm
# -----------------------------------------------------------------------------
# Dictionnaire province pour les villes de destination
province_dict <- in_both_clean %>%
  mutate(province = GEO_clean %>% str_extract("(?<=, ).*$") %>% str_trim()) %>%
  select(name_cancensus, province)

predict_data <- coords %>%
  filter(name_cancensus != "Wood Buffalo") %>%
  left_join(province_dict, by = "name_cancensus") %>%
  mutate(
    dist_km = distHaversine(c(wb_lon, wb_lat), cbind(longitude, latitude)) / 1000,
    border  = as.integer(province != "Alberta")
  ) %>%
  filter(dist_km > 0) %>%
  mutate(
    flux_predicted  = exp(predict(gravity_model,
                                  newdata = data.frame(dist_km  = dist_km,
                                                       border   = border,
                                                       pop_dest = Population))),
    total_predicted = sum(flux_predicted, na.rm = TRUE),
    share_fm        = flux_predicted / total_predicted
  )

sum(predict_data$share_fm, na.rm = TRUE)  # doit être 1 ✓

shares_fm <- predict_data %>% select(name_cancensus, share_fm)

# Validation externe : parts prédites vs. parts observées (2016/2017)
shares_fm_obs <- data1 %>%
  filter(GEO == "Wood Buffalo (CA), Alberta",
         REF_DATE == "2016/2017",
         !grepl("Area outside", `Geography.of.destination`)) %>%
  mutate(
    name_cancensus = `Geography.of.destination` %>%
      str_remove("\\s*\\(C[MA]+\\)") %>% str_remove(",.*") %>% str_trim(),
    total_out    = sum(VALUE, na.rm = TRUE),
    share_fm_obs = VALUE / total_out
  ) %>%
  select(name_cancensus, share_fm_obs)

validation <- predict_data %>%
  select(name_cancensus, share_fm) %>%
  left_join(shares_fm_obs, by = "name_cancensus") %>%
  filter(!is.na(share_fm_obs))

val_model <- lm(share_fm_obs ~ share_fm, data = validation)
summary(val_model)
# Résultats de la validation externe (N = 149) :
# β = 3.638 (p<0.001) -- les parts prédites sont corrélées aux parts observées
# R² = 0.60 ✓ -- le modèle de gravité explique 60% de la variance des flux réels
# L'intercept négatif (-0.019) reflète le fait que les parts prédites sont
# plus concentrées sur Edmonton que les parts observées -- attendu car Edmonton
# est la grande ville albertaine la plus proche, surpondérée par le modèle
cor(validation$share_fm, validation$share_fm_obs, use = "complete.obs")
# r = 0.774 ✓ -- corrélation forte entre parts prédites et observées

val_model_no_edm <- lm(share_fm_obs ~ share_fm,
                       data = validation %>% filter(name_cancensus != "Edmonton"))
summary(val_model_no_edm)

validation <- validation %>%
  mutate(label = ifelse(share_fm_obs > 0.05 | share_fm > 0.04, name_cancensus, NA))

r_val   <- round(cor(validation$share_fm, validation$share_fm_obs,
                     use = "complete.obs"), 3)
r2_val  <- round(summary(val_model)$r.squared, 2)
n_val   <- nrow(validation)
ann_lab <- sprintf("italic(r) == %.3f~~italic(R)^2 == %.2f~~italic(N) == %d",
                   r_val, r2_val, n_val)

p_validation <- ggplot(validation, aes(x = share_fm, y = share_fm_obs)) +
  geom_abline(slope = 1, intercept = 0,
              linetype = "dashed", color = "gray60", linewidth = 0.6) +
  geom_smooth(method = "lm", se = TRUE,
              color = "#2166ac", fill = "#bdc9e1", linewidth = 0.85) +
  geom_point(color = "#4d4d4d", size = 2, alpha = 0.65) +
  geom_text_repel(aes(label = label),
                  size = 2.9, color = "#333333",
                  box.padding = 0.35, max.overlaps = 12,
                  segment.color = "gray70", segment.size = 0.3) +
  annotate("text", x = Inf, y = -Inf,
           label = ann_lab, parse = TRUE,
           hjust = 1.05, vjust = -0.6,
           size = 3.2, color = "gray30") +
  scale_x_continuous(labels = scales::percent_format(accuracy = 0.1),
                     expand  = expansion(mult = c(0.02, 0.05))) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 0.1),
                     expand  = expansion(mult = c(0.05, 0.05))) +
  labs(
    x = "Gravity-predicted migration share",
    y = "Observed migration share (2016/2017)"
  ) +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        axis.title        = element_text(size = 10))

ggsave("figures/validation_fort_mcmurray.pdf", plot = p_validation,
       width = 6.5, height = 5)


# -----------------------------------------------------------------------------
# V.4 Construction des panels d'exposition (WL baseline + WL+FM robustesse)
# -----------------------------------------------------------------------------
# panel_instrument_wl  : exposition Williams Lake 2017 uniquement (analyse principale)
# panel_instrument_full : expositions WL + FM (vérifications de robustesse PARTIE IX)

panel_instrument_wl <- expand.grid(
  name_cancensus = unique(in_both_clean$name_cancensus),
  REF_DATE       = unique(data2_two$REF_DATE)
) %>%
  as_tibble() %>%
  left_join(shares_wl, by = "name_cancensus") %>%
  mutate(
    Post2017  = as.integer(REF_DATE >= 2017),
    Post2018  = as.integer(REF_DATE >= 2018),
    share_wl  = ifelse(name_cancensus == "Williams Lake", 0, share_wl),
    share_wl  = ifelse(name_cancensus %in% c("Campbellton", "Hawkesbury"), 0, share_wl),
    Z_hat_wl  = replace_na(share_wl, 0) * Post2017,
    Z_hat_wl18 = replace_na(share_wl, 0) * Post2018   # Section 7.5 : timing t>=2018
  )

panel_instrument_full <- expand.grid(
  name_cancensus = unique(predict_data$name_cancensus),
  REF_DATE       = unique(data2_two$REF_DATE)
) %>%
  as_tibble() %>%
  left_join(shares_wl, by = "name_cancensus") %>%
  left_join(shares_fm,  by = "name_cancensus") %>%
  mutate(
    Post2016 = as.integer(REF_DATE >= 2016),
    Post2017 = as.integer(REF_DATE >= 2017),
    share_wl = ifelse(name_cancensus == "Williams Lake", 0, share_wl),
    share_fm = ifelse(name_cancensus %in% c("Campbellton", "Hawkesbury"), 0, share_fm),
    Z_hat_wl = replace_na(share_wl, 0) * Post2017,
    Z_hat_fm = replace_na(share_fm, 0) * Post2016,
    Z_hat    = Z_hat_wl + Z_hat_fm
  )


# =============================================================================
# PARTIE VI -- CONSTRUCTION DES PANELS
# =============================================================================

# -----------------------------------------------------------------------------
# VI.1 Fonctions de merge et cylindrage
# -----------------------------------------------------------------------------
# merge_loyers() joint un panel d'exposition avec un dataset de loyers.
# cylindrer() exclut les villes avec ≥ 1 NA dans log_rent (sauf Leamington,
# qui a des NA en 2023-2024 mais une exposition non-nulle utile à l'identification).

name_bridge <- in_both_clean %>% select(name_cancensus, GEO_cmhc)

merge_loyers <- function(panel_instrument, data2_type) {
  rent_biprov <- data2_type %>%
    filter(grepl("Campbellton|Hawkesbury", GEO)) %>%
    mutate(name_cancensus = ifelse(grepl("Campbellton", GEO),
                                   "Campbellton", "Hawkesbury")) %>%
    group_by(name_cancensus, REF_DATE) %>%
    summarise(rent = mean(VALUE, na.rm = TRUE), .groups = "drop") %>%
    mutate(rent = ifelse(is.nan(rent), NA, rent))

  panel_instrument %>%
    left_join(name_bridge, by = "name_cancensus") %>%
    left_join(data2_type %>% select(REF_DATE, GEO, VALUE) %>%
                rename(GEO_cmhc = GEO, rent = VALUE),
              by = c("GEO_cmhc", "REF_DATE")) %>%
    left_join(rent_biprov, by = c("name_cancensus", "REF_DATE")) %>%
    mutate(rent = ifelse(is.na(rent.x), rent.y, rent.x)) %>%
    select(-rent.x, -rent.y) %>%
    filter(!name_cancensus %in% c("Wasaga Beach", "Whitehorse")) %>%
    filter(REF_DATE >= 2008, REF_DATE <= 2024) %>%
    mutate(rent     = ifelse(rent == 0, NA, rent),
           log_rent = log(rent))
}

cylindrer <- function(panel) {
  villes_na <- panel %>%
    group_by(name_cancensus) %>%
    summarise(n_NA   = sum(is.na(log_rent)),
              mean_Z = mean(Z_hat_wl, na.rm = TRUE)) %>%
    filter(n_NA > 0)

  panel %>%
    filter(!name_cancensus %in% villes_na$name_cancensus[
      villes_na$name_cancensus != "Leamington"])
}


# -----------------------------------------------------------------------------
# VI.2 Panels baseline et robustesse FM
# -----------------------------------------------------------------------------
panel_balanced_wl   <- cylindrer(merge_loyers(panel_instrument_wl,   data2_two))
panel_balanced_full <- cylindrer(merge_loyers(panel_instrument_full, data2_two))

cat("Panel baseline (WL, 2ch, apt3+) :",
    length(unique(panel_balanced_wl$name_cancensus)), "villes ×",
    length(unique(panel_balanced_wl$REF_DATE)), "années =",
    nrow(panel_balanced_wl), "observations\n")
# → 131 villes × 17 années = 2 227 lignes panel
# feols drope automatiquement les 2 NAs Leamington 2023-2024 lors de l'estimation
# → N estimation = 2 225, cohérent avec toutes les tables du manuscrit
# Note : on passe de 145 villes (échantillon analytique) à 131 après cylindrage
# Les 14 villes exclues ont ≥1 NA dans log_rent sur 2008-2024, toutes avec
# Z_hat_wl ≈ 0 → leur exclusion n'affecte pas l'identification
cat("NA dans log_rent :", sum(is.na(panel_balanced_wl$log_rent)),
    "(Leamington 2023-2024)\n")
# → 2 NA restants, tous concentrés sur Leamington 2023-2024 (post-traitement)
# → Confirme N = 2 227 - 2 = 2 225 dans toutes les régressions feols

# Variables province et bc (utilisées PARTIE VII et PARTIE IX)
panel_balanced_wl <- panel_balanced_wl %>%
  mutate(
    province = GEO_cmhc %>% str_extract("[^,]+$") %>% str_trim(),
    bc       = as.integer(province == "British Columbia")
  )

panel_balanced_full <- panel_balanced_full %>%
  mutate(
    year_rel_fm = REF_DATE - 2016,
    share_fm0   = replace_na(share_fm, 0)
  )

# Sous-échantillon BC uniquement
panel_bc <- panel_balanced_wl %>%
  filter(province == "British Columbia")


# -----------------------------------------------------------------------------
# VI.3 Panels d'hétérogénéité -- 4 structures × 4 types d'unité
# -----------------------------------------------------------------------------
# Toutes utilisent panel_instrument_wl (exposition Williams Lake 2017 uniquement).
# Nommage : panel_sX_uY où X = numéro structure, Y = type d'unité.

panel_s1_bachelor <- cylindrer(merge_loyers(panel_instrument_wl, data2_s1_bachelor))
panel_s1_one      <- cylindrer(merge_loyers(panel_instrument_wl, data2_s1_one))
panel_s1_two      <- panel_balanced_wl   # déjà créé = baseline
panel_s1_three    <- cylindrer(merge_loyers(panel_instrument_wl, data2_s1_three))

panel_s2_bachelor <- cylindrer(merge_loyers(panel_instrument_wl, data2_s2_bachelor))
panel_s2_one      <- cylindrer(merge_loyers(panel_instrument_wl, data2_s2_one))
panel_s2_two      <- cylindrer(merge_loyers(panel_instrument_wl, data2_s2_two))
panel_s2_three    <- cylindrer(merge_loyers(panel_instrument_wl, data2_s2_three))

panel_s3_bachelor <- cylindrer(merge_loyers(panel_instrument_wl, data2_s3_bachelor))
panel_s3_one      <- cylindrer(merge_loyers(panel_instrument_wl, data2_s3_one))
panel_s3_two      <- cylindrer(merge_loyers(panel_instrument_wl, data2_s3_two))
panel_s3_three    <- cylindrer(merge_loyers(panel_instrument_wl, data2_s3_three))

panel_s4_bachelor <- cylindrer(merge_loyers(panel_instrument_wl, data2_s4_bachelor))
panel_s4_one      <- cylindrer(merge_loyers(panel_instrument_wl, data2_s4_one))
panel_s4_two      <- cylindrer(merge_loyers(panel_instrument_wl, data2_s4_two))
panel_s4_three    <- cylindrer(merge_loyers(panel_instrument_wl, data2_s4_three))


# =============================================================================
# PARTIE VII -- VÉRIFICATIONS D'IDENTIFICATION
# =============================================================================
# Les vérifications d'identification précèdent intentionnellement les résultats
# principaux (PARTIE VIII) : on valide la stratégie avant de présenter les
# estimations. C'est la pratique standard depuis Angrist & Pischke (2009).

# Variable d'année relative au feu (nécessaire pour les event studies)
panel_balanced_wl <- panel_balanced_wl %>%
  mutate(
    year_rel_wl = REF_DATE - 2017,
    share_wl0   = replace_na(share_wl, 0)
  )


# -----------------------------------------------------------------------------
# VII.1 Event study correct -- test des pre-trends (es_wl_correct)
# -----------------------------------------------------------------------------
# POURQUOI : L'event study correct interagit les années relatives avec share_wl
# (part pré-feu FIXE), pas avec Z_hat_wl. Cette distinction est essentielle :
# Z_hat_wl est mécaniquement nul avant 2017 (par construction), donc interagir
# Z_hat_wl avec des indicatrices d'année produirait des zéros pré-traitement
# qui NE sont PAS informatifs sur les parallel trends.
# En revanche, share_wl × year_indicator est non-trivial avant 2017 et permet
# de tester si les villes plus exposées avaient des trajectoires de loyer
# différentes AVANT le feu.

es_wl_correct <- feols(
  log_rent ~ i(year_rel_wl, share_wl0, ref = -1) |
    name_cancensus + REF_DATE,
  data    = panel_balanced_wl,
  cluster = ~name_cancensus
)

iplot(es_wl_correct,
      main = "Pre-trends Test : Williams Lake Exposure × Year (95% CI shown)",
      xlab = "Years relative to 2017 wildfire",
      ylab = "Coefficient on pre-fire migration share (share_wl)",
      ci.width = 0.25)   # largeur des crochets IC


# -----------------------------------------------------------------------------
# VII.2 Dynamique post-traitement (es_wl_post_dynamics)
# -----------------------------------------------------------------------------
# NOTE : les coefficients avant k=0 sont mécaniquement nuls car Z_hat_wl = 0
# avant 2017. Ce graphique illustre la dynamique post-traitement seulement.
# Ne pas interpréter les zéros pré-traitement comme une validation
# des parallel trends -- utiliser es_wl_correct pour ce test.

es_wl_post_dynamics <- feols(
  log_rent ~ i(year_rel_wl, Z_hat_wl, ref = -1) |
    name_cancensus + REF_DATE,
  data    = panel_balanced_wl,
  cluster = ~name_cancensus
)

iplot(es_wl_post_dynamics,
      main = "Post-Treatment Dynamics : BC 2017 - Williams Lake (95% CI shown)",
      xlab = "Years relative to wildfire (0 = 2017)",
      ylab = "Coefficient on Z_WL (Williams Lake exposure index)",
      ci.width = 0.25)
# Effet positif et croissant à partir de k=0, persistant jusqu'en 2024 ✓


# -----------------------------------------------------------------------------
# VII.3 Event study Fort McMurray (robustesse de l'exposition FM)
# -----------------------------------------------------------------------------
es_fm <- feols(log_rent ~ i(year_rel_fm, Z_hat_fm, ref = -1) |
                 name_cancensus + REF_DATE,
               data = panel_balanced_full, cluster = ~name_cancensus)

iplot(es_fm,
      main = "Event Study - Fort McMurray 2016 (95% CI shown)",
      xlab = "Years relative to wildfire (0 = 2016)",
      ylab = "Coefficient on Z_FM (Fort McMurray exposure index)",
      ci.width = 0.25)
# Aucun effet post-traitement détectable ✓ -- déplacement temporaire (~4 mois)

range(panel_balanced_full$year_rel_fm)
# doit donner -8 et +8 (2008 à 2024)


# -----------------------------------------------------------------------------
# VII.3b Pre-trends test Fort McMurray (es_fm_correct)
# -----------------------------------------------------------------------------
# POURQUOI : identique à VII.1 pour Williams Lake. On interagit share_fm0
# (part gravitaire pré-feu FIXE, non-nulle avant 2016) avec les indicatrices
# d'année. Cela permet de tester si les villes plus exposées à Fort McMurray
# avaient des trajectoires de loyer différentes AVANT le feu.
# Contrairement à es_fm (qui utilise Z_hat_fm = share_fm × Post2016,
# mécaniquement nul avant 2016), es_fm_correct génère une vraie variation
# pré-traitement et constitue un test valide des parallel trends.

es_fm_correct <- feols(
  log_rent ~ i(year_rel_fm, share_fm0, ref = -1) |
    name_cancensus + REF_DATE,
  data    = panel_balanced_full,
  cluster = ~name_cancensus
)

iplot(es_fm_correct,
      main = "Pre-trends Test : Fort McMurray Exposure × Year (95% CI shown)",
      xlab = "Years relative to 2016 wildfire",
      ylab = "Coefficient on gravity-predicted pre-fire share (share_fm)",
      ci.width = 0.25)
# Si les parallel trends tiennent : coefficients proches de 0 et non significatifs
# avant k = 0, puis éventuellement un effet post-traitement (ou non, ici attendu nul)


# -----------------------------------------------------------------------------
# VII.4 Placebo tests -- validité de l'identification
# -----------------------------------------------------------------------------
# POURQUOI : Si l'indice d'exposition ne capte pas de tendances de loyer
# préexistantes, il ne doit pas "prédire" des hausses de loyer sur des
# pseudo-années de traitement pré-2017. Des placebos significatifs
# indiquent une violation des parallel trends.
# Référence : Angrist & Pischke (Mostly Harmless Econometrics, ch.5)

placebo_years <- 2010:2016

# --- Placebo 1 : spécification baseline (FE ville + année) ---
placebo_baseline <- map_dfr(placebo_years, function(yr) {
  panel_balanced_wl %>%
    mutate(Z_placebo = replace_na(share_wl, 0) * as.integer(REF_DATE >= yr)) %>%
    feols(log_rent ~ Z_placebo | name_cancensus + REF_DATE,
          data = ., cluster = ~name_cancensus) %>%
    broom::tidy() %>%
    filter(term == "Z_placebo") %>%
    mutate(placebo_year = yr, spec = "Baseline (City+Year FE)")
})
# Résultats placebo baseline :
# 2010 : β = 0.217 (p=0.057) | 2011 : β = 0.304** | 2012 : β = 0.386***
# 2013 : β = 0.445*** | 2014 : β = 0.551*** | 2015 : β = 0.669*** | 2016 : β = 0.760***
# → PROBLÈME : tous significatifs et croissants vers 2017
# → Les villes exposées (BC) avaient des loyers qui augmentaient plus vite
#   que les non-exposées AVANT le feu → violation des parallel trends
# → Cause : confounding provincial -- les villes exposées sont quasi-toutes en BC
#   qui a connu une hausse des loyers supérieure à la moyenne canadienne 2010-2017

# --- Placebo 2 : avec contrôle de la tendance BC linéaire ---
placebo_bc_trend <- map_dfr(placebo_years, function(yr) {
  panel_balanced_wl %>%
    mutate(Z_placebo = replace_na(share_wl, 0) * as.integer(REF_DATE >= yr)) %>%
    feols(log_rent ~ Z_placebo + bc:REF_DATE | name_cancensus + REF_DATE,
          data = ., cluster = ~name_cancensus) %>%
    broom::tidy() %>%
    filter(term == "Z_placebo") %>%
    mutate(placebo_year = yr, spec = "BC linear trend")
})
# Résultats placebo avec tendance BC :
# 2010 : β = -0.422 (p=0.065) | 2011-2016 : tous non significatifs (p > 0.10)
# → Les placebo deviennent non significatifs une fois la tendance BC contrôlée ✓
# → MAIS le coefficient principal rf_bc_trend devient également non significatif
#   (β = 0.261, p=0.22) → la tendance BC et l'exposition sont trop colinéaires

# --- Placebo 3 : avec tendances linéaires ville (spécification la plus exigeante) ---
placebo_trends <- map_dfr(placebo_years, function(yr) {
  panel_balanced_wl %>%
    mutate(Z_placebo = replace_na(share_wl, 0) * as.integer(REF_DATE >= yr)) %>%
    feols(log_rent ~ Z_placebo | name_cancensus + REF_DATE +
            name_cancensus[[REF_DATE]],
          data = ., cluster = ~name_cancensus) %>%
    broom::tidy() %>%
    filter(term == "Z_placebo") %>%
    mutate(placebo_year = yr, spec = "City linear trends")
})
# Résultats placebo avec tendances ville :
# 2010-2014 : β négatifs et significatifs (-0.38 à -0.66)
# 2015 : β ≈ 0 (n.s.)
# 2016 : β = +0.365** -- significatif et positif
# → Pattern en forme de U inversé -- les tendances linéaires ville sur-ajustent
#   et créent des artefacts. Spécification trop agressive.

# --- Placebo 4 : échantillon BC uniquement (24 villes) ---
# panel_bc défini en PARTIE VI
placebo_bc_only <- map_dfr(placebo_years, function(yr) {
  panel_bc %>%
    mutate(Z_placebo = replace_na(share_wl, 0) * as.integer(REF_DATE >= yr)) %>%
    feols(log_rent ~ Z_placebo | name_cancensus + REF_DATE,
          data = ., cluster = ~name_cancensus) %>%
    broom::tidy() %>%
    filter(term == "Z_placebo") %>%
    mutate(placebo_year = yr, spec = "BC only")
})
# Résultats placebo BC uniquement :
# Tous non significatifs (p > 0.55 partout) ✓
# → Au sein de la BC, les villes exposées et non-exposées avaient des trajectoires
#   similaires avant 2017 -- le confounding est bien provincial, pas intra-provincial

placebo_all <- bind_rows(placebo_baseline, placebo_bc_trend,
                         placebo_trends,   placebo_bc_only)

placebo_summary <- placebo_all %>%
  select(placebo_year, estimate, std.error, p.value, spec) %>%
  mutate(
    stars = case_when(
      p.value < 0.01 ~ "***",
      p.value < 0.05 ~ "**",
      p.value < 0.10 ~ "*",
      TRUE           ~ ""
    ),
    cell = sprintf("%.3f%s\n(%.3f)", estimate, stars, std.error)
  )


# =============================================================================
# PARTIE VIII -- RÉSULTATS PRINCIPAUX
# =============================================================================
# POURQUOI CET ORDRE : les vérifications d'identification (PARTIE VII) ont
# établi les conditions dans lesquelles les résultats suivants sont
# interprétables. L'ordre reflète la logique de Angrist & Pischke : montrer
# d'abord que l'identification est consistente avec les données, puis présenter
# les estimations.

# Fonction d'estimation (forme réduite standard, réutilisée dans PARTIE X)
estimer_rf <- function(panel) {
  feols(log_rent ~ Z_hat_wl | name_cancensus + REF_DATE,
        data = panel, cluster = ~name_cancensus)
}

# Fonction d'estimation multi-spécifications (progression des contrôles)
# Retourne une liste nommée de 5 spécifications pour Table 1
estimate_main_specs <- function(panel) {
  list(
    ols      = feols(log_rent ~ Z_hat_wl,
                     data = panel, cluster = ~name_cancensus),
    fe_year  = feols(log_rent ~ Z_hat_wl | REF_DATE,
                     data = panel, cluster = ~name_cancensus),
    fe_city  = feols(log_rent ~ Z_hat_wl | name_cancensus,
                     data = panel, cluster = ~name_cancensus),
    baseline = feols(log_rent ~ Z_hat_wl | name_cancensus + REF_DATE,
                     data = panel, cluster = ~name_cancensus),
    trends   = feols(log_rent ~ Z_hat_wl | name_cancensus + REF_DATE +
                       name_cancensus[[REF_DATE]],
                     data = panel, cluster = ~name_cancensus)
  )
}

main_specs <- estimate_main_specs(panel_balanced_wl)

rf_ols      <- main_specs$ols
rf_fe_year  <- main_specs$fe_year
rf_fe_city  <- main_specs$fe_city
rf_baseline <- main_specs$baseline
rf_trends   <- main_specs$trends

# --- Résultats obtenus ---
# Col  Spécification            β        Sig.
# (1)  OLS pur                  3.838    ***   ← biais important sans contrôles
# (2)  + FE année               2.488    **    ← tendances nationales absorbées
# (3)  + FE ville               3.397    ***   ← niveaux fixes absorbés
# (4)  + FE ville + année       0.818    ***   ← spécification principale ✓
# (5)  + tendances ville        0.594    ***   ← effet non drivé par tendances ✓
#
# INTERPRÉTATION (reduced-form uniquement) :
# Consistent with a wildfire-induced housing demand shock, cities with stronger
# pre-fire migration links to Williams Lake experienced higher rent growth
# after 2017. The estimate is consistent with, but does not establish,
# a causal effect of displacement on rents given the pre-trend concerns
# documented in PARTIE VII.
# La chute de β de 3.84 (OLS) à 0.82 (FE ville+année) confirme que les
# spécifications naïves capturaient un biais d'endogénéité.
# La stabilité entre (4) et (5) est le résultat clé pour la thèse.

# ─── Vérification des R² et Within R² rapportés dans Table 1 (tab:main_result) ─
# Le manuscrit hard-code R² = (0.055, 0.322, 0.662, 0.940, 0.975) et
# Within R² = (---, 0.032, 0.065, 0.021, 0.007). On les imprime ici pour
# pouvoir auditer ces valeurs à chaque exécution.
cat("=== Vérification R^2 / Within R^2 -- Table 1 (tab:main_result) ===\n")
cat(sprintf("(1) OLS         : R2 = %.3f\n",
            fitstat(rf_ols, "r2")$r2))
cat(sprintf("(2) Year FE     : R2 = %.3f  Within R2 = %.3f\n",
            fitstat(rf_fe_year,  "r2")$r2,  fitstat(rf_fe_year,  "wr2")$wr2))
cat(sprintf("(3) City FE     : R2 = %.3f  Within R2 = %.3f\n",
            fitstat(rf_fe_city,  "r2")$r2,  fitstat(rf_fe_city,  "wr2")$wr2))
cat(sprintf("(4) Baseline    : R2 = %.3f  Within R2 = %.3f\n",
            fitstat(rf_baseline, "r2")$r2,  fitstat(rf_baseline, "wr2")$wr2))
cat(sprintf("(5) City trends : R2 = %.3f  Within R2 = %.3f\n",
            fitstat(rf_trends,   "r2")$r2,  fitstat(rf_trends,   "wr2")$wr2))


# =============================================================================
# PARTIE IX -- VÉRIFICATIONS DE ROBUSTESSE
# =============================================================================

# -----------------------------------------------------------------------------
# IX.1 Expositions alternatives (FM)
# -----------------------------------------------------------------------------
rf_fm       <- feols(log_rent ~ Z_hat_fm | name_cancensus + REF_DATE,
                     data = panel_balanced_full, cluster = ~name_cancensus)
rf_combined <- feols(log_rent ~ Z_hat_wl + Z_hat_fm | name_cancensus + REF_DATE,
                     data = panel_balanced_full, cluster = ~name_cancensus)
rf_agg      <- feols(log_rent ~ Z_hat | name_cancensus + REF_DATE,
                     data = panel_balanced_full, cluster = ~name_cancensus)
rf_unbal    <- feols(log_rent ~ Z_hat_wl | name_cancensus + REF_DATE,
                     data = merge_loyers(panel_instrument_wl, data2_two),
                     cluster = ~name_cancensus)
# Résultats Table 2 -- Robustesse des expositions :
# (1) Baseline Z_hat_wl          : β = +0.818*** -- résultat principal
# (2) Z_hat_fm seul               : β = +0.055  (n.s.) -- FM non significatif
#     → retour rapide des évacués Fort McMurray (~4 mois), marchés albertains élastiques
# (3) Z_hat_wl + Z_hat_fm         : β_wl = +0.835*** / β_fm = -0.670 (n.s.)
#     → WL reste robuste, FM toujours non significatif
# (4) Z_hat agrégé (WL + FM)      : β = +0.637*** -- dilution par FM non significatif
# (5) Panel non-cylindré           : β_wl = +0.828*** -- cylindrage sans conséquence


# -----------------------------------------------------------------------------
# IX.2 Tendances spécifiques et échantillon BC
# -----------------------------------------------------------------------------
# POURQUOI : La principale menace à l'identification est que les villes exposées
# sont géographiquement concentrées en BC, qui connaissait une tendance haussière
# des loyers indépendamment du feu. Ces spécifications adressent directement
# cette préoccupation.

# (A) Tendance linéaire spécifique à la Colombie-Britannique
# Tests whether the main result is simply driven by BC's particular rent dynamics
rf_bc_trend <- feols(
  log_rent ~ Z_hat_wl | name_cancensus + REF_DATE + bc[[REF_DATE]],
  data    = panel_balanced_wl,
  cluster = ~name_cancensus
)

# (B) Tendances linéaires province-spécifiques (spécification plus générale)
# Allows each province to have its own linear rent trajectory
rf_province_trends <- feols(
  log_rent ~ Z_hat_wl | name_cancensus + REF_DATE + province[[REF_DATE]],
  data    = panel_balanced_wl,
  cluster = ~name_cancensus
)

# (C) Échantillon Colombie-Britannique uniquement
# Compares exposed and less-exposed cities within BC only, reducing the
# risk of confounding from comparing BC to other provinces.
# panel_bc défini en PARTIE VI
rf_bc_only <- feols(
  log_rent ~ Z_hat_wl | name_cancensus + REF_DATE,
  data    = panel_bc,
  cluster = ~name_cancensus
)
# Note: smaller sample (~24 cities) → less precise, but informative for
# identification validity.

# (E) Province×Year fixed effects -- TEST CENTRAL D'IDENTIFICATION
# POURQUOI : Ce test compare uniquement des villes DANS LA MÊME PROVINCE ET LA MÊME
# ANNÉE, éliminant directement le confondant "dynamisme BC" qui affecte la
# spécification baseline. C'est le test le plus exigeant disponible avec ces données.
# INTERPRÉTATION : Si significatif → effets WL existent au sein des provinces.
# Si non significatif → le résultat baseline est entièrement drivé par des
# différences interprovinciates, ce qui est la critique centrale de l'Avocat du Diable.
# Note technique (fixest) : province^REF_DATE crée une FE pour chaque cellule
# province×année, absorbant toute tendance commune à toutes les villes d'une province
# la même année. Les cellules singletons (1 seule ville dans province-année) sont
# automatiquement droppées par feols.
rf_prov_year <- feols(
  log_rent ~ Z_hat_wl | name_cancensus + province^REF_DATE,
  data    = panel_balanced_wl,
  cluster = ~name_cancensus
)
cat(sprintf("Province×Year FE : beta = %.3f (se = %.3f) N = %d (singletons droppés)\n",
            coef(rf_prov_year)["Z_hat_wl"],
            se(rf_prov_year)["Z_hat_wl"],
            nobs(rf_prov_year)))
# Résultat attendu : β ≈ 0.164 (s.e. ≈ 0.260, p > 0.10), N ≈ 2191
# DIAGNOSTIC CENTRAL : résultat non significatif → confondant BC provincial non absorbé
# dans la spécification baseline.

# (F) Region×Year fixed effects (West/East)
# Moins saturé que province×année : compare villes au sein de la même macro-région
# (West = BC/AB/SK/MB, East = ON/QC/Maritimes) et de la même année.
# Moins de cellules vides que province×année → plus de puissance.
panel_balanced_wl <- panel_balanced_wl %>%
  mutate(
    region     = ifelse(
      province %in% c("British Columbia", "Alberta", "Saskatchewan", "Manitoba"),
      "West", "East"
    ),
    Z_hat_wl18 = share_wl0 * as.integer(REF_DATE >= 2018)
  )

rf_region_year <- feols(
  log_rent ~ Z_hat_wl | name_cancensus + region^REF_DATE,
  data    = panel_balanced_wl,
  cluster = ~name_cancensus
)
cat(sprintf("Region×Year FE   : beta = %.3f (se = %.3f) N = %d\n",
            coef(rf_region_year)["Z_hat_wl"],
            se(rf_region_year)["Z_hat_wl"],
            nobs(rf_region_year)))
# Résultat attendu : β ≈ 0.957*** (s.e. ≈ 0.308) -- West/East ne neutralise pas
# entièrement le confondant BC, mais preserve une identification géographiquement plus restreinte.

# (D) Timing t>=2018 : redéfinition de l'indicateur post-traitement (Section 7.5)
# Le feu de Williams Lake survient en juillet 2017 mais les loyers SCHL sont
# mesurés en octobre. L'indicateur 1[t>=2017] inclut donc une observation
# partiellement traitée. On attend une légère hausse du coefficient (le k=0
# de l'event study est faible et non significatif) -- ce que confirme β=0.873***.
rf_t2018 <- feols(
  log_rent ~ Z_hat_wl18 | name_cancensus + REF_DATE,
  data    = panel_balanced_wl,
  cluster = ~name_cancensus
)

# =============================================================================
# DIAGNOSTIC FINAL ET IMPLICATIONS POUR LA THÈSE
# =============================================================================
#
# RÉSUMÉ DU PROBLÈME D'IDENTIFICATION :
#
# Spécification                    β_wl    Sig.   Placebo
# Baseline (City+Year FE)          0.818   ***    problématique (tendance BC)
# + Tendance BC linéaire           0.261   n.s.   non significatifs
# + Tendances ville                0.594   ***    pattern aberrant (sur-ajustement)
# BC uniquement (24 villes)        0.144   n.s.   non significatifs
#
# DIAGNOSTIC :
# Le résultat principal (β = 0.818***) capture deux effets superposés :
# (1) L'effet reduced-form du choc BC 2017 sur les loyers des villes d'accueil
# (2) La tendance haussière spécifique aux marchés locatifs de BC 2010-2017
# Ces deux effets sont difficiles à séparer car les villes exposées sont
# quasi-exclusivement en BC.
#
# FORMULATION POUR LA THÈSE :
# "Our main estimate should be interpreted with caution as exposed cities are
#  predominantly located in British Columbia, which experienced above-average
#  rent growth over the 2010-2017 period independently of the wildfire.
#  While city-specific time trends partially absorb this confounding
#  (Column 5, β = 0.594***), placebo tests suggest residual bias cannot
#  be fully ruled out. Our estimates therefore represent reduced-form evidence
#  consistent with a wildfire-induced housing demand shock, and may constitute
#  an upper bound of the causal effect of climate-induced displacement on
#  rental prices."
#
# Cette formulation est académiquement honnête et montre la maturité
# méthodologique que les jurys de M2 valorisent.


# =============================================================================
# PARTIE IX.3  Hard Checks -- Canal politique Vancouver + Trimming outliers
# =============================================================================
# POURQUOI :
# Vancouver (s_WL = 0.174) est la 2e ville la plus exposée ET a subi deux
# chocs de politique immobilière quasi-simultanés indépendants du feu :
#   - BC Foreign Buyers Tax (août 2016) : taxation des acheteurs étrangers
#     → report de la demande de l'achat vers la location
#   - Speculation & Vacancy Tax (annoncée 2018)
# Si Vancouver drive π via ce canal politique et non via la migration,
# l'interprétation Williams Lake-spécifique serait affaiblie. Ces trois tests ciblent ce canal alternatif.

# Hard Check 1 : Exclure Vancouver du panel
panel_drop_van <- panel_balanced_wl %>%
  filter(name_cancensus != "Vancouver")

rf_drop_van <- feols(
  log_rent ~ Z_hat_wl | name_cancensus + REF_DATE,
  data    = panel_drop_van,
  cluster = ~name_cancensus
)

# Hard Check 2 : Contrôler par Vancouver × Post2017
# Absorbe tout choc post-2017 spécifique à Vancouver (politique ou autre)
# sans retirer la ville. Si π_ZWL reste stable, le résultat n'est pas drivé
# par Vancouver via un canal non-migratoire.
panel_balanced_wl <- panel_balanced_wl %>%
  mutate(
    d_van    = as.integer(name_cancensus == "Vancouver"),
    van_post = d_van * Post2017
  )

rf_van_post <- feols(
  log_rent ~ Z_hat_wl + van_post | name_cancensus + REF_DATE,
  data    = panel_balanced_wl,
  cluster = ~name_cancensus
)

# Soft Check (trimming) : Exclure Kamloops ET Vancouver
# Ces deux villes (s_WL = 0.212 et 0.174) représentent ~40% du poids
# d'identification. Si π reste significatif sur les 41 villes restantes
# (s_WL ≤ 0.10), l'effet n'est pas mécanique sur 2 points de levier.
panel_drop_kam_van <- panel_balanced_wl %>%
  filter(!name_cancensus %in% c("Vancouver", "Kamloops"))

rf_drop_kam_van <- feols(
  log_rent ~ Z_hat_wl | name_cancensus + REF_DATE,
  data    = panel_drop_kam_van,
  cluster = ~name_cancensus
)

cat("--- Hard Checks IX.3 ---\n")
cat(sprintf("Hard Check 1 (drop Vancouver)    : beta = %.3f (se = %.3f) | N = %d\n",
            coef(rf_drop_van)["Z_hat_wl"],
            se(rf_drop_van)["Z_hat_wl"],
            nobs(rf_drop_van)))
cat(sprintf("Hard Check 2 (Vancouver x Post)  : beta_Z = %.3f (se = %.3f) | beta_VP = %.3f\n",
            coef(rf_van_post)["Z_hat_wl"],
            se(rf_van_post)["Z_hat_wl"],
            coef(rf_van_post)["van_post"]))
cat(sprintf("Soft Check   (drop Kam+Van)       : beta = %.3f (se = %.3f) | N = %d\n",
            coef(rf_drop_kam_van)["Z_hat_wl"],
            se(rf_drop_kam_van)["Z_hat_wl"],
            nobs(rf_drop_kam_van)))

# =============================================================================
# PARTIE IX.4bis  Échantillon tronqué 2019 (Robustesse COVID)
# =============================================================================
# POURQUOI : Les coefficients de l'event study continuent à croître jusqu'en 2023-2024.
# La pandémie COVID-19 (2020-2022) a provoqué des migrations interrégionales massives
# au Canada, notamment vers les villes secondaires de CB (Kamloops, Kelowna, Prince George).
# Ces villes sont précisément les plus exposées à Williams Lake. Si les coefficients
# post-2017 sont concentrés dans la période 2020-2024, cela suggère un confondant COVID
# plutôt qu'un effet persistant du feu de 2017.
# TEST : En tronquant l'échantillon à 2019 (avant COVID), on isole la fenêtre 2017-2019
# où seul l'effet WL peut être actif.

panel_pre_covid <- panel_balanced_wl %>%
  filter(REF_DATE <= 2019)

rf_pre_covid <- feols(
  log_rent ~ Z_hat_wl | name_cancensus + REF_DATE,
  data    = panel_pre_covid,
  cluster = ~name_cancensus
)
cat(sprintf("Pre-COVID (2008-2019): beta = %.3f (se = %.3f) N_obs = %d, N_cities = %d\n",
            coef(rf_pre_covid)["Z_hat_wl"],
            se(rf_pre_covid)["Z_hat_wl"],
            nobs(rf_pre_covid),
            length(unique(panel_pre_covid$name_cancensus))))
# LECTURE :
# Si β_pre_covid > 0 et significatif → l'effet existe avant COVID → plus crédible
# Si β_pre_covid ≈ 0 → l'effet est concentré en 2020-2024 → confondant COVID probable


# =============================================================================
# PARTIE IX.4ter  Tendances estimées uniquement sur la période pré-traitement
# =============================================================================
# POURQUOI : La spécification avec tendances linéaires par ville (eq. 3) estime
# les pentes sur 2008-2024, y compris les années post-traitement. Si le traitement
# a un effet persistant, la tendance estimée pour les villes exposées est
# "contaminée" par l'effet traitement lui-même, biaisant π̂ vers le bas.
# SOLUTION ANGRIST-PISCHKE : estimer les tendances uniquement sur 2008-2016
# (pré-traitement), puis détrendre le panel entier avec ces pentes exogènes.
# Cette spécification impose une restriction plus faible sur le contrefactuel.

pre_period_data <- panel_balanced_wl %>% filter(REF_DATE < 2017)

trend_fits_pre <- pre_period_data %>%
  group_by(name_cancensus) %>%
  summarise(
    trend_slope = tryCatch(
      coef(lm(log_rent ~ REF_DATE, data = pick(everything())))["REF_DATE"],
      error = function(e) NA_real_
    ),
    .groups = "drop"
  )

panel_detrended <- panel_balanced_wl %>%
  left_join(trend_fits_pre, by = "name_cancensus") %>%
  mutate(
    log_rent_dt = log_rent - trend_slope * (REF_DATE - 2016)
  )

rf_pretrend_only <- feols(
  log_rent_dt ~ Z_hat_wl | name_cancensus + REF_DATE,
  data    = panel_detrended %>% filter(!is.na(log_rent_dt)),
  cluster = ~name_cancensus
)
cat(sprintf("Tendances pré-2017 seulement : beta = %.3f (se = %.3f) N = %d\n",
            coef(rf_pretrend_only)["Z_hat_wl"],
            se(rf_pretrend_only)["Z_hat_wl"],
            nobs(rf_pretrend_only)))
# LECTURE :
# Si similaire à rf_trends (0.594) → la contamination des tendances par le traitement est mineure.
# Si sensiblement différent → la spécification rf_trends absorbe de l'effet traitement.


# =============================================================================
# PARTIE IX.5  Placebo Quesnel -- Soft Check (données déjà disponibles)
# =============================================================================
# POURQUOI : Quesnel est une ville BC de taille comparable à Williams Lake
# (~10 000 hab.), forestière, géographiquement proche, qui n'a PAS brûlé en 2017.
# Si Z_quesnel = share_quesnel × Post2017 prédit les loyers comme Z_WL,
# alors s_WL proxie la "connectivité BC générale" et non le choc spécifique.
# Si Z_quesnel est nul → l'identification est bien spécifique au feu WL.

dataset_quesnel_clean <- data1 %>%
  filter(GEO == "Quesnel (CA), British Columbia",
         REF_DATE == "2016/2017",
         !grepl("Area outside", `Geography.of.destination`)) %>%
  mutate(
    total_out    = sum(VALUE, na.rm = TRUE),
    share_quesnel = VALUE / total_out
  )

shares_quesnel <- dataset_quesnel_clean %>%
  mutate(name_cancensus = `Geography.of.destination` %>%
           str_remove("\\s*\\(C[MA]+\\)") %>% str_remove(",.*") %>% str_trim()) %>%
  select(name_cancensus, share_quesnel) %>%
  distinct(name_cancensus, .keep_all = TRUE)

panel_quesnel <- panel_balanced_wl %>%
  left_join(shares_quesnel, by = "name_cancensus") %>%
  mutate(
    share_quesnel = replace_na(share_quesnel, 0),
    Z_quesnel     = share_quesnel * Post2017
  )

rf_placebo_quesnel <- feols(
  log_rent ~ Z_quesnel | name_cancensus + REF_DATE,
  data    = panel_quesnel,
  cluster = ~name_cancensus
)

cat(sprintf("Placebo Quesnel : beta = %.3f (se = %.3f) p = %.3f\n",
            coef(rf_placebo_quesnel)["Z_quesnel"],
            se(rf_placebo_quesnel)["Z_quesnel"],
            coeftable(rf_placebo_quesnel)["Z_quesnel", "Pr(>|t|)"]))
# PROBLÈME RÉSIDUEL : Quesnel est géographiquement et économiquement proche de WL.
# Un résultat non-nul pour Quesnel peut refléter une connectivité BC générale
# (les deux villes partagent les mêmes destinations), pas l'absence de spécificité WL.
# → Solution : placebo avec une ville GÉOGRAPHIQUEMENT DISTANTE (voir IX.5bis ci-dessous).


# =============================================================================
# PARTIE IX.5bis  Placebo distant -- Greater Sudbury (Ontario)
# =============================================================================
# POURQUOI : Quesnel est trop proche de WL pour constituer un placebo propre.
# Il faut une ville de profil économique similaire (forestière/minière, ~10-30k hab.)
# mais GÉOGRAPHIQUEMENT DISTANTE de la BC, de sorte que sa connectivité migratoire
# ne soit pas corrélée avec celle de Williams Lake.
#
# CHOIX : Greater Sudbury (CMA), Ontario
#   - Ville minière (nickel), ~166 000 habitants (beaucoup plus grande que WL)
#   - À ~5 000 km de Williams Lake → réseau migratoire orthogonal à WL
#   - N'a pas subi de catastrophe naturelle majeure en 2017
#   - Présente dans la table 17-10-0141-01 → shares disponibles
#
# Si Z_sudbury = share_sudbury × Post2017 prédit les loyers →
#   l'effet BC 2017 capture une régularité canadienne générale (toutes villes connectées
#   à un pôle économique distinct connaissent le même schéma → non spécifique WL)
# Si Z_sudbury est nul →
#   l'identification est géographiquement spécifique (réseau BC)

dataset_sudbury_clean <- data1 %>%
  filter(GEO == "Greater Sudbury (CMA), Ontario",
         REF_DATE == "2016/2017",
         !grepl("Area outside", `Geography.of.destination`)) %>%
  mutate(
    total_out      = sum(VALUE, na.rm = TRUE),
    share_sudbury  = VALUE / total_out
  )

shares_sudbury <- dataset_sudbury_clean %>%
  mutate(name_cancensus = `Geography.of.destination` %>%
           str_remove("\\s*\\(C[MA]+\\)") %>% str_remove(",.*") %>% str_trim()) %>%
  select(name_cancensus, share_sudbury) %>%
  distinct(name_cancensus, .keep_all = TRUE)

panel_sudbury <- panel_balanced_wl %>%
  left_join(shares_sudbury, by = "name_cancensus") %>%
  mutate(
    share_sudbury = replace_na(share_sudbury, 0),
    Z_sudbury     = share_sudbury * Post2017
  )

rf_placebo_sudbury <- feols(
  log_rent ~ Z_sudbury | name_cancensus + REF_DATE,
  data    = panel_sudbury,
  cluster = ~name_cancensus
)
cat(sprintf("Placebo Sudbury (distant) : beta = %.3f (se = %.3f) p = %.3f\n",
            coef(rf_placebo_sudbury)["Z_sudbury"],
            se(rf_placebo_sudbury)["Z_sudbury"],
            coeftable(rf_placebo_sudbury)["Z_sudbury", "Pr(>|t|)"]))
# LECTURE :
# Si p > 0.10 → le résultat WL est géographiquement spécifique (BC réseau).
# Si p < 0.10 → une connectivité migratoire quelconque prédit les loyers →
#              l'effet n'est pas propre à WL (remet en cause l'interprétation).


# =============================================================================
# PARTIE IX.6  Gravity shares pour WL -- Soft Check (données déjà disponibles)
# =============================================================================
# POURQUOI : Les parts observées s_WL pourraient refléter des tendances
# migratoires pré-existantes corrélées aux loyers. En utilisant des parts
# prédites par gravité (distance + population uniquement), on s'assure que
# la variation est purement géographique et prédéterminée.
# Même logique que l'exposition Fort McMurray (cf. PARTIE V.2).

wl_lat <- 52.1418
wl_lon <- -122.1418

gravity_wl <- coords %>%
  filter(name_cancensus != "Williams Lake") %>%
  left_join(province_dict, by = "name_cancensus") %>%
  rowwise() %>%
  mutate(
    dist_km = distHaversine(c(wl_lon, wl_lat),
                            c(longitude, latitude)) / 1000,
    border  = as.integer(province != "British Columbia")
  ) %>%
  ungroup() %>%
  filter(dist_km > 0) %>%
  mutate(
    flux_pred   = exp(predict(gravity_model,
                              newdata = data.frame(dist_km  = dist_km,
                                                   border   = border,
                                                   pop_dest = Population))),
    total_pred  = sum(flux_pred, na.rm = TRUE),
    share_grav  = flux_pred / total_pred
  ) %>%
  select(name_cancensus, share_grav)

panel_grav <- panel_balanced_wl %>%
  left_join(gravity_wl, by = "name_cancensus") %>%
  mutate(
    share_grav = replace_na(share_grav, 0),
    Z_grav     = share_grav * Post2017
  )

rf_gravity_wl <- feols(
  log_rent ~ Z_grav | name_cancensus + REF_DATE,
  data    = panel_grav,
  cluster = ~name_cancensus
)

cat(sprintf("Gravity WL      : beta = %.3f (se = %.3f) p = %.3f\n",
            coef(rf_gravity_wl)["Z_grav"],
            se(rf_gravity_wl)["Z_grav"],
            coeftable(rf_gravity_wl)["Z_grav", "Pr(>|t|)"]))

# -----------------------------------------------------------------------------
# PPML gravity shares pour WL — robustesse Santos Silva & Tenreyro (2006)
# gravity_ppml est estimé en PARTIE V.2 (fepois sur flux FM + WL ensemble).
# Les prédictions PPML sont en niveaux (pas de exp() nécessaire).
# -----------------------------------------------------------------------------
gravity_wl_ppml <- coords %>%
  filter(name_cancensus != "Williams Lake") %>%
  left_join(province_dict, by = "name_cancensus") %>%
  rowwise() %>%
  mutate(
    dist_km = distHaversine(c(wl_lon, wl_lat),
                            c(longitude, latitude)) / 1000,
    border  = as.integer(province != "British Columbia")
  ) %>%
  ungroup() %>%
  filter(dist_km > 0) %>%
  mutate(
    flux_pred_ppml  = predict(gravity_ppml,
                              newdata = data.frame(dist_km  = dist_km,
                                                   border   = border,
                                                   pop_dest = Population)),
    total_pred_ppml = sum(flux_pred_ppml, na.rm = TRUE),
    share_grav_ppml = flux_pred_ppml / total_pred_ppml
  ) %>%
  select(name_cancensus, share_grav_ppml)

panel_grav_ppml <- panel_balanced_wl %>%
  left_join(gravity_wl_ppml, by = "name_cancensus") %>%
  mutate(
    share_grav_ppml = replace_na(share_grav_ppml, 0),
    Z_grav_ppml     = share_grav_ppml * Post2017
  )

rf_gravity_ppml_wl <- feols(
  log_rent ~ Z_grav_ppml | name_cancensus + REF_DATE,
  data    = panel_grav_ppml,
  cluster = ~name_cancensus
)

cat(sprintf("Gravity WL PPML : beta = %.3f (se = %.3f) p = %.3f\n",
            coef(rf_gravity_ppml_wl)["Z_grav_ppml"],
            se(rf_gravity_ppml_wl)["Z_grav_ppml"],
            coeftable(rf_gravity_ppml_wl)["Z_grav_ppml", "Pr(>|t|)"]))

# =============================================================================
# PARTIE IX.4  Code préparé -- données externes requises (non exécutable)
# =============================================================================

# ── M2-A : ForestShare × Post2017 (données NAICS 11+21 StatCan) ─────────────
# Une fois forest_employment_city.csv chargé (1 ligne par ville, colonne
# forest_share = emplois NAICS 11+21 / emplois totaux) :
#
# forest_data <- read.csv("raw/forest_share/forest_employment_city.csv") %>%
#   rename(name_cancensus = city)
#
# panel_forest <- panel_balanced_wl %>%
#   left_join(forest_data, by = "name_cancensus") %>%
#   mutate(forest_share = replace_na(forest_share, 0),
#          forest_post  = forest_share * Post2017)
#
# rf_forest <- feols(
#   log_rent ~ Z_hat_wl + forest_post | name_cancensus + REF_DATE,
#   data    = panel_forest,
#   cluster = ~name_cancensus
# )
# Lecture : si coef(rf_forest)["Z_hat_wl"] ≈ rf_baseline → canal migration domine.

# ── M2-B : Placebo WL' -- Quesnel (même table 17-10-0141, autre origine) ────
# Construire share_quesnel depuis la Table 17-10-0141 en filtrant sur
# Geography.of.origin == "Quesnel, British Columbia" pour la période 2015/2016.
#
# share_quesnel <- dataset_wl_raw %>%          # même table, filtre origine = Quesnel
#   filter(grepl("Quesnel", Geography.of.origin),
#          REF_DATE %in% c("2015/2016")) %>%
#   group_by(name_cancensus) %>%
#   summarise(share_placebo = VALUE / sum(VALUE, na.rm = TRUE), .groups = "drop")
#
# panel_placebo <- panel_balanced_wl %>%
#   left_join(share_quesnel, by = "name_cancensus") %>%
#   mutate(share_placebo = replace_na(share_placebo, 0),
#          Z_placebo     = share_placebo * Post2017)
#
# rf_placebo_quesnel <- feols(
#   log_rent ~ Z_placebo | name_cancensus + REF_DATE,
#   data    = panel_placebo,
#   cluster = ~name_cancensus
# )
# Lecture : si p > 0.10 → l'identification est spécifique à Williams Lake.
# Si p < 0.10 → s_WL proxie la connectivité BC générale.

# ── M2-C : Gravity shares pour WL (pas de données externes -- uses coords) ──
# Williams Lake coords (hardcodées, source : Google Maps)
# wl_lat <- 52.1418; wl_lon <- -122.1418
#
# gravity_wl <- coords %>%
#   filter(name_cancensus != "Williams Lake") %>%
#   mutate(
#     dist_km    = geosphere::distHaversine(
#                    cbind(wl_lon, wl_lat),
#                    cbind(longitude, latitude)) / 1000,
#     grav_raw   = Population / dist_km^1.5,
#     share_grav = grav_raw / sum(grav_raw, na.rm = TRUE)
#   ) %>%
#   select(name_cancensus, share_grav)
#
# panel_grav <- panel_balanced_wl %>%
#   left_join(gravity_wl, by = "name_cancensus") %>%
#   mutate(share_grav = replace_na(share_grav, 0),
#          Z_grav     = share_grav * Post2017)
#
# rf_gravity <- feols(
#   log_rent ~ Z_grav | name_cancensus + REF_DATE,
#   data    = panel_grav,
#   cluster = ~name_cancensus
# )

# ── M2-D : Taux de vacance comme variable dépendante (CMHC Table 34-10-0127) ─
# Une fois vacancy_rate_city.csv chargé (même structure que data2_two) :
#
# data_vacancy <- read.csv("raw/34100127-eng/34100127.csv") %>%
#   filter(grepl("2 bedroom", `Type.of.unit`),
#          grepl("row and apartment|apartment.*3 units", `Type.of.structure`,
#                ignore.case = TRUE)) %>%
#   select(REF_DATE, GEO, VALUE) %>%
#   rename(vacancy = VALUE) %>%
#   mutate(REF_DATE = as.integer(REF_DATE))
#
# panel_vacancy <- merge_loyers_generic(panel_instrument_wl, data_vacancy) %>%
#   cylindrer()
#
# rf_vacancy <- feols(
#   vacancy ~ Z_hat_wl | name_cancensus + REF_DATE,
#   data    = panel_vacancy,
#   cluster = ~name_cancensus
# )
# Lecture : si coef < 0 et significatif → exposition réduit la vacance =
# validation indépendante du mécanisme (pression demande, pas offre).

# =============================================================================
# PARTIE IX.7 -- TESTS DE ROBUSTESSE COMPLÉMENTAIRES (BLINDAGE CAUSAL)
# =============================================================================
# POURQUOI : la PARTIE VII traite les placebos en temps (années fictives
# 2010-2016) et la PARTIE IX.3-IX.6 traite les hard checks Vancouver et les
# soft checks Quesnel/gravité. Cette section ajoute quatre vérifications
# additionnelles centrales dans la pratique Angrist-Pischke :
#   (i)   Leave-one-out itératif sur les villes les plus exposées : prouve
#         que l'identification n'est pas concentrée sur 1-2 leverage points.
#   (ii)  Spécifications alternatives de clustering (province, bidimensionnel)
#         et de pondération (WLS population, WLS loyer 2008).
#   (iii) Hétérogénéité géographique : effet par taille de ville et par
#         région (West vs East), au-delà de l'hétérogénéité structurelle 4x4.
#   (iv)  Inférence par permutation (randomization inference) : p-value
#         calibrée empiriquement sans hypothèse asymptotique.
# Toutes les estimations s'appuient sur panel_balanced_wl (131 villes × 17 ans).

# -----------------------------------------------------------------------------
# IX.7.1  Leave-one-out itératif sur les villes les plus exposées
# -----------------------------------------------------------------------------
# Logique : on retire séquentiellement les k villes les plus exposées
# (k = 1, 2, 3, 5, 10) et on ré-estime la spécification baseline. Si le
# coefficient ne survit qu'avec Kamloops et Vancouver, drop des deux
# l'effondre. S'il survit, la variation identifiante est distribuée sur
# l'ensemble du support de s_WL,d.

exposed_ranked <- shares_wl %>%
  filter(share_wl > 0) %>%
  arrange(desc(share_wl)) %>%
  pull(name_cancensus)

ranks_to_drop <- c(1, 2, 3, 5, 10)

loo_models <- lapply(ranks_to_drop, function(k) {
  drop_cities <- exposed_ranked[seq_len(k)]
  feols(log_rent ~ Z_hat_wl | name_cancensus + REF_DATE,
        data    = panel_balanced_wl %>%
                    filter(!name_cancensus %in% drop_cities),
        cluster = ~name_cancensus)
})
names(loo_models) <- paste0("Drop top ", ranks_to_drop)

cat("\n=== LEAVE-ONE-OUT PAR RANG D'EXPOSITION (IX.7.1) ===\n")
loo_summary <- data.frame(
  k_dropped = ranks_to_drop,
  coef      = vapply(loo_models, function(m) unname(coef(m)["Z_hat_wl"]),
                     numeric(1)),
  se        = vapply(loo_models,
                     function(m) sqrt(vcov(m)["Z_hat_wl", "Z_hat_wl"]),
                     numeric(1)),
  n_cities  = vapply(loo_models,
                     function(m) length(unique(m$fixef_id[[1]])),
                     numeric(1))
)
print(loo_summary)


# -----------------------------------------------------------------------------
# IX.7.2  Spécifications alternatives : clustering et pondération
# -----------------------------------------------------------------------------
# La spécification baseline cluster les SE au niveau ville (131 clusters).
# Trois alternatives stressent l'inférence et la pondération :
#   (a) Clustering province (10 clusters) : capture la corrélation spatiale
#       intra-provinciale, pertinent si le choc BC se diffuse régionalement.
#   (b) Clustering bidimensionnel ville × année : corrige une éventuelle
#       corrélation transversale entre villes la même année (chocs macro).
#   (c) WLS par population : poids = taille du marché locatif, donne plus
#       d'importance aux grandes villes (réduit le bruit des petits marchés).
#   (d) WLS par loyer 2008 : poids = niveau de loyer initial, contrôle pour
#       les différences structurelles de prix à la base.

rf_cluster_prov <- feols(
  log_rent ~ Z_hat_wl | name_cancensus + REF_DATE,
  data = panel_balanced_wl, cluster = ~province
)

rf_cluster_twoway <- feols(
  log_rent ~ Z_hat_wl | name_cancensus + REF_DATE,
  data = panel_balanced_wl, cluster = ~name_cancensus + REF_DATE
)

panel_wls_pop <- panel_balanced_wl %>%
  left_join(coords %>% select(name_cancensus, Population),
            by = "name_cancensus") %>%
  filter(!is.na(Population), Population > 0)

rf_wls_pop <- feols(
  log_rent ~ Z_hat_wl | name_cancensus + REF_DATE,
  data    = panel_wls_pop,
  cluster = ~name_cancensus,
  weights = ~Population
)

rent_2008 <- panel_balanced_wl %>%
  filter(REF_DATE == 2008) %>%
  select(name_cancensus, rent_2008 = rent)

panel_wls_rent <- panel_balanced_wl %>%
  left_join(rent_2008, by = "name_cancensus") %>%
  filter(!is.na(rent_2008), rent_2008 > 0)

rf_wls_rent <- feols(
  log_rent ~ Z_hat_wl | name_cancensus + REF_DATE,
  data    = panel_wls_rent,
  cluster = ~name_cancensus,
  weights = ~rent_2008
)

cat("\n=== CLUSTERING ET PONDÉRATION (IX.7.2) ===\n")
print(etable(rf_baseline, rf_cluster_prov, rf_cluster_twoway,
             rf_wls_pop, rf_wls_rent,
             headers = c("Baseline", "Cl. province", "Cl. two-way",
                         "WLS Pop", "WLS Rent 2008"),
             digits  = 3))


# -----------------------------------------------------------------------------
# IX.7.3  Hétérogénéité par taille de ville et par région
# -----------------------------------------------------------------------------
# Découpe taille : small (<50k), medium (50-200k), large (>=200k)
# Découpe région : West (BC, AB, SK, MB) vs East (ON, QC, Atlantique)
# Logique : si la rigidité de l'offre amplifie la réponse au choc de demande,
# l'effet doit être plus marqué dans les petites villes (stock thin, foncier
# rare) et dans l'Ouest (province BC, supply-constrained selon Saiz 2010).

panel_het <- panel_balanced_wl %>%
  left_join(coords %>% select(name_cancensus, Population),
            by = "name_cancensus") %>%
  mutate(
    size_cat = case_when(
      is.na(Population)   ~ NA_character_,
      Population <  50000 ~ "small",
      Population < 200000 ~ "medium",
      TRUE                ~ "large"
    ),
    region = ifelse(
      province %in% c("British Columbia", "Alberta",
                      "Saskatchewan",      "Manitoba"),
      "West", "East"
    )
  )

rf_size_small  <- feols(log_rent ~ Z_hat_wl | name_cancensus + REF_DATE,
                        data    = panel_het %>%
                                  filter(size_cat == "small"),
                        cluster = ~name_cancensus)
rf_size_medium <- feols(log_rent ~ Z_hat_wl | name_cancensus + REF_DATE,
                        data    = panel_het %>%
                                  filter(size_cat == "medium"),
                        cluster = ~name_cancensus)
rf_size_large  <- feols(log_rent ~ Z_hat_wl | name_cancensus + REF_DATE,
                        data    = panel_het %>%
                                  filter(size_cat == "large"),
                        cluster = ~name_cancensus)

rf_region_west <- feols(log_rent ~ Z_hat_wl | name_cancensus + REF_DATE,
                        data    = panel_het %>% filter(region == "West"),
                        cluster = ~name_cancensus)
rf_region_east <- feols(log_rent ~ Z_hat_wl | name_cancensus + REF_DATE,
                        data    = panel_het %>% filter(region == "East"),
                        cluster = ~name_cancensus)

cat("\n=== HÉTÉROGÉNÉITÉ TAILLE × RÉGION (IX.7.3) ===\n")
print(etable(rf_size_small, rf_size_medium, rf_size_large,
             rf_region_west, rf_region_east,
             headers = c("Small <50k", "Medium 50-200k", "Large >=200k",
                         "West", "East"),
             digits  = 3))


# -----------------------------------------------------------------------------
# IX.7.3bis  Hétérogénéité supply-side : marchés tendus vs. relâchés (Saiz 2010)
# -----------------------------------------------------------------------------
# Prédiction théorique (Saiz 2010) : la réponse en loyers à un choc de demande
# est 2-3× plus forte dans les marchés où l'offre est inélastique.
# Proxy : loyer moyen pré-2017 (2008-2016) au-dessus/en-dessous de la médiane
# → les marchés historiquement chers = plus contraints en offre (terrain rare,
#   zoning strict, géographie physique limitante)
# Spécification : log R_dt = π₁ Z_dt + π₂ (Z_dt × Tight_d) + γ_d + λ_t + ε_dt
# Prédiction Saiz : π₂ > 0

# 1. Construire l'indicateur de tension pré-traitement
pre2017_rent_summary <- panel_balanced_wl %>%
  filter(REF_DATE <= 2016) %>%
  group_by(name_cancensus) %>%
  summarise(mean_pre_rent = mean(rent, na.rm = TRUE), .groups = "drop") %>%
  mutate(tight_market = as.integer(mean_pre_rent >= median(mean_pre_rent, na.rm = TRUE)))

median_pre_rent <- median(pre2017_rent_summary$mean_pre_rent, na.rm = TRUE)
n_tight <- sum(pre2017_rent_summary$tight_market)
n_slack <- sum(1L - pre2017_rent_summary$tight_market)

cat(sprintf("\n=== SUPPLY-SIDE HETEROGENEITY (IX.7.3bis) ===\n"))
cat(sprintf("Médiane loyer pré-2017 : $%.0f\n", median_pre_rent))
cat(sprintf("Marchés tendus (tight, above median) : n=%d villes\n", n_tight))
cat(sprintf("Marchés relâchés (slack, below median) : n=%d villes\n", n_slack))

# 2. Joindre au panel et créer terme d'interaction
panel_supply_het <- panel_balanced_wl %>%
  left_join(pre2017_rent_summary %>% select(name_cancensus, mean_pre_rent, tight_market),
            by = "name_cancensus") %>%
  mutate(Z_tight = Z_hat_wl * tight_market)

# 3. Régression avec interaction formelle
rf_supply_interaction <- feols(
  log_rent ~ Z_hat_wl + Z_tight | name_cancensus + REF_DATE,
  data    = panel_supply_het,
  cluster = ~name_cancensus
)

# 4. Régressions séparées par sous-marché (complément interprétatif)
rf_supply_tight <- feols(
  log_rent ~ Z_hat_wl | name_cancensus + REF_DATE,
  data    = panel_supply_het %>% filter(tight_market == 1L),
  cluster = ~name_cancensus
)
rf_supply_slack <- feols(
  log_rent ~ Z_hat_wl | name_cancensus + REF_DATE,
  data    = panel_supply_het %>% filter(tight_market == 0L),
  cluster = ~name_cancensus
)

print(etable(rf_supply_interaction, rf_supply_tight, rf_supply_slack,
             headers = c("Interaction", "Tight (above median)", "Slack (below median)"),
             digits  = 3))

# 5. Export LaTeX
tex_supply_het <- etable(
  rf_supply_interaction, rf_supply_tight, rf_supply_slack,
  headers = c("Interaction", "Tight markets", "Slack markets"),
  digits  = 3,
  tex     = TRUE,
  notes   = paste0(
    "\\textit{Notes:} Tight markets: cities with mean 2008--2016 rent ",
    "above the sample median (\\$", sprintf("%.0f", median_pre_rent), "). ",
    "Column~(1): pooled specification with interaction term ",
    "$Z_{WL,dt}\\times\\mathrm{Tight}_d$. ",
    "Columns~(2)--(3): separate regressions on tight and slack subsamples. ",
    "All specifications include city and year fixed effects. ",
    "Standard errors clustered at the city level. ",
    "* $p<0.10$, ** $p<0.05$, *** $p<0.01$."
  )
)
# Encapsuler dans un environnement table LaTeX propre
supply_het_full <- paste0(
  "\\begin{table}[htbp]\n",
  "  \\centering\\small\n",
  "  \\caption{Supply-Side Heterogeneity: Tight vs.\\ Slack Rental Markets}\n",
  "  \\label{tab:supply_het}\n",
  tex_supply_het,
  "\\end{table}\n"
)
writeLines(supply_het_full, "tables/table_supply_het.tex")
cat("→ table_supply_het.tex exportée\n")

# -----------------------------------------------------------------------------
# IX.7.3ter  Coefficient plot combiné — city size / region / supply constraints
# Remplace les références aux deux tables (tab:heterogeneity_size_region et
# tab:supply_het) par une figure unique (fig:heterogeneity_all_dims) dans le
# texte principal. Les tables restent en annexe pour la précision numérique.
# -----------------------------------------------------------------------------
all_dims <- data.frame(
  label     = c("Small cities", "Medium cities", "Large cities",
                "West", "East",
                "Tight supply", "Slack supply"),
  dimension = c(rep("City size", 3), rep("Region", 2), rep("Supply", 2)),
  coef      = c(coef(rf_size_small)["Z_hat_wl"],
                coef(rf_size_medium)["Z_hat_wl"],
                coef(rf_size_large)["Z_hat_wl"],
                coef(rf_region_west)["Z_hat_wl"],
                coef(rf_region_east)["Z_hat_wl"],
                coef(rf_supply_tight)["Z_hat_wl"],
                coef(rf_supply_slack)["Z_hat_wl"]),
  se        = c(sqrt(vcov(rf_size_small)["Z_hat_wl","Z_hat_wl"]),
                sqrt(vcov(rf_size_medium)["Z_hat_wl","Z_hat_wl"]),
                sqrt(vcov(rf_size_large)["Z_hat_wl","Z_hat_wl"]),
                sqrt(vcov(rf_region_west)["Z_hat_wl","Z_hat_wl"]),
                sqrt(vcov(rf_region_east)["Z_hat_wl","Z_hat_wl"]),
                sqrt(vcov(rf_supply_tight)["Z_hat_wl","Z_hat_wl"]),
                sqrt(vcov(rf_supply_slack)["Z_hat_wl","Z_hat_wl"])),
  stringsAsFactors = FALSE
)
all_dims$label     <- factor(all_dims$label,
                              levels = rev(c("Small cities","Medium cities","Large cities",
                                             "West","East","Tight supply","Slack supply")))
all_dims$dimension <- factor(all_dims$dimension,
                              levels = c("City size","Region","Supply"))

baseline_coef <- unname(coef(rf_baseline)["Z_hat_wl"])

p_all_dims <- ggplot(all_dims,
                     aes(x = coef, y = label)) +
  geom_vline(xintercept = baseline_coef,
             linetype = "dashed", colour = "grey50", linewidth = 0.4) +
  geom_vline(xintercept = 0, colour = "black", linewidth = 0.3) +
  geom_errorbarh(aes(xmin = coef - 1.96 * se,
                     xmax = coef + 1.96 * se),
                 height = 0.25, linewidth = 0.5, colour = "#444444") +
  geom_point(size = 2.5, colour = "#444444") +
  facet_grid(dimension ~ ., scales = "free_y", space = "free_y") +
  labs(x = "Coefficient on WL exposure index", y = NULL,
       caption = "Dashed line: baseline TWFE estimate (0.818). Error bars: 95% CI.") +
  theme_classic(base_size = 10) +
  theme(panel.grid.major.x = element_line(colour = "grey90", linewidth = 0.3),
        axis.line.y        = element_blank(),
        axis.ticks.y       = element_blank(),
        strip.background   = element_blank(),
        strip.text.y       = element_text(face = "bold", size = 9, angle = 0,
                                          hjust = 0))

ggsave("figures/fig_heterogeneity_all_dims.pdf", p_all_dims,
       width = 7, height = 4.5)
ggsave("figures/fig_heterogeneity_all_dims.png", p_all_dims,
       width = 7, height = 4.5, dpi = 300)
cat("fig_heterogeneity_all_dims exported.\n")

# -----------------------------------------------------------------------------
# IX.7.4  Inférence par permutation (randomization inference, in-space placebo)
# -----------------------------------------------------------------------------
# Logique : sous l'hypothèse nulle "aucun effet causal de l'exposition à WL",
# l'attribution des shares aux villes est arbitraire. On permute donc
# aléatoirement les shares pré-feu entre les 131 villes et on ré-estime le
# coefficient sur chaque permutation. La distribution empirique des
# coefficients sous H0 donne une p-value calibrée sans hypothèse asymptotique
# (Fisher 1935, Young 2019 QJE). Si le coef baseline tombe dans une queue
# extrême de cette distribution, H0 est rejetée.

set.seed(20260514)
n_perm <- 500

shares_unique <- panel_balanced_wl %>%
  distinct(name_cancensus, share_wl)

permutation_coefs <- vapply(seq_len(n_perm), function(i) {
  perm <- shares_unique %>%
    mutate(share_perm = sample(share_wl)) %>%
    select(name_cancensus, share_perm)
  panel_p <- panel_balanced_wl %>%
    left_join(perm, by = "name_cancensus") %>%
    mutate(Z_perm = share_perm * Post2017)
  m <- feols(log_rent ~ Z_perm | name_cancensus + REF_DATE,
             data = panel_p)
  unname(coef(m)["Z_perm"])
}, numeric(1))

p_two_sided <- mean(abs(permutation_coefs) >= abs(baseline_coef))
p_one_sided <- mean(permutation_coefs       >= baseline_coef)

cat("\n=== RANDOMIZATION INFERENCE (IX.7.4) ===\n")
cat(sprintf("Coefficient baseline       : %.4f\n",        baseline_coef))
cat(sprintf("Distribution placebo (mean): %.4f (sd=%.4f)\n",
            mean(permutation_coefs), sd(permutation_coefs)))
cat(sprintf("Quantiles 2.5%% / 97.5%%     : %.4f / %.4f\n",
            quantile(permutation_coefs, 0.025),
            quantile(permutation_coefs, 0.975)))
cat(sprintf("p-value bilatérale         : %.4f\n", p_two_sided))
cat(sprintf("p-value unilatérale        : %.4f\n", p_one_sided))


# -----------------------------------------------------------------------------
# IX.7.5bis  Poids de Rotemberg -- Diagnostique d'influence shift-share
# -----------------------------------------------------------------------------
# POURQUOI (Goldsmith-Pinkham, Sorkin & Swift 2020, section III.B) :
# Dans un design shift-share avec un seul choc agrégé (le feu de WL), le coefficient
# OLS peut s'écrire comme une moyenne pondérée par les poids de Rotemberg :
#   π̂ = Σ_d w_d^R × α̂_d
# où α̂_d est l'estimateur just-identified pour la ville d seule et w_d^R est son poids.
# Les poids de Rotemberg mesurent la contribution de chaque ville à l'identification.
#
# CAS SIMPLE (panel cylindré, même Post pour toutes les villes, FE ville+année) :
# Dans la spécification réduite log_rent_dt = π × Z_dt + ε_dt,
# après déméanage ville×année, la contribution de la ville d est proportionnelle à
# la covariance de s_d × Post_t avec log_rent_dt, divisée par la variance totale
# de Z_dt. Avec un seul choc binaire Post, le poids se simplifie à :
#   w_d^R = s_d^2 / Σ_{d'} s_{d'}^2
# (vrai si même nombre d'années pre/post pour toutes les villes -- panel cylindré ici).
#
# NOTE : Cette formule suppose l'exogénéité des parts s_d (condition Goldsmith 2020).
# Elle ne dépend pas de la structure des erreurs ni du clustering.

rotemberg <- shares_wl %>%
  replace_na(list(share_wl = 0)) %>%
  mutate(
    s2           = share_wl^2,
    w_rotemberg  = s2 / sum(s2, na.rm = TRUE),
    w_pct        = w_rotemberg * 100
  ) %>%
  arrange(desc(w_rotemberg)) %>%
  filter(share_wl > 0)

cat("\n=== POIDS DE ROTEMBERG (IX.7.5bis) ===\n")
cat(sprintf("%-20s %8s %10s %12s\n",
            "Ville", "s_WL (%)", "w_R (%)", "w_R cumulé"))
cumul <- 0
for (i in seq_len(min(nrow(rotemberg), 15))) {
  cumul <- cumul + rotemberg$w_pct[i]
  cat(sprintf("%-20s %8.2f %10.2f %12.2f\n",
              rotemberg$name_cancensus[i],
              rotemberg$share_wl[i] * 100,
              rotemberg$w_pct[i],
              cumul))
}
cat(sprintf("\nTop 5 villes : %.1f%% du poids total d'identification\n",
            sum(rotemberg$w_pct[1:5])))
cat(sprintf("Top 2 (Kamloops + Vancouver) : %.1f%%\n",
            sum(rotemberg$w_pct[1:2])))

# Export tableau Rotemberg pour appendice LaTeX
rotemberg_top <- rotemberg %>% slice_head(n = 10)

tex_rotemberg <- paste0(
  "\\begin{table}[htbp]\n",
  "  \\caption{Rotemberg Weights by Pre-Fire Exposure Rank}\n",
  "  \\label{tab:rotemberg}\n",
  "  \\centering\\small\n",
  "  \\begin{threeparttable}\n",
  "  \\begin{tabular}{lrrrr}\n",
  "    \\toprule\n",
  "    City & $s_{WL,d}$ (\\%) & $s_{WL,d}^2$ & $w_d^R$ (\\%) & Cumulative (\\%) \\\\\n",
  "    \\midrule\n"
)
cumul_tex <- 0
for (i in seq_len(nrow(rotemberg_top))) {
  cumul_tex <- cumul_tex + rotemberg_top$w_pct[i]
  tex_rotemberg <- paste0(tex_rotemberg,
    sprintf("    %s & %.2f & %.4f & %.2f & %.2f \\\\\n",
            rotemberg_top$name_cancensus[i],
            rotemberg_top$share_wl[i] * 100,
            rotemberg_top$s2[i],
            rotemberg_top$w_pct[i],
            cumul_tex))
}
tex_rotemberg <- paste0(tex_rotemberg,
  "    \\addlinespace\n",
  sprintf("    \\textit{All other cities (N=%d)} & %.2f & --- & %.2f & 100.00 \\\\\n",
          nrow(rotemberg) - 10,
          sum(rotemberg$share_wl[11:nrow(rotemberg)], na.rm = TRUE) * 100,
          sum(rotemberg$w_pct[11:nrow(rotemberg)], na.rm = TRUE)),
  "    \\bottomrule\n",
  "  \\end{tabular}\n",
  "  \\begin{tablenotes}[flushleft]\\footnotesize\n",
  "    \\item \\textit{Notes:} Rotemberg weights $w_d^R = s_{WL,d}^2 / \\sum_{d'} s_{WL,d'}^2$",
  " (\\citealt{goldsmith2020}, Section~III.B). In this balanced panel with a single",
  " binary Post indicator, these weights equal the squared pre-fire migration share",
  " normalised to sum to one across the 43 exposed cities. Each weight measures",
  " city~$d$'s contribution to the overall identification of $\\hat{\\pi}$.",
  " Non-exposed cities ($s_{WL,d}=0$) have zero Rotemberg weight and serve as the",
  " comparison group, not as identifying variation.\n",
  "  \\end{tablenotes}\n",
  "  \\end{threeparttable}\n",
  "\\end{table}\n"
)
writeLines(tex_rotemberg, "tables/table_rotemberg.tex")
cat("table_rotemberg.tex exporté\n")


# -----------------------------------------------------------------------------
# IX.7.5  Export LaTeX des trois tableaux complémentaires
# -----------------------------------------------------------------------------
etable(loo_models$`Drop top 1`,  loo_models$`Drop top 2`,
       loo_models$`Drop top 3`,  loo_models$`Drop top 5`,
       loo_models$`Drop top 10`,
       headers = c("Drop top 1", "Drop top 2", "Drop top 3",
                   "Drop top 5", "Drop top 10"),
       fitstat = ~ n + r2 + wr2 + rmse,
       title   = "Robustness: Leave-one-out by exposure rank",
       label   = "tab:robustness_loo",
       digits  = 3, tex = TRUE,
       file    = "tables/table_robustness_loo.tex", replace = TRUE)

etable(rf_baseline, rf_cluster_prov, rf_cluster_twoway,
       rf_wls_pop,  rf_wls_rent,
       headers = c("Baseline", "Cluster province", "Cluster two-way",
                   "WLS Pop", "WLS Rent 2008"),
       title   = "Robustness: Alternative standard errors and weighting",
       label   = "tab:robustness_clustering",
       digits  = 3, tex = TRUE,
       file    = "tables/table_robustness_clustering.tex", replace = TRUE)

etable(rf_size_small, rf_size_medium, rf_size_large,
       rf_region_west, rf_region_east,
       headers = c("Small (<50k)", "Medium (50-200k)", "Large ($\\geq$200k)",
                   "West", "East"),
       title   = "Heterogeneity by city size and region",
       label   = "tab:heterogeneity_size_region",
       digits  = 3, tex = TRUE,
       file    = "tables/table_heterogeneity_size_region.tex", replace = TRUE)


# =============================================================================
# PARTIE X -- ANALYSE D'HÉTÉROGÉNÉITÉ (4×4)
# =============================================================================
# POURQUOI : L'hétérogénéité cross-structure valide que le résultat principal
# n'est pas un artefact du choix de la structure de référence.
# L'hétérogénéité cross-type d'unité révèle des différences de profil de
# demande entre déplacés (studios vs. grandes unités).
# 16 régressions : 4 structures × 4 types d'unité.
# Toutes utilisent Z_hat_wl (BC 2017) avec FE ville + année.

# Structure S1 : Apartment structures of three units and over (BASELINE)
rf_s1_bachelor <- estimer_rf(panel_s1_bachelor)
rf_s1_one      <- estimer_rf(panel_s1_one)
rf_s1_two      <- rf_baseline                   # = Table 1
rf_s1_three    <- estimer_rf(panel_s1_three)

# Structure S2 : Row and apartment structures of three units and over
rf_s2_bachelor <- estimer_rf(panel_s2_bachelor)
rf_s2_one      <- estimer_rf(panel_s2_one)
rf_s2_two      <- estimer_rf(panel_s2_two)
rf_s2_three    <- estimer_rf(panel_s2_three)

# Structure S3 : Row structures of three units and over
# → S3 bachelor non estimable (2 villes après cylindrage, singletons FE)
#   Vérifié avec tryCatch -- cellule marquée "---" dans le tableau LaTeX
rf_s3_bachelor <- NULL  # non estimable
rf_s3_one      <- estimer_rf(panel_s3_one)    # 14 villes -- interpréter avec prudence
rf_s3_two      <- estimer_rf(panel_s3_two)    # 53 villes ✓
rf_s3_three    <- estimer_rf(panel_s3_three)  # 47 villes ✓

# Structure S4 : Apartment structures of six units and over
rf_s4_bachelor <- estimer_rf(panel_s4_bachelor)
rf_s4_one      <- estimer_rf(panel_s4_one)
rf_s4_two      <- estimer_rf(panel_s4_two)
rf_s4_three    <- estimer_rf(panel_s4_three)

# Récapitulatif des coefficients (NA pour S3 bachelor)
safe_coef <- function(model) {
  if (is.null(model)) return(NA_real_)
  stats::coef(model)["Z_hat_wl"]
}

cat("\n=== MATRICE D'HÉTÉROGÉNÉITÉ 4×4 (β_wl) ===\n")
cat(sprintf("%-50s %8s %8s %8s %8s\n",
            "Structure", "Bachelor", "1ch", "2ch", "3ch"))
cat(sprintf("%-50s %8.3f %8.3f %8.3f %8.3f\n",
            "Apt 3+ units (BASELINE)",
            safe_coef(rf_s1_bachelor), safe_coef(rf_s1_one),
            safe_coef(rf_s1_two),      safe_coef(rf_s1_three)))
cat(sprintf("%-50s %8.3f %8.3f %8.3f %8.3f\n",
            "Row & apt 3+ units",
            safe_coef(rf_s2_bachelor), safe_coef(rf_s2_one),
            safe_coef(rf_s2_two),      safe_coef(rf_s2_three)))
cat(sprintf("%-50s %8s %8.3f %8.3f %8.3f\n",
            "Row structures 3+ units",
            "---",
            safe_coef(rf_s3_one), safe_coef(rf_s3_two), safe_coef(rf_s3_three)))
cat(sprintf("%-50s %8.3f %8.3f %8.3f %8.3f\n",
            "Apt 6+ units",
            safe_coef(rf_s4_bachelor), safe_coef(rf_s4_one),
            safe_coef(rf_s4_two),      safe_coef(rf_s4_three)))
# Résultats de la matrice 4×4 :
#
# Structure                    Bachelor   1ch     2ch     3ch
# S1 Apt 3+ units (BASELINE)   0.916***  0.753***  0.818***  0.864***
# S2 Row & apt 3+ units        0.883***  0.760***  0.783***  0.608**
# S3 Row 3+ units              ---       -4.288*** 0.424**   0.401
# S4 Apt 6+ units              0.787***  0.730***  0.798***  0.792***
#
# INTERPRÉTATION :
# Message principal : l'effet BC 2017 est robuste à travers S1, S2 et S4.
# S3 (row structures) : segment marginal, faible couverture géographique,
# résultats bruités -- à ne pas interpréter causalement.
# S1 (baseline) : effet plus fort pour les studios (0.916) que pour le 1 chambre
# (0.753) -- cohérent avec un profil de déplacés incluant des jeunes actifs.


# =============================================================================
# PARTIE XI -- CALIBRATION ÉCONOMIQUE
# =============================================================================
# POURQUOI : Convertir les log-points en dollars courants permet de situer
# l'effet dans le contexte des loyers canadiens et de comparer avec la
# littérature (Saiz 2007 : +1% population → +1% loyer).

beta_rf_wl <- coef(rf_baseline)["Z_hat_wl"]
loyer_2015 <- mean(panel_balanced_wl$rent[panel_balanced_wl$REF_DATE == 2015],
                   na.rm = TRUE)
cat("Loyer de référence 2015 :", round(loyer_2015, 0), "$\n")
# → Loyer moyen 2015 = 848$ (loyer 2 chambres, apt 3+ unités, moyenne des 131 villes)

panel_balanced_wl %>%
  filter(REF_DATE == 2017) %>%
  arrange(desc(Z_hat_wl)) %>%
  head(5) %>%
  mutate(
    effet_pct       = beta_rf_wl * Z_hat_wl * 100,
    effet_monetaire = loyer_2015 * (exp(beta_rf_wl * Z_hat_wl) - 1)
  ) %>%
  select(name_cancensus, share_wl, effet_pct, effet_monetaire)
# Top 5 villes les plus exposées en 2017 :
# Kamloops      : share_wl = 21.2% → effet +17.4% → +161$/mois
# Vancouver     : share_wl = 17.4% → effet +14.2% → +129$/mois
# Prince George : share_wl =  9.97% → effet  +8.2% →  +72$/mois
# Kelowna       : share_wl =  6.91% → effet  +5.7% →  +49$/mois
# Quesnel       : share_wl =  6.75% → effet  +5.5% →  +48$/mois
# → Les villes BC proches de Williams Lake absorbent la majorité des déplacés
# → L'effet monétaire est économiquement significatif pour Kamloops et Vancouver
# → Cohérent avec Saiz (2007) : +1% population → +1% loyer


# =============================================================================
# PARTIE XI.bis -- TABLES 4 ET 5 : CALIBRATION ÉCONOMIQUE ET MARKET TIGHTNESS
# =============================================================================
# POURQUOI : Ces deux tables n'étaient pas générées dans le pipeline R.
# Elles étaient construites directement dans un LLM, ce qui rompt la
# reproductibilité. Ce bloc les intègre dans le script en s'appuyant sur les
# objets déjà disponibles (rf_baseline, rf_trends, loyer_2015, shares_wl)
# et en hard-codant les données CMHC Rental Market Survey 2017 (snapshot
# statique d'octobre 2017, source : CMHC 2017 British Columbia Rental Market
# Report, private apartment structures ≥ 3 units, two-bedroom segment).

# -----------------------------------------------------------------------------
# XI.bis.1  Paramètres communs
# -----------------------------------------------------------------------------
evacuees_total <- 11000   # résidents sous ordre d'évacuation obligatoire (juillet 2017)
q_base         <- 0.25    # fraction de déplacés générant une demande locative (scénario intermédiaire)
h              <- 2       # personnes par ménage locataire (hypothèse de calibration)

beta_baseline <- coef(rf_baseline)["Z_hat_wl"]   # 0.818 -- spécification col. (4)
beta_trend    <- coef(rf_trends)["Z_hat_wl"]     # 0.594 -- spécification col. (5) avec tendances ville
loyer_ref     <- loyer_2015                       # 848$ -- loyer moyen 2015 (2 chambres, apt 3+ unités)

# -----------------------------------------------------------------------------
# XI.bis.2  Table 4 -- Calibration économique des magnitudes reduced-form
# -----------------------------------------------------------------------------
# Colonnes :
#   s_WL,d          : part migratoire pré-feu (depuis shares_wl)
#   displaced       : 11 000 × s_WL,d  (ménages potentiellement alloués à la ville d)
#   households      : H_d = q × 11 000 × s_WL,d / h  (ménages locataires potentiels, scénario q=0.25)
#   pi_s            : π̂ × s_WL,d  (différentiel log-loyer prédit, baseline)
#   delta_pct       : (exp(π̂ × s_WL,d) − 1) × 100  (effet en %, baseline)
#   delta_dollar    : Δ% × 848$  (effet en dollars courants, baseline)
#   delta_pct_trend : même calcul avec π̂ = 0.594 (city-trend)

top5 <- shares_wl %>%
  arrange(desc(share_wl)) %>%
  slice(1:5)

table4 <- top5 %>%
  mutate(
    displaced        = evacuees_total * share_wl,
    households       = q_base * evacuees_total * share_wl / h,
    pi_s             = beta_baseline * share_wl,
    delta_pct        = (exp(beta_baseline * share_wl) - 1) * 100,
    delta_dollar     = loyer_ref * (exp(beta_baseline * share_wl) - 1),
    delta_pct_trend  = (exp(beta_trend * share_wl) - 1) * 100
  )

print(table4)
# Vérification attendue (arrondi) :
# Kamloops      s=0.212  disp=2332  H=292  π̂s=0.173  Δ%=18.9  Δ$=160  Δ%trend=13.5
# Vancouver     s=0.174  disp=1914  H=239  π̂s=0.142  Δ%=15.3  Δ$=130  Δ%trend=10.9
# Prince George s=0.100  disp=1097  H=137  π̂s=0.082  Δ%= 8.5  Δ$= 72  Δ%trend= 6.1
# Kelowna       s=0.069  disp= 760  H= 95  π̂s=0.056  Δ%= 5.8  Δ$= 49  Δ%trend= 4.2
# Quesnel       s=0.068  disp= 743  H= 93  π̂s=0.055  Δ%= 5.7  Δ$= 48  Δ%trend= 4.1

# -----------------------------------------------------------------------------
# XI.bis.3  Table 5 -- Market tightness : demande potentielle vs unités vacantes
# -----------------------------------------------------------------------------
# Source CMHC Rental Market Survey 2017 :
#   CMHC (2017). British Columbia Rental Market Report, octobre 2017.
#   Segment : appartements privés ≥ 3 unités, deux chambres.
#   URL : https://www.cmhc-schl.gc.ca/professionals/housing-markets-data-and-research/
#          housing-research/rental-market-reports/british-columbia
# Ces données sont statiques (snapshot d'octobre 2017) → hard-coded justifié.

rms2017 <- data.frame(
  name_cancensus  = c("Kamloops", "Vancouver", "Prince George", "Kelowna", "Quesnel"),
  units_2br       = c(1499,       26375,        1505,            2341,       318),
  vacancy_rate    = c(0.011,      0.010,         0.030,           0.002,      0.030)
)

table5 <- rms2017 %>%
  left_join(
    table4 %>% select(name_cancensus, share_wl, households),
    by = "name_cancensus"
  ) %>%
  mutate(
    implied_vacant   = units_2br * vacancy_rate,
    # Ménages potentiels sous trois scénarios de q (h=2 dans chaque cas)
    H_q10  = 0.10 * evacuees_total * share_wl / h,
    H_q25  = 0.25 * evacuees_total * share_wl / h,   # = households de table4
    H_q50  = 0.50 * evacuees_total * share_wl / h,
    # Ratios H/V
    HV_q10 = H_q10 / implied_vacant,
    HV_q25 = H_q25 / implied_vacant,
    HV_q50 = H_q50 / implied_vacant
  )

print(table5 %>% select(name_cancensus, units_2br, vacancy_rate,
                         implied_vacant, H_q25, HV_q10, HV_q25, HV_q50))
# Vérification attendue (arrondi) :
# Kamloops      1499  1.1%  16.5  291.5  7.1  17.7  35.4
# Vancouver    26375  1.0% 263.8  239.3  0.4   0.9   1.8
# Prince George 1505  3.0%  45.2  137.1  1.2   3.0   6.1
# Kelowna       2341  0.2%   4.7   95.0  8.1  20.3  40.6
# Quesnel        318  3.0%   9.5   92.8  3.9   9.7  19.5


# =============================================================================
# PARTIE XII -- STATISTIQUES DESCRIPTIVES
# =============================================================================
# POURQUOI : Les statistiques descriptives permettent au lecteur de situer
# l'échantillon et de vérifier la plausibilité des ordres de grandeur.

# Panel propre : uniquement les observations utilisées dans l'estimation
# feols drope automatiquement les 2 NA de Leamington 2023-2024 → on fait pareil
panel_for_stats <- panel_balanced_wl %>%
  filter(!is.na(log_rent))

cat("Observations utilisées dans l'estimation :", nrow(panel_for_stats), "\n")
# → 2225 ✓ (cohérent avec le N reporté par feols)

# Panel A -- Statistiques globales sur l'ensemble du panel (2008-2024)
stats_global <- panel_for_stats %>%
  summarise(
    across(c(rent, log_rent, share_wl),
           list(mean   = ~mean(., na.rm = TRUE),
                median = ~median(., na.rm = TRUE),
                sd     = ~sd(., na.rm = TRUE),
                min    = ~min(., na.rm = TRUE),
                max    = ~max(., na.rm = TRUE)),
           .names = "{.col}_{.fn}")
  )

# Variation within/between (sur panel_for_stats)
stats_within <- panel_for_stats %>%
  group_by(name_cancensus) %>%
  mutate(rent_demeaned = rent - mean(rent, na.rm = TRUE)) %>%
  ungroup() %>%
  summarise(
    sd_within  = sd(rent_demeaned, na.rm = TRUE),
    sd_between = sd(tapply(rent, name_cancensus, mean, na.rm = TRUE))
  )

# Panel B -- Comparaison exposed vs non-exposed × pre/post 2017
stats_panel_b <- panel_for_stats %>%
  mutate(exposed = share_wl > 0,
         period  = ifelse(REF_DATE < 2017, "Pre-2017", "Post-2017")) %>%
  group_by(exposed, period) %>%
  summarise(
    rent_mean = mean(rent, na.rm = TRUE),
    rent_sd   = sd(rent, na.rm = TRUE),
    n_villes  = n_distinct(name_cancensus),
    n_obs     = n(),
    .groups   = "drop"
  )

# Affichage pour vérification
print(stats_global)
print(stats_within)
print(stats_panel_b)


# =============================================================================
# PARTIE XII.bis -- FIGURES DESCRIPTIVES POUR LA SECTION DONNÉES
# =============================================================================
# POURQUOI : Quatre figures destinées à la Section 3 (Data) du manuscrit, à
# insérer avant les tableaux de régression. Thème minimaliste académique :
# pas de quadrillage, axes nus, palette discrète noir/gris/blanc. Légendes
# techniques uniquement, toute interprétation économique reste dans le texte.

# -----------------------------------------------------------------------------
# Thème commun aux quatre figures
# -----------------------------------------------------------------------------
thesis_theme <- theme_classic(base_size = 11) +
  theme(
    panel.grid       = element_blank(),
    axis.line        = element_line(color = "black", linewidth = 0.4),
    axis.ticks       = element_line(color = "black", linewidth = 0.3),
    legend.position  = "bottom",
    legend.title     = element_blank(),
    legend.key       = element_rect(fill = "white", color = NA),
    plot.title       = element_blank(),
    plot.background  = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    strip.background = element_blank()
  )

# -----------------------------------------------------------------------------
# FIGURE D1 -- Évolution annuelle du loyer moyen, exposées vs non-exposées
# -----------------------------------------------------------------------------
# Moyenne par groupe et par année du loyer nominal des appartements de
# deux chambres dans les structures de trois unités et plus, sur le panel
# cylindré 131 villes x 2008-2024. Ligne verticale en 2017 = année du feu.

fig_d1_data <- panel_balanced_wl %>%
  mutate(group = ifelse(share_wl > 0,
                        "Exposed (N = 43)",
                        "Non-exposed (N = 88)")) %>%
  group_by(REF_DATE, group) %>%
  summarise(mean_rent = mean(rent, na.rm = TRUE), .groups = "drop")

label_d1 <- fig_d1_data %>% filter(REF_DATE == max(REF_DATE))

fig_d1 <- ggplot(fig_d1_data,
                 aes(x = REF_DATE, y = mean_rent,
                     linetype = group, shape = group, color = group)) +
  geom_line(linewidth = 0.6) +
  geom_point(size = 1.8, fill = "white") +
  scale_color_manual(values = c("Exposed (N = 43)"    = "#E69F00",
                                "Non-exposed (N = 88)" = "#999999")) +
  geom_vline(xintercept = 2016.75, linetype = "dotted",
             color = "black", linewidth = 0.4) +
  annotate("text", x = 2016.9,
           y = max(fig_d1_data$mean_rent) * 0.98,
           label = "Williams Lake fire (July 2017)",
           hjust = 0, size = 3.0) +
  geom_text(data = label_d1, aes(label = group),
            nudge_x = 0.3, size = 2.8, hjust = 0) +
  scale_x_continuous(breaks = seq(2008, 2024, 2),
                     expand = expansion(mult = c(0.02, 0.22))) +
  scale_shape_manual(values = c(16, 1)) +
  scale_linetype_manual(values = c("solid", "dashed")) +
  labs(x = "Year",
       y = "Mean monthly rent (CAD)") +
  thesis_theme +
  theme(legend.position = "none")

ggsave("figures/fig_rent_evolution.pdf", fig_d1,
       width = 7, height = 4.2)

# -----------------------------------------------------------------------------
# FIGURE D2 -- Distribution du log loyer, pré-2017 vs post-2017, par groupe
# -----------------------------------------------------------------------------
# Boxplots côte à côte : log loyer agrégé sur la période pré-2017 et la
# période post-2017, séparé entre les villes exposées et non-exposées.
# Permet de visualiser à la fois le niveau et la dispersion par groupe.

fig_d2_data <- panel_balanced_wl %>%
  mutate(group  = ifelse(share_wl > 0, "Exposed", "Non-exposed"),
         period = ifelse(REF_DATE >= 2017, "Post-2017", "Pre-2017"),
         period = factor(period, levels = c("Pre-2017", "Post-2017")))

fig_d2 <- ggplot(fig_d2_data,
                 aes(x = period, y = log_rent, fill = group)) +
  geom_boxplot(width = 0.55,
               outlier.size = 0.6, outlier.shape = 21,
               color = "black", linewidth = 0.3) +
  scale_fill_manual(values = c("grey45", "white")) +
  annotate("text", x = 0.78, y = 7.7, label = "Exposed",
           size = 3, fontface = "italic", color = "grey20") +
  annotate("text", x = 1.22, y = 7.7, label = "Non-exposed",
           size = 3, fontface = "italic", color = "grey20") +
  labs(x = NULL, y = "Log monthly rent (Canadian dollars)") +
  thesis_theme +
  theme(legend.position = "none")

ggsave("figures/fig_rent_distribution.pdf", fig_d2,
       width = 7, height = 4.2)

# -----------------------------------------------------------------------------
# FIGURE D3 -- Distribution de l'intensité de traitement parmi les exposées
# -----------------------------------------------------------------------------
# Histogramme des parts migratoires pré-feu s_WL,d pour les 43 villes ayant
# reçu au moins un sortant de Williams Lake sur la vague 2016/2017. Met en
# évidence la concentration extrême sur les premières destinations.

fig_d3_data <- shares_wl %>%
  filter(share_wl > 0) %>%
  mutate(share_pct = share_wl * 100)

top3_d3 <- fig_d3_data %>% slice_max(share_pct, n = 3)

fig_d3 <- ggplot(fig_d3_data, aes(x = share_pct)) +
  geom_histogram(binwidth = 1,
                 fill = "grey55", color = "black",
                 linewidth = 0.25, boundary = 0) +
  geom_text(data = top3_d3,
            aes(y = 1.3, label = name_cancensus),
            size = 2.8, angle = 0, hjust = 0.5) +
  scale_x_continuous(breaks = seq(0, 22, 2),
                     expand = expansion(mult = c(0.01, 0.02))) +
  scale_y_continuous(breaks = seq(0, 40, 5),
                     expand = expansion(mult = c(0, 0.05))) +
  labs(x = expression("Pre-fire migration share " * s["WL,d"] *
                      " (% of Williams Lake out-migrants)"),
       y = "Number of exposed cities (N = 43)") +
  thesis_theme

ggsave("figures/fig_treatment_intensity.pdf", fig_d3,
       width = 7, height = 4.2)

# -----------------------------------------------------------------------------
# FIGURE D4 -- Carte spatiale de l'intensité d'exposition
# -----------------------------------------------------------------------------
# Carte par points (longitude/latitude, coord_quickmap), taille du point
# proportionnelle à s_WL,d, forme distinguant les villes exposées des
# non-exposées. Williams Lake marqué par une étoile.

fig_d4_data <- panel_balanced_wl %>%
  distinct(name_cancensus, share_wl) %>%
  left_join(coords, by = "name_cancensus") %>%
  filter(!is.na(longitude), !is.na(latitude)) %>%
  mutate(exposure = ifelse(share_wl > 0, "Exposed", "Non-exposed"),
         size_var = pmax(share_wl, 0.001))

fig_d4 <- ggplot(fig_d4_data,
                 aes(x = longitude, y = latitude)) +
  geom_point(aes(size = size_var, shape = exposure,
                 fill = exposure, color = exposure),
             alpha = 0.75, stroke = 0.4) +
  scale_fill_manual(values = c("Exposed" = "#E69F00", "Non-exposed" = "#999999"),
                    name = NULL) +
  scale_color_manual(values = c("Exposed" = "#E69F00", "Non-exposed" = "#999999"),
                     name = NULL) +
  scale_size_continuous(range = c(0.7, 8),
                        breaks = c(0.01, 0.05, 0.10, 0.20),
                        labels = c("1%", "5%", "10%", "20%"),
                        name = "Pre-fire share\nfrom Williams Lake (%)") +
  scale_shape_manual(values = c("Exposed" = 21, "Non-exposed" = 4),
                     name = NULL,
                     guide = guide_legend(override.aes =
                                          list(size = 3,
                                               fill = c("#E69F00", NA),
                                               color = c("#E69F00", "#999999")))) +
  annotate("point", x = -122.14, y = 52.13,
           shape = 8, size = 3.2, stroke = 0.8, color = "#0072B2") +
  annotate("text", x = -122.14, y = 52.7,
           label = "Williams Lake", size = 2.9) +
  geom_text_repel(data = fig_d4_data %>% slice_max(size_var, n = 5),
                  aes(label = name_cancensus),
                  size = 2.5, box.padding = 0.3,
                  min.segment.length = 0, segment.color = "gray60",
                  max.overlaps = 15) +
  labs(x = "Longitude", y = "Latitude") +
  coord_quickmap() +
  thesis_theme +
  theme(legend.position = "right",
        legend.box = "vertical",
        legend.title = element_text(size = 9))

ggsave("figures/fig_spatial_exposure.pdf", fig_d4,
       width = 8, height = 5)


# =============================================================================
# PARTIE XII.5 -- CARTE D'EXPOSITION SPATIALE (figure principale du mémoire)
# =============================================================================
# Cette section reproduit figures/fig_spatial_exposure_qgis.pdf, la carte
# utilisée dans le corps du mémoire (Figure 1). Elle lit cities_exposure_arcgis.csv
# (exporté depuis export_arcgis.R) et ne dépend d'aucun objet R antérieur.
# Packages requis : rnaturalearth, rnaturalearthdata (chargés en PARTIE I).

{
  pts_map <- read.csv("cities_exposure_arcgis.csv", stringsAsFactors = FALSE)

  exposed_map     <- pts_map |> filter(exposure_label == "Exposed")
  non_exposed_map <- pts_map |> filter(exposure_label == "Non-exposed")
  origin_map      <- pts_map |> filter(exposure_label == "Origin (Williams Lake)")

  canada_prov_map <- tryCatch(
    sf::read_sf("https://raw.githubusercontent.com/nvkelso/natural-earth-vector/master/geojson/ne_50m_admin_1_states_provinces.geojson") |>
      filter(admin == "Canada"),
    error = function(e) NULL
  )

  canada_map <- ne_countries(country = "Canada",
                             scale = "medium", returnclass = "sf")
  usa_map    <- ne_countries(country = "United States of America",
                             scale = "medium", returnclass = "sf")

  pal_map <- c("0-1%" = "#fdcc8a", "1-3%" = "#fc8d59",
               "3-7%" = "#e34a33", ">7%"  = "#b30000")

  exposed_map <- exposed_map |>
    mutate(bin = cut(share_wl_pct,
                     breaks = c(0, 1, 3, 7, 100),
                     labels = c("0-1%", "1-3%", "3-7%", ">7%"),
                     include.lowest = TRUE))

  p_map <- ggplot() +
    geom_sf(data = canada_map, fill = "#f5f5f0", colour = NA) +
    { if (!is.null(canada_prov_map))
        geom_sf(data = canada_prov_map, fill = NA,
                colour = "grey65", linewidth = 0.22)
      else
        geom_sf(data = canada_map, fill = NA,
                colour = "grey60", linewidth = 0.3) } +
    geom_sf(data = usa_map, fill = "#ebebeb",
            colour = "grey70", linewidth = 0.2) +
    geom_point(data = non_exposed_map,
               aes(x = longitude, y = latitude),
               shape = 21, size = 1.6,
               fill = "grey80", colour = "grey50", stroke = 0.3, alpha = 0.7) +
    geom_point(data = exposed_map,
               aes(x = longitude, y = latitude,
                   fill = bin, size = share_wl_pct),
               shape = 21, colour = "white", stroke = 0.35, alpha = 0.9) +
    geom_point(data = origin_map,
               aes(x = longitude, y = latitude),
               shape = 23, size = 5, fill = "#2c7bb6",
               colour = "white", stroke = 0.8) +
    geom_text_repel(data = origin_map,
                    aes(x = longitude, y = latitude, label = name_cancensus),
                    size = 2.5, fontface = "bold", colour = "#2c7bb6",
                    nudge_y = 1.2, segment.colour = "grey50",
                    segment.size = 0.3) +
    scale_fill_manual(values = pal_map,
                      name = "Pre-fire share\n(% of WL out-migrants)") +
    scale_size_continuous(range = c(2, 9), guide = "none") +
    coord_sf(xlim = c(-140, -52), ylim = c(42, 62), expand = FALSE) +
    labs(
      title    = "Pre-fire Williams Lake Out-Migration Shares",
      subtitle = "Share of 2016/2017 Williams Lake out-migrants received by each city",
      caption  = paste0(
        "Notes: Filled circles are the 43 cities with positive pre-fire exposure (s > 0). ",
        "Hollow grey circles are the 88 non-exposed cities (s = 0). ",
        "Diamond marks Williams Lake (origin). ",
        "Exposure shares drawn from Statistics Canada internal migration data, 2016/2017.")
    ) +
    theme_void(base_size = 10) +
    theme(
      plot.title      = element_text(face = "bold", size = 11, hjust = 0,
                                     margin = margin(b = 3)),
      plot.subtitle   = element_text(size = 8.5, hjust = 0, colour = "grey30",
                                     margin = margin(b = 6)),
      plot.caption    = element_text(size = 6.5, hjust = 0, colour = "grey40",
                                     margin = margin(t = 6)),
      legend.position      = c(0.02, 0.25),
      legend.justification = c(0, 0.5),
      legend.title    = element_text(size = 7.5, face = "bold"),
      legend.text     = element_text(size = 7),
      legend.key.size = unit(0.45, "cm"),
      plot.background = element_rect(fill = "white", colour = NA),
      plot.margin     = margin(8, 8, 8, 8)
    )

  ggsave("figures/fig_spatial_exposure_qgis.pdf",
         p_map, width = 7, height = 4.8, device = "pdf")
  ggsave("figures/fig_spatial_exposure_qgis.png",
         p_map, width = 7, height = 4.8, dpi = 300)
  cat("  fig_spatial_exposure_qgis.pdf exporté.\n")
}


# =============================================================================
# PARTIE XIII -- EXPORT TABLEAUX ET FIGURES
# =============================================================================
# POURQUOI : Centraliser tous les exports en fin de script garantit que les
# fichiers sont générés après que tous les objets sont définis, et facilite
# la maintenance (un seul endroit à modifier pour changer un chemin de sortie).

# Labels communs pour etable()
coef_map <- c(
  "Z_hat_wl" = "$Z_{WL}$ (BC 2017)",
  "Z_hat_fm" = "$Z_{FM}$ (Fort McMurray 2016)",
  "Z_hat"    = "$\\hat{Z}$ (agrégé)"
)
stats_map <- c("n" = "Observations", "r2" = "$R^2$", "wr2" = "Within $R^2$")


# -----------------------------------------------------------------------------
# XIII.1 Table 1 -- Résultat principal (5 spécifications progressives)
# -----------------------------------------------------------------------------
tex_t1 <- etable(
  rf_ols, rf_fe_year, rf_fe_city, rf_baseline, rf_trends,
  dict        = c(coef_map, stats_map),
  depvar      = FALSE,
  fixef.group = list("City FE"          = "name_cancensus",
                     "Year FE"          = "REF_DATE",
                     "City trends"      = "name_cancensus[[REF_DATE]]"),
  fitstat     = ~ n + r2 + wr2 + rmse,
  digits      = 3, digits.stats = 3, se.below = TRUE,
  signif.code = c("***" = 0.01, "**" = 0.05, "*" = 0.10),
  tex         = TRUE, style.tex = style.tex("aer"),
  title       = "Effect of BC 2017 Wildfire Displacement on Canadian Rental Prices",
  label       = "tab:main_result",
  notes       = paste(
    "Dependent variable: log monthly rent, 2-bedroom apartments,",
    "buildings $\\geq$ 3 units (CMHC). $Z_{WL}$: pre-fire observed migration",
    "shares from Williams Lake $\\times$ Post$_{2017}$",
    "(Goldsmith-Pinkham, Sorkin \\& Swift 2020).",
    "Column (5) includes city-specific linear time trends.",
    "The results should be interpreted as reduced-form evidence consistent",
    "with a wildfire-induced housing demand shock.",
    "Standard errors clustered at the city level. Period: 2008--2024.",
    "*** $p<0.01$, ** $p<0.05$, * $p<0.10$."
  )
)

tex_t1 <- sub(
  "& \\(1\\)\\s*& \\(2\\)\\s*& \\(3\\)\\s*& \\(4\\)\\s*& \\(5\\)\\\\\\\\",
  "& (1) OLS & (2) Year FE & (3) City FE & (4) City+Year FE & (5) City trends \\\\\\\\",
  tex_t1
)
cat(tex_t1)
writeLines(tex_t1, "tables/table1_main.tex")


# -----------------------------------------------------------------------------
# XIII.2 Table 2 -- Robustesse des expositions (5 colonnes)
# -----------------------------------------------------------------------------
tex_t2 <- etable(
  rf_baseline, rf_fm, rf_combined, rf_agg, rf_unbal,
  dict        = c(coef_map, stats_map),
  depvar      = FALSE,
  fixef.group = list("City FE" = "name_cancensus", "Year FE" = "REF_DATE"),
  fitstat     = ~ n + r2 + wr2 + rmse,
  digits      = 3, digits.stats = 3, se.below = TRUE,
  signif.code = c("***" = 0.01, "**" = 0.05, "*" = 0.10),
  tex         = TRUE, style.tex = style.tex("aer"),
  title       = "Robustness: Alternative Exposure Specifications",
  label       = "tab:robustness_instrument",
  notes       = paste(
    "Dependent variable: log monthly rent, 2-bedroom apartments,",
    "buildings $\\geq$ 3 units (CMHC).",
    "$Z_{WL}$: Williams Lake exposure index (observed pre-fire shares $\\times$ Post$_{2017}$).",
    "$Z_{FM}$: Fort McMurray exposure index (gravity-predicted shares $\\times$ Post$_{2016}$).",
    "Standard errors clustered at the city level. Period: 2008--2024.",
    "Column (1) uses the main panel of 131 cities ($N=2{,}225$). Columns (2)--(4) use the panel augmented with Fort McMurray destination cities ($N=2{,}242$). Column (5) uses the unbalanced panel ($N=2{,}382$).",
    "*** $p<0.01$, ** $p<0.05$, * $p<0.10$."
  )
)

tex_t2 <- sub(
  "& \\(1\\)\\s*& \\(2\\)\\s*& \\(3\\)\\s*& \\(4\\)\\s*& \\(5\\)\\\\\\\\",
  "& (1) Baseline & (2) FM only & (3) WL+FM & (4) Aggregated & (5) Unbalanced \\\\\\\\",
  tex_t2
)
cat(tex_t2)
writeLines(tex_t2, "tables/table2_robustness.tex")


# -----------------------------------------------------------------------------
# XIII.2bis  Checks Angrist (COVID, pre-period trends, placebos distants)
# -----------------------------------------------------------------------------
# La table Province×Année est produite définitivement en PARTIE XV (avec
# le stopifnot N=2191 et les noms de variables consolidés). L'export ci-dessous
# couvre uniquement les checks supplémentaires (table_angrist_checks.tex).

# Export tables pour les nouveaux checks Angrist
etable(rf_baseline, rf_pre_covid, rf_pretrend_only,
       rf_placebo_sudbury, rf_placebo_quesnel,
       dict        = c(coef_map, stats_map),
       depvar      = FALSE,
       fixef.group = list("City FE" = "name_cancensus", "Year FE" = "REF_DATE"),
       fitstat     = ~ n + r2 + wr2 + rmse,
       digits      = 3, digits.stats = 3, se.below = TRUE,
       signif.code = c("***" = 0.01, "**" = 0.05, "*" = 0.10),
       headers     = c("Baseline", "Pre-COVID ($\\leq$ 2019)", "Pre-period trends",
                       "Placebo Sudbury", "Placebo Quesnel"),
       tex         = TRUE, style.tex = style.tex("aer"),
       title       = "Robustness: COVID Truncation, Pre-Period Trends, and Distant Placebo",
       label       = "tab:angrist_checks",
       notes       = paste(
         "Dependent variable: log monthly rent, two-bedroom apartments, buildings $\\geq 3$ units (CMHC).",
         "Column~(1): baseline TWFE (2008--2024).",
         "Column~(2): sample truncated to 2019 to exclude COVID-era mobility shocks.",
         "Column~(3): log rent detrended using city-specific linear trends estimated on",
         "the pre-treatment period (2008--2016) only; city and year FEs included on detrended outcome.",
         "Column~(4): placebo origin Greater Sudbury (Ontario), ~5,000 km from Williams Lake;",
         "pre-fire 2016/2017 migration shares $\\times$ Post$_{2017}$.",
         "Column~(5): placebo origin Quesnel (British Columbia), geographically proximate to WL.",
         "Standard errors clustered at the city level.",
         "*** $p<0.01$, ** $p<0.05$, * $p<0.10$."
       ),
       file    = "tables/table_angrist_checks.tex", replace = TRUE)


# -----------------------------------------------------------------------------
# XIII.3 Table de robustesse identification (7 colonnes)
# -----------------------------------------------------------------------------
tex_robust <- etable(
  rf_baseline, rf_trends, rf_bc_trend, rf_province_trends,
  rf_bc_only,  rf_unbal,  rf_combined, rf_t2018,
  dict        = c(coef_map, stats_map),
  depvar      = FALSE,
  fixef.group = list(
    "City FE"             = "name_cancensus",
    "Year FE"             = "REF_DATE",
    "City trends"         = "name_cancensus[[REF_DATE]]",
    "BC linear trend"     = "bc[[REF_DATE]]",
    "Province trends"     = "province[[REF_DATE]]"
  ),
  fitstat     = ~ n + r2 + wr2 + rmse,
  digits      = 3, digits.stats = 3, se.below = TRUE,
  signif.code = c("***" = 0.01, "**" = 0.05, "*" = 0.10),
  tex         = TRUE, style.tex = style.tex("aer"),
  title       = "Robustness: Identification Checks",
  label       = "tab:robustness_identification",
  notes       = paste(
    "Dependent variable: log monthly rent, 2-bedroom apartments,",
    "buildings $\\geq$ 3 units (CMHC). $Z_{WL}$: pre-fire migration",
    "shares $\\times$ Post$_{2017}$ (columns 1--7);",
    "$Z_{WL,18}$: pre-fire shares $\\times$ Post$_{2018}$ (column 8).",
    "Column (8) addresses the partial-year overlap between the July 2017",
    "wildfire and the October 2017 CMHC rent measurement.",
    "Standard errors clustered at the city level. Period: 2008--2024.",
    "*** $p<0.01$, ** $p<0.05$, * $p<0.10$."
  )
)

tex_robust <- sub(
  "& \\(1\\)\\s*& \\(2\\)\\s*& \\(3\\)\\s*& \\(4\\)\\s*& \\(5\\)\\s*& \\(6\\)\\s*& \\(7\\)\\s*& \\(8\\)\\\\\\\\",
  "& (1) Baseline & (2) City trends & (3) BC trend & (4) Prov.\\ trends & (5) BC only & (6) Unbalanced & (7) WL+FM & (8) $t\\\\geq 2018$ \\\\\\\\",
  tex_robust
)
cat(tex_robust)
writeLines(tex_robust, "tables/table_robustness_identification.tex")


# -----------------------------------------------------------------------------
# XIII.4 Table 3 -- Hétérogénéité 4×4, orientation recommandée
# -----------------------------------------------------------------------------
# Révision demandée : unit types en lignes, structures en colonnes.
# Chaque cellule contient β et (SE) sur deux lignes.
# La cellule baseline (2 bedrooms, Apt >= 3 units) est marquée avec †.
# La cellule Row >= 3 units / Bachelor est marquée "---" (non estimable).

all_models_by_structure <- list(
  "Apt $\\geq 3$ units"        = list(Bachelor = rf_s1_bachelor,
                                      `1 bedroom` = rf_s1_one,
                                      `2 bedrooms` = rf_s1_two,
                                      `3 bedrooms` = rf_s1_three),
  "Row \\& apt $\\geq 3$ units" = list(Bachelor = rf_s2_bachelor,
                                      `1 bedroom` = rf_s2_one,
                                      `2 bedrooms` = rf_s2_two,
                                      `3 bedrooms` = rf_s2_three),
  "Row $\\geq 3$ units"        = list(Bachelor = rf_s3_bachelor,
                                      `1 bedroom` = rf_s3_one,
                                      `2 bedrooms` = rf_s3_two,
                                      `3 bedrooms` = rf_s3_three),
  "Apt $\\geq 6$ units"        = list(Bachelor = rf_s4_bachelor,
                                      `1 bedroom` = rf_s4_one,
                                      `2 bedrooms` = rf_s4_two,
                                      `3 bedrooms` = rf_s4_three)
)

fmt_cell <- function(model, baseline = FALSE) {
  if (is.null(model)) return(list(coef = "\\multicolumn{1}{c}{---}", se = ""))
  b     <- stats::coef(model)["Z_hat_wl"]
  se_v  <- fixest::se(model)["Z_hat_wl"]
  p_v   <- fixest::pvalue(model)["Z_hat_wl"]
  stars <- ifelse(p_v < 0.01, "$^{***}$",
                  ifelse(p_v < 0.05, "$^{**}$",
                         ifelse(p_v < 0.10, "$^{*}$", "")))
  dag   <- ifelse(baseline, "$^{\\dag}$", "")
  list(
    coef = sprintf("%.3f%s%s", b, stars, dag),
    se   = sprintf("(%.3f)", se_v)
  )
}

unit_labels <- c("Bachelor", "1 bedroom", "2 bedrooms", "3 bedrooms")
structure_labels <- names(all_models_by_structure)

body <- ""
for (u in unit_labels) {
  coef_cells <- sapply(structure_labels, function(s) {
    is_baseline <- (u == "2 bedrooms" && s == "Apt $\\geq 3$ units")
    fmt_cell(all_models_by_structure[[s]][[u]], baseline = is_baseline)$coef
  })
  se_cells <- sapply(structure_labels, function(s) {
    fmt_cell(all_models_by_structure[[s]][[u]])$se
  })

  body <- paste0(body,
                 u, " & ",
                 paste(coef_cells, collapse = " & "), " \\\\\n",
                 " & ", paste(se_cells, collapse = " & "), " \\\\\n")
  if (u != tail(unit_labels, 1)) body <- paste0(body, "\\addlinespace\n")
}

tex_t3 <- paste0(
  "\\begin{table}[htbp]\n",
  "  \\caption{Heterogeneity by Unit Type and Housing Structure: ",
  "Post-2017 Williams Lake Exposure}\n",
  "  \\label{tab:heterogeneity_4x4}\n",
  "  \\bigskip\n",
  "  \\centering\n",
  "  \\begin{tabular}{lcccc}\n",
  "    \\toprule\n",
  "    & \\multicolumn{4}{c}{\\textit{Housing structure}} \\\\\n",
  "    \\cmidrule(lr){2-5}\n",
  "    Unit type & Apt $\\geq 3$ units & Row \\& apt $\\geq 3$ units & ",
  "Row $\\geq 3$ units & Apt $\\geq 6$ units \\\\\n",
  "    \\midrule\n",
  body,
  "    \\\\\n",
  "    City FE & \\multicolumn{4}{c}{$\\checkmark$} \\\\\n",
  "    Year FE & \\multicolumn{4}{c}{$\\checkmark$} \\\\\n",
  "    \\bottomrule\n",
  "  \\end{tabular}\n",
  "  \\par\\raggedright\\smallskip\n",
  "  \\footnotesize\n",
  "  Each cell reports the coefficient on $Z_{WL}$ from a separate regression ",
  "of log monthly rent on the post-2017 Williams Lake exposure index with ",
  "city and year fixed effects. Clustered standard errors at the city level ",
  "in parentheses. $^{\\dag}$~Baseline specification (2 bedrooms, apt ",
  "$\\geq$ 3 units), identical to Table~\\ref{tab:main_result}. ",
  "\\texttt{---}: cell not estimable due to insufficient observations after ",
  "balancing (row structures, bachelor units: 2 cities). ",
  "Row structures, 1-bedroom units: 14 cities, interpret with caution. ",
  "Period: 2008--2024. $^{***}$~$p<0.01$, $^{**}$~$p<0.05$, ",
  "$^{*}$~$p<0.10$.\n",
  "\\end{table}\n"
)

cat(tex_t3)
writeLines(tex_t3, "tables/table3_heterogeneity.tex")


# -----------------------------------------------------------------------------
# XIII.4bis  Table 4 -- Calibration économique (generated from table4 object)
# -----------------------------------------------------------------------------

fmt_pct  <- function(x) sprintf("%.1f", x)
fmt_dol  <- function(x) sprintf("\\$%d",  as.integer(round(x)))
fmt_shr  <- function(x) sprintf("%.3f", x)
fmt_disp <- function(x) format(round(x), big.mark = ",", trim = TRUE)

rows_t4 <- apply(table4, 1, function(r) {
  city <- r["name_cancensus"]
  paste0(
    "    ", city,
    " & ", fmt_shr(as.numeric(r["share_wl"])),
    " & ", fmt_disp(as.numeric(r["displaced"])),
    " & ", round(as.numeric(r["households"])),
    " & ", fmt_shr(as.numeric(r["pi_s"])),
    " & ", fmt_pct(as.numeric(r["delta_pct"])),
    " & ", fmt_dol(as.numeric(r["delta_dollar"])),
    " & ", fmt_pct(as.numeric(r["delta_pct_trend"])),
    " \\\\\n"
  )
})

tex_t4 <- paste0(
  "\\begin{table}[htbp]\n",
  "  \\caption{Economic Calibration of the Reduced-Form Magnitudes}\n",
  "  \\label{tab:calibration}\n",
  "  \\bigskip\n  \\centering\n",
  "  \\begin{tabular}{lrrrrrrr}\n",
  "    \\toprule\n",
  "    City & $s_{WL,d}$ & $11{,}000{\\times}s_{WL,d}$ & Households$^{a}$",
  " & $\\hat{\\pi}s_{WL,d}$ & $\\Delta$\\%rent$^{b}$ & $\\Delta$\\$rent$^{c}$",
  " & $\\Delta$\\%trend$^{d}$ \\\\\n",
  "    \\midrule\n",
  paste(rows_t4, collapse = ""),
  "    \\bottomrule\n",
  "  \\end{tabular}\n",
  "  \\par\\raggedright\\smallskip\n",
  "  \\footnotesize\n",
  "  Top-five most exposed cities in the main sample. ",
  "$^{a}$~Potential renter households under the illustrative scenario $(q,h)=(0.25,2)$: ",
  "$\\hat{H}_d = q\\times11{,}000\\times s_{WL,d}/h$. ",
  "$^{b}$~$\\Delta\\%=(\\exp(\\hat{\\pi}s_{WL,d})-1)\\times100$ using $\\hat{\\pi}=0.818$. ",
  "$^{c}$~$\\Delta\\$=\\Delta\\%\\times\\$848$, the 2015 CMHC average monthly rent used for scale. ",
  "$^{d}$~Same calculation using the city-trend coefficient $\\hat{\\pi}=0.594$ from column~(5) of ",
  "Table~\\ref{tab:main_result}. Household counts are scenarios, not estimates of realised ",
  "post-fire relocations.\n",
  "\\end{table}\n"
)

cat(tex_t4)
writeLines(tex_t4, "tables/table4_calibration.tex")


# -----------------------------------------------------------------------------
# XIII.4ter  Table 5 -- Market tightness (generated from table5 object)
# -----------------------------------------------------------------------------

fmt_vac  <- function(x) sprintf("%.1f\\%%", x * 100)
fmt_num1 <- function(x) sprintf("%.1f", x)

rows_t5 <- apply(table5, 1, function(r) {
  city <- r["name_cancensus"]
  paste0(
    "    ", city,
    " & ", format(as.integer(r["units_2br"]), big.mark = ",", trim = TRUE),
    " & ", fmt_vac(as.numeric(r["vacancy_rate"])),
    " & ", fmt_num1(as.numeric(r["implied_vacant"])),
    " & ", fmt_num1(as.numeric(r["H_q25"])),
    " & ", fmt_num1(as.numeric(r["HV_q10"])),
    " & ", fmt_num1(as.numeric(r["HV_q25"])),
    " & ", fmt_num1(as.numeric(r["HV_q50"])),
    " \\\\\n"
  )
})

tex_t5 <- paste0(
  "\\begin{table}[htbp]\n",
  "  \\caption{Potential Rental Demand Relative to Vacant Two-Bedroom Apartments in 2017}\n",
  "  \\label{tab:vacancy_calibration}\n",
  "  \\bigskip\n  \\centering\n",
  "  \\begin{tabular}{lrrrrrrr}\n",
  "    \\toprule\n",
  "    City & 2017 2BR units & 2017 2BR vacancy & Implied vacant units",
  " & Households$^{a}$ & H/V, $q{=}.10$ & H/V, $q{=}.25$ & H/V, $q{=}.50$ \\\\\n",
  "    \\midrule\n",
  paste(rows_t5, collapse = ""),
  "    \\bottomrule\n",
  "  \\end{tabular}\n",
  "  \\par\\raggedright\\smallskip\n",
  "  \\footnotesize\n",
  "  CMHC Rental Market Survey data refer to private apartment structures with at least three ",
  "rental units. Two-bedroom units and vacancy rates are taken from the 2017 British Columbia ",
  "Rental Market Report tables, measured in October 2017. Implied vacant units equal ",
  "two-bedroom units multiplied by the vacancy rate. ",
  "$^{a}$~Potential renter households under the intermediate scenario $(q,h)=(0.25,2)$ ",
  "from Table~\\ref{tab:calibration}. $H/V$ denotes potential renter households divided by ",
  "implied vacant two-bedroom units. This table is a market-tightness plausibility exercise, ",
  "not an estimate of realised migration or a test of the reduced-form identification strategy.\n",
  "\\end{table}\n"
)

cat(tex_t5)
writeLines(tex_t5, "tables/table5_tightness.tex")


# -----------------------------------------------------------------------------
# XIII.5 Statistiques descriptives
# -----------------------------------------------------------------------------
panel_for_stats <- panel_for_stats %>%
  mutate(Exposed = ifelse(share_wl > 0, "Exposed", "Non-Exposed"))

datasummary(
  (`Rent ($)` = rent) + (`Log Rent` = log_rent) + (`Exposure share (Share WL)` = share_wl) ~
    1 * (Mean + SD + Min + Max) + Exposed * (Mean + SD),
  data  = panel_for_stats,
  output = "tables/table1_summary_statistics.tex",
  title  = "Summary Statistics: Full Sample and by Exposure Status",
  notes  = paste("Notes: 'Exposed' refers to cities with a strictly positive",
                 "pre-fire migration share from Williams Lake.",
                 "The sample consists of 2,225 observations (131 cities, 2008--2024,",
                 "excluding 2 missing values for Leamington in 2023--2024)."),
  label  = "tab:summary_stats"
)

# Export Panel A/B au format thèse
tex_desc <- paste0(
  "\\begin{table}[htbp]\n",
  "  \\caption*{\\textbf{Table 1: Descriptive Statistics}}\n",
  "  \\centering\\scriptsize\n",
  "  \\begin{tabular}{lcccccc}\n",
  "    \\toprule\n",
  "    Variable & Obs. & Mean & Median & Std. Dev. & Min & Max \\\\\n",
  "    \\midrule\n",
  sprintf("    Monthly rent (\\$) & 2,225 & %.0f & %.0f & %.0f & %.0f & %.0f \\\\\n",
          stats_global$rent_mean, stats_global$rent_median,
          stats_global$rent_sd,   stats_global$rent_min, stats_global$rent_max),
  sprintf("    Log monthly rent & 2,225 & %.2f & %.2f & %.3f & %.2f & %.2f \\\\\n",
          stats_global$log_rent_mean, stats_global$log_rent_median,
          stats_global$log_rent_sd,   stats_global$log_rent_min, stats_global$log_rent_max),
  sprintf("    $s_{WL}$ (migration share) & 2,225 & %.3f & %.3f & %.3f & %.3f & %.3f \\\\\n",
          stats_global$share_wl_mean, stats_global$share_wl_median,
          stats_global$share_wl_sd,   stats_global$share_wl_min, stats_global$share_wl_max),
  "    \\addlinespace\n",
  sprintf("    \\multicolumn{7}{l}{\\textit{Within-city std. dev.: \\$%.0f \\quad Between-city std. dev.: \\$%.0f}} \\\\\n",
          stats_within$sd_within, stats_within$sd_between),
  "    \\multicolumn{7}{l}{\\textit{43 cities: $s_{WL,d}>0$ (exposed) \\quad 88 cities: $s_{WL,d}=0$ (non-exposed)}} \\\\\n",
  "    \\bottomrule\n",
  "  \\end{tabular}\n\n",
  "  \\smallskip\n",
  "  \\textit{Panel B: Average monthly rent (\\$) by exposure and period}\n",
  "  \\smallskip\n\n",
  "  \\begin{tabularx}{\\textwidth}{Xcccc}\n",
  "    \\toprule\n",
  "    & \\multicolumn{2}{c}{Non-exposed} & \\multicolumn{2}{c}{Exposed} \\\\\n",
  "    \\cmidrule(lr){2-3}\\cmidrule(lr){4-5}\n",
  "    & Pre-2017 & Post-2017 & Pre-2017 & Post-2017 \\\\\n",
  "    \\midrule\n",
  sprintf("    Mean rent (\\$) & %.0f & %.0f & %.0f & %.0f \\\\\n",
          stats_panel_b$rent_mean[stats_panel_b$exposed==FALSE & stats_panel_b$period=="Pre-2017"],
          stats_panel_b$rent_mean[stats_panel_b$exposed==FALSE & stats_panel_b$period=="Post-2017"],
          stats_panel_b$rent_mean[stats_panel_b$exposed==TRUE  & stats_panel_b$period=="Pre-2017"],
          stats_panel_b$rent_mean[stats_panel_b$exposed==TRUE  & stats_panel_b$period=="Post-2017"]),
  sprintf("    Std. dev. & (%.0f) & (%.0f) & (%.0f) & (%.0f) \\\\\n",
          stats_panel_b$rent_sd[stats_panel_b$exposed==FALSE & stats_panel_b$period=="Pre-2017"],
          stats_panel_b$rent_sd[stats_panel_b$exposed==FALSE & stats_panel_b$period=="Post-2017"],
          stats_panel_b$rent_sd[stats_panel_b$exposed==TRUE  & stats_panel_b$period=="Pre-2017"],
          stats_panel_b$rent_sd[stats_panel_b$exposed==TRUE  & stats_panel_b$period=="Post-2017"]),
  "    $\\Delta$ Post--Pre (\\$) & \\multicolumn{2}{c}{$+228$ $(+31\\%)$} & \\multicolumn{2}{c}{$+324$ $(+36\\%)$} \\\\\n",
  "    Cities & \\multicolumn{2}{c}{88} & \\multicolumn{2}{c}{43} \\\\\n",
  "    Observations & 792 & 704 & 387 & 342 \\\\\n",
  "    \\bottomrule\n",
  "  \\end{tabularx}\n",
  "\\end{table}\n"
)

writeLines(tex_desc, "tables/table_desc_stats.tex")
cat("table_desc_stats.tex exporte\n")


# -----------------------------------------------------------------------------
# XIII.6 Figures event studies et placebos
# -----------------------------------------------------------------------------
# Convention graphique commune à tous les event studies :
#   - Ruban gris très léger (alpha = 0.10) : IC 95% en arrière-plan
#   - Barres d'erreur verticales avec crochets (geom_errorbar) : IC 95% au premier plan
#   - Ligne + points : coefficients estimés
#   - Pointillé vertical à x = -0.5 : séparation pré/post traitement
#   - Tirets horizontaux à y = 0 : hypothèse nulle
#   - Points vides (shape = 1) : périodes non-estimées (zéros mécaniques)
#   - Points pleins (shape = 16) : périodes estimées
# Tous les IC sont à 95% (z = 1.96), erreurs-standard clusterisées par ville.

# Fonction utilitaire : extraire coefficients + SE + IC depuis un objet feols
extract_es <- function(model) {
  data.frame(
    estimate = stats::coef(model),
    se       = fixest::se(model)
  ) %>%
    mutate(
      ci_low  = estimate - 1.96 * se,
      ci_high = estimate + 1.96 * se,
      k       = as.numeric(gsub(".*::([-0-9]+).*", "\\1", rownames(.)))
    )
}

# Thème commun pour tous les event studies
theme_es <- theme_bw(base_size = 12) +
  theme(
    panel.grid.minor  = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.title        = element_text(face = "bold", size = 12),
    plot.subtitle     = element_text(color = "gray40", size = 10),
    plot.caption      = element_text(color = "gray40", size = 8, hjust = 0),
    axis.title        = element_text(size = 11)
  )

# Couleurs
COL_WL  <- "#2166ac"   # bleu -- Williams Lake
COL_PRE <- "#d6604d"   # rouge -- pre-trends test
COL_FM  <- "#4dac26"   # vert -- Fort McMurray


# ── Figure 1 : Pre-trends test (es_wl_correct) ──────────────────────────────
# Tous les coefficients sont ESTIMÉS (share_wl × year, ref = -1).
# On rajoute manuellement le point de référence k = -1 (normalisé à 0, SE = 0).

es_correct_plot <- bind_rows(
  extract_es(es_wl_correct),
  data.frame(estimate = 0, se = 0, ci_low = 0, ci_high = 0, k = -1)
) %>%
  arrange(k) %>%
  mutate(reference = (k == -1))   # TRUE uniquement pour le point de référence

ref_y_wl <- min(es_correct_plot$ci_low, na.rm = TRUE) * 0.85

p_es_wl_pre <- ggplot(es_correct_plot, aes(x = k, y = estimate)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.5) +
  geom_vline(xintercept = -0.5, linetype = "dotted", color = "gray60", linewidth = 0.5) +
  geom_ribbon(aes(ymin = ci_low, ymax = ci_high),
              alpha = 0.10, fill = COL_PRE, na.rm = TRUE) +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high),
                width = 0.25, linewidth = 0.65, color = COL_PRE, na.rm = TRUE) +
  geom_line(color = COL_PRE, linewidth = 0.75, na.rm = TRUE) +
  geom_point(aes(shape = reference), color = COL_PRE, size = 2.8, fill = "white") +
  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 22), guide = "none") +
  annotate("text", x = -1, y = ref_y_wl,
           label = "Reference (k = -1)", size = 2.8, color = "gray40") +
  scale_x_continuous(breaks = min(es_correct_plot$k):max(es_correct_plot$k)) +
  labs(
    x = "Years relative to 2017 wildfire",
    y = "Estimated coefficient (log monthly rent)"
  ) +
  theme_es +
  theme(legend.position = "none")
ggsave("figures/event_study_wl_pretrends.pdf", p_es_wl_pre, width = 8, height = 5)
ggsave("figures/event_study_wl_pretrends.png", p_es_wl_pre, width = 8, height = 5, dpi = 300)


# ── Figure 2 : Post-treatment dynamics (es_wl_post_dynamics) ────────────────
# Avant 2017 : Z_hat_wl = 0 par construction → coefficients mécaniquement nuls,
# pas d'IC estimable. Ces points sont affichés en cercles vides pour signaler
# qu'ils ne constituent PAS un test de parallel trends.

es_post_plot <- bind_rows(
  # Périodes pré-traitement : zéros mécaniques, pas d'IC
  data.frame(
    estimate  = 0,
    se        = NA_real_,
    ci_low    = NA_real_,
    ci_high   = NA_real_,
    k         = -9:-2,
    estimated = FALSE
  ),
  # Référence k = -1 (normalisée à 0)
  data.frame(estimate = 0, se = 0, ci_low = 0, ci_high = 0,
             k = -1, estimated = FALSE),
  # Périodes post-traitement : coefficients estimés avec IC
  extract_es(es_wl_post_dynamics) %>% mutate(estimated = TRUE)
) %>% arrange(k)

p_es_wl_post <- ggplot(es_post_plot, aes(x = k, y = estimate)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.5) +
  geom_vline(xintercept = -0.5, linetype = "dotted", color = "gray60", linewidth = 0.5) +
  # IC 95% -- ruban léger (uniquement post-traitement car na.rm = TRUE)
  geom_ribbon(aes(ymin = ci_low, ymax = ci_high),
              alpha = 0.10, fill = COL_WL, na.rm = TRUE) +
  # IC 95% -- crochets (uniquement sur les points estimés)
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high),
                width = 0.25, linewidth = 0.65, color = COL_WL, na.rm = TRUE) +
  # Ligne de connexion (na.rm pour éviter les ruptures sur NA)
  geom_line(color = COL_WL, linewidth = 0.75, na.rm = TRUE) +
  # Points : cercles vides = zéros mécaniques, cercles pleins = estimés
  geom_point(aes(shape = estimated, fill = estimated),
             color = COL_WL, size = 2.8) +
  scale_shape_manual(values = c(`FALSE` = 21, `TRUE` = 16),
                     labels = c("Mechanical zero (no CI)", "Estimated"),
                     name   = NULL) +
  scale_fill_manual(values  = c(`FALSE` = "white", `TRUE` = COL_WL),
                    guide   = "none") +
  scale_x_continuous(breaks = -9:7) +
  labs(
    title    = "Post-Treatment Dynamics: BC 2017 Wildfire (Williams Lake)",
    subtitle = "Pre-period: mechanically zero (Z₂⁰¹⁷ = 0 ∀ k < 0) -- see pre-trends test for parallel trends assessment",
    x        = "Years relative to wildfire (0 = 2017)",
    y        = expression(paste("Coefficient on ", Z[WL])),
    caption  = "95% CI on post-treatment periods. Clustered SE at city level. Open circles: no CI (mechanical zeros)."
  ) +
  theme_es +
  theme(legend.position = "bottom",
        legend.key.size  = unit(0.5, "cm"))
ggsave("figures/event_study_wl_post_dynamics.pdf", p_es_wl_post, width = 8, height = 5)
ggsave("figures/event_study_wl_post_dynamics.png", p_es_wl_post, width = 8, height = 5, dpi = 300)


# ── Figure 3a : Pre-trends test Fort McMurray ───────────────────────────────
es_fm_correct_plot <- bind_rows(
  extract_es(es_fm_correct),
  data.frame(estimate = 0, se = 0, ci_low = 0, ci_high = 0, k = -1)
) %>%
  arrange(k) %>%
  mutate(reference = (k == -1))

ref_y_fm <- min(es_fm_correct_plot$ci_low, na.rm = TRUE) * 0.85

p_es_fm_pre <- ggplot(es_fm_correct_plot, aes(x = k, y = estimate)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.5) +
  geom_vline(xintercept = -0.5, linetype = "dotted", color = "gray60", linewidth = 0.5) +
  geom_ribbon(aes(ymin = ci_low, ymax = ci_high),
              alpha = 0.10, fill = COL_FM, na.rm = TRUE) +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high),
                width = 0.25, linewidth = 0.65, color = COL_FM, na.rm = TRUE) +
  geom_line(color = COL_FM, linewidth = 0.75, na.rm = TRUE) +
  geom_point(aes(shape = reference), color = COL_FM, size = 2.8, fill = "white") +
  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 22), guide = "none") +
  annotate("text", x = -1, y = ref_y_fm,
           label = "Reference (k = -1)", size = 2.8, color = "gray40") +
  scale_x_continuous(breaks = min(es_fm_correct_plot$k):max(es_fm_correct_plot$k)) +
  labs(
    x = "Years relative to 2016 wildfire",
    y = "Estimated coefficient (log monthly rent)"
  ) +
  theme_es +
  theme(legend.position = "none")
ggsave("figures/event_study_fm_pretrends.pdf", p_es_fm_pre, width = 8, height = 5)
ggsave("figures/event_study_fm_pretrends.png", p_es_fm_pre, width = 8, height = 5, dpi = 300)


# ── Figure 3b : Post-treatment dynamics Fort McMurray ───────────────────────
es_fm_plot <- bind_rows(
  data.frame(
    estimate  = 0,
    se        = NA_real_,
    ci_low    = NA_real_,
    ci_high   = NA_real_,
    k         = -8:-2,
    estimated = FALSE
  ),
  data.frame(estimate = 0, se = 0, ci_low = 0, ci_high = 0,
             k = -1, estimated = FALSE),
  extract_es(es_fm) %>% mutate(estimated = TRUE)
) %>% arrange(k)

p_es_fm_post <- ggplot(es_fm_plot, aes(x = k, y = estimate)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.5) +
  geom_vline(xintercept = -0.5, linetype = "dotted", color = "gray60", linewidth = 0.5) +
  geom_ribbon(aes(ymin = ci_low, ymax = ci_high),
              alpha = 0.10, fill = COL_FM, na.rm = TRUE) +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high),
                width = 0.25, linewidth = 0.65, color = COL_FM, na.rm = TRUE) +
  geom_line(color = COL_FM, linewidth = 0.75, na.rm = TRUE) +
  geom_point(aes(shape = estimated, fill = estimated),
             color = COL_FM, size = 2.8) +
  scale_shape_manual(values = c(`FALSE` = 21, `TRUE` = 16),
                     labels = c("Mechanical zero", "Estimated"),
                     name   = NULL) +
  scale_fill_manual(values = c(`FALSE` = "white", `TRUE` = COL_FM),
                    guide  = "none") +
  scale_x_continuous(breaks = -8:8) +
  labs(
    title    = "Post-Treatment Dynamics: Fort McMurray 2016 (Wood Buffalo)",
    subtitle = "Pre-period: mechanically zero (Z_FM = 0 ∀ k < 0) -- no detectable post-treatment effect",
    x        = "Years relative to wildfire (0 = 2016)",
    y        = expression(paste("Coefficient on ", hat(Z)[FM])),
    caption  = "95% CI on post-treatment periods. Clustered SE at city level."
  ) +
  theme_es +
  theme(legend.position = "bottom",
        legend.key.size  = unit(0.5, "cm"))
ggsave("figures/event_study_fm_clean.pdf", p_es_fm_post, width = 8, height = 5)
ggsave("figures/event_study_fm_clean.png", p_es_fm_post, width = 8, height = 5, dpi = 300)


# ── Figure 4 : Placebo tests (4 spécifications) ─────────────────────────────
# Pour les placebos, on utilise geom_errorbar seul (pas de ruban) :
# les 7 pseudo-années sont des points discrets, le ruban n'est pas adapté.
# Chaque panneau = une spécification. L'axe y est libre (scales = "free_y").

placebo_plot_data <- placebo_all %>%
  mutate(
    ci_low  = estimate - 1.96 * std.error,
    ci_high = estimate + 1.96 * std.error,
    # Renommer les spécifications pour l'affichage dans les panneaux
    spec_label = factor(spec,
      levels = c("Baseline (City+Year FE)", "BC linear trend",
                 "City linear trends",     "BC only"),
      labels = c("(1) Baseline\n(City + Year fixed effects)",
                 "(2) British Columbia linear trend",
                 "(3) City-specific trends",
                 "(4) British Columbia only"))
  )

# Palettes de couleurs pour les 4 panneaux
pal_placebo <- c("#2166ac", "#d6604d", "#4dac26", "#762a83")

p_placebo <- ggplot(placebo_plot_data,
                    aes(x = placebo_year, y = estimate, color = spec_label)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray45", linewidth = 0.5) +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high),
                width = 0.25, linewidth = 0.65) +
  geom_line(linewidth = 0.75) +
  geom_point(size = 2.8, shape = 16) +
  facet_wrap(~spec_label, scales = "free_y", ncol = 2) +
  scale_x_continuous(breaks = 2010:2016,
                     labels = c("2010","2011","2012","2013","2014","2015","2016")) +
  scale_color_manual(values = pal_placebo, guide = "none") +
  labs(
    x = "Pseudo-treatment year",
    y = "Coefficient"
  ) +
  theme_es +
  theme(strip.text = element_text(face = "bold", size = 10),
        strip.background = element_rect(fill = "gray95"))
ggsave("figures/placebo_tests.pdf", p_placebo, width = 10, height = 6)
ggsave("figures/placebo_tests.png", p_placebo, width = 10, height = 6, dpi = 300)

cat("Tous les tableaux et figures ont ete exportes.\n")


# -----------------------------------------------------------------------------
# XIII.7  Scatter : s_{WL,d} vs post-2017 rent differential (reduced-form)
# -----------------------------------------------------------------------------
# Visual argument for section 4.3.
# Strategy: remove year FEs from log_rent (demean by year across all cities),
# then compute city-level mean of demeaned rent pre-2017 and post-2017.
# The y-axis is the difference (post minus pre), the x-axis is share_wl.
# Only the 43 exposed cities (share_wl > 0) are plotted.

year_fe_scatter <- panel_balanced_wl %>%
  group_by(REF_DATE) %>%
  summarise(year_mean = mean(log_rent, na.rm = TRUE), .groups = "drop")

df_scatter <- panel_balanced_wl %>%
  left_join(year_fe_scatter, by = "REF_DATE") %>%
  mutate(log_rent_dm = log_rent - year_mean) %>%
  group_by(name_cancensus, share_wl) %>%
  summarise(
    pre  = mean(log_rent_dm[REF_DATE < 2017],  na.rm = TRUE),
    post = mean(log_rent_dm[REF_DATE >= 2017], na.rm = TRUE),
    diff = post - pre,
    .groups = "drop"
  ) %>%
  filter(share_wl > 0)

labels_scatter <- c("Kamloops", "Vancouver", "Prince George", "Kelowna", "Quesnel")

p_scatter <- ggplot(df_scatter, aes(x = share_wl, y = diff)) +
  geom_point(color = "steelblue", size = 2.5, alpha = 0.8) +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE, color = "black",
              fill = "grey80", linewidth = 0.8) +
  geom_label_repel(
    data        = filter(df_scatter, name_cancensus %in% labels_scatter),
    aes(label   = name_cancensus),
    size        = 3.2,
    box.padding = 0.4,
    min.segment.length = 0
  ) +
  labs(
    x = expression(paste("Pre-fire migration share  ", italic(s)[WL*","*d])),
    y = "Post-2017 minus pre-2017 log rent\n(year fixed effects removed)"
  ) +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank())

ggsave("figures/scatter_reducedform.pdf", plot = p_scatter, width = 7, height = 5)
cat("scatter_reducedform.pdf exported.\n")


# -----------------------------------------------------------------------------
# XIII.8 Table Hard Checks -- Vancouver policy confound + Trimming
# -----------------------------------------------------------------------------
coef_map_hc <- c(
  "Z_hat_wl" = "$Z_{WL}$ (BC 2017)",
  "van_post"  = "Vancouver $\\times$ Post$_{2017}$"
)

tex_hc <- etable(
  rf_baseline, rf_drop_van, rf_van_post, rf_drop_kam_van,
  dict        = c(coef_map_hc, stats_map),
  depvar      = FALSE,
  fixef.group = list("City FE" = "name_cancensus", "Year FE" = "REF_DATE"),
  fitstat     = ~ n + r2 + wr2 + rmse,
  digits      = 3, digits.stats = 3, se.below = TRUE,
  signif.code = c("***" = 0.01, "**" = 0.05, "*" = 0.10),
  tex         = TRUE, style.tex = style.tex("aer"),
  title       = "Robustness: Hard Checks on Policy Confounding and Exposure Specificity",
  label       = "tab:hard_checks",
  notes       = paste(
    "Dependent variable: log monthly rent, 2-bedroom apartments,",
    "buildings $\\geq$ 3 units (CMHC). $Z_{WL}$: pre-fire migration shares from",
    "Williams Lake $\\times$ Post$_{2017}$. Baseline replicates column (4) of",
    "Table~\\ref{tab:main_result}.",
    "Column (2): Vancouver excluded from sample (drops BC housing policy confound).",
    "Column (3): adds a Vancouver $\\times$ Post$_{2017}$ dummy absorbing any",
    "post-2017 shock specific to Vancouver (BC Foreign Buyers Tax 2016,",
    "Speculation \\& Vacancy Tax 2018), without removing the city.",
    "Column (4): both Kamloops ($s_{WL}=0.212$) and Vancouver ($s_{WL}=0.174$)",
    "excluded --- the two highest-exposure cities representing",
    "$\\approx$40\\% of total identification weight.",
    "Standard errors clustered at the city level. Period: 2008--2024.",
    "*** $p<0.01$, ** $p<0.05$, * $p<0.10$."
  )
)

tex_hc <- sub(
  "& \\(1\\)\\s*& \\(2\\)\\s*& \\(3\\)\\s*& \\(4\\)\\\\\\\\",
  paste("& (1) Baseline & (2) Drop Vancouver",
        "& (3) Van. $\\\\times$ Post & (4) Drop Kam.+Van. \\\\\\\\"),
  tex_hc
)
writeLines(tex_hc, "tables/table_hard_checks.tex")
cat("table_hard_checks.tex exported.\n")


# -----------------------------------------------------------------------------
# XIII.9 Table Soft Checks -- Gravity WL + Placebo Quesnel
# -----------------------------------------------------------------------------
coef_map_sc <- c(
  "Z_hat_wl"    = "$Z_{WL}$ observed shares",
  "Z_grav"      = "$Z_{WL}^{grav}$ gravity OLS",
  "Z_grav_ppml" = "$Z_{WL}^{grav,ppml}$ gravity PPML",
  "Z_quesnel"   = "$Z_{Quesnel}$ placebo shares"
)

tex_sc <- etable(
  rf_baseline, rf_gravity_wl, rf_gravity_ppml_wl, rf_placebo_quesnel,
  dict        = c(coef_map_sc, stats_map),
  depvar      = FALSE,
  fixef.group = list("City FE" = "name_cancensus", "Year FE" = "REF_DATE"),
  fitstat     = ~ n + r2 + wr2 + rmse,
  digits      = 3, digits.stats = 3, se.below = TRUE,
  signif.code = c("***" = 0.01, "**" = 0.05, "*" = 0.10),
  tex         = TRUE, style.tex = style.tex("aer"),
  title       = "Robustness: Soft Checks on Exposure Definition and Specificity",
  label       = "tab:soft_checks",
  notes       = paste(
    "Dependent variable: log monthly rent, 2-bedroom apartments,",
    "buildings $\\geq$ 3 units (CMHC).",
    "Column (1): baseline with observed pre-fire migration shares.",
    "Column (2): exposure index replaced by OLS log-linear gravity-predicted shares",
    "$\\hat{s}_{WL,d}^{grav} \\propto \\text{Pop}_d / \\text{dist}(WL,d)^{\\hat{\\beta}_1}$.",
    "Column (3): same shares predicted by Poisson Pseudo-Maximum Likelihood",
    "\\citep{santossilva2006}, which handles zero flows and corrects heteroskedasticity",
    "in multiplicative gravity models.",
    "Column (4): placebo exposure index from Quesnel (BC) pre-fire migration",
    "shares $\\times$ Post$_{2017}$ --- Quesnel did not burn in 2017.",
    "Standard errors clustered at the city level. Period: 2008--2024.",
    "*** $p<0.01$, ** $p<0.05$, * $p<0.10$."
  )
)

tex_sc <- sub(
  "& \\(1\\)\\s*& \\(2\\)\\s*& \\(3\\)\\s*& \\(4\\)\\\\\\\\",
  "& (1) Baseline & (2) Gravity OLS & (3) Gravity PPML & (4) Placebo Quesnel \\\\\\\\",
  tex_sc
)
writeLines(tex_sc, "tables/table_soft_checks.tex")
cat("table_soft_checks.tex exported.\n")

# =============================================================================
# PARTIE XIV -- INFÉRENCE HONNÊTE (RAMBACHAN-ROTH)   [MIS EN COMMENTAIRE]
# Pour réactiver : sélectionner tout le bloc et décommenter (Ctrl+Shift+C).
# Temps d'exécution estimé : 3-5 minutes (bisection = 8 appels HonestDiD).
# =============================================================================
#
# suppressPackageStartupMessages(library(HonestDiD))
#
# # --- XIV.1  Extraire betahat et sigma depuis es_wl_correct -------------------
# all_coefs_rr <- coef(es_wl_correct)
# all_vcov_rr  <- vcov(es_wl_correct)
#
# coef_k_rr    <- as.numeric(gsub(".*::(-?[0-9]+):.*", "\\1", names(all_coefs_rr)))
# pre_order_rr  <- which(coef_k_rr < -1)[order(coef_k_rr[coef_k_rr < -1])]
# post_order_rr <- which(coef_k_rr >= 0)[order(coef_k_rr[coef_k_rr >= 0])]
# numPre_rr  <- length(pre_order_rr)
# numPost_rr <- length(post_order_rr)
#
# betahat_rr <- c(all_coefs_rr[pre_order_rr], all_coefs_rr[post_order_rr])
# sigma_rr   <- all_vcov_rr[c(pre_order_rr, post_order_rr),
#                            c(pre_order_rr, post_order_rr)]
# l_vec_rr <- c(0, rep(1/7, 7))
#
# # --- XIV.2  M-bar=0 : estimateur ponctuel et IC classique --------------------
# sigma_post_rr <- sigma_rr[(numPre_rr+1):(numPre_rr+numPost_rr),
#                            (numPre_rr+1):(numPre_rr+numPost_rr)]
# est_m0_rr <- as.numeric(l_vec_rr %*% betahat_rr[(numPre_rr+1):(numPre_rr+numPost_rr)])
# se_m0_rr  <- sqrt(as.numeric(t(l_vec_rr) %*% sigma_post_rr %*% l_vec_rr))
# lb_m0_rr  <- est_m0_rr - 1.96 * se_m0_rr
# ub_m0_rr  <- est_m0_rr + 1.96 * se_m0_rr
#
# # --- XIV.3  HonestDiD C-LF pour M-bar in {0.5, 1.0, 1.5} -------------------
# rr_res_rr <- HonestDiD::createSensitivityResults_relativeMagnitudes(
#   betahat        = betahat_rr,   sigma          = sigma_rr,
#   numPrePeriods  = numPre_rr,    numPostPeriods = numPost_rr,
#   Mbarvec        = c(0.5, 1.0, 1.5),
#   l_vec          = l_vec_rr,     alpha          = 0.05)
#
# rr_table_rr <- rbind(
#   data.frame(Mbar = 0.0, lb = lb_m0_rr,        ub = ub_m0_rr),
#   data.frame(Mbar = 0.5, lb = rr_res_rr$lb[1], ub = rr_res_rr$ub[1]),
#   data.frame(Mbar = 1.0, lb = rr_res_rr$lb[2], ub = rr_res_rr$ub[2]),
#   data.frame(Mbar = 1.5, lb = rr_res_rr$lb[3], ub = rr_res_rr$ub[3]))
#
# # --- XIV.4  Valeur de rupture M* par bisection (5 itérations) ---------------
# call_hd_rr <- function(Mbar) {
#   res <- HonestDiD::createSensitivityResults_relativeMagnitudes(
#     betahat        = betahat_rr,  sigma          = sigma_rr,
#     numPrePeriods  = numPre_rr,   numPostPeriods = numPost_rr,
#     Mbarvec        = Mbar,        l_vec          = l_vec_rr,  alpha = 0.05)
#   res$lb
# }
# lo_rr <- 0; hi_rr <- 0.5
# for (i_rr in 1:5) {
#   mid_rr <- (lo_rr + hi_rr) / 2
#   if (call_hd_rr(mid_rr) > 0) lo_rr <- mid_rr else hi_rr <- mid_rr
# }
# mstar_rr <- hi_rr
#
# # --- XIV.5  Export table_rambachan_roth.tex ----------------------------------
# fmt_rr <- function(Mbar, lb, ub) {
#   if (abs(Mbar) < 1e-9)
#     sprintf("      $%.1f$ & $%.3f$ & $[%.3f,\\ %.3f]$ \\\\", Mbar, (lb+ub)/2, lb, ub)
#   else
#     sprintf("      $%.1f$ & ---     & $[%.3f,\\ %.3f]$ \\\\", Mbar, lb, ub)
# }
# rows_rr <- mapply(fmt_rr, rr_table_rr$Mbar, rr_table_rr$lb, rr_table_rr$ub)
# tex_rr <- c(
#   "\\begin{table}[htbp]",
#   "   \\caption{\\label{tab:rambachan_roth} Sensitivity to Parallel-Trend Violations: Rambachan--Roth Confidence Sets}",
#   "   \\bigskip", "   \\centering\\small", "   \\begin{threeparttable}",
#   "      \\begin{tabular}{lcc}", "         \\toprule",
#   "         $\\bar{M}$ & Point estimate & 95\\% confidence set \\\\",
#   "         \\midrule", rows_rr, "         \\bottomrule",
#   "      \\end{tabular}", "      \\begin{tablenotes}\\footnotesize",
#   "         \\item \\textit{Notes:} C-LF method of \\citealt{rambachan2023}.",
#   paste0("           Breakdown value $\\bar{M}^{*}=",
#          sprintf("%.2f", mstar_rr), "$. Panel: 131 cities, 2008--2024."),
#   "      \\end{tablenotes}", "   \\end{threeparttable}", "\\end{table}")
# writeLines(tex_rr, "tables/table_rambachan_roth.tex")
# saveRDS(list(rr_table = rr_table_rr, mstar = mstar_rr,
#              betahat = betahat_rr, sigma = sigma_rr),
#         "processed/rr_results.rds")
# cat("Rambachan-Roth exporté.\n")

# =============================================================================
# PARTIE XV -- ABSORPTION RÉGIONALE (Section 7.4) + TIMING t>=2018 (Section 7.5)
# =============================================================================
# Trois spécifications qui testent si l'effet survit à l'absorption provinciale :
#   (2) Province×year FE  : λ_{p(d)t}  -- contrôle provincial maximal
#   (3) Région×year FE    : λ_{r(d)t}  -- Ouest vs Est
#   (4) Within-West       : sous-éch. BC/AB/SK/MB, year FE communs
# Plus : t>=2018 (section 7.5) -- redéfinition de l'indicateur post-traitement.
# Output : table_prov_year.tex
# =============================================================================

# -- XV.1  Sous-panel Ouest ----------------------------------------------------
# region et Z_hat_wl18 sont définis en PARTIE IX.2 (F) et disponibles ici.

panel_west_xv <- panel_balanced_wl %>% filter(region == "West")

# -- XV.2  Régressions ---------------------------------------------------------
rf_xv_base       <- feols(log_rent ~ Z_hat_wl   | name_cancensus + REF_DATE,
                          data = panel_balanced_wl, cluster = ~name_cancensus)

rf_xv_prov_year  <- feols(log_rent ~ Z_hat_wl   | name_cancensus + province^REF_DATE,
                          data = panel_balanced_wl, cluster = ~name_cancensus)

rf_xv_region_year <- feols(log_rent ~ Z_hat_wl  | name_cancensus + region^REF_DATE,
                           data = panel_balanced_wl, cluster = ~name_cancensus)

rf_xv_within_west <- feols(log_rent ~ Z_hat_wl  | name_cancensus + REF_DATE,
                           data = panel_west_xv,   cluster = ~name_cancensus)

rf_xv_t2018      <- feols(log_rent ~ Z_hat_wl18 | name_cancensus + REF_DATE,
                          data = panel_balanced_wl, cluster = ~name_cancensus)

# Vérification explicite : Province×Year FE doit droper exactement 34 singletons
# (cohérent avec la note de bas de la Table prov_year du manuscrit, N=2,191).
stopifnot(nobs(rf_xv_prov_year) == 2191L)
cat(sprintf("Province-Year FE: %d singleton obs dropped (expected 34)\n",
            2225L - nobs(rf_xv_prov_year)))

cat("\n=== PARTIE XV : ABSORPTION RÉGIONALE ===\n")
print(etable(rf_xv_base, rf_xv_prov_year, rf_xv_region_year, rf_xv_within_west,
             headers = c("(1) Baseline", "(2) Prov×Year", "(3) Region×Year",
                         "(4) Within-West"),
             digits  = 3))

cat(sprintf("\nTiming t>=2018 : %.4f (se=%.4f)\n",
            coef(rf_xv_t2018)["Z_hat_wl18"], se(rf_xv_t2018)["Z_hat_wl18"]))

# -- XV.3  Export table_prov_year.tex ------------------------------------------
coef_map_xv  <- c(Z_hat_wl   = "$Z_{WL}$")
stats_map_xv <- c(n = "Observations", r2 = "R$^{2}$", wr2 = "Within R$^{2}$")

tex_py <- etable(
  rf_xv_base, rf_xv_prov_year, rf_xv_region_year, rf_xv_within_west,
  dict        = c(coef_map_xv, stats_map_xv),
  depvar      = FALSE,
  fixef.group = list(
    "City FE"          = "name_cancensus",
    "Year FE"          = "REF_DATE",
    "Province×Year FE" = "province^REF_DATE",
    "Region×Year FE"   = "region^REF_DATE"
  ),
  fitstat     = ~ n + r2 + wr2 + rmse,
  digits      = 3, digits.stats = 3, se.below = TRUE,
  signif.code = c("***" = 0.01, "**" = 0.05, "*" = 0.10),
  tex         = TRUE, style.tex = style.tex("aer"),
  title       = "Robustness: Province-Year and Regional Absorption",
  label       = "tab:prov_year",
  notes       = paste(
    "Dependent variable: log monthly rent, 2-bedroom apartments,",
    "buildings $\\geq 3$ units (CMHC). $Z_{WL}$: shift-share exposure",
    "(pre-fire migration shares $\\times$ Post$_{2017}$). Column (1):",
    "baseline with city and year fixed effects. Column (2): province-by-year",
    "fixed effects replace year FE, retaining only within-province cross-city",
    "variation. Column (3): region-by-year fixed effects (West = BC, AB, SK, MB",
    "vs.~Rest of Canada). Column (4): within-Western subsample",
    "($N_{cities}=51$, four Western provinces). Province-year FE removes",
    "34 singleton observations. Standard errors clustered at the city level.",
    "*** $p<0.01$, ** $p<0.05$, * $p<0.10$."
  )
)

tex_py <- sub(
  "& \\(1\\)\\s*& \\(2\\)\\s*& \\(3\\)\\s*& \\(4\\)\\\\\\\\",
  "& (1) Baseline & (2) Prov.$\\\\times$Year & (3) Region$\\\\times$Year & (4) Within-West \\\\\\\\",
  tex_py
)

# Clarify the year-effects row in the exported LaTeX table:
# common year FE are present in columns (1) and (4);
# in columns (2) and (3), they are absorbed by province-by-year
# and region-by-year fixed effects, respectively.
tex_py <- gsub(
  "(?m)^\\s*Year FE[^\\n]*\\\\\\\\\\s*$",
  "Year FE & $\\\\checkmark$ &  &  & $\\\\checkmark$ \\\\\\\\",
  tex_py,
  perl = TRUE
)

writeLines(tex_py, "tables/table_prov_year.tex")
cat("table_prov_year.tex exporté.\n")

# =============================================================================
# PARTIE XVI -- VALIDATION DU CANAL MIGRATOIRE (Section 5)
# =============================================================================
# Objectif :
#   (i) montrer que les parts pré-feu s_WL,d prédisent les arrivées observées
#       depuis Williams Lake après la vague 2016/2017 ;
#   (ii) distinguer ce canal d'une simple connectivité générale de la BC via
#       un placebo Quesnel ;
#   (iii) produire les figures et le tableau directement insérables dans Overleaf.
#
# Convention temporelle :
#   wave_rel = -1 : 2016/2017, strictement pré-feu ;
#   wave_rel =  0 : 2017/2018, première vague post-feu ;
#   wave_rel =  1 : 2018/2019 ;
#   wave_rel =  2 : 2019/2020 ;
#   wave_rel =  3 : 2020/2021.
#
# Remarque importante :
#   Le tableau Statistique Canada peut enregistrer les flux nuls soit comme
#   zéro explicite, soit par absence de ligne selon la version téléchargée.
#   Le code complète le panel destination × vague et remplace les absences par 0.
#   Vérifier le diagnostic `missing_before_zero_*` avant rédaction finale.

waves_channel <- c("2016/2017", "2017/2018", "2018/2019",
                   "2019/2020", "2020/2021")

wave_dictionary <- tibble::tibble(
  REF_DATE_mig = waves_channel,
  wave_rel     = c(-1L, 0L, 1L, 2L, 3L),
  post_flow    = as.integer(c(FALSE, TRUE, TRUE, TRUE, TRUE)),
  rent_year_concurrent = c(2016L, 2017L, 2018L, 2019L, 2020L),
  rent_year_lagged     = c(2017L, 2018L, 2019L, 2020L, 2021L)
)

balanced_destinations <- panel_balanced_wl %>%
  distinct(name_cancensus, share_wl0)

build_origin_arrival_panel <- function(origin_regex, share_df, share_name,
                                       flow_name = "arrivals") {
  raw_flows <- data1 %>%
    filter(str_detect(GEO, origin_regex),
           REF_DATE %in% waves_channel,
           !grepl("Area outside", `Geography.of.destination`)) %>%
    mutate(
      name_cancensus = `Geography.of.destination` %>%
        str_remove("\\s*\\(C[MA]+\\)") %>%
        str_remove(",.*") %>%
        str_trim(),
      REF_DATE_mig = REF_DATE
    ) %>%
    group_by(name_cancensus, REF_DATE_mig) %>%
    summarise(arrivals_raw = sum(VALUE, na.rm = TRUE), .groups = "drop")

  grid <- tidyr::expand_grid(
    name_cancensus = balanced_destinations$name_cancensus,
    REF_DATE_mig   = waves_channel
  )

  missing_before_zero <- grid %>%
    left_join(raw_flows, by = c("name_cancensus", "REF_DATE_mig")) %>%
    summarise(n_missing_before_zero = sum(is.na(arrivals_raw))) %>%
    pull(n_missing_before_zero)

  panel <- grid %>%
    left_join(raw_flows, by = c("name_cancensus", "REF_DATE_mig")) %>%
    left_join(wave_dictionary, by = "REF_DATE_mig") %>%
    left_join(share_df, by = "name_cancensus") %>%
    mutate(
      arrivals_raw = replace_na(arrivals_raw, 0),
      !!flow_name  := arrivals_raw,
      !!share_name := replace_na(.data[[share_name]], 0)
    ) %>%
    select(-arrivals_raw)

  list(panel = panel, missing_before_zero = missing_before_zero)
}

# --- XVI.1 Williams Lake: arrivées et amplification post-feu ------------------
shares_wl_channel <- shares_wl %>%
  rename(share_wl_channel = share_wl)

wl_channel_obj <- build_origin_arrival_panel(
  origin_regex = "^Williams Lake \\(CA\\), British Columbia$",
  share_df     = shares_wl_channel,
  share_name   = "share_wl_channel",
  flow_name    = "arrivals_wl"
)

panel_mig_wl <- wl_channel_obj$panel %>%
  mutate(
    share_wl_channel = ifelse(name_cancensus == "Williams Lake", 0,
                              share_wl_channel),
    Z_channel_wl = share_wl_channel * post_flow
  )

missing_before_zero_wl <- wl_channel_obj$missing_before_zero
cat(sprintf("Channel validation WL: rows completed as zero before estimation = %d\n",
            missing_before_zero_wl))

# Average post-vs-pre relation in arrivals
fs_wl_avg <- feols(
  arrivals_wl ~ Z_channel_wl | name_cancensus + REF_DATE_mig,
  data    = panel_mig_wl,
  cluster = ~name_cancensus
)

# Dynamic relation: one coefficient per post wave, relative to 2016/2017
fs_wl_event <- feols(
  arrivals_wl ~ i(wave_rel, share_wl_channel, ref = -1) |
    name_cancensus + REF_DATE_mig,
  data    = panel_mig_wl,
  cluster = ~name_cancensus
)

# Count-data robustness: PPML with two-way fixed effects.
# This is supportive only; the paper's main channel diagnostic remains the
# OLS count specification for transparent interpretation.
fs_wl_ppml <- fepois(
  arrivals_wl ~ Z_channel_wl | name_cancensus + REF_DATE_mig,
  data    = panel_mig_wl,
  cluster = ~name_cancensus
)

# Descriptive table by exposure groups: zero-share cities + terciles among positives
wl_group_map <- panel_mig_wl %>%
  distinct(name_cancensus, share_wl_channel) %>%
  arrange(share_wl_channel) %>%
  mutate(exposure_group = "No pre-fire WL link")

wl_pos_groups <- wl_group_map %>%
  filter(share_wl_channel > 0) %>%
  mutate(exposure_group = paste0("Positive-share tercile ",
                                 dplyr::ntile(share_wl_channel, 3)))

wl_group_map <- wl_group_map %>%
  filter(share_wl_channel <= 0) %>%
  bind_rows(wl_pos_groups) %>%
  mutate(exposure_group = factor(
    exposure_group,
    levels = c("No pre-fire WL link",
               "Positive-share tercile 1",
               "Positive-share tercile 2",
               "Positive-share tercile 3")
  ))

channel_desc_wl <- panel_mig_wl %>%
  left_join(wl_group_map, by = c("name_cancensus", "share_wl_channel")) %>%
  group_by(exposure_group, post_flow) %>%
  summarise(
    destinations = n_distinct(name_cancensus),
    total_arrivals = sum(arrivals_wl, na.rm = TRUE),
    mean_arrivals_per_city_wave = mean(arrivals_wl, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(period = ifelse(post_flow == 1, "Post-fire waves", "Pre-fire wave")) %>%
  select(exposure_group, period, destinations, total_arrivals,
         mean_arrivals_per_city_wave)

print(channel_desc_wl)

# --- XVI.2 Placebo Quesnel: même exercice ------------------------------------
shares_quesnel_channel <- data1 %>%
  filter(GEO == "Quesnel (CA), British Columbia",
         REF_DATE == "2016/2017",
         !grepl("Area outside", `Geography.of.destination`)) %>%
  mutate(
    total_out = sum(VALUE, na.rm = TRUE),
    share_quesnel_channel = VALUE / total_out,
    name_cancensus = `Geography.of.destination` %>%
      str_remove("\\s*\\(C[MA]+\\)") %>%
      str_remove(",.*") %>%
      str_trim()
  ) %>%
  select(name_cancensus, share_quesnel_channel) %>%
  distinct(name_cancensus, .keep_all = TRUE)

q_channel_obj <- build_origin_arrival_panel(
  origin_regex = "^Quesnel \\(CA\\), British Columbia$",
  share_df     = shares_quesnel_channel,
  share_name   = "share_quesnel_channel",
  flow_name    = "arrivals_quesnel"
)

panel_mig_q <- q_channel_obj$panel %>%
  mutate(
    share_quesnel_channel = ifelse(name_cancensus == "Quesnel", 0,
                                   share_quesnel_channel),
    Z_channel_q = share_quesnel_channel * post_flow
  )

missing_before_zero_q <- q_channel_obj$missing_before_zero
cat(sprintf("Channel validation Quesnel: rows completed as zero before estimation = %d\n",
            missing_before_zero_q))

fs_q_avg <- feols(
  arrivals_quesnel ~ Z_channel_q | name_cancensus + REF_DATE_mig,
  data    = panel_mig_q,
  cluster = ~name_cancensus
)

fs_q_event <- feols(
  arrivals_quesnel ~ i(wave_rel, share_quesnel_channel, ref = -1) |
    name_cancensus + REF_DATE_mig,
  data    = panel_mig_q,
  cluster = ~name_cancensus
)

# --- XVI.3 Exports: figures du canal -----------------------------------------
# Note: iplot() produit un graphique base R (pas ggplot2) → pdf()/dev.off()
# est ici le seul dispositif d'export correct (ggsave() ne capture pas base R).
pdf("figures/event_study_migration_wl.pdf", width = 7.2, height = 4.8)
iplot(fs_wl_event,
      main = "",
      xlab = "Migration wave relative to 2016/2017",
      ylab = "Coefficient on pre-fire Williams Lake share",
      ci.width = 0.25)
dev.off()

pdf("figures/event_study_migration_quesnel.pdf", width = 7.2, height = 4.8)
iplot(fs_q_event,
      main = "",
      xlab = "Migration wave relative to 2016/2017",
      ylab = "Coefficient on pre-fire Quesnel share",
      ci.width = 0.25)
dev.off()

channel_plot_data <- panel_mig_wl %>%
  left_join(wl_group_map, by = c("name_cancensus", "share_wl_channel")) %>%
  group_by(exposure_group, REF_DATE_mig) %>%
  summarise(mean_arrivals = mean(arrivals_wl, na.rm = TRUE), .groups = "drop")

p_channel_groups <- ggplot(channel_plot_data,
                           aes(x = REF_DATE_mig, y = mean_arrivals,
                               group = exposure_group,
                               linetype = exposure_group,
                               shape = exposure_group)) +
  geom_line(linewidth = 0.55) +
  geom_point(size = 1.8) +
  geom_vline(xintercept = "2017/2018", linetype = "dotted",
             linewidth = 0.4) +
  labs(x = "Migration wave",
       y = "Mean arrivals from Williams Lake per destination",
       linetype = NULL,
       shape = NULL) +
  thesis_theme +
  theme(axis.text.x = element_text(angle = 25, hjust = 1))

ggsave("figures/fig_migration_flows_by_exposure.pdf", p_channel_groups,
       width = 7.4, height = 4.6)

# --- XVI.4 Export: table de validation du canal ------------------------------
coef_map_channel <- c(
  "Z_channel_wl" = "$s_{WL,d} \\times \\mathrm{PostFlow}$",
  "Z_channel_q"  = "$s_{Q,d} \\times \\mathrm{PostFlow}$"
)

tex_channel <- etable(
  fs_wl_avg, fs_q_avg,
  dict        = coef_map_channel,
  depvar      = FALSE,
  fixef.group = list(
    "Destination FE" = "name_cancensus",
    "Wave FE"        = "REF_DATE_mig"
  ),
  fitstat     = ~ n + r2 + wr2 + rmse,
  digits      = 3, digits.stats = 3, se.below = TRUE,
  signif.code = c("***" = 0.01, "**" = 0.05, "*" = 0.10),
  tex         = TRUE, style.tex = style.tex("aer"),
  title       = "Validation of the Migration Channel: Post-Fire Arrivals",
  label       = "tab:channel_validation",
  notes       = paste(
    "Dependent variable: observed annual migration arrivals from the origin city",
    "to each destination in the balanced rental-market sample.",
    "Column (1) uses Williams Lake as origin. Column (2) uses Quesnel as a",
    "placebo origin. The pre-fire wave is 2016/2017; post-fire waves are",
    "2017/2018 through 2020/2021. Missing destination-wave rows are completed",
    "as zeros after the diagnostic check reported in the replication code.",
    "Standard errors clustered at the destination-city level."
  )
)

tex_channel <- sub(
  "& \\(1\\)\\s*& \\(2\\)\\\\\\\\",
  "& (1) Williams Lake & (2) Placebo Quesnel \\\\\\\\",
  tex_channel
)

writeLines(tex_channel, "tables/table_channel_validation.tex")
cat("table_channel_validation.tex exported.\n")


# =============================================================================
# PARTIE XVII -- 2SLS EXPLORATOIRE SUR FENÊTRE LIMITÉE DE FLUX MIGRATOIRES
# =============================================================================
# Objectif :
#   Relier l'exposition prédéterminée à un flux réalisé, puis les flux réalisés
#   aux loyers, sur la fenêtre où les données de migrations existent.
#
# Lecture :
#   Cet exercice n'est PAS le cœur du papier. Il est informatif sur le mécanisme,
#   mais il hérite des mêmes préoccupations de tendances régionales que le
#   reduced form et il souffre d'un problème d'alignement temporel entre les
#   vagues de migration (juillet-juin) et le loyer CMHC mesuré à l'automne.
#
# Deux alignements sont reportés :
#   (A) concurrent : wave 2017/2018 -> rent year 2017 ;
#   (B) lagged     : wave 2017/2018 -> rent year 2018.
# Le second est préférable conceptuellement pour que la vague de migration
# précède davantage la mesure de loyer, mais il ne résout pas parfaitement
# la granularité annuelle.

build_iv_panel <- function(alignment = c("concurrent", "lagged")) {
  alignment <- match.arg(alignment)

  year_col <- ifelse(alignment == "concurrent",
                     "rent_year_concurrent", "rent_year_lagged")

  mig_iv <- panel_mig_wl %>%
    mutate(REF_DATE = .data[[year_col]]) %>%
    select(name_cancensus, REF_DATE, arrivals_wl, Z_channel_wl,
           share_wl_channel, post_flow, wave_rel, REF_DATE_mig)

  panel_balanced_wl %>%
    select(name_cancensus, REF_DATE, log_rent, province, bc, share_wl0) %>%
    inner_join(mig_iv, by = c("name_cancensus", "REF_DATE")) %>%
    filter(!is.na(log_rent))
}

panel_iv_concurrent <- build_iv_panel("concurrent")
panel_iv_lagged     <- build_iv_panel("lagged")

iv_wl_concurrent <- feols(
  log_rent ~ 1 | name_cancensus + REF_DATE | arrivals_wl ~ Z_channel_wl,
  data    = panel_iv_concurrent,
  cluster = ~name_cancensus
)

iv_wl_lagged <- feols(
  log_rent ~ 1 | name_cancensus + REF_DATE | arrivals_wl ~ Z_channel_wl,
  data    = panel_iv_lagged,
  cluster = ~name_cancensus
)

# Reduced form and first stage on the same two IV windows.
rf_iv_concurrent <- feols(
  log_rent ~ Z_channel_wl | name_cancensus + REF_DATE,
  data    = panel_iv_concurrent,
  cluster = ~name_cancensus
)

fs_iv_concurrent <- feols(
  arrivals_wl ~ Z_channel_wl | name_cancensus + REF_DATE,
  data    = panel_iv_concurrent,
  cluster = ~name_cancensus
)

rf_iv_lagged <- feols(
  log_rent ~ Z_channel_wl | name_cancensus + REF_DATE,
  data    = panel_iv_lagged,
  cluster = ~name_cancensus
)

fs_iv_lagged <- feols(
  arrivals_wl ~ Z_channel_wl | name_cancensus + REF_DATE,
  data    = panel_iv_lagged,
  cluster = ~name_cancensus
)

cat("\n=== Exploratory FE-2SLS, concurrent alignment ===\n")
print(summary(iv_wl_concurrent, stage = 1:2))
cat("\n=== Exploratory FE-2SLS, lagged alignment ===\n")
print(summary(iv_wl_lagged, stage = 1:2))

# Table compacte pour l'annexe : first stage, reduced form, 2SLS
coef_map_iv <- c(
  "Z_channel_wl" = "$s_{WL,d} \\times \\mathrm{PostFlow}$",
  "fit_arrivals_wl" = "Observed arrivals from Williams Lake"
)

tex_iv <- etable(
  fs_iv_concurrent, rf_iv_concurrent, iv_wl_concurrent,
  fs_iv_lagged,     rf_iv_lagged,     iv_wl_lagged,
  dict        = coef_map_iv,
  depvar      = FALSE,
  fitstat     = ~ n + r2 + wr2 + rmse + ivf1,
  digits      = 3, digits.stats = 3, se.below = TRUE,
  signif.code = c("***" = 0.01, "**" = 0.05, "*" = 0.10),
  tex         = TRUE, style.tex = style.tex("aer"),
  title       = "Exploratory Limited-Window Migration First Stage and FE-2SLS",
  label       = "tab:iv_limited_window",
  notes       = paste(
    "Columns (1)--(3) align migration wave 2017/2018 with rent year 2017.",
    "Columns (4)--(6) align the same wave with rent year 2018.",
    "Columns (1) and (4) are first stages; columns (2) and (5) are reduced forms;",
    "columns (3) and (6) are FE-2SLS estimates in which observed Williams Lake",
    "arrivals are instrumented by pre-fire Williams Lake shares interacted with",
    "the post-flow indicator. All specifications include destination-city and",
    "rent-year fixed effects. Standard errors clustered at the destination-city level.",
    "These estimates are exploratory and are not used as the paper's main causal",
    "estimand because temporal alignment is imperfect and the pre-trend concerns",
    "documented in the main design remain relevant."
  )
)

writeLines(tex_iv, "tables/table_iv_limited_window.tex")
cat("table_iv_limited_window.tex exported.\n")


# =============================================================================
# COEFFICIENT PLOTS (remplacent les tables pour publication)
# Les tables restent exportées au-dessus pour consultation interne.
# =============================================================================

# -----------------------------------------------------------------------------
# FIGURE : Leave-one-out coefficient plot (remplace tab:robustness_loo)
# -----------------------------------------------------------------------------
loo_plot_data <- data.frame(
  label  = c("Baseline\n(no drop)", "Drop top 1\n(Kamloops)",
             "Drop top 2", "Drop top 3", "Drop top 5", "Drop top 10"),
  coef   = c(coef(rf_baseline)["Z_hat_wl"],
             vapply(loo_models, function(m) unname(coef(m)["Z_hat_wl"]),
                    numeric(1))),
  se     = c(sqrt(vcov(rf_baseline)["Z_hat_wl", "Z_hat_wl"]),
             vapply(loo_models,
                    function(m) sqrt(vcov(m)["Z_hat_wl", "Z_hat_wl"]),
                    numeric(1))),
  is_baseline = c(TRUE, FALSE, FALSE, FALSE, FALSE, FALSE)
)
loo_plot_data$label <- factor(loo_plot_data$label,
                               levels = rev(loo_plot_data$label))

p_loo <- ggplot(loo_plot_data,
                aes(x = coef, y = label, colour = is_baseline)) +
  geom_vline(xintercept = coef(rf_baseline)["Z_hat_wl"],
             linetype = "dashed", colour = "grey50", linewidth = 0.4) +
  geom_vline(xintercept = 0, colour = "black", linewidth = 0.3) +
  geom_errorbarh(aes(xmin = coef - 1.96 * se,
                     xmax = coef + 1.96 * se),
                 height = 0.2, linewidth = 0.5) +
  geom_point(size = 2.5) +
  scale_colour_manual(values = c("TRUE" = "#000000", "FALSE" = "#444444"),
                      guide = "none") +
  labs(x = "Coefficient on WL exposure index",
       y = NULL,
       title = NULL) +
  theme_classic(base_size = 11) +
  theme(panel.grid.major.x = element_line(colour = "grey90", linewidth = 0.3),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank())

ggsave("figures/fig_loo_coefplot.pdf", p_loo,
       width = 6, height = 3.5, device = pdf)
ggsave("figures/fig_loo_coefplot.png", p_loo,
       width = 6, height = 3.5, dpi = 300)
cat("fig_loo_coefplot exported.\n")


# -----------------------------------------------------------------------------
# FIGURE : Heterogeneity 4x4 coefficient plot (remplace tab:heterogeneity_4x4)
# -----------------------------------------------------------------------------
unit_map <- c("Bachelor" = "Bachelor",
              "1 bedroom" = "1 bedroom",
              "2 bedrooms" = "2 bedrooms",
              "3 bedrooms" = "3 bedrooms")

struct_map <- c(
  "Apt $\\geq 3$ units"         = "Apt >=3 units",
  "Row \\& apt $\\geq 3$ units" = "Row & apt >=3",
  "Row $\\geq 3$ units"         = "Row >=3 units",
  "Apt $\\geq 6$ units"         = "Apt >=6 units"
)

het_rows <- list()
for (s in names(all_models_by_structure)) {
  for (u in names(unit_map)) {
    m <- all_models_by_structure[[s]][[u]]
    if (is.null(m)) next
    b  <- unname(coef(m)["Z_hat_wl"])
    se_v <- sqrt(vcov(m)["Z_hat_wl", "Z_hat_wl"])
    het_rows[[length(het_rows) + 1]] <- data.frame(
      structure = struct_map[s],
      unit      = unit_map[u],
      coef      = b,
      se        = se_v,
      is_baseline = (s == "Apt $\\geq 3$ units" && u == "2 bedrooms"),
      stringsAsFactors = FALSE
    )
  }
}
het_data <- do.call(rbind, het_rows)
het_data$unit      <- factor(het_data$unit,
                              levels = rev(c("Bachelor", "1 bedroom",
                                             "2 bedrooms", "3 bedrooms")))
het_data$structure <- factor(het_data$structure,
                              levels = struct_map)

p_het <- ggplot(het_data,
                aes(x = coef, y = unit,
                    colour = is_baseline, shape = is_baseline)) +
  geom_vline(xintercept = coef(rf_baseline)["Z_hat_wl"],
             linetype = "dashed", colour = "grey50", linewidth = 0.4) +
  geom_vline(xintercept = 0, colour = "black", linewidth = 0.3) +
  geom_errorbarh(aes(xmin = coef - 1.96 * se,
                     xmax = coef + 1.96 * se),
                 height = 0.25, linewidth = 0.5) +
  geom_point(size = 2.5) +
  scale_colour_manual(values = c("TRUE" = "#000000", "FALSE" = "#444444"),
                      guide = "none") +
  scale_shape_manual(values = c("TRUE" = 18, "FALSE" = 16), guide = "none") +
  facet_wrap(~ structure, nrow = 1) +
  labs(x = "Coefficient on WL exposure index", y = NULL) +
  theme_classic(base_size = 10) +
  theme(panel.grid.major.x = element_line(colour = "grey90", linewidth = 0.3),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 9))

ggsave("figures/fig_heterogeneity_coefplot.pdf", p_het,
       width = 9, height = 3.2, device = "pdf")
ggsave("figures/fig_heterogeneity_coefplot.png", p_het,
       width = 9, height = 3.2, dpi = 300)
cat("fig_heterogeneity_coefplot exported.\n")

# =============================================================================
# PARTIE XVIII -- DIAGNOSTIC POIDS TWFE (de Chaisemartin & D'Haultfoeuille 2020)
# =============================================================================
#
# CONTEXTE :
# Le post LinkedIn de Clément de Chaisemartin (mai 2026) attire l'attention sur
# les poids potentiellement négatifs dans les estimateurs TWFE avec traitement
# CONTINU (extension Callaway-Goodman-Bacon-Sant'Anna / CGBSA, 2024).
#
# QUESTION :
# Le TWFE de la spécification baseline (rf_baseline) assigne-t-il des poids
# négatifs à certaines ATTs ? Si oui, quelle proportion ?
#
# PERTINENCE POUR CE DESIGN :
# - Z_hat_wl = share_wl × Post2017 est un traitement CONTINU en intensité
# - Choc SIMULTANÉ : toutes les villes passent de D=0 à D=share_wl en 2017
# - La préoccupation Goodman-Bacon (staggered timing) est adressée dans le
#   footnote §4. La préoccupation CGBSA (hétérogénéité des ATTs par intensité)
#   est structurellement limitée ici car l'intensité est fixe dans le temps
#   (= 0 avant 2017, = share_wl après) → pas de variation intra-ville de D.
#
# HYPOTHÈSE A PRIORI :
# Fraction de poids négatifs faible ou nulle, pour deux raisons :
# (1) Choc simultané → pas d'"early adopters" qui servent de contrôle pour
#     les "late adopters" — le canal principal des poids négatifs de GB (2021).
# (2) Intensité fixe post-traitement → la variation exploitée est purement
#     cross-sectionelle (Kamloops vs Vancouver vs non-exposé) ; dans ce cas,
#     les poids du TWFE sont proportionnels à (D_d - D_bar)² × T_post, qui
#     sont toujours positifs.
# =============================================================================

library(TwoWayFEWeights)

# Données : retirer les 2 NA Leamington 2023-2024
df_twfe_diag <- panel_balanced_wl %>%
  filter(!is.na(log_rent), !is.na(Z_hat_wl))

cat("\n=============================================================\n")
cat("PARTIE XVIII — DIAGNOSTIC POIDS TWFE (TwoWayFEWeights)\n")
cat("=============================================================\n")
cat(sprintf("Panel : %d villes × %d années = %d observations\n",
            length(unique(df_twfe_diag$name_cancensus)),
            length(unique(df_twfe_diag$REF_DATE)),
            nrow(df_twfe_diag)))
cat(sprintf("Z_hat_wl : min = %.4f, max = %.4f, moy = %.4f\n",
            min(df_twfe_diag$Z_hat_wl),
            max(df_twfe_diag$Z_hat_wl),
            mean(df_twfe_diag$Z_hat_wl)))
cat(sprintf("Villes exposées (Z_hat_wl > 0 au moins 1 année) : %d\n",
            df_twfe_diag %>%
              filter(Z_hat_wl > 0) %>%
              distinct(name_cancensus) %>%
              nrow()))

# --- Diagnostic principal ---
# type = "feTR" : TWFE avec "treatment restriction" (spécification standard)
# C'est l'équivalent de la décomposition de de Chaisemartin & D'Haultfoeuille
# (2020, AER) pour traitements continus.
cat("\n--- Calcul des poids TWFE (type = feTR) ---\n")
twfe_w <- tryCatch(
  twowayfeweights(
    data = df_twfe_diag,
    Y    = "log_rent",
    G    = "name_cancensus",
    T    = "REF_DATE",
    D    = "Z_hat_wl",
    type = "feTR"
  ),
  error = function(e) {
    cat("Erreur :", conditionMessage(e), "\n"); NULL
  }
)

# --- Extraction et interprétation ---
if (!is.null(twfe_w)) {

  # Affichage du résumé natif du package
  cat("\n--- Résumé du package TwoWayFEWeights ---\n")
  print(summary(twfe_w))

  # Calcul manuel à partir du vecteur de poids
  w <- twfe_w$weights
  n_tot <- length(w)
  n_pos <- sum(w > 0, na.rm = TRUE)
  n_neg <- sum(w < 0, na.rm = TRUE)
  sum_pos <- sum(w[w > 0], na.rm = TRUE)
  sum_neg <- sum(w[w < 0], na.rm = TRUE)   # valeur absolue = part "perturbatrice"

  cat("\n=== RÉSULTATS NUMÉRIQUES ===\n")
  cat(sprintf("Nombre d'ATTs estimés         : %d\n",    n_tot))
  cat(sprintf("ATTs avec poids positif       : %d (%.1f%%)\n",
              n_pos, 100 * n_pos / n_tot))
  cat(sprintf("ATTs avec poids négatif       : %d (%.1f%%)\n",
              n_neg, 100 * n_neg / n_tot))
  cat(sprintf("Somme des poids positifs      : %.4f\n",  sum_pos))
  cat(sprintf("Somme des poids négatifs (val abs) : %.4f\n", abs(sum_neg)))
  cat(sprintf("Coefficient TWFE (rf_baseline): %.4f\n",
              coef(rf_baseline)["Z_hat_wl"]))
  cat(sprintf("Minimum de poids             : %.6f\n",  min(w, na.rm = TRUE)))
  cat(sprintf("Maximum de poids             : %.6f\n",  max(w, na.rm = TRUE)))

  # Seuils d'interprétation (de Chaisemartin & D'Haultfoeuille 2020 recommandent
  # de comparer |sum(poids négatifs)| à la magnitude du coefficient)
  ratio_neg <- abs(sum_neg) / abs(coef(rf_baseline)["Z_hat_wl"])

  cat("\n=== INTERPRÉTATION ===\n")
  if (n_neg == 0) {
    cat(">> AUCUN poids négatif détecté.\n")
    cat("   Le TWFE est une moyenne pondérée d'ATTs avec poids 100% positifs.\n")
    cat("   La critique CGBSA (traitements continus) ne s'applique pas ici.\n")
    cat("   Explication : choc simultané + intensité fixe post-2017 → les poids\n")
    cat("   sont proportionnels à (D_d - D_bar)^2 × 1[post], tous ≥ 0.\n")
    cat("   APPORT POUR LA SOUTENANCE : confirme le footnote §4 EMPIRIQUEMENT.\n")
    cat("   Valeur : forte (1 ligne en réponse à la question sur CGBSA).\n")
  } else if (ratio_neg < 0.05) {
    cat(sprintf(">> Poids négatifs TRÈS FAIBLES (|sum_neg| / |beta| = %.3f < 5%%).\n",
                ratio_neg))
    cat("   Négligeable économiquement. Le TWFE reste interprétable.\n")
    cat("   APPORT POUR LA SOUTENANCE : robustesse confirmée.\n")
  } else if (ratio_neg < 0.20) {
    cat(sprintf(">> Poids négatifs MODÉRÉS (|sum_neg| / |beta| = %.3f).\n", ratio_neg))
    cat("   La non-linéarité des effets par intensité d'exposition est possible.\n")
    cat("   Recommandation : mentionner comme limitation, non centrale.\n")
  } else {
    cat(sprintf(">> Poids négatifs IMPORTANTS (|sum_neg| / |beta| = %.3f).\n", ratio_neg))
    cat("   Le TWFE mixe des effets positifs et négatifs de façon problématique.\n")
    cat("   Recommandation : citer CGBSA et envisager un estimateur robuste.\n")
  }
}

cat("\n=== CONCLUSION GÉNÉRALE ===\n")
cat("Ce diagnostic répond à la question soulevée par le post LinkedIn de\n")
cat("de Chaisemartin (mai 2026) sur les poids TWFE avec traitement continu.\n")
cat("Il est DISTINCT de la critique Goodman-Bacon (staggered adoption) déjà\n")
cat("couverte dans le footnote §4 du manuscrit.\n")
cat("Valeur pour le mémoire : footnote additionnelle de 2 lignes, si les\n")
cat("résultats confirment l'absence de poids négatifs.\n")
cat("Valeur pour la soutenance : réponse directe si le jury cite CGBSA.\n")
cat("Valeur pour un referee AER/RES : robustesse standard attendue (à ajouter\n")
cat("dans une soumission future avec le package twowayfeweights cité).\n")
