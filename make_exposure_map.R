library(ggplot2)
library(sf)
library(dplyr)
library(scales)
library(ggrepel)
library(rnaturalearth)
library(rnaturalearthdata)

setwd("~/Desktop/MASTER THESIS")

# ── Data ──────────────────────────────────────────────────────────────────────
pts <- read.csv("cities_exposure_arcgis.csv", stringsAsFactors = FALSE)

exposed     <- pts |> filter(exposure_label == "Exposed")
non_exposed <- pts |> filter(exposure_label == "Non-exposed")
origin      <- pts |> filter(exposure_label == "Origin (Williams Lake)")

# ── Basemap ───────────────────────────────────────────────────────────────────
# Try to get provinces; fall back to country outline if unavailable
canada_prov <- tryCatch(
  sf::read_sf("https://raw.githubusercontent.com/nvkelso/natural-earth-vector/master/geojson/ne_50m_admin_1_states_provinces.geojson") |>
    filter(admin == "Canada"),
  error = function(e) NULL
)

canada <- ne_countries(country = "Canada",   scale = "medium", returnclass = "sf")
usa    <- ne_countries(country = "United States of America",
                       scale = "medium", returnclass = "sf")

# ── Palette ───────────────────────────────────────────────────────────────────
# Four exposure bins (in percent)
exposed <- exposed |>
  mutate(bin = cut(share_wl_pct,
                   breaks = c(0, 1, 3, 7, 100),
                   labels = c("0–1%", "1–3%", "3–7%", ">7%"),
                   include.lowest = TRUE))

pal <- c("0–1%" = "#fdcc8a", "1–3%" = "#fc8d59",
         "3–7%" = "#e34a33", ">7%"  = "#b30000")

# ── Plot ──────────────────────────────────────────────────────────────────────
p <- ggplot() +

  # Canada base fill
  geom_sf(data = canada, fill = "#f5f5f0", colour = NA) +

  # Provincial borders if available
  { if (!is.null(canada_prov))
      geom_sf(data = canada_prov, fill = NA, colour = "grey65", linewidth = 0.22)
    else
      geom_sf(data = canada, fill = NA, colour = "grey60", linewidth = 0.3) } +

  # US border (light)
  geom_sf(data = usa, fill = "#ebebeb", colour = "grey70", linewidth = 0.2) +

  # Non-exposed cities
  geom_point(data = non_exposed,
             aes(x = longitude, y = latitude),
             shape = 21, size = 1.6,
             fill = "grey80", colour = "grey50", stroke = 0.3,
             alpha = 0.7) +

  # Exposed cities — size + colour by share
  geom_point(data = exposed,
             aes(x = longitude, y = latitude,
                 fill = bin, size = share_wl_pct),
             shape = 21, colour = "white", stroke = 0.35, alpha = 0.9) +

  # Williams Lake origin — distinctive marker
  geom_point(data = origin,
             aes(x = longitude, y = latitude),
             shape = 23, size = 5, fill = "#2c7bb6",
             colour = "white", stroke = 0.8) +

  geom_text_repel(data = origin,
                  aes(x = longitude, y = latitude, label = name_cancensus),
                  size = 2.5, fontface = "bold", colour = "#2c7bb6",
                  nudge_y = 1.2, segment.colour = "grey50",
                  segment.size = 0.3) +

  # Scales
  scale_fill_manual(values = pal, name = "Pre-fire share\n(% of WL out-migrants)") +
  scale_size_continuous(range = c(2, 9), guide = "none") +

  # Bounding box: Canada west + centre (where the action is)
  coord_sf(xlim = c(-140, -52), ylim = c(42, 62), expand = FALSE) +

  # Labels
  labs(
    title    = "Pre-fire Williams Lake Out-Migration Shares",
    subtitle = "Share of 2016/2017 Williams Lake out-migrants received by each city",
    caption  = paste0("Notes: Filled circles are the 43 cities with positive pre-fire exposure (s > 0). ",
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
    legend.position = c(0.02, 0.25),
    legend.justification = c(0, 0.5),
    legend.title    = element_text(size = 7.5, face = "bold"),
    legend.text     = element_text(size = 7),
    legend.key.size = unit(0.45, "cm"),
    plot.background = element_rect(fill = "white", colour = NA),
    plot.margin     = margin(8, 8, 8, 8)
  )

# ── Export ────────────────────────────────────────────────────────────────────
ggsave("figures/fig_spatial_exposure_qgis.pdf",
       p, width = 7, height = 4.8, device = "pdf")

ggsave("figures/fig_spatial_exposure_qgis.png",
       p, width = 7, height = 4.8, dpi = 300)

cat("Map saved to figures/fig_spatial_exposure_qgis.pdf\n")
