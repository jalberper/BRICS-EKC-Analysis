################################################################################
# EKC + Instituciones (WDI + V-Dem) para BRICS (1990–2024)
# PRE-ARBITRAJE Q1: Cointegración + Dependencia cruzada + ECM-MG/CS-ARDL + Robustez
#
# Variables institucionales:
#   - estabilidad_politica_vdem
#   - gobernanza_vdem
#
# MAIN SPECS (recomendación Q1):
#   (A) CS-ECM-MG (ARDL(1,1,1) por país + cross-sectional averages) para M2/M3
#   (B) ECM-MG (sin CS averages) para M1 como comparación
#
# ROBUSTEZ:
#   - PCCE (plm::pcce) para M2/M3 (Pesaran 2006)
#   - FE (within) + Driscoll-Kraay para M1–M3
#   - FE en primeras diferencias (Plan B, validez I(1) sin cointegración fuerte)
#
# José Alberto Pérez Gómez
# @jalberper
# Diciembre 2025
# 
# 
################################################################################

# =============================================================================
# BLOQUE 0) PAQUETES
# =============================================================================
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(readr)
  library(stringr)
  library(lubridate)
  library(ggplot2)
  
  library(wbstats)
  library(jsonlite)
  library(httr)
  
  library(plm)
  library(lmtest)
  library(sandwich)
  
  library(modelsummary)
  library(flextable)
  library(officer)
  
  library(MASS)   # mvrnorm para bootstrap paramétrico de turning points
})

# =============================================================================
# BLOQUE 1) RUTAS + FREEZE
# =============================================================================
root_dir   <- getwd()
dir_raw    <- file.path(root_dir, "data_raw")
dir_clean  <- file.path(root_dir, "data_clean")
dir_tables <- file.path(root_dir, "tables")
dir_figs   <- file.path(root_dir, "figures")

dir.create(dir_raw,   showWarnings = FALSE, recursive = TRUE)
dir.create(dir_clean, showWarnings = FALSE, recursive = TRUE)
dir.create(dir_tables, showWarnings = FALSE, recursive = TRUE)
dir.create(dir_figs,   showWarnings = FALSE, recursive = TRUE)

freeze_date <- Sys.Date()

save_freeze <- function(df, basename){
  csv_path <- file.path(dir_raw, paste0(basename, "_", freeze_date, ".csv"))
  rds_path <- file.path(dir_raw, paste0(basename, "_", freeze_date, ".rds"))
  readr::write_csv(df, csv_path)
  saveRDS(df, rds_path)
  message("Freeze guardado: ", basename)
  invisible(list(csv = csv_path, rds = rds_path))
}

# =============================================================================
# BLOQUE 2) PARÁMETROS + BRICS + INDICADORES
# =============================================================================
brics_iso3 <- c("BRA","RUS","IND","CHN","ZAF")
start_year <- 1990
end_year   <- 2024

# WDI indicators (consistentes con tu trabajo)
# PIB pc PPP constante 2017 int$: NY.GDP.PCAP.PP.KD
# CO2 pc (AR5): EN.GHG.CO2.PC.CE.AR5
# Energía pc (kgoe): EG.USE.PCAP.KG.OE
# Controles: industria% PIB, urbanización
wdi_indicators <- c(
  GDP_pc       = "NY.GDP.PCAP.PP.KD",
  CO2_pc       = "EN.GHG.CO2.PC.CE.AR5",
  energy_pc    = "EG.USE.PCAP.KG.OE",
  industry_pib = "NV.IND.TOTL.ZS",
  urban        = "SP.URB.TOTL.IN.ZS"
)

# =============================================================================
# BLOQUE 3) DESCARGA WDI (freeze)
# =============================================================================
# Indicadores WDI
wdi_indicators <- c(
  GDP_pc       = "NY.GDP.PCAP.PP.KD",
  CO2_pc       = "EN.GHG.CO2.PC.CE.AR5",
  energy_pc    = "EG.USE.PCAP.KG.OE",
  industry_pib = "NV.IND.TOTL.ZS",
  urban        = "SP.URB.TOTL.IN.ZS"
)

# Descarga WDI
wdi_df <- wb_data(
  indicator   = wdi_indicators,
  country     = brics_iso3,
  start_date  = start_year,
  end_date    = end_year,
  return_wide = TRUE
)

# Verifica que SÍ se creó
stopifnot(exists("wdi_df"))

# Limpieza mínima y segura
wdi_df <- wdi_df %>%
  rename(year = date) %>%
  mutate(year = as.integer(year)) %>%
  dplyr::select(dplyr::any_of(c(
    "country", "iso3c", "year",
    names(wdi_indicators)
  )))

# Guarda cache
saveRDS(wdi_df, wdi_cache_rds)

# Revisión rápida
glimpse(wdi_df)

# Cobertura rápida
cov_wdi <- wdi_df %>%
  group_by(iso3c) %>%
  summarise(
    min_year = min(year, na.rm = TRUE),
    max_year = max(year, na.rm = TRUE),
    na_GDP = sum(is.na(GDP_pc)),
    na_CO2 = sum(is.na(CO2_pc)),
    na_EN  = sum(is.na(energy_pc)),
    .groups = "drop"
  )
print(cov_wdi)

# =============================================================================
# BLOQUE 4) IMF WEO PPPPC PARA RUSIA 2022–2024 (freeze) - SOLO SI WDI GDP_pc ES NA
# =============================================================================
# Nota metodológica: IMF PPPPC no es idéntico a WDI NY.GDP.PCAP.PP.KD,
# pero puede usarse como "proxy de extensión" para 2022–2024, documentado y con robustez.
imf_cache_rds <- file.path(dir_raw, paste0("IMF_PPPPC_RUS_", freeze_date, ".rds"))

fetch_imf_datamapper <- function(indicator = "PPPPC", country = "RUS"){
  url <- paste0("https://www.imf.org/external/datamapper/api/v1/", indicator, "/", country)
  resp <- httr::GET(url)
  httr::stop_for_status(resp)
  txt <- httr::content(resp, as = "text", encoding = "UTF-8")
  jsonlite::fromJSON(txt, simplifyVector = FALSE)
}

if (file.exists(imf_cache_rds)) {
  imf_weo <- readRDS(imf_cache_rds)
  message("IMF cargado desde freeze: ", imf_cache_rds)
} else {
  imf_json <- fetch_imf_datamapper("PPPPC","RUS")
  weo_vals <- imf_json$values$PPPPC$RUS
  
  imf_weo <- tibble::tibble(
    iso3c = "RUS",
    year  = as.integer(names(weo_vals)),
    GDP_pc_PPP_imf = as.numeric(unlist(weo_vals))
  ) %>%
    filter(year %in% 2022:2024)
  
  saveRDS(imf_weo, imf_cache_rds)
  save_freeze(imf_weo, "IMF_PPPPC_RUS")
}
# =============================================================================
# BLOQUE 5) OWID PARA COMPLETAR ENERGÍA/CO2 (freeze) + CONVERSIÓN A kgoe
# =============================================================================
# - Energía OWID como kWh/persona -> convertir a kgoe/persona para comparabilidad.
# Conversión:
#   1 toe = 11,630 kWh = 1,000 kgoe  =>  kgoe = kWh * (1000/11630)

kwh_to_kgoe <- function(kwh){ kwh * (1000/11630) }

# Función genérica para bajar un dataset de OWID Grapher por slug y devolverlo como tibble
owid_grapher_csv <- function(slug) {
  url <- paste0("https://ourworldindata.org/grapher/", slug, ".csv")
  readr::read_csv(url, show_col_types = FALSE)
}

owid_energy_cache <- file.path(dir_raw, paste0("OWID_energy_pc_kwh_", freeze_date, ".rds"))
owid_co2_cache    <- file.path(dir_raw, paste0("OWID_co2_pc_", freeze_date, ".rds"))

# Slugs CORRECTOS
slug_energy <- "per-capita-energy-use"     # kWh/person (OWID Grapher)
slug_co2    <- "co-emissions-per-capita"   # tCO2/person (OWID Grapher)

# ----------------------------
# Energía OWID
# ----------------------------
if (file.exists(owid_energy_cache)) {
  df_energy_owid <- readRDS(owid_energy_cache)
  message("OWID energía cargado desde freeze: ", owid_energy_cache)
} else {
  df_energy_owid <- owid_grapher_csv(slug_energy)
  saveRDS(df_energy_owid, owid_energy_cache)
  save_freeze(df_energy_owid, "OWID_energy_use_per_person")
}

# ----------------------------
# CO2 OWID
# ----------------------------
if (file.exists(owid_co2_cache)) {
  df_co2_owid <- readRDS(owid_co2_cache)
  message("OWID CO2 cargado desde freeze: ", owid_co2_cache)
} else {
  df_co2_owid <- owid_grapher_csv(slug_co2)
  saveRDS(df_co2_owid, owid_co2_cache)
  save_freeze(df_co2_owid, "OWID_co2_emissions_per_capita")
}

# =============================================================================
# Normalización OWID -> iso3c/year y conversión energía a kgoe
# =============================================================================
# OWID Grapher estándar: Entity | Code | Year | <value>
# La columna de valores (4ta) puede cambiar de nombre, por eso la renombramos robustamente.

# Energía (kWh/persona)
val_energy <- names(df_energy_owid)[4]
energy_owid_brics <- df_energy_owid %>%
  dplyr::rename(country_owid = Entity, iso3c = Code, year = Year) %>%
  dplyr::rename(energy_pc_owid_kwh = !!val_energy) %>%
  dplyr::mutate(year = as.integer(year)) %>%
  dplyr::filter(iso3c %in% brics_iso3, year >= start_year, year <= end_year) %>%
  dplyr::transmute(
    iso3c, year,
    energy_pc_owid_kwh,
    energy_pc_owid_kgoe = kwh_to_kgoe(energy_pc_owid_kwh)
  )

# CO2 (tCO2/persona)
val_co2 <- names(df_co2_owid)[4]
co2_owid_brics <- df_co2_owid %>%
  dplyr::rename(country_owid = Entity, iso3c = Code, year = Year) %>%
  dplyr::rename(CO2_pc_owid = !!val_co2) %>%
  dplyr::mutate(year = as.integer(year)) %>%
  dplyr::filter(iso3c %in% brics_iso3, year >= start_year, year <= end_year) %>%
  dplyr::transmute(
    iso3c, year,
    CO2_pc_owid
  )
dplyr::glimpse(energy_owid_brics)
dplyr::glimpse(co2_owid_brics)

# =============================================================================
# BLOQUE 6) ENSAMBLE MACRO (WDI + IMF + OWID) + PROVENANCE
# =============================================================================
datos <- wdi_df %>%
  left_join(imf_weo, by = c("iso3c","year")) %>%
  left_join(energy_owid_brics, by = c("iso3c","year")) %>%
  left_join(co2_owid_brics,    by = c("iso3c","year")) %>%
  mutate(
    # GDP: llenar solo donde WDI falta y IMF existe
    GDP_pc_filled = if_else(is.na(GDP_pc) & !is.na(GDP_pc_PPP_imf), GDP_pc_PPP_imf, GDP_pc),
    
    # Energía: llenar solo donde WDI falta y OWID(kgoe) existe
    energy_pc_filled = if_else(is.na(energy_pc) & !is.na(energy_pc_owid_kgoe), energy_pc_owid_kgoe, energy_pc),
    
    # CO2: llenar solo donde WDI falta y OWID existe (si aplica)
    CO2_pc_filled = if_else(is.na(CO2_pc) & !is.na(CO2_pc_owid), CO2_pc_owid, CO2_pc)
  )

provenance <- datos %>%
  transmute(
    iso3c, year,
    src_GDP = case_when(
      !is.na(GDP_pc) ~ "WDI:NY.GDP.PCAP.PP.KD",
      is.na(GDP_pc) & !is.na(GDP_pc_PPP_imf) ~ "IMF:PPPPC (proxy extension)",
      TRUE ~ NA_character_
    ),
    src_energy = case_when(
      !is.na(energy_pc) ~ "WDI:EG.USE.PCAP.KG.OE",
      is.na(energy_pc) & !is.na(energy_pc_owid_kgoe) ~ "OWID:energy-use-per-person (kWh->kgoe)",
      TRUE ~ NA_character_
    ),
    src_CO2 = case_when(
      !is.na(CO2_pc) ~ "WDI:EN.GHG.CO2.PC.CE.AR5",
      is.na(CO2_pc) & !is.na(CO2_pc_owid) ~ "OWID:co-emissions-per-capita (proxy extension)",
      TRUE ~ NA_character_
    )
  )

save_freeze(provenance, "Provenance_BRICS")

saveRDS(datos, file.path(dir_clean, paste0("macro_WDI_IMF_OWID_", freeze_date, ".rds")))
readr::write_csv(datos, file.path(dir_clean, paste0("macro_WDI_IMF_OWID_", freeze_date, ".csv")))

# =============================================================================
# BLOQUE 7) PANEL MACRO BASE + TRANSFORMACIONES (EKC)
# =============================================================================
panel_brics <- datos %>%
  filter(iso3c %in% brics_iso3) %>%
  filter(
    !is.na(GDP_pc_filled), GDP_pc_filled > 0,
    !is.na(CO2_pc_filled), CO2_pc_filled > 0,
    !is.na(energy_pc_filled), energy_pc_filled > 0,
    !is.na(industry_pib), industry_pib >= 0, industry_pib <= 100,
    !is.na(urban), urban >= 0, urban <= 100
  ) %>%
  mutate(
    ln_GDPpc   = log(GDP_pc_filled),
    ln_CO2pc   = log(CO2_pc_filled),
    ln_energy_pc = log(energy_pc_filled)
  )

# Centrado global (pool) para VIF y estabilidad numérica:
panel_brics <- panel_brics %>%
  mutate(
    ln_GDPpc_centered = ln_GDPpc - mean(ln_GDPpc, na.rm = TRUE),
    ln_GDPpc_sq       = ln_GDPpc_centered^2
  )

# =============================================================================
# BLOQUE 8) V-DEM (v15) + CONSTRUCCIÓN DE ÍNDICES (TUS VARIABLES EXACTAS)
# =============================================================================
# Fuente recomendada:
# - Usar vdemdata (R) y congelar en RDS para reproducibilidad.
# =============================================================================

# --- Paths freeze ---
vdem_path_rds <- file.path(dir_raw, paste0("V-Dem-CY-Full+Others-v15_", freeze_date, ".rds"))
vdem_path_csv <- file.path(dir_raw, "V-Dem-CY-Full+Others-v15.csv")  # opcional

# --- 8.1 Cargar o crear freeze (RDS) ---
if (file.exists(vdem_path_rds)) {
  
  vdem_raw <- readRDS(vdem_path_rds)
  message("V-Dem cargado desde freeze local (RDS): ", vdem_path_rds)
  
} else {
  
  # Instalar/cargar vdemdata (GitHub oficial)
  if (!requireNamespace("vdemdata", quietly = TRUE)) {
    if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
    remotes::install_github("vdeminstitute/vdemdata")
  }
  library(vdemdata)
  
  # Dataset completo Country-Year (Full+Others)
  vdem_raw <- vdemdata::vdem
  
  # Congelar en RDS
  saveRDS(vdem_raw, vdem_path_rds)
  message("V-Dem v15 congelado en RDS: ", vdem_path_rds)
  
  # (Opcional) exportar CSV local para auditoría humana (archivo grande)
  if (!file.exists(vdem_path_csv)) {
    if (!requireNamespace("data.table", quietly = TRUE)) install.packages("data.table")
    data.table::fwrite(vdem_raw, vdem_path_csv)
    message("V-Dem v15 exportado a CSV local: ", vdem_path_csv)
  }
  
  # (Opcional) bitácora freeze si la estás usando:
  # save_freeze(vdem_raw, "VDEM_CY_FullOthers_v15")
}

# Chequeo duro
stopifnot(is.data.frame(vdem_raw), nrow(vdem_raw) > 0)

# --- 8.2 Armonización de llaves (ISO3) y validaciones mínimas ---
# vdemdata usa 'country_text_id' como ISO3 (ej. BRA, CHN, IND, RUS, ZAF).
# Estandarizamos a 'iso3c' para merge con WDI/OWID.
if ("country_text_id" %in% names(vdem_raw)) {
  vdem_raw <- vdem_raw %>%
    dplyr::rename(iso3c = country_text_id)
}

if (!("iso3c" %in% names(vdem_raw))) {
  stop("V-Dem no tiene columna ISO3 ('iso3c' o 'country_text_id'). Revisa estructura de vdemdata::vdem.")
}
if (!("year" %in% names(vdem_raw))) {
  stop("V-Dem no tiene columna 'year'. Revisa estructura de vdemdata::vdem.")
}

# Asegurar tipo de year antes de filtrar
vdem_raw <- vdem_raw %>%
  dplyr::mutate(year = as.integer(year))

# --- 8.3 Filtro BRICS + ventana temporal ---
vdem_brics <- vdem_raw %>%
  dplyr::filter(
    iso3c %in% brics_iso3,
    year >= start_year, year <= end_year
  ) %>%
  dplyr::select(iso3c, year, dplyr::everything())

# --- 8.4 Variables exactas ---
# --- VARIABLES EXACTAS (V-Dem v15) ---
# Estabilidad política:
# - v2x_clphy   : Physical violence index
# - v2svstterr : State threats – territorial

# Gobernanza:
# - v2x_rule   : Rule of law
# - v2castate  : State capacity
# - v2x_corr   : Political corruption (invertida)

var_clphy   <- "v2x_clphy"
var_stterr  <- "v2svstterr"
var_rulelaw <- "v2x_rule"
var_castate <- "v2castate"
var_corr    <- "v2x_corr"

# Utilidad: promedio con mínimo k observado
row_mean_min_k <- function(df, k_min) {
  x  <- as.matrix(df)
  ok <- rowSums(!is.na(x))
  out <- rowMeans(x, na.rm = TRUE)
  out[ok < k_min] <- NA_real_
  out
}

# Chequeo de variables obligatorias
must_exist <- function(v){
  v <- v[!is.na(v)]
  missing <- setdiff(v, names(vdem_brics))
  if (length(missing) > 0) stop("Faltan variables V-Dem: ", paste(missing, collapse = ", "))
  invisible(TRUE)
}
must_exist(c(var_clphy, var_stterr, var_rulelaw, var_castate, var_corr))

# --- 8.5 Construcción de índices (z-scores pooled) ---
vdem_indices <- vdem_brics %>%
  dplyr::transmute(
    iso3c, year,
    clphy_raw   = .data[[var_clphy]],
    stterr_raw  = .data[[var_stterr]],
    rol_raw     = .data[[var_rulelaw]],
    castate_raw = .data[[var_castate]],
    corr_raw    = .data[[var_corr]],
    corr_good   = -1 * corr_raw
  ) %>%
  dplyr::mutate(
    # z-scores pooled (1990–2024 BRICS). Mantiene comparabilidad temporal.
    z_clphy   = as.numeric(scale(clphy_raw)),
    z_stterr  = as.numeric(scale(stterr_raw)),
    z_rol     = as.numeric(scale(rol_raw)),
    z_castate = as.numeric(scale(castate_raw)),
    z_corrctl = as.numeric(scale(corr_good))
  ) %>%
  dplyr::mutate(
    estabilidad_politica_vdem = row_mean_min_k(dplyr::select(., z_clphy, z_stterr), k_min = 2),
    gobernanza_vdem           = row_mean_min_k(dplyr::select(., z_rol, z_castate, z_corrctl), k_min = 3)
  ) %>%
  dplyr::select(iso3c, year, estabilidad_politica_vdem, gobernanza_vdem)

# --- 8.6 Chequeo de cobertura (sanity check) ---
coverage_check <- vdem_indices %>%
  dplyr::group_by(iso3c) %>%
  dplyr::summarise(
    min_year = min(year, na.rm = TRUE),
    max_year = max(year, na.rm = TRUE),
    n_years = dplyr::n(),
    na_estabilidad = sum(is.na(estabilidad_politica_vdem)),
    na_gobernanza  = sum(is.na(gobernanza_vdem)),
    .groups = "drop"
  )

print(coverage_check)

# =============================================================================
# BLOQUE 9) PANEL FINAL (macro + V-Dem) + INTERACCIONES (within-country)
# =============================================================================
# Objetivos:
# (i) Join seguro (sin duplicados / sin many-to-many).
# (ii) Validar que el merge NO perdió años por tipos de dato.
# (iii) Construir interacciones con interpretación estrictamente within-country.
# =============================================================================

# ----------------------------
# 9.0 Chequeos duros de llaves
# ----------------------------
# Asegurar tipos de datos consistentes
panel_brics <- panel_brics %>% dplyr::mutate(year = as.integer(year))
vdem_indices <- vdem_indices %>% dplyr::mutate(year = as.integer(year))

# Unicidad de llaves
stopifnot(
  nrow(panel_brics) == nrow(dplyr::distinct(panel_brics, iso3c, year))
)
stopifnot(
  nrow(vdem_indices) == nrow(dplyr::distinct(vdem_indices, iso3c, year))
)

# Cobertura esperada de V-Dem sobre el panel macro:
# (V-Dem ya está completo 1990–2024, así que cualquier pérdida sería por merge/macro)
anti_missing <- dplyr::anti_join(panel_brics, vdem_indices, by = c("iso3c", "year"))
if (nrow(anti_missing) > 0) {
  message("ALERTA: Hay filas macro sin match en V-Dem (revisar llaves year/iso3c). Ejemplos:")
  print(head(anti_missing, 10))
  stop("Merge macro–V-Dem incompleto. Corrige antes de continuar.")
}

# ----------------------------
# 9.1 Merge final macro + V-Dem
# ----------------------------
panel_brics2 <- panel_brics %>%
  dplyr::left_join(vdem_indices, by = c("iso3c","year")) %>%
  dplyr::arrange(iso3c, year)

# Chequeo: después del join NO debe haber NA en instituciones
na_inst <- panel_brics2 %>%
  dplyr::summarise(
    na_estab = sum(is.na(estabilidad_politica_vdem)),
    na_gov   = sum(is.na(gobernanza_vdem))
  )
print(na_inst)
stopifnot(na_inst$na_estab == 0, na_inst$na_gov == 0)

# ----------------------------
# 9.2 Centrados para interacciones (within-country)
# ----------------------------
# Nota conceptual para el paper:
# - Estas interacciones identifican si CAMBIOS temporales dentro de cada país
#   en estabilidad/gobernanza moderan el efecto del PIB/energía.
# - No interpretarlas como “niveles institucionales entre países”.

panel_brics2 <- panel_brics2 %>%
  dplyr::group_by(iso3c) %>%
  dplyr::mutate(
    ln_energy_w = ln_energy_pc - mean(ln_energy_pc, na.rm = TRUE),
    gov_w       = gobernanza_vdem - mean(gobernanza_vdem, na.rm = TRUE),
    stab_w      = estabilidad_politica_vdem - mean(estabilidad_politica_vdem, na.rm = TRUE),
    
    # Interacciones:
    logpib_estabilidad = ln_GDPpc_centered * stab_w,
    energia_gobierno   = ln_energy_w * gov_w
  ) %>%
  dplyr::ungroup()

# ----------------------------
# 9.3 Auditoría de balance (si afirmas “balanced”, que sea real)
# ----------------------------
bal <- panel_brics2 %>%
  dplyr::count(iso3c, name = "n_years") %>%
  dplyr::summarise(min = min(n_years), max = max(n_years))

print(bal)
print(range(panel_brics2$year, na.rm = TRUE))

# Auditoría adicional: años faltantes por país (si quieres documentarlo)
missing_years <- panel_brics2 %>%
  dplyr::group_by(iso3c) %>%
  dplyr::summarise(
    min_year = min(year, na.rm = TRUE),
    max_year = max(year, na.rm = TRUE),
    n_years  = dplyr::n(),
    .groups  = "drop"
  )
print(missing_years)

na_core <- panel_brics2 %>%
  dplyr::group_by(iso3c) %>%
  dplyr::summarise(
    na_ln_co2    = sum(is.na(ln_CO2pc)),
    na_ln_gdp    = sum(is.na(ln_GDPpc_centered)),
    na_ln_gdp2   = sum(is.na(ln_GDPpc_sq)),
    na_ln_energy = sum(is.na(ln_energy_pc)),
    .groups = "drop"
  )
print(na_core)

panel_brics2 %>%
  dplyr::summarise(
    sd_logpib_estab = sd(logpib_estabilidad, na.rm = TRUE),
    sd_energia_gob  = sd(energia_gobierno, na.rm = TRUE)
  ) %>% print()
colnames(panel_brics2)
# ----------------------------
# 9.4 Guardado
# ----------------------------
saveRDS(panel_brics2, file.path(dir_clean, paste0("panel_BRICS_macro_vdem_", freeze_date, ".rds")))
readr::write_csv(panel_brics2, file.path(dir_clean, paste0("panel_BRICS_macro_vdem_", freeze_date, ".csv")))

# =============================================================================
# =============================================================================
# BLOQUE 10) OBJETOS PANEL (plm) + PREP DINÁMICA (lags/diffs) PARA CS-ECM/PMG
# =============================================================================
# Este bloque:
# - Declara pdata.frame (plm)
# - Asegura orden y estructura de panel
# - Crea lags y primeras diferencias (coherente con I(1)/I(0))
# - Deja un dataset "dinámico" listo para:
#     (i) ECM-MG / ARDL-PMG
#     (ii) CCE/CS-ECM con promedios transversales (si aplica)
# =============================================================================

# ----------------------------
# 10.1 Panel plm
# ----------------------------
panel_brics2 <- panel_brics2 %>%
  dplyr::arrange(iso3c, year) %>%
  dplyr::mutate(year = as.integer(year))

pdata <- plm::pdata.frame(panel_brics2, index = c("iso3c", "year"), drop.index = FALSE, row.names = TRUE)

# Sanity: panel balanceado y rango temporal
bal_plm <- plm::pdim(pdata)
print(bal_plm)  # N, T, balanced?

# ----------------------------
# 10.2 Variables necesarias para dinámica
# ----------------------------
# Convención:
# - Variables I(1): ln_CO2pc, ln_GDPpc_centered, ln_GDPpc_sq, ln_energy_pc
# - Variables I(0): estabilidad_politica_vdem, gobernanza_vdem (y sus centrados within)

# Lags en niveles (t-1)
pdata$L_ln_CO2pc   <- plm::lag(pdata$ln_CO2pc, 1)
pdata$L_ln_GDPpc   <- plm::lag(pdata$ln_GDPpc_centered, 1)
pdata$L_ln_GDPpc2  <- plm::lag(pdata$ln_GDPpc_sq, 1)
pdata$L_ln_energy  <- plm::lag(pdata$ln_energy_pc, 1)

pdata$L_stab       <- plm::lag(pdata$estabilidad_politica_vdem, 1)
pdata$L_gov        <- plm::lag(pdata$gobernanza_vdem, 1)

pdata$L_logpib_estabilidad <- plm::lag(pdata$logpib_estabilidad, 1)
pdata$L_energia_gobierno   <- plm::lag(pdata$energia_gobierno, 1)

# Primeras diferencias (Δ) – forma correcta en plm
pdata$D_ln_CO2pc   <- pdata$ln_CO2pc - plm::lag(pdata$ln_CO2pc, 1)
pdata$D_ln_GDPpc   <- pdata$ln_GDPpc_centered - plm::lag(pdata$ln_GDPpc_centered, 1)
pdata$D_ln_GDPpc2  <- pdata$ln_GDPpc_sq - plm::lag(pdata$ln_GDPpc_sq, 1)
pdata$D_ln_energy  <- pdata$ln_energy_pc - plm::lag(pdata$ln_energy_pc, 1)

pdata$D_stab       <- pdata$estabilidad_politica_vdem - plm::lag(pdata$estabilidad_politica_vdem, 1)
pdata$D_gov        <- pdata$gobernanza_vdem - plm::lag(pdata$gobernanza_vdem, 1)

pdata$D_logpib_estabilidad <- pdata$logpib_estabilidad - plm::lag(pdata$logpib_estabilidad, 1)
pdata$D_energia_gobierno   <- pdata$energia_gobierno - plm::lag(pdata$energia_gobierno, 1)

summary(pdata$D_ln_CO2pc)
summary(pdata$D_ln_GDPpc)

# Instituciones (I(0)) – cambios (Δ) por si decides incluir short-run institucional
pdata$D_stab <- pdata$estabilidad_politica_vdem - plm::lag(pdata$estabilidad_politica_vdem, 1)
pdata$D_gov  <- pdata$gobernanza_vdem           - plm::lag(pdata$gobernanza_vdem, 1)

# Interacciones (ya definidas en niveles within-country) – cambios (Δ)
pdata$D_logpib_estabilidad <- pdata$logpib_estabilidad - plm::lag(pdata$logpib_estabilidad, 1)
pdata$D_energia_gobierno   <- pdata$energia_gobierno   - plm::lag(pdata$energia_gobierno, 1)

# ----------------------------
# 10.3 Dataset "dinámico" para estimación (quita la primera observación por país)
# ----------------------------
dyn_vars <- c(
  "iso3c", "year",
  "ln_CO2pc", "ln_GDPpc_centered", "ln_GDPpc_sq", "ln_energy_pc",
  "estabilidad_politica_vdem", "gobernanza_vdem",
  "logpib_estabilidad", "energia_gobierno",
  "L_ln_CO2pc", "L_ln_GDPpc", "L_ln_GDPpc2", "L_ln_energy",
  "L_stab", "L_gov", "L_logpib_estabilidad", "L_energia_gobierno",
  "D_ln_CO2pc", "D_ln_GDPpc", "D_ln_GDPpc2", "D_ln_energy",
  "D_stab", "D_gov", "D_logpib_estabilidad", "D_energia_gobierno"
)

pdata_dyn <- as.data.frame(pdata)[, dyn_vars]
pdata_dyn <- pdata_dyn %>% tidyr::drop_na(D_ln_CO2pc, L_ln_CO2pc)  # elimina t=1990 por país (lags/diffs)

# Chequeo final de tamaño esperado: 5 países * (35-1)=170 obs
message("Obs esperadas (balanced dinámico): ", length(unique(pdata_dyn$iso3c)) * (length(unique(panel_brics2$year)) - 1))
message("Obs reales en pdata_dyn: ", nrow(pdata_dyn))
stopifnot(nrow(pdata_dyn) == length(unique(pdata_dyn$iso3c)) * (length(unique(panel_brics2$year)) - 1))

# Guardar opcional para inspección
saveRDS(pdata_dyn, file.path(dir_clean, paste0("pdata_dyn_", freeze_date, ".rds")))

stopifnot(sum(is.na(pdata$D_stab)) == length(unique(pdata$iso3c)))
stopifnot(sum(is.na(pdata$D_gov))  == length(unique(pdata$iso3c)))

# =============================================================================
# BLOQUE 11) FE (within) + DRISCOLL–KRAAY (ROBUSTEZ, NO MAIN)
# =============================================================================
# - Estos FE en niveles se reportan SOLO como robustez descriptiva.
# - La especificación principal debe ser dinámica (CS-ECM/CCE/PMG) por I(1)+cointegración.
# =============================================================================

# Asegurar paquetes
if (!requireNamespace("plm", quietly = TRUE)) stop("Falta paquete plm")
if (!requireNamespace("sandwich", quietly = TRUE)) stop("Falta paquete sandwich")
if (!requireNamespace("lmtest", quietly = TRUE)) stop("Falta paquete lmtest")
if (!requireNamespace("readr", quietly = TRUE)) stop("Falta paquete readr")
if (!requireNamespace("tibble", quietly = TRUE)) stop("Falta paquete tibble")
if (!requireNamespace("dplyr", quietly = TRUE)) stop("Falta paquete dplyr")

# ----------------------------
# 11.1 Modelos FE (within)
# ----------------------------
M1_fe <- plm::plm(
  ln_CO2pc ~ ln_GDPpc_centered + ln_GDPpc_sq + ln_energy_pc,
  data = pdata, model = "within", effect = "individual"
)

M2_fe <- plm::plm(
  ln_CO2pc ~ ln_GDPpc_centered + ln_GDPpc_sq + ln_energy_pc +
    estabilidad_politica_vdem + gobernanza_vdem,
  data = pdata, model = "within", effect = "individual"
)

M3_fe <- plm::plm(
  ln_CO2pc ~ ln_GDPpc_centered + ln_GDPpc_sq + ln_energy_pc +
    estabilidad_politica_vdem + gobernanza_vdem +
    logpib_estabilidad + energia_gobierno,
  data = pdata, model = "within", effect = "individual"
)

# ----------------------------
# 11.2 Errores Driscoll–Kraay (SCC)
# ----------------------------
# Regla práctica para maxlag en DK: usar orden ~ floor(T^(1/3)).
# Aquí T=35 => 35^(1/3)=3.27 => 3.
T_panel <- length(unique(panel_brics2$year))
maxlag_default <- max(1, floor(T_panel^(1/3)))

vcovDK <- function(mod, maxlag = maxlag_default){
  plm::vcovSCC(mod, type = "HC1", maxlag = maxlag)
}

# Coeficientes con SE DK
M1_fe_DK <- lmtest::coeftest(M1_fe, vcov. = vcovDK(M1_fe))
M2_fe_DK <- lmtest::coeftest(M2_fe, vcov. = vcovDK(M2_fe))
M3_fe_DK <- lmtest::coeftest(M3_fe, vcov. = vcovDK(M3_fe))

# Guardar resultados DK en CSV (tidy)
tidy_coeftest <- function(ct, model_name){
  out <- data.frame(
    term      = rownames(ct),
    estimate  = unname(ct[,1]),
    std_error = unname(ct[,2]),
    statistic = unname(ct[,3]),
    p_value   = unname(ct[,4]),
    model     = model_name,
    row.names = NULL
  )
  out
}

dk_table <- dplyr::bind_rows(
  tidy_coeftest(M1_fe_DK, "M1_FE_DK"),
  tidy_coeftest(M2_fe_DK, "M2_FE_DK"),
  tidy_coeftest(M3_fe_DK, "M3_FE_DK")
)

readr::write_csv(dk_table, file.path(dir_tables, "Table_FE_DK_Coefficients.csv"))

# ----------------------------
# 11.3 Pesaran CD: ojo con interpretación bajo FE
# ----------------------------
# 
# - Reportar CD sobre RESIDUALES del FE (equivalente a test sobre el error común restante).
# - Con N=5, el poder es limitado; aun así es un check estándar.

cd_M1 <- plm::pcdtest(M1_fe, test = "cd")
cd_M2 <- plm::pcdtest(M2_fe, test = "cd")
cd_M3 <- plm::pcdtest(M3_fe, test = "cd")

cd_table <- tibble::tibble(
  model = c("M1_FE","M2_FE","M3_FE"),
  statistic = c(unname(cd_M1$statistic), unname(cd_M2$statistic), unname(cd_M3$statistic)),
  p_value   = c(cd_M1$p.value, cd_M2$p.value, cd_M3$p.value),
  N = rep(plm::pdim(pdata)$nT$n, 3),
  T = rep(plm::pdim(pdata)$nT$T, 3),
  note = c(
    "FE residual CD test (levels). DK corrects SE only; coefficients may remain biased under strong common factors.",
    "FE residual CD test (levels). Prefer CCE/CS-ECM as main if CD significant.",
    "FE residual CD test (levels). Prefer CCE/CS-ECM as main if CD significant."
  )
)

readr::write_csv(cd_table, file.path(dir_tables, "Table_PesaranCD_FE.csv"))

# (Opcional) imprimir en consola
print(cd_table)
message("DK maxlag usado: ", maxlag_default)

# =============================================================================
# 13.2 — Panel Unit Root Tests (IPS) 
#(Appendix A1, Panel B)
# Objetivo:
# - Generar Panel B de la Table A1 (Unit Roots) con pruebas IPS.
# - Exportar CSV a /tables para que el bloque de exportación a Word lo detecte.
# =============================================================================

# --- Seguridad: crear carpeta si no existe ---
if (!dir.exists(dir_tables)) dir.create(dir_tables, recursive = TRUE)

# --- Helper (opcional): pretty_term fallback ---
if (!exists("pretty_term")) {
  pretty_term <- function(x) x
}

# --- Variables a testear (ajusta si lo necesitas) ---
ips_vars <- c(
  "ln_CO2pc",
  "ln_GDPpc_centered",
  "ln_GDPpc_sq",
  "ln_energy_pc",
  "estabilidad_politica_vdem",
  "gobernanza_vdem",
  "logpib_estabilidad",
  "energia_gobierno"
)

# --- Chequeo de existencia en pdata ---
missing_ips <- setdiff(ips_vars, names(pdata))
if (length(missing_ips) > 0) {
  warning("IPS: faltan variables en pdata y se omitirán: ", paste(missing_ips, collapse = ", "))
  ips_vars <- setdiff(ips_vars, missing_ips)
}

# Función robusta para correr IPS y extraer estadístico + p-value
# -----------------------------------------------------------------------------
run_ips_one <- function(v, exo = "intercept", lags = "AIC") {
  fml <- stats::as.formula(paste0(v, " ~ 1"))
  
  out <- tryCatch(
    plm::purtest(
      fml,
      data = pdata,
      test = "ips",
      exo  = exo,
      lags = lags
    ),
    error = function(e) e
  )
  
  if (inherits(out, "error")) {
    return(tibble::tibble(
      variable  = v,
      statistic = NA_real_,
      p_value   = NA_real_,
      exo       = exo,
      lags      = as.character(lags),
      status    = paste("ERROR:", out$message)
    ))
  }
  
  stat <- NA_real_
  pval <- NA_real_
  
  # --- purtest devuelve un objeto "purtest"
  # --- out$statistic es un objeto "htest" y AHÍ está el p-value
  if (!is.null(out$statistic) && inherits(out$statistic, "htest")) {
    # estadístico
    if (!is.null(out$statistic$statistic) && is.numeric(out$statistic$statistic)) {
      stat <- as.numeric(out$statistic$statistic)[1]
    }
    # p-value
    if (!is.null(out$statistic$p.value) && is.numeric(out$statistic$p.value)) {
      pval <- as.numeric(out$statistic$p.value)[1]
    }
  } else {
    # Fallbacks por si cambian estructuras
    if (!is.null(out$statistic) && is.numeric(out$statistic)) {
      stat <- as.numeric(out$statistic)[1]
    }
    if (!is.null(out$p.value) && is.numeric(out$p.value)) {
      pval <- as.numeric(out$p.value)[1]
    }
    if (is.list(out$statistic) && !is.null(out$statistic$statistic)) {
      stat <- as.numeric(out$statistic$statistic)[1]
    }
    if (is.list(out$statistic) && !is.null(out$statistic$p.value)) {
      pval <- as.numeric(out$statistic$p.value)[1]
    }
  }
  
  tibble::tibble(
    variable  = v,
    statistic = stat,
    p_value   = pval,
    exo       = exo,
    lags      = as.character(lags),
    status    = "OK"
  )
}

# -----------------------------------------------------------------------------
# Correr IPS y exportar
# -----------------------------------------------------------------------------
ips_tbl <- dplyr::bind_rows(
  lapply(ips_vars, run_ips_one, exo = "intercept", lags = "AIC")
) %>%
  dplyr::mutate(
    statistic = round(statistic, 3),
    p_value   = round(p_value, 3),
    variable  = pretty_term(variable)
  ) %>%
  dplyr::rename(
    Variable  = variable,
    Statistic = statistic,
    `p-value` = p_value,
    Exo       = exo,
    Lags      = lags,
    Status    = status
  )

ips_path <- file.path(dir_tables, "IPS_UnitRoot_Tests.csv")
readr::write_csv(ips_tbl, ips_path)

print(ips_tbl)
message("OK: IPS exportado a: ", ips_path)

list.files(dir_tables, pattern = "IPS", full.names = TRUE)

# =============================================================================
# BLOQUE 12) CCE / PCCE (Pesaran 2006) – MANEJO DE DEPENDENCIA CRUZADA (CD)
# =============================================================================
# Objetivo:
# - Si Pesaran CD rechaza independencia (M2/M3), DK NO basta (corrige SE, no sesgo).
# - CCE/PCCE controla factores comunes no observados (common correlated effects).
#
# Nota:
# - La implementación y argumentos exactos dependen del paquete que provee pcce().
# - Aquí forzamos una implementación "Pooled CCE" (PCCE), consistente con Pesaran (2006).
# =============================================================================

# ---- Chequeo: que exista pcce en tu sesión ----
if (!exists("pcce")) {
  stop("No existe la función pcce() en tu entorno. Carga el paquete correspondiente (p.ej., 'plm' si aplica, o el paquete donde la tengas).")
}

# =============================================================================
# 12.1 PCCE pooled (recomendado para CD)
# =============================================================================

M2_pcce <- pcce(
  ln_CO2pc ~ ln_GDPpc_centered + ln_GDPpc_sq + ln_energy_pc +
    estabilidad_politica_vdem + gobernanza_vdem,
  data = pdata,
  model = "pooled"
)

M3_pcce <- pcce(
  ln_CO2pc ~ ln_GDPpc_centered + ln_GDPpc_sq + ln_energy_pc +
    estabilidad_politica_vdem + gobernanza_vdem +
    logpib_estabilidad + energia_gobierno,
  data = pdata,
  model = "pooled"
)

# Resúmenes (para log)
print(summary(M2_pcce))
print(summary(M3_pcce))

# =============================================================================
# 12.2 (Opcional) CCE manual como chequeo de consistencia (NO reemplaza pcce)
# =============================================================================
# - Crear promedios transversales por año de Y y X, incluirlos en FE.
# - Esto se conoce como CCEP "manual" si incluyes Ybar y Xbar.
# - OJO: esto es inferior a una implementación pcce/cce formal cuando existe.
# - Este bloque crea M2_pcce y M3_pcce, y (opcional) CCE manual para chequeo.
# - Al final guarda lo que exista sin romper el script.
# =============================================================================

stopifnot(exists("pcce"))

# Helper: guardar coeficientes de forma uniforme
tidy_any <- function(obj, model_name){
  s <- summary(obj)
  # Intento 1: coef matrix estándar
  if (!is.null(s$coefficients)) {
    cm <- as.data.frame(s$coefficients)
    cm$term <- rownames(cm)
    rownames(cm) <- NULL
    names(cm)[1:4] <- c("estimate","std_error","statistic","p_value")
    cm$model <- model_name
    return(cm[, c("model","term","estimate","std_error","statistic","p_value")])
  }
  # Fallback: si summary no trae 'coefficients', usar coef()
  cm <- data.frame(
    model = model_name,
    term = names(coef(obj)),
    estimate = as.numeric(coef(obj)),
    std_error = NA_real_,
    statistic = NA_real_,
    p_value = NA_real_
  )
  cm
}

# ----------------------------
# 12.1 CCEMG (MAIN): M1, M2, M3
# ----------------------------
M1_ccemg <- pcce(
  ln_CO2pc ~ ln_GDPpc_centered + ln_GDPpc_sq + ln_energy_pc,
  data  = pdata,
  model = "mg"
)

M2_ccemg <- pcce(
  ln_CO2pc ~ ln_GDPpc_centered + ln_GDPpc_sq + ln_energy_pc +
    estabilidad_politica_vdem + gobernanza_vdem,
  data  = pdata,
  model = "mg"
)

M3_ccemg <- pcce(
  ln_CO2pc ~ ln_GDPpc_centered + ln_GDPpc_sq + ln_energy_pc +
    estabilidad_politica_vdem + gobernanza_vdem +
    logpib_estabilidad + energia_gobierno,
  data  = pdata,
  model = "mg"
)

message("OK: CCEMG estimado para M1–M3.")
print(summary(M1_ccemg))
print(summary(M2_ccemg))
print(summary(M3_ccemg))

# ----------------------------
# 12.2 PCCE pooled (sensibilidad): M2 y M3 (opcional)
# ----------------------------
M2_pcce <- pcce(
  ln_CO2pc ~ ln_GDPpc_centered + ln_GDPpc_sq + ln_energy_pc +
    estabilidad_politica_vdem + gobernanza_vdem,
  data  = pdata,
  model = "p"
)

M3_pcce <- pcce(
  ln_CO2pc ~ ln_GDPpc_centered + ln_GDPpc_sq + ln_energy_pc +
    estabilidad_politica_vdem + gobernanza_vdem +
    logpib_estabilidad + energia_gobierno,
  data  = pdata,
  model = "p"
)

message("OK: PCCE pooled estimado para M2–M3 (sensibilidad).")
print(summary(M2_pcce))
print(summary(M3_pcce))

# ----------------------------
# 12.3 Tidy outputs para tablas (CCE-MG + PCCE pooled)
# ----------------------------
tidy_pcce_safe <- function(obj, model_name){
  s <- tryCatch(summary(obj), error = function(e) NULL)
  
  # 1) Caso: summary trae coefficients (matriz/data.frame, con 2–4+ columnas)
  if (!is.null(s) && !is.null(s$coefficients)) {
    cm <- s$coefficients
    
    if (is.list(cm) && !is.matrix(cm) && !is.data.frame(cm)) {
      cm <- tryCatch(as.data.frame(cm), error = function(e) NULL)
    }
    
    if (is.matrix(cm) || is.data.frame(cm)) {
      cm <- as.data.frame(cm)
      cm$term <- rownames(cm)
      rownames(cm) <- NULL
      
      num_cols <- setdiff(names(cm), "term")
      k <- length(num_cols)
      
      res <- data.frame(
        model     = model_name,
        term      = cm$term,
        estimate  = NA_real_,
        std_error = NA_real_,
        statistic = NA_real_,
        p_value   = NA_real_,
        stringsAsFactors = FALSE
      )
      
      if (k >= 1) res$estimate  <- as.numeric(cm[[ num_cols[1] ]])
      if (k >= 2) res$std_error <- as.numeric(cm[[ num_cols[2] ]])
      if (k >= 3) res$statistic <- as.numeric(cm[[ num_cols[3] ]])
      if (k >= 4) res$p_value   <- as.numeric(cm[[ num_cols[4] ]])
      
      return(res)
    }
  }
  
  # 2) Fallback: coef() + (si existe) vcov()
  b <- tryCatch(coef(obj), error = function(e) NULL)
  if (is.null(b)) stop("No pude extraer coef() del objeto pcce().")
  
  res <- data.frame(
    model     = model_name,
    term      = names(b),
    estimate  = as.numeric(b),
    std_error = NA_real_,
    statistic = NA_real_,
    p_value   = NA_real_,
    stringsAsFactors = FALSE
  )
  
  V <- tryCatch(vcov(obj), error = function(e) NULL)
  if (!is.null(V)) {
    se <- sqrt(diag(V))
    if (!is.null(names(se))) se <- se[res$term]
    res$std_error <- as.numeric(se)
    res$statistic <- res$estimate / res$std_error
    res$p_value   <- 2 * stats::pnorm(abs(res$statistic), lower.tail = FALSE)
  }
  
  res
}

# Crear tablas tidy (ESTO ES LO QUE TE FALTABA)
ccemg_table <- dplyr::bind_rows(
  tidy_pcce_safe(M1_ccemg, "M1_CCEMG"),
  tidy_pcce_safe(M2_ccemg, "M2_CCEMG"),
  tidy_pcce_safe(M3_ccemg, "M3_CCEMG")
)

pcce_table <- dplyr::bind_rows(
  tidy_pcce_safe(M2_pcce, "M2_PCCE"),
  tidy_pcce_safe(M3_pcce, "M3_PCCE")
)

# Guardar CSV (para auditoría y para Word)
readr::write_csv(ccemg_table, file.path(dir_tables, "Table_CCE_MG_Coefficients.csv"))
readr::write_csv(pcce_table,  file.path(dir_tables, "Table_CCE_Pooled_Coefficients.csv"))

message("OK: ccemg_table y pcce_table creadas y exportadas a /tables.")

# ----------------------------
# 12.4 Guardado seguro de modelos (RDS)
# ----------------------------
models_to_save <- list(
  M1_ccemg = M1_ccemg,
  M2_ccemg = M2_ccemg,
  M3_ccemg = M3_ccemg,
  M2_pcce  = M2_pcce,
  M3_pcce  = M3_pcce
)

saveRDS(models_to_save, file.path(dir_clean, paste0("models_CCE_", freeze_date, ".rds")))
message("Modelos CCE guardados en: ", file.path(dir_clean, paste0("models_CCE_", freeze_date, ".rds")))

# =============================================================================
# =============================================================================
# BLOQUE 13) PLAN B: FE EN PRIMERAS DIFERENCIAS (VALIDEZ I(1), PIERDE LR)
# =============================================================================
# - Este bloque sirve como "Plan B" si  objeta niveles con I(1).
# - En FD la interpretación es de CORTO PLAZO (Δ), no long-run.
# - FD elimina efectos fijos de país por construcción; aún así usamos plm para consistencia.
# =============================================================================

panel_diff <- panel_brics2 %>%
  dplyr::group_by(iso3c) %>%
  dplyr::arrange(year, .by_group = TRUE) %>%
  dplyr::mutate(
    # Diferencias primeras (usar niveles originales: más defendible)
    d_ln_CO2   = ln_CO2pc - dplyr::lag(ln_CO2pc),
    d_ln_GDP   = ln_GDPpc_centered - dplyr::lag(ln_GDPpc_centered),
    d_ln_GDP2  = ln_GDPpc_sq - dplyr::lag(ln_GDPpc_sq),
    d_ln_EN    = ln_energy_pc - dplyr::lag(ln_energy_pc),
    
    d_stab     = estabilidad_politica_vdem - dplyr::lag(estabilidad_politica_vdem),
    d_gov      = gobernanza_vdem - dplyr::lag(gobernanza_vdem),
    
    # Interacciones: diferenciar las ya construidas en niveles
    d_logpib_stab = logpib_estabilidad - dplyr::lag(logpib_estabilidad),
    d_en_gov      = energia_gobierno   - dplyr::lag(energia_gobierno)
  ) %>%
  dplyr::ungroup() %>%
  tidyr::drop_na(d_ln_CO2, d_ln_GDP, d_ln_GDP2, d_ln_EN)

pdata_d <- plm::pdata.frame(panel_diff, index = c("iso3c","year"))

# ----------------------------
# 13.1 Modelos FD
# ----------------------------
M1_fd <- plm::plm(
  d_ln_CO2 ~ d_ln_GDP + d_ln_GDP2 + d_ln_EN,
  data = pdata_d, model = "pooling"
)

M2_fd <- plm::plm(
  d_ln_CO2 ~ d_ln_GDP + d_ln_GDP2 + d_ln_EN + d_stab + d_gov,
  data = pdata_d, model = "pooling"
)

M3_fd <- plm::plm(
  d_ln_CO2 ~ d_ln_GDP + d_ln_GDP2 + d_ln_EN + d_stab + d_gov + d_logpib_stab + d_en_gov,
  data = pdata_d, model = "pooling"
)

# Nota: en FD, "pooling" es suficiente porque ya eliminaste FE con Δ.
# Si prefieres mantener "within" por consistencia, no cambia mucho aquí; pero pooling es más limpio.

# ----------------------------
# 13.2 SE Driscoll–Kraay en FD (para reporte consistente)
# ----------------------------
T_panel <- length(unique(panel_brics2$year))
maxlag_default <- max(1, floor(T_panel^(1/3)))  # ~3 para T=35

vcovDK_d <- function(mod, maxlag = maxlag_default){
  plm::vcovSCC(mod, type = "HC1", maxlag = maxlag)
}

M1_fd_DK <- lmtest::coeftest(M1_fd, vcov. = vcovDK_d(M1_fd))
M2_fd_DK <- lmtest::coeftest(M2_fd, vcov. = vcovDK_d(M2_fd))
M3_fd_DK <- lmtest::coeftest(M3_fd, vcov. = vcovDK_d(M3_fd))

# Tidy para CSV (auditoría)
tidy_coeftest <- function(ct, model_name){
  data.frame(
    term      = rownames(ct),
    estimate  = unname(ct[,1]),
    std_error = unname(ct[,2]),
    statistic = unname(ct[,3]),
    p_value   = unname(ct[,4]),
    model     = model_name,
    row.names = NULL
  )
}

fd_dk_table <- dplyr::bind_rows(
  tidy_coeftest(M1_fd_DK, "M1_FD_DK"),
  tidy_coeftest(M2_fd_DK, "M2_FD_DK"),
  tidy_coeftest(M3_fd_DK, "M3_FD_DK")
)

readr::write_csv(fd_dk_table, file.path(dir_tables, "Table_FD_DK_Coefficients.csv"))

# ----------------------------
# 13.3 Auditoría de tamaño (debe ser 170 si es balanced)
# ----------------------------
expected_fd <- length(unique(panel_brics2$iso3c)) * (length(unique(panel_brics2$year)) - 1)
message("Obs esperadas FD: ", expected_fd)
message("Obs reales FD: ", nrow(panel_diff))
stopifnot(nrow(panel_diff) == expected_fd)

# ----------------------------
# 13.4 Guardar modelos FD para el bloque final (Word)
# ----------------------------
saveRDS(
  list(M1_fd=M1_fd, M2_fd=M2_fd, M3_fd=M3_fd,
       M1_fd_DK=M1_fd_DK, M2_fd_DK=M2_fd_DK, M3_fd_DK=M3_fd_DK),
  file.path(dir_clean, paste0("models_FD_", freeze_date, ".rds"))
)
# =============================================================================
# BLOQUE 14) MAIN: ECM-MG / CS-ECM-MG (ARDL(1,1,1)) POR PAÍS (MANUAL + TRANSPARENTE)
# =============================================================================
# Para cada país i:
#   Δy_t = c_i + ψ_i*Δy_{t-1} + φ_i*y_{t-1} + Σ b_{k,i}*x_{k,t-1} + Σ γ_{k,i}*Δx_{k,t} + u_t
#
# Long-run (por país): θ_{k,i} = - b_{k,i} / φ_i      (requiere φ_i < 0)
#
# CS-ECM (Chudik & Pesaran style augmentation):
# - Agrega promedios transversales por año en niveles y en diferencias:
#   cs_y, cs_x  y  d_cs_y, d_cs_x
# =============================================================================

# ----------------------------
# Helpers
# ----------------------------
make_cs_averages <- function(df, y, xvars, t = "year"){
  df %>%
    dplyr::group_by(.data[[t]]) %>%
    dplyr::summarise(
      cs_y = mean(.data[[y]], na.rm = TRUE),
      dplyr::across(dplyr::all_of(xvars), ~ mean(.x, na.rm = TRUE), .names = "cs_{.col}"),
      .groups = "drop"
    ) %>%
    dplyr::rename(year = .data[[t]])   # asegura columna year
}

make_diffs <- function(df, vars, id = "iso3c", t = "year"){
  df %>%
    dplyr::group_by(.data[[id]]) %>%
    dplyr::arrange(.data[[t]], .by_group = TRUE) %>%
    dplyr::mutate(dplyr::across(dplyr::all_of(vars), ~ .x - dplyr::lag(.x), .names = "d_{.col}")) %>%
    dplyr::ungroup()
}

estimate_ecm_country <- function(df_i, y, xvars, add_cs = FALSE, cs_lvl = NULL, cs_dif = NULL,
                                 include_dy_l1 = TRUE){
  df_i <- df_i %>% dplyr::arrange(year)
  
  df_i2 <- make_diffs(df_i, vars = c(y, xvars), id = "iso3c", t = "year") %>%
    dplyr::mutate(
      y_l1  = dplyr::lag(.data[[y]], 1),
      dy_l1 = dplyr::lag(.data[[paste0("d_", y)]], 1)
    )
  
  for (x in xvars) {
    df_i2[[paste0(x, "_l1")]] <- dplyr::lag(df_i2[[x]], 1)
  }
  
  if (add_cs) {
    if (is.null(cs_lvl) || is.null(cs_dif)) stop("CS-ECM requiere cs_lvl y cs_dif.")
    df_i2 <- df_i2 %>%
      dplyr::left_join(cs_lvl, by = "year") %>%
      dplyr::left_join(cs_dif, by = "year")
  }
  
  rhs_terms <- c(
    "y_l1",
    paste0(xvars, "_l1"),
    paste0("d_", xvars)
  )
  
  if (include_dy_l1) rhs_terms <- c("dy_l1", rhs_terms)
  
  if (add_cs) {
    rhs_terms <- c(rhs_terms,
                   "cs_y", paste0("cs_", xvars),
                   "d_cs_y", paste0("d_cs_", xvars))
  }
  
  fml <- stats::as.formula(paste0("d_", y, " ~ ", paste(rhs_terms, collapse = " + ")))
  
  est_df <- df_i2 %>% tidyr::drop_na(dplyr::all_of(c(paste0("d_", y), rhs_terms)))
  if (nrow(est_df) < 12) return(NULL)
  
  fit <- stats::lm(fml, data = est_df)
  
  coefs <- stats::coef(fit)
  if (!("y_l1" %in% names(coefs))) return(NULL)
  
  phi <- unname(coefs["y_l1"])
  b_l1 <- sapply(xvars, function(x){
    nm <- paste0(x, "_l1")
    if (nm %in% names(coefs)) unname(coefs[nm]) else NA_real_
  })
  lr <- -b_l1 / phi
  
  list(fit = fit, phi = phi, b_l1 = b_l1, lr = lr, n = nrow(est_df))
}  
  

mg_from_country_models <- function(res_list, xvars, require_phi_negative = TRUE){
  ok <- purrr::keep(res_list, ~ !is.null(.x))
  
  if (length(ok) == 0) return(NULL)
  
  # Extrae
  countries <- names(ok)
  phi_i <- sapply(ok, `[[`, "phi")
  n_i   <- sapply(ok, `[[`, "n")
  
  # Opcional: exigir convergencia (phi<0)
  if (require_phi_negative) {
    keep_idx <- which(phi_i < 0 & is.finite(phi_i))
    ok <- ok[keep_idx]
    if (length(ok) == 0) return(NULL)
    countries <- names(ok)
    phi_i <- sapply(ok, `[[`, "phi")
    n_i   <- sapply(ok, `[[`, "n")
  }
  
  lr_mat <- do.call(rbind, lapply(ok, function(x) x$lr))
  colnames(lr_mat) <- paste0("LR_", xvars)
  
  lr_mg  <- colMeans(lr_mat, na.rm = TRUE)
  phi_mg <- mean(phi_i, na.rm = TRUE)
  
  country_table <- tibble::tibble(
    iso3c = countries,
    phi   = as.numeric(phi_i),
    n     = as.integer(n_i)
  ) %>% dplyr::bind_cols(as.data.frame(lr_mat))
  
  list(
    country_table = country_table,
    lr_mg = lr_mg,
    phi_mg = phi_mg
  )
}

# ----------------------------
# Variables base (niveles) para ECM
# ----------------------------
y_var <- "ln_CO2pc"

x_M1 <- c("ln_GDPpc_centered","ln_GDPpc_sq","ln_energy_pc")
x_M2 <- c("ln_GDPpc_centered","ln_GDPpc_sq","ln_energy_pc",
          "estabilidad_politica_vdem","gobernanza_vdem")
x_M3 <- c("ln_GDPpc_centered","ln_GDPpc_sq","ln_energy_pc",
          "estabilidad_politica_vdem","gobernanza_vdem",
          "logpib_estabilidad","energia_gobierno")

# ----------------------------
# CS averages (niveles) y diferencias para CS-ECM
# ----------------------------
cs_lvl_M2 <- make_cs_averages(panel_brics2, y = y_var, xvars = x_M2, t = "year")
cs_lvl_M3 <- make_cs_averages(panel_brics2, y = y_var, xvars = x_M3, t = "year")

cs_dif_M2 <- cs_lvl_M2 %>%
  dplyr::arrange(year) %>%
  dplyr::mutate(
    d_cs_y = cs_y - dplyr::lag(cs_y)
  ) %>%
  dplyr::mutate(
    dplyr::across(dplyr::all_of(paste0("cs_", x_M2)),
                  ~ .x - dplyr::lag(.x),
                  .names = "d_{.col}")
  ) %>%
  dplyr::select(year, d_cs_y, dplyr::all_of(paste0("d_cs_", x_M2)))

cs_dif_M3 <- cs_lvl_M3 %>%
  dplyr::arrange(year) %>%
  dplyr::mutate(
    d_cs_y = cs_y - dplyr::lag(cs_y)
  ) %>%
  dplyr::mutate(
    dplyr::across(dplyr::all_of(paste0("cs_", x_M3)),
                  ~ .x - dplyr::lag(.x),
                  .names = "d_{.col}")
  ) %>%
  dplyr::select(year, d_cs_y, dplyr::all_of(paste0("d_cs_", x_M3)))

# ----------------------------
# Estimar ECM por país (MG)
# ----------------------------
split_cty <- split(panel_brics2, panel_brics2$iso3c)

# M1: ECM-MG (sin CS)
ecm_M1 <- lapply(split_cty, estimate_ecm_country, y = y_var, xvars = x_M1, add_cs = FALSE)
mg_M1  <- mg_from_country_models(ecm_M1, x_M1, require_phi_negative = TRUE)

# M2/M3: CS-ECM-MG (MAIN para CD)
ecm_M2_cs <- lapply(split_cty, estimate_ecm_country, y = y_var, xvars = x_M2,
                    add_cs = TRUE, cs_lvl = cs_lvl_M2, cs_dif = cs_dif_M2)
mg_M2_cs  <- mg_from_country_models(ecm_M2_cs, x_M2, require_phi_negative = TRUE)

ecm_M3_cs <- lapply(split_cty, estimate_ecm_country, y = y_var, xvars = x_M3,
                    add_cs = TRUE, cs_lvl = cs_lvl_M3, cs_dif = cs_dif_M3)
mg_M3_cs  <- mg_from_country_models(ecm_M3_cs, x_M3, require_phi_negative = TRUE)

# ----------------------------
# Chequeos críticos
# ----------------------------
phi_check <- tibble::tibble(
  model    = c("M1_ECM_MG","M2_CS_ECM_MG","M3_CS_ECM_MG"),
  phi_mean = c(mg_M1$phi_mg, mg_M2_cs$phi_mg, mg_M3_cs$phi_mg)
)
print(phi_check)

# φ por país (para apéndice / diagnóstico)
country_phi <- dplyr::bind_rows(
  mg_M1$country_table    %>% dplyr::mutate(model="M1_ECM"),
  mg_M2_cs$country_table %>% dplyr::mutate(model="M2_CS_ECM"),
  mg_M3_cs$country_table %>% dplyr::mutate(model="M3_CS_ECM")
)

# Long-run MG (coeficientes de largo plazo)
lr_mg_table <- dplyr::bind_rows(
  tibble::tibble(term = names(mg_M1$lr_mg),    estimate = as.numeric(mg_M1$lr_mg),    model = "M1_ECM_MG"),
  tibble::tibble(term = names(mg_M2_cs$lr_mg), estimate = as.numeric(mg_M2_cs$lr_mg), model = "M2_CS_ECM_MG"),
  tibble::tibble(term = names(mg_M3_cs$lr_mg), estimate = as.numeric(mg_M3_cs$lr_mg), model = "M3_CS_ECM_MG")
)

# (ii) Sensibilidad a “CS augmentation”

#Con N=5, la CS-augmentation a veces “absorbe demasiado”. Necesitas verificar 
# que los coeficientes LR tengan signos razonables y que no estén explotando 
# por φ cercano a 0.
# =============================================================================
# =============================================================================
# BLOQUE 14B) DIAGNÓSTICOS ECM (CONVERGENCIA, OVERSHOOTING, EXCLUSIONES)
# =============================================================================

# Países incluidos por modelo (ya filtrados por phi<0 en mg_from_country_models)
incl_M1 <- mg_M1$country_table$iso3c
incl_M2 <- mg_M2_cs$country_table$iso3c
incl_M3 <- mg_M3_cs$country_table$iso3c

excl_M1 <- setdiff(names(split_cty), incl_M1)
excl_M2 <- setdiff(names(split_cty), incl_M2)
excl_M3 <- setdiff(names(split_cty), incl_M3)

exclusions <- tibble::tibble(
  model = c("M1_ECM_MG","M2_CS_ECM_MG","M3_CS_ECM_MG"),
  included = c(paste(incl_M1, collapse=", "), paste(incl_M2, collapse=", "), paste(incl_M3, collapse=", ")),
  excluded = c(paste(excl_M1, collapse=", "), paste(excl_M2, collapse=", "), paste(excl_M3, collapse=", "))
)
print(exclusions)
readr::write_csv(exclusions, file.path(dir_tables, "Table_ECM_exclusions_phi_rule.csv"))

# Diagnóstico phi (overshooting si phi < -1 en anual)
phi_diag <- country_phi %>%
  dplyr::transmute(
    model, iso3c, phi,
    overshoot = (phi < -1),
    abs_phi = abs(phi)
  ) %>%
  dplyr::arrange(model, dplyr::desc(abs_phi))

print(phi_diag)
readr::write_csv(phi_diag, file.path(dir_tables, "Table_ECM_phi_diagnostics.csv"))

phi_summary <- phi_diag %>%
  dplyr::group_by(model) %>%
  dplyr::summarise(
    n = dplyr::n(),
    phi_min = min(phi, na.rm=TRUE),
    phi_max = max(phi, na.rm=TRUE),
    share_overshoot = mean(overshoot, na.rm=TRUE),
    .groups="drop"
  )
print(phi_summary)
readr::write_csv(phi_summary, file.path(dir_tables, "Table_ECM_phi_summary.csv"))

# =============================================================================
# BLOQUE 14C) SENSIBILIDAD: CS-ECM-MG SIN Δy_{t-1} (para estabilizar φ)
# =============================================================================
ecm_M2_cs_nody <- lapply(split_cty, estimate_ecm_country, y = y_var, xvars = x_M2,
                         add_cs = TRUE, cs_lvl = cs_lvl_M2, cs_dif = cs_dif_M2,
                         include_dy_l1 = FALSE)
mg_M2_cs_nody  <- mg_from_country_models(ecm_M2_cs_nody, x_M2, require_phi_negative = TRUE)

ecm_M3_cs_nody <- lapply(split_cty, estimate_ecm_country, y = y_var, xvars = x_M3,
                         add_cs = TRUE, cs_lvl = cs_lvl_M3, cs_dif = cs_dif_M3,
                         include_dy_l1 = FALSE)
mg_M3_cs_nody  <- mg_from_country_models(ecm_M3_cs_nody, x_M3, require_phi_negative = TRUE)

phi_check_nody <- tibble::tibble(
  model    = c("M2_CS_ECM_MG_noDY","M3_CS_ECM_MG_noDY"),
  phi_mean = c(mg_M2_cs_nody$phi_mg, mg_M3_cs_nody$phi_mg),
  n_countries = c(nrow(mg_M2_cs_nody$country_table), nrow(mg_M3_cs_nody$country_table))
)
print(phi_check_nody)

# Guardar resultados para tablas finales
saveRDS(list(mg_M2_cs_nody=mg_M2_cs_nody, mg_M3_cs_nody=mg_M3_cs_nody,
             ecm_M2_cs_nody=ecm_M2_cs_nody, ecm_M3_cs_nody=ecm_M3_cs_nody),
        file.path(dir_clean, paste0("models_ECM_MG_noDY_", freeze_date, ".rds")))

# =============================================================================
# BLOQUE 14D) DIAGNÓSTICOS PARA noDY (debe ir después de phi_check_nody)
# =============================================================================

# Construir country_phi para noDY
country_phi_nody <- dplyr::bind_rows(
  mg_M2_cs_nody$country_table %>% dplyr::mutate(model="M2_CS_ECM_noDY"),
  mg_M3_cs_nody$country_table %>% dplyr::mutate(model="M3_CS_ECM_noDY")
)

# Exclusiones (debería ser vacío en ambos)
incl_M2_nody <- mg_M2_cs_nody$country_table$iso3c
incl_M3_nody <- mg_M3_cs_nody$country_table$iso3c

excl_M2_nody <- setdiff(names(split_cty), incl_M2_nody)
excl_M3_nody <- setdiff(names(split_cty), incl_M3_nody)

exclusions_nody <- tibble::tibble(
  model = c("M2_CS_ECM_noDY","M3_CS_ECM_noDY"),
  included = c(paste(incl_M2_nody, collapse=", "), paste(incl_M3_nody, collapse=", ")),
  excluded = c(paste(excl_M2_nody, collapse=", "), paste(excl_M3_nody, collapse=", "))
)
print(exclusions_nody)
readr::write_csv(exclusions_nody, file.path(dir_tables, "Table_ECM_exclusions_noDY.csv"))

# Diagnóstico phi (overshooting)
phi_diag_nody <- country_phi_nody %>%
  dplyr::transmute(
    model, iso3c, phi,
    overshoot = (phi < -1),
    abs_phi = abs(phi)
  ) %>%
  dplyr::arrange(model, dplyr::desc(abs_phi))

print(phi_diag_nody)
readr::write_csv(phi_diag_nody, file.path(dir_tables, "Table_ECM_phi_diagnostics_noDY.csv"))

phi_summary_nody <- phi_diag_nody %>%
  dplyr::group_by(model) %>%
  dplyr::summarise(
    n = dplyr::n(),
    phi_min = min(phi, na.rm=TRUE),
    phi_max = max(phi, na.rm=TRUE),
    share_overshoot = mean(overshoot, na.rm=TRUE),
    .groups="drop"
  )
print(phi_summary_nody)
readr::write_csv(phi_summary_nody, file.path(dir_tables, "Table_ECM_phi_summary_noDY.csv"))

# =============================================================================
# BLOQUE 14E) LONG-RUN MG PARA noDY (MAIN)
# =============================================================================

lr_mg_nody <- dplyr::bind_rows(
  tibble::tibble(term = names(mg_M2_cs_nody$lr_mg),
                 estimate = as.numeric(mg_M2_cs_nody$lr_mg),
                 model = "M2_CS_ECM_noDY"),
  tibble::tibble(term = names(mg_M3_cs_nody$lr_mg),
                 estimate = as.numeric(mg_M3_cs_nody$lr_mg),
                 model = "M3_CS_ECM_noDY")
)

print(lr_mg_nody)
readr::write_csv(lr_mg_nody, file.path(dir_tables, "Table_LongRun_MG_noDY.csv"))

# Guardar objetos noDY
saveRDS(list(ecm_M2_cs_nody=ecm_M2_cs_nody, mg_M2_cs_nody=mg_M2_cs_nody,
             ecm_M3_cs_nody=ecm_M3_cs_nody, mg_M3_cs_nody=mg_M3_cs_nody),
        file.path(dir_clean, paste0("models_ECM_MG_noDY_", freeze_date, ".rds")))

# =============================================================================
# BLOQUE 14F) EXTENDED ESTABLE: M3a y M3b (una interacción a la vez)
# =============================================================================

x_M3a <- c("ln_GDPpc_centered","ln_GDPpc_sq","ln_energy_pc",
           "estabilidad_politica_vdem","gobernanza_vdem",
           "logpib_estabilidad")

x_M3b <- c("ln_GDPpc_centered","ln_GDPpc_sq","ln_energy_pc",
           "estabilidad_politica_vdem","gobernanza_vdem",
           "energia_gobierno")

cs_lvl_M3a <- make_cs_averages(panel_brics2, y = y_var, xvars = x_M3a, t = "year")
cs_lvl_M3b <- make_cs_averages(panel_brics2, y = y_var, xvars = x_M3b, t = "year")

cs_dif_M3a <- cs_lvl_M3a %>%
  dplyr::arrange(year) %>%
  dplyr::mutate(d_cs_y = cs_y - dplyr::lag(cs_y)) %>%
  dplyr::mutate(dplyr::across(dplyr::all_of(paste0("cs_", x_M3a)),
                              ~ .x - dplyr::lag(.x), .names = "d_{.col}")) %>%
  dplyr::select(year, d_cs_y, dplyr::all_of(paste0("d_cs_", x_M3a)))

cs_dif_M3b <- cs_lvl_M3b %>%
  dplyr::arrange(year) %>%
  dplyr::mutate(d_cs_y = cs_y - dplyr::lag(cs_y)) %>%
  dplyr::mutate(dplyr::across(dplyr::all_of(paste0("cs_", x_M3b)),
                              ~ .x - dplyr::lag(.x), .names = "d_{.col}")) %>%
  dplyr::select(year, d_cs_y, dplyr::all_of(paste0("d_cs_", x_M3b)))

ecm_M3a_cs_nody <- lapply(split_cty, estimate_ecm_country, y = y_var, xvars = x_M3a,
                          add_cs = TRUE, cs_lvl = cs_lvl_M3a, cs_dif = cs_dif_M3a,
                          include_dy_l1 = FALSE)
mg_M3a_cs_nody  <- mg_from_country_models(ecm_M3a_cs_nody, x_M3a, require_phi_negative = TRUE)

ecm_M3b_cs_nody <- lapply(split_cty, estimate_ecm_country, y = y_var, xvars = x_M3b,
                          add_cs = TRUE, cs_lvl = cs_lvl_M3b, cs_dif = cs_dif_M3b,
                          include_dy_l1 = FALSE)
mg_M3b_cs_nody  <- mg_from_country_models(ecm_M3b_cs_nody, x_M3b, require_phi_negative = TRUE)

phi_check_M3ab <- tibble::tibble(
  model = c("M3a_CS_ECM_noDY","M3b_CS_ECM_noDY"),
  phi_mean = c(mg_M3a_cs_nody$phi_mg, mg_M3b_cs_nody$phi_mg),
  n_countries = c(nrow(mg_M3a_cs_nody$country_table), nrow(mg_M3b_cs_nody$country_table))
)
print(phi_check_M3ab)

# LR tables
lr_M3ab <- dplyr::bind_rows(
  tibble::tibble(term = names(mg_M3a_cs_nody$lr_mg), estimate = as.numeric(mg_M3a_cs_nody$lr_mg), model="M3a_CS_ECM_noDY"),
  tibble::tibble(term = names(mg_M3b_cs_nody$lr_mg), estimate = as.numeric(mg_M3b_cs_nody$lr_mg), model="M3b_CS_ECM_noDY")
)
print(lr_M3ab)
readr::write_csv(lr_M3ab, file.path(dir_tables, "Table_LongRun_MG_M3ab_noDY.csv"))

# ----------------------------
# A. Validez del turning point
# ----------------------------
tp_validity <- function(tp_gdp, b2, gdp_range){
  # gdp_range = c(min, max)
  if (length(gdp_range) != 2) return(FALSE)
  if (any(!is.finite(gdp_range))) return(FALSE)
  if (!is.finite(tp_gdp)) return(FALSE)
  if (is.na(b2)) return(FALSE)
  
  # Dentro del soporte observado + curvatura definida (b2 != 0)
  (tp_gdp >= gdp_range[1]) && (tp_gdp <= gdp_range[2]) && (b2 != 0)
}

# ----------------------------
# B. Clasificación de forma 
# ----------------------------
shape_classifier <- function(b1, b2){
  if (is.na(b1) || is.na(b2) || !is.finite(b1) || !is.finite(b2)) return(NA_character_)
  if (b2 < 0 && b1 > 0) return("EKC (inverted-U): b1>0, b2<0")
  if (b2 > 0 && b1 < 0) return("U-shape (recarbonization): b1<0, b2>0")
  if (b2 < 0 && b1 < 0) return("Monotonic decreasing (concave): b1<0, b2<0")
  if (b2 > 0 && b1 > 0) return("Monotonic increasing (convex): b1>0, b2>0")
  return("Undefined / borderline")
}

# ----------------------------
# C. Turning point desde coeficientes (tu fórmula)
# ----------------------------
turning_point_from_b <- function(b1, b2, mean_lnGDP){
  # TP* en x centrado = -b1/(2*b2); luego descentrar y exp()
  # Manejo de casos degenerados
  if (is.na(b1) || is.na(b2) || !is.finite(b1) || !is.finite(b2) || b2 == 0) return(NA_real_)
  tp_x  <- -b1 / (2*b2)
  ln_tp <- tp_x + mean_lnGDP
  exp(ln_tp)
}

# ----------------------------
# D. Rangos de soporte (panel vs por país)
# ----------------------------
gdp_support_panel <- function(df, gdp_var = "GDP_pc_filled"){
  rng <- range(df[[gdp_var]], na.rm = TRUE)
  if (any(!is.finite(rng))) return(c(NA_real_, NA_real_))
  rng
}

gdp_support_by_country <- function(df, id = "iso3c", gdp_var = "GDP_pc_filled"){
  df %>%
    dplyr::group_by(.data[[id]]) %>%
    dplyr::summarise(
      gdp_min = min(.data[[gdp_var]], na.rm = TRUE),
      gdp_max = max(.data[[gdp_var]], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      gdp_min = ifelse(is.finite(gdp_min), gdp_min, NA_real_),
      gdp_max = ifelse(is.finite(gdp_max), gdp_max, NA_real_)
    )
}

# =============================================================================
# BLOQUE 15) TURNING POINTS + SOPORTE + BOOTSTRAP MG (ilustrativo)
# =============================================================================
# NOTA METODOLÓGICA (para el paper):
# - Los turning points basados en FE estático en niveles son ILUSTRATIVOS,
#   porque el modelo main es CS-ECM. En arbitraje, el TP relevante debe derivarse
#   de coeficientes de largo plazo (LR) del CS-ECM.
# - Aquí reportamos:
#   (A) TP "ilustrativo" por modelos FE (M1–M3), con criterio dentro del soporte.
#   (B) TP por país (MG) con bootstrap paramétrico MVN para IC, también ilustrativo.
# - Criterio EKC (inverted-U) en niveles:
#     b2 < 0 y TP dentro del rango observado de PIB pc.
#   (si b2 > 0 => U-shape / recarbonización; TP no aplica para EKC estándar).

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(MASS)
})

# ----------------------------
# 15.1 Helpers robustos
# ----------------------------

# Devuelve un vector PIBpc observado para soporte
get_gdp_pc_series <- function(df){
  if ("GDP_pc_filled" %in% names(df)) return(df$GDP_pc_filled)
  if ("GDP_pc" %in% names(df)) return(df$GDP_pc)
  # fallback: si no hay GDP_pc pero sí ln_GDPpc
  if ("ln_GDPpc" %in% names(df)) return(exp(df$ln_GDPpc))
  stop("No encuentro GDP_pc_filled ni GDP_pc ni ln_GDPpc para construir el soporte del PIB pc.")
}

# Turning point en PIB pc:
# modelo: y = ... + b1*x + b2*x^2, con x = ln_GDPpc_centered
# TP en x: -b1/(2*b2)
# volver a lnGDP: lnGDP_TP = x_TP + mean(ln_GDPpc)
# PIBpc_TP = exp(lnGDP_TP)
turning_point_from_b <- function(b1, b2, mean_lnGDP){
  if (is.na(b1) || is.na(b2)) return(NA_real_)
  if (abs(b2) < 1e-10) return(NA_real_)        # evita explosiones numéricas
  tp_x  <- -b1 / (2*b2)
  ln_tp <- tp_x + mean_lnGDP
  exp(ln_tp)
}

# Evalúa criterio EKC dentro del soporte
tp_validity <- function(tp, b2, gdp_range){
  is.finite(tp) && !is.na(tp) &&
    (tp >= gdp_range[1]) && (tp <= gdp_range[2]) &&
    (!is.na(b2) && b2 < 0)
}

# ----------------------------
# 15.2 Soporte observado y media lnGDP (pool)
# ----------------------------
gdp_series <- get_gdp_pc_series(panel_brics2)
gdp_range  <- range(gdp_series, na.rm = TRUE)

if (!("ln_GDPpc" %in% names(panel_brics2))) {
  # si no existe, intenta reconstruir desde GDP_pc
  if ("GDP_pc" %in% names(panel_brics2)) {
    panel_brics2 <- panel_brics2 %>% mutate(ln_GDPpc = log(GDP_pc))
  } else if ("GDP_pc_filled" %in% names(panel_brics2)) {
    panel_brics2 <- panel_brics2 %>% mutate(ln_GDPpc = log(GDP_pc_filled))
  } else {
    stop("No existe ln_GDPpc y no puedo reconstruirla (falta GDP_pc).")
  }
}

mean_lnGDP_pool <- mean(panel_brics2$ln_GDPpc, na.rm = TRUE)

# ----------------------------
# 15.3 (A) TP panel: FE estático (ilustrativo) para M1–M3
# ----------------------------
tp_panel_models <- function(mod, mod_name, mean_lnGDP, gdp_range){
  b <- tryCatch(coef(mod), error = function(e) NULL)
  if (is.null(b)) return(NULL)
  
  req <- c("ln_GDPpc_centered","ln_GDPpc_sq")
  if (!all(req %in% names(b))) return(NULL)
  
  b1 <- unname(b["ln_GDPpc_centered"])
  b2 <- unname(b["ln_GDPpc_sq"])
  
  # TP en GDPpc (nivel) + TP en lnGDP centrado (para auditoría)
  tp_lnGDP_centered <- ifelse(is.na(b1) || is.na(b2) || b2 == 0, NA_real_, -b1/(2*b2))
  tp_gdp <- turning_point_from_b(b1, b2, mean_lnGDP)
  
  tibble::tibble(
    model = mod_name,
    b1 = b1,
    b2 = b2,
    tp_lnGDP_centered = tp_lnGDP_centered,
    tp_gdp_pc = tp_gdp,
    within_support = tp_validity(tp_gdp, b2, gdp_range),
    ekc_strict = (is.finite(b1) && is.finite(b2) && b1 > 0 && b2 < 0),
    u_strict   = (is.finite(b1) && is.finite(b2) && b1 < 0 && b2 > 0),
    shape = shape_classifier(b1, b2)
  )
}

mean_lnGDP_pool <- mean(panel_brics2$ln_GDPpc, na.rm = TRUE)
gdp_range_panel <- gdp_support_panel(panel_brics2, gdp_var = "GDP_pc_filled")

tp_panel_table <- dplyr::bind_rows(
  if (exists("M1_fe")) tp_panel_models(M1_fe, "M1_FE", mean_lnGDP_pool, gdp_range_panel) else NULL,
  if (exists("M2_fe")) tp_panel_models(M2_fe, "M2_FE", mean_lnGDP_pool, gdp_range_panel) else NULL,
  if (exists("M3_fe")) tp_panel_models(M3_fe, "M3_FE", mean_lnGDP_pool, gdp_range_panel) else NULL
)


readr::write_csv(tp_panel_table, file.path(dir_tables, "Table_TurningPoint_panel_FE_illustrative.csv"))
print(tp_panel_table)

# ----------------------------
# 15.4 (B) TP por país (MG): bootstrap paramétrico MVN (ilustrativo)
# ----------------------------
# Regresión por país (estática) para TP + IC.
# Importante: Si b2 >= 0, NO hay EKC inverted-U -> TP no aplica.
# Aun así reportamos TP_hat pero marcamos within_support=FALSE.
# =============================================================================
# BLOQUE 15B) TURNING POINTS POR PAÍS (BOOTSTRAP) + SOPORTE POR PAÍS (FIX DEFINITIVO)
# =============================================================================

# ----------------------------
# Helpers
# ----------------------------

# Devuelve serie PIBpc observada para soporte
get_gdp_pc_series <- function(df){
  if ("GDP_pc_filled" %in% names(df)) return(df$GDP_pc_filled)
  if ("GDP_pc" %in% names(df)) return(df$GDP_pc)
  if ("ln_GDPpc" %in% names(df)) return(exp(df$ln_GDPpc))
  stop("No encuentro GDP_pc_filled ni GDP_pc ni ln_GDPpc para construir soporte del PIB pc.")
}

# Turning point en PIBpc:
# x = ln_GDPpc_centered, TP en x: -b1/(2*b2)
# volver a lnGDP: lnGDP_TP = x_TP + mean(ln_GDPpc)  => PIBpc_TP = exp(lnGDP_TP)
turning_point_from_b <- function(b1, b2, mean_lnGDP){
  if (is.na(b1) || is.na(b2)) return(NA_real_)
  if (!is.finite(b1) || !is.finite(b2)) return(NA_real_)
  if (abs(b2) < 1e-10) return(NA_real_)  # evita explosión numérica
  tp_x  <- -b1 / (2*b2)
  ln_tp <- tp_x + mean_lnGDP
  exp(ln_tp)
}

# Soporte PIBpc por país
gdp_support_by_country <- function(df, id = "iso3c", gdp_var = "GDP_pc_filled"){
  gdp_series <- if (gdp_var %in% names(df)) df[[gdp_var]] else get_gdp_pc_series(df)
  df2 <- df
  df2$.__gdp__ <- gdp_series
  
  df2 %>%
    dplyr::group_by(.data[[id]]) %>%
    dplyr::summarise(
      gdp_min = min(.__gdp__, na.rm = TRUE),
      gdp_max = max(.__gdp__, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::rename(iso3c = .data[[id]])
}

# ----------------------------
# Bootstrap TP por país
# ----------------------------
boot_tp_country <- function(df_i, B = 800){
  needed <- c("ln_CO2pc","ln_GDPpc_centered","ln_GDPpc_sq","ln_energy_pc")
  miss <- setdiff(needed, names(df_i))
  if (length(miss) > 0) {
    stop("Faltan variables en df_i para TP bootstrap: ", paste(miss, collapse=", "))
  }
  
  # Asegurar ln_GDPpc para mean_lnGDP_i (si no existe, reconstruir)
  if (!("ln_GDPpc" %in% names(df_i))) {
    if ("GDP_pc_filled" %in% names(df_i)) {
      df_i <- df_i %>% dplyr::mutate(ln_GDPpc = log(GDP_pc_filled))
    } else if ("GDP_pc" %in% names(df_i)) {
      df_i <- df_i %>% dplyr::mutate(ln_GDPpc = log(GDP_pc))
    } else {
      stop("No existe ln_GDPpc y no puedo reconstruirla (falta GDP_pc/GDP_pc_filled).")
    }
  }
  
  fit <- stats::lm(ln_CO2pc ~ ln_GDPpc_centered + ln_GDPpc_sq + ln_energy_pc, data = df_i)
  
  # Coeficientes
  bhat_all <- stats::coef(fit)
  if (!all(c("ln_GDPpc_centered","ln_GDPpc_sq") %in% names(bhat_all))) {
    return(tibble::tibble(
      tp_hat = NA_real_, tp_p025 = NA_real_, tp_p975 = NA_real_,
      b1 = NA_real_, b2 = NA_real_, n = nrow(df_i),
      note = "Missing b1/b2 in fit"
    ))
  }
  
  bhat <- bhat_all[c("ln_GDPpc_centered","ln_GDPpc_sq")]
  V <- tryCatch(
    stats::vcov(fit)[c("ln_GDPpc_centered","ln_GDPpc_sq"), c("ln_GDPpc_centered","ln_GDPpc_sq")],
    error = function(e) NULL
  )
  
  if (any(!is.finite(bhat)) || is.null(V) || any(!is.finite(V))) {
    return(tibble::tibble(
      tp_hat = NA_real_, tp_p025 = NA_real_, tp_p975 = NA_real_,
      b1 = unname(bhat[1]), b2 = unname(bhat[2]), n = nrow(df_i),
      note = "Non-finite coefficients/variance"
    ))
  }
  
  mean_lnGDP_i <- mean(df_i$ln_GDPpc, na.rm = TRUE)
  tp_hat <- turning_point_from_b(bhat[1], bhat[2], mean_lnGDP = mean_lnGDP_i)
  
  draws <- tryCatch(
    MASS::mvrnorm(n = B, mu = as.numeric(bhat), Sigma = V),
    error = function(e) NULL
  )
  if (is.null(draws)) {
    return(tibble::tibble(
      tp_hat = tp_hat, tp_p025 = NA_real_, tp_p975 = NA_real_,
      b1 = unname(bhat[1]), b2 = unname(bhat[2]), n = nrow(df_i),
      note = "MVN draw failed"
    ))
  }
  colnames(draws) <- c("b1","b2")
  
  tp_draw <- apply(draws, 1, function(z){
    turning_point_from_b(z[1], z[2], mean_lnGDP = mean_lnGDP_i)
  })
  tp_draw <- tp_draw[is.finite(tp_draw)]
  
  if (length(tp_draw) < 30) {
    return(tibble::tibble(
      tp_hat = tp_hat, tp_p025 = NA_real_, tp_p975 = NA_real_,
      b1 = unname(bhat[1]), b2 = unname(bhat[2]), n = nrow(df_i),
      note = "Too few finite TP draws"
    ))
  }
  
  tibble::tibble(
    tp_hat  = tp_hat,
    tp_p025 = as.numeric(stats::quantile(tp_draw, 0.025, na.rm = TRUE)),
    tp_p975 = as.numeric(stats::quantile(tp_draw, 0.975, na.rm = TRUE)),
    b1 = unname(bhat[1]),
    b2 = unname(bhat[2]),
    n  = nrow(df_i),
    note = NA_character_
  )
}

# ----------------------------
# Construir tp_by_country DESDE CERO (no depende del estado previo)
# ----------------------------
tp_by_country <- panel_brics2 %>%
  dplyr::group_by(iso3c) %>%
  dplyr::group_modify(~ boot_tp_country(.x, B = 800)) %>%
  dplyr::ungroup()

# ----------------------------
# Soporte PIBpc por país + bandera within_support_country
# ----------------------------
gdp_country_support <- gdp_support_by_country(
  panel_brics2,
  id = "iso3c",
  gdp_var = "GDP_pc_filled"
)

# Chequeo rápido
stopifnot(all(c("iso3c","gdp_min","gdp_max") %in% names(gdp_country_support)))

tp_by_country <- tp_by_country %>%
  dplyr::left_join(gdp_country_support, by = "iso3c")

# Verificación inmediata
stopifnot(all(c("gdp_min","gdp_max") %in% names(tp_by_country)))

tp_by_country <- tp_by_country %>%
  dplyr::mutate(
    tol = 0.02 * (gdp_max - gdp_min),
    within_support_country =
      is.finite(tp_hat) &
      (tp_hat >= (gdp_min - tol)) &
      (tp_hat <= (gdp_max + tol)) &
      is.finite(b2) & (b2 < 0)
  ) %>%
  dplyr::select(-tol)

# ----------------------------
# Output
# ----------------------------
print(
  tp_by_country %>%
    dplyr::select(
      iso3c, tp_hat, gdp_min, gdp_max, b2,
      within_support_country, note
    )
)

readr::write_csv(tp_by_country, file.path(dir_tables, "Table_TurningPoints_byCountry_bootstrap.csv"))
message("OK: Table_TurningPoints_byCountry_bootstrap.csv guardada en: ", dir_tables)

tp_by_country <- tp_by_country %>%
  dplyr::mutate(
    within_support_panel =
      is.finite(tp_hat) &
      (tp_hat >= min(gdp_min, na.rm = TRUE)) &
      (tp_hat <= max(gdp_max, na.rm = TRUE)) &
      is.finite(b2) & (b2 < 0)
  )
print(
  tp_by_country %>%
    dplyr::select(
      iso3c, tp_hat, gdp_min, gdp_max, b2,
      within_support_panel, within_support_country, note
    )
)

# ----------------------------
# 15.5 (Opcional) Turning point "long-run" desde CS-ECM (si quieres)
# ----------------------------
# OJO: Solo es válido si tu especificación LR realmente es cuadrática en lnGDP_centered.
# Para CS-ECM, el TP debería usar los LR betas del modelo main (mg_M2_cs_nody o mg_M3a...).
# En tus resultados actuales LR b2 suele ser >0, así que EKC TP probablemente NO aplica.
# Lo dejamos como plantilla por si en otra especificación obtienes b2<0.

calc_tp_from_lr <- function(lr_b1, lr_b2, mean_lnGDP){
  turning_point_from_b(lr_b1, lr_b2, mean_lnGDP)
}

# Ejemplo (solo si tienes LR para ln_GDPpc_centered y ln_GDPpc_sq):
if (exists("mg_M2_cs_nody") &&
    all(c("LR_ln_GDPpc_centered","LR_ln_GDPpc_sq") %in% names(mg_M2_cs_nody$lr_mg))) {
  
  lr_b1 <- unname(mg_M2_cs_nody$lr_mg["LR_ln_GDPpc_centered"])
  lr_b2 <- unname(mg_M2_cs_nody$lr_mg["LR_ln_GDPpc_sq"])
  
  tp_lr <- calc_tp_from_lr(lr_b1, lr_b2, mean_lnGDP_pool)
  tp_lr_valid <- tp_validity(tp_lr, lr_b2, gdp_range)
  
  tp_lr_row <- tibble::tibble(
    model = "M2_CS_ECM_noDY_LR",
    lr_b1 = lr_b1,
    lr_b2 = lr_b2,
    tp_gdp_pc = tp_lr,
    within_support = tp_lr_valid,
    shape = case_when(
      is.na(lr_b2) ~ NA_character_,
      lr_b2 < 0 ~ "Inverted-U candidate (EKC)",
      lr_b2 > 0 ~ "U-shape / recarbonization",
      TRUE ~ "Undefined"
    )
  )
  
  readr::write_csv(tp_lr_row, file.path(dir_tables, "Table_TurningPoint_LR_CS_ECM.csv"))
  print(tp_lr_row)
}

message("BLOQUE 15 listo. Archivos guardados en: ", dir_tables)

# =============================================================================
# BLOQUE: TABLAS WORD (Journal-ready)
# =============================================================================
# =============================================================================
# BLOQUE: TABLES 1–3 (Main text) — Word, journal-ready
# - Table 1: Descriptive statistics
# - Table 2: Main results (CCE-MG, M1–M3) con asteriscos
# - Table 3: Turning point implied by CCE-MG (baseline EKC: M2) con notas
#
# Requisitos:
# - Usa tus objetos/archivos ya creados:
#     * panel_brics2  (data frame)  [si no existe, intenta cargar el último RDS]
#     * tables/Table_CCE_MG_Coefficients.csv  (ccemg_table)
# - Exporta 3 .docx a /tables
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(officer)
  library(flextable)
})

# ----------------------------
# 0) Rutas 
# ----------------------------
if (!exists("dir_tables")) dir_tables <- file.path(getwd(), "tables")
if (!exists("dir_clean"))  dir_clean  <- file.path(getwd(), "data_clean")
if (!dir.exists(dir_tables)) dir.create(dir_tables, recursive = TRUE, showWarnings = FALSE)

# ----------------------------
# 1) Helpers
# ----------------------------
fmt_num <- function(x, d = 3){
  ifelse(is.na(x), NA_character_, formatC(x, format = "f", digits = d))
}

stars_from_p <- function(p){
  # Vectorizado
  out <- rep("", length(p))
  out[!is.na(p) & p < 0.10] <- "*"
  out[!is.na(p) & p < 0.05] <- "**"
  out[!is.na(p) & p < 0.01] <- "***"
  out
}

to_cell <- function(est, se, p, d = 3){
  # Vectorizado; evita el error "condition has length > 1"
  est_s <- fmt_num(est, d)
  se_s  <- fmt_num(se,  d)
  st    <- stars_from_p(p)
  ifelse(
    is.na(est),
    NA_character_,
    ifelse(
      is.na(se),
      paste0(est_s, st),
      paste0(est_s, st, " (", se_s, ")")
    )
  )
}

make_wide_cells <- function(df_long, model_order){
  # df_long: columnas model, term, estimate, std_error, p_value (p_value puede ser NA)
  df_long %>%
    mutate(
      model = factor(model, levels = model_order),
      cell  = to_cell(estimate, std_error, p_value, d = 3)
    ) %>%
    select(model, term, cell) %>%
    tidyr::pivot_wider(names_from = model, values_from = cell) %>%
    arrange(term)
}

pretty_term <- function(x){
  # Etiquetas journal-friendly (ajusta si quieres)
  x %>%
    str_replace_all("^LR_", "") %>%
    str_replace_all("ln_GDPpc_centered", "ln GDP pc (centered)") %>%
    str_replace_all("ln_GDPpc_sq",       "ln GDP pc squared") %>%
    str_replace_all("ln_energy_pc",      "ln Energy pc") %>%
    str_replace_all("estabilidad_politica_vdem", "Political stability (V-Dem)") %>%
    str_replace_all("gobernanza_vdem",           "Governance (V-Dem)") %>%
    str_replace_all("logpib_estabilidad",        "ln GDP × Stability") %>%
    str_replace_all("energia_gobierno",          "ln Energy × Governance")
}

add_footer_lastrow <- function(ft, footer_text){
  # Inserta una última fila "Source/Notes/Created..." dentro de la tabla
  stopifnot(is.character(footer_text), length(footer_text) == 1)
  
  cn <- colnames(ft$body$dataset)
  footer_vals <- as.list(setNames(rep("", length(cn)), cn))
  footer_vals[[cn[1]]] <- footer_text
  
  ft2 <- flextable::add_body_row(ft, values = footer_vals, top = FALSE)
  ft2 <- flextable::merge_at(ft2, i = nrow(ft2$body$dataset), j = 1:length(cn))
  ft2 <- flextable::align(ft2, i = nrow(ft2$body$dataset), align = "left", part = "body")
  ft2 <- flextable::fontsize(ft2, i = nrow(ft2$body$dataset), size = 9, part = "body")
  ft2 <- flextable::italic(ft2, i = nrow(ft2$body$dataset), part = "body")
  ft2
}

journal_ft <- function(df){
  flextable(df) %>%
    theme_booktabs() %>%
    autofit()
}

# ----------------------------
# 2) Asegurar panel_brics2 (para Table 1 y Table 3 mean lnGDP)
# ----------------------------
if (!exists("panel_brics2")) {
  # intenta cargar el último panel guardado en dir_clean
  rds_candidates <- list.files(dir_clean, pattern = "^panel_BRICS_macro_vdem_.*\\.rds$", full.names = TRUE)
  if (length(rds_candidates) == 0) {
    stop("No existe 'panel_brics2' en memoria y no encontré panel_BRICS_macro_vdem_*.rds en /data_clean.")
  }
  rds_candidates <- rds_candidates[order(file.info(rds_candidates)$mtime, decreasing = TRUE)]
  panel_brics2 <- readRDS(rds_candidates[1])
  message("Cargué panel_brics2 desde: ", rds_candidates[1])
}

# reconstruir ln_GDPpc si hiciera falta (para mean_lnGDP_pool)
if (!("ln_GDPpc" %in% names(panel_brics2))) {
  if ("GDP_pc_filled" %in% names(panel_brics2)) {
    panel_brics2 <- panel_brics2 %>% mutate(ln_GDPpc = log(GDP_pc_filled))
  } else if ("GDP_pc" %in% names(panel_brics2)) {
    panel_brics2 <- panel_brics2 %>% mutate(ln_GDPpc = log(GDP_pc))
  } else {
    stop("No encuentro ln_GDPpc ni GDP_pc_filled/GDP_pc para reconstruir ln_GDPpc.")
  }
}

# ----------------------------
# TABLE 1 — Descriptive Statistics (FIX)
# ----------------------------
desc_vars <- c(
  "ln_CO2pc",
  "ln_GDPpc",
  "ln_energy_pc",
  "estabilidad_politica_vdem",
  "gobernanza_vdem"
)

missing_desc <- setdiff(desc_vars, names(panel_brics2))
if (length(missing_desc) > 0) {
  stop("Faltan variables para Table 1: ", paste(missing_desc, collapse = ", "))
}

tab1 <- panel_brics2 %>%
  dplyr::select(dplyr::all_of(desc_vars)) %>%
  tidyr::pivot_longer(
    cols = everything(),
    names_to = "variable",
    values_to = "value"
  ) %>%
  dplyr::group_by(variable) %>%
  dplyr::summarise(
    Mean = mean(value, na.rm = TRUE),
    SD   = sd(value, na.rm = TRUE),
    Min  = min(value, na.rm = TRUE),
    Max  = max(value, na.rm = TRUE),
    Obs  = sum(!is.na(value)),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    variable = pretty_term(variable),
    Mean = fmt_num(Mean, 3),
    SD   = fmt_num(SD, 3),
    Min  = fmt_num(Min, 3),
    Max  = fmt_num(Max, 3)
  ) %>%
  dplyr::rename(Variable = variable)

ft1 <- journal_ft(tab1)
ft1 <- add_footer_lastrow(
  ft1,
  paste0(
    "Source: Authors’ calculations using World Development Indicators (WDI) and V-Dem v15. ",
    "Notes: Descriptive statistics for BRICS countries over 1990–2024. ",
    "Created: ", format(Sys.time(), "%Y-%m-%d %H:%M"), "."
  )
)

doc1 <- read_docx() %>%
  body_add_par(
    "Table 1. Descriptive Statistics (BRICS, 1990–2024)",
    style = "heading 1"
  ) %>%
  body_add_flextable(ft1)

print(doc1, target = file.path(dir_tables, "Table_1_Descriptive_Statistics.docx"))

# =============================================================================
# TABLE 2 — Main Results (CCE-MG, M1–M3) — Word 
# Requiere: Table_CCE_MG_Coefficients.csv ya creado en /tables
# Output: tables/Table_2_Main_Results_CCE_MG.docx
# =============================================================================

# 0) Checks de rutas/insumos
# ----------------------------
if (!exists("dir_tables")) stop("No existe dir_tables. Define dir_tables <- file.path(getwd(),'tables').")
if (!dir.exists(dir_tables)) dir.create(dir_tables, recursive = TRUE)

ccemg_csv <- file.path(dir_tables, "Table_CCE_MG_Coefficients.csv")
if (!file.exists(ccemg_csv)) stop("No encuentro: ", ccemg_csv, "\nPrimero genera/guarda Table_CCE_MG_Coefficients.csv.")

# ----------------------------
# 1) Helpers: formato, asteriscos, etiquetas
# ----------------------------
fmt_num <- function(x, d = 3){
  ifelse(is.na(x) | !is.finite(x), NA_character_, formatC(x, format = "f", digits = d))
}

sig_stars <- function(p){
  ifelse(is.na(p), "",
         ifelse(p < 0.01, "***",
                ifelse(p < 0.05, "**",
                       ifelse(p < 0.10, "*", ""))))
}

to_cell <- function(est, se, p, d = 3){
  # vectorizado (evita error: "condition has length > 1")
  est_chr <- fmt_num(est, d)
  se_chr  <- fmt_num(se,  d)
  stars   <- sig_stars(p)
  out <- ifelse(is.na(est_chr), NA_character_, paste0(est_chr, stars, "\n(", se_chr, ")"))
  out
}

pretty_term <- function(x){
  # Ajusta nombres a estilo paper (edítalo a tu gusto)
  x <- as.character(x)
  x <- str_replace_all(x, "^LR_", "")
  dplyr::recode(
    x,
    "ln_GDPpc_centered" = "ln(GDPpc) (centered)",
    "ln_GDPpc_sq"       = "[ln(GDPpc) (centered)]^2",
    "ln_energy_pc"      = "ln(Energy use per capita)",
    "estabilidad_politica_vdem" = "Political stability (V-Dem index)",
    "gobernanza_vdem"           = "Governance (V-Dem index)",
    "logpib_estabilidad"        = "ln(GDPpc) × Stability (within)",
    "energia_gobierno"          = "ln(Energy) × Governance (within)",
    .default = x
  )
}

make_wide_cells <- function(df_long, model_order){
  df_long %>%
    dplyr::mutate(
      model = factor(model, levels = model_order),
      cell  = to_cell(estimate, std_error, p_value, d = 3)
    ) %>%
    dplyr::select(model, term, cell) %>%
    tidyr::pivot_wider(names_from = model, values_from = cell) %>%
    dplyr::arrange(term)
}

# ----------------------------
# 2) Cargar y estandarizar ccemg_table (desde CSV)
# ----------------------------
ccemg_long <- readr::read_csv(ccemg_csv, show_col_types = FALSE)

# Normaliza nombres posibles
# Esperado: model, term, estimate, std_error, statistic, p_value
need_cols <- c("model","term","estimate","std_error","p_value")
miss_cols <- setdiff(need_cols, names(ccemg_long))
if (length(miss_cols) > 0) stop("Faltan columnas en CCE-MG CSV: ", paste(miss_cols, collapse = ", "))

ccemg_long <- ccemg_long %>%
  dplyr::mutate(
    estimate  = suppressWarnings(as.numeric(estimate)),
    std_error = suppressWarnings(as.numeric(std_error)),
    p_value   = suppressWarnings(as.numeric(p_value)),
    term      = as.character(term),
    model     = as.character(model)
  )

# Mantén SOLO M1–M3 en orden
model_order <- c("M1_CCEMG","M2_CCEMG","M3_CCEMG")
ccemg_long <- ccemg_long %>%
  dplyr::filter(model %in% model_order)

# ----------------------------
# 3) Tabla en formato ancho (journal-style: coef + (SE) con asteriscos)
# ----------------------------
tab2_wide <- ccemg_long %>%
  dplyr::mutate(term = pretty_term(term)) %>%
  make_wide_cells(model_order = model_order) %>%
  dplyr::rename(
    Term = term,
    `M1` = M1_CCEMG,
    `M2` = M2_CCEMG,
    `M3` = M3_CCEMG
  )

# ----------------------------
# 4) Flextable (Word) + Footer: Source/Notes + leyenda asteriscos
# ----------------------------
ft2 <- flextable::flextable(tab2_wide)

ft2 <- ft2 %>%
  flextable::autofit() %>%
  flextable::fontsize(size = 10, part = "all") %>%
  flextable::fontsize(size = 11, part = "header") %>%
  flextable::bold(part = "header") %>%
  flextable::align(align = "left",  j = 1, part = "all") %>%
  flextable::align(align = "center", j = 2:ncol(tab2_wide), part = "all") %>%
  flextable::valign(valign = "top", part = "body") %>%
  flextable::line_spacing(space = 1.0, part = "all")

# (Opcional) ancho uniforme de modelos (ajusta a tu gusto)
# ft2 <- flextable::width(ft2, j = 1, width = 2.4)
# ft2 <- flextable::width(ft2, j = 2:ncol(tab2_wide), width = 1.6)

# Nota/leyenda como “última fila” visual (pie de tabla)
note_lines <- c(
  "Source: Authors’ calculations based on WDI (World Bank), V-Dem v15, IMF Datamapper (PPPPC, proxy extension for RUS 2022–2024 when needed), and OWID (proxy extensions for missing energy/CO2 when needed).",
  "Notes: Long-run coefficients computed as θ_k^LR = -β_k/φ for each country, then averaged (CCE Mean Group). Each cell shows coefficient with standard error in parentheses. Short-run dynamics available in Appendix. Significance: *** p<0.01, ** p<0.05, * p<0.10."
)

# añadimos como líneas debajo de la tabla (más estándar Q1 que “fila” dentro del cuerpo)
# Si INSISTES en que sea fila dentro de la tabla, dímelo y lo adapto.
ft2 <- flextable::add_footer_lines(ft2, values = note_lines)
ft2 <- flextable::fontsize(ft2, part = "footer", size = 9)
ft2 <- flextable::align(ft2, part = "footer", align = "left")

# ----------------------------
# 5) Exportar Word (FIX)
# ----------------------------
doc <- officer::read_docx()

doc <- doc %>%
  officer::body_add_par(
    "Table 2. Main Results (CCE-MG, BRICS 1990–2024)",
    style = "heading 1"
  ) %>%
  flextable::body_add_flextable(ft2)

out_path <- file.path(dir_tables, "Table_2_Main_Results_CCE_MG.docx")
print(doc, target = out_path)

message("OK: Table 2 exportada a -> ", out_path)

# ----------------------------
# ----------------------------
# 5) TABLE 3 — Turning Point implied by CCE-MG (baseline EKC: M2)
# ----------------------------
# Turning point: x = ln_GDPpc_centered, y = ... + b1*x + b2*x^2
# TP in x: -b1/(2*b2); then lnGDP_TP = x_TP + mean(ln_GDPpc_pool); GDPpc_TP = exp(lnGDP_TP)

mean_lnGDP_pool <- mean(panel_brics2$ln_GDPpc, na.rm = TRUE)

# Extraer b1, b2 y p-values desde ccemg_long (M2_CCEMG)
get_coef <- function(df, term_name){
  row <- df %>% dplyr::filter(model == "M2_CCEMG", term == term_name)
  if (nrow(row) != 1) return(list(est = NA_real_, se = NA_real_, p = NA_real_))
  list(est = row$estimate[1], se = row$std_error[1], p = row$p_value[1])
}

b1 <- get_coef(ccemg_long, "ln_GDPpc_centered")
b2 <- get_coef(ccemg_long, "ln_GDPpc_sq")

tp_gdp <- NA_real_
tp_x   <- NA_real_
shape  <- NA_character_

if (is.finite(b1$est) && is.finite(b2$est) && abs(b2$est) > 1e-10) {
  tp_x  <- -b1$est / (2*b2$est)
  ln_tp <- tp_x + mean_lnGDP_pool
  tp_gdp <- exp(ln_tp)
  
  if (b2$est < 0 && b1$est > 0) shape <- "Inverted-U (EKC)"
  if (b2$est > 0 && b1$est < 0) shape <- "U-shape (recarbonization)"
  if (b2$est < 0 && b1$est < 0) shape <- "Monotone decreasing (concave)"
  if (b2$est > 0 && b1$est > 0) shape <- "Monotone increasing (convex)"
}

# Soporte observado del PIB pc (pool)
gdp_series <- if ("GDP_pc_filled" %in% names(panel_brics2)) panel_brics2$GDP_pc_filled else exp(panel_brics2$ln_GDPpc)
gdp_min <- min(gdp_series, na.rm = TRUE)
gdp_max <- max(gdp_series, na.rm = TRUE)

within_support <- is.finite(tp_gdp) &&
  (tp_gdp >= gdp_min) && (tp_gdp <= gdp_max) &&
  is.finite(b2$est) && (b2$est < 0)

# (A) Fila pooled (CCE-MG, M2)
tab3 <- tibble::tibble(
  Model = "CCE-MG (M2 baseline)",
  `β1: ln GDP pc (centered)` = to_cell(b1$est, b1$se, b1$p, d = 3),
  `β2: ln GDP pc squared`    = to_cell(b2$est, b2$se, b2$p, d = 3),
  `Turning point (GDP pc)`   = ifelse(is.finite(tp_gdp), fmt_num(tp_gdp, 2), NA_character_),
  Shape                      = shape,
  `Within observed support`  = ifelse(isTRUE(within_support), "Yes", "No")
)

# (B) Filas por país (debajo de la pooled row), usando tp_by_country
if (!exists("tp_by_country")) stop("No existe tp_by_country. Ejecuta primero el bloque de Turning points por país.")

country_names <- c(
  "BRA" = "Brazil",
  "RUS" = "Russia",
  "IND" = "India",
  "CHN" = "China",
  "ZAF" = "South Africa"
)

# Asegurar iso3c (por si tp_by_country usa otro nombre)
if (!("iso3c" %in% names(tp_by_country))) {
  if ("Unit" %in% names(tp_by_country)) {
    tp_by_country <- tp_by_country %>% dplyr::rename(iso3c = Unit)
  } else if ("country" %in% names(tp_by_country)) {
    tp_by_country <- tp_by_country %>% dplyr::rename(iso3c = country)
  } else {
    stop("tp_by_country no tiene columna iso3c (ni Unit/country). Revisa names(tp_by_country).")
  }
}

tab3_countries <- tp_by_country %>%
  dplyr::filter(iso3c %in% names(country_names)) %>%
  dplyr::mutate(
    Model = country_names[iso3c],
    `β1: ln GDP pc (centered)` = to_cell(b1, NA_real_, NA_real_, d = 3),
    `β2: ln GDP pc squared`    = to_cell(b2, NA_real_, NA_real_, d = 3),
    `Turning point (GDP pc)`   = ifelse(is.finite(tp_hat), fmt_num(tp_hat, 2), NA_character_),
    Shape = dplyr::case_when(
      is.na(b1) | is.na(b2) ~ NA_character_,
      b2 < 0 & b1 > 0 ~ "Inverted-U (EKC)",
      b2 > 0 & b1 < 0 ~ "U-shape (recarbonization)",
      b2 < 0 & b1 < 0 ~ "Monotone decreasing (concave)",
      b2 > 0 & b1 > 0 ~ "Monotone increasing (convex)",
      TRUE ~ "Undefined"
    ),
    `Within observed support` = dplyr::if_else(isTRUE(within_support_country), "Yes", "No")
  ) %>%
  # Ordenar ANTES de eliminar iso3c
  dplyr::arrange(factor(iso3c, levels = c("BRA","RUS","IND","CHN","ZAF"))) %>%
  dplyr::select(
    Model,
    `β1: ln GDP pc (centered)`,
    `β2: ln GDP pc squared`,
    `Turning point (GDP pc)`,
    Shape,
    `Within observed support`
  )

tab3 <- dplyr::bind_rows(tab3, tab3_countries)

# (Opcional) CSV replicable
readr::write_csv(tab3, file.path(dir_tables, "Table_3_Turning_Point_CCE_MG_M2_withCountries.csv"))

ft3 <- journal_ft(tab3)
ft3 <- add_footer_lastrow(
  ft3,
  paste0(
    "Source: Authors’ calculations based on CCE-MG (M2) and country-specific baseline EKC fits. ",
    "Notes: Pooled turning point computed as exp( -β1/(2β2) + mean(ln GDP pc) ) using CCE-MG (M2) coefficients and pooled mean ln(GDPpc). ",
    "Country rows use tp_by_country (baseline EKC by country; see Table_TurningPoints_byCountry_bootstrap.csv for bootstrap CIs). ",
    "Within support: pooled row refers to pooled observed GDPpc range; country rows refer to country-specific observed GDPpc support (±2% tolerance). ",
    "Significance: *** p<0.01, ** p<0.05, * p<0.10 (pooled row only). ",
    "Created: ", format(Sys.time(), "%Y-%m-%d %H:%M"), "."
  )
)

doc3 <- officer::read_docx()
doc3 <- doc3 %>%
  officer::body_add_par("Table 3. Turning Point Implied by CCE-MG (Baseline EKC: M2)", style = "heading 1") %>%
  flextable::body_add_flextable(ft3)

print(doc3, target = file.path(dir_tables, "Table_3_Turning_Point_CCE_MG_M2.docx"))

message("Listo. Exportadas Table 1–3 a: ", dir_tables)

# =============================================================================
# BLOQUE: Table 4 + Appendices (A1–A3, B1) en Word
# =============================================================================
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(readr)
  library(flextable)
  library(officer)
})

# -----------------------------------------------------------------------------
# 0) FIX: conflictos con select() (MASS::select rompe dplyr::select)
# -----------------------------------------------------------------------------
# Usa dplyr::select SIEMPRE.

# -----------------------------------------------------------------------------
# 1) Helpers: stars, formatting, labels, style + footer
# -----------------------------------------------------------------------------
sig_stars <- function(p){
  dplyr::case_when(
    is.na(p) ~ "",
    p < 0.01 ~ "***",
    p < 0.05 ~ "**",
    p < 0.10 ~ "*",
    TRUE     ~ ""
  )
}
fmt_num <- function(x, d = 3){
  ifelse(is.na(x), "", formatC(x, format="f", digits=d))
}
to_cell <- function(est, se, p, d = 3){
  if (is.na(est)) return("")
  stars <- sig_stars(p)
  est_s <- fmt_num(est, d)
  if (is.na(se)) return(paste0(est_s, stars))
  se_s  <- fmt_num(se, d)
  paste0(est_s, stars, "\n(", se_s, ")")
}
term_labels <- c(
  "(Intercept)" = "Constant",
  "ln_GDPpc_centered" = "ln GDPpc (centered)",
  "ln_GDPpc_sq" = "[ln GDPpc (centered)]\u00B2",
  "ln_energy_pc" = "ln Energy use (pc)",
  "estabilidad_politica_vdem" = "Political stability (V-Dem index)",
  "gobernanza_vdem" = "Governance (V-Dem index)",
  "logpib_estabilidad" = "ln GDPpc \u00D7 Political stability",
  "energia_gobierno" = "ln Energy \u00D7 Governance"
)
apply_term_labels <- function(df){
  df %>% mutate(term = dplyr::recode(term, !!!term_labels, .default = term))
}
style_journal_ft <- function(ft){
  ft %>%
    fontsize(size = 10, part = "all") %>%
    font(fontname = "Times New Roman", part = "all") %>%
    align(align = "center", part = "header") %>%
    align(align = "left", part = "body", j = 1) %>%
    valign(valign = "top", part = "all") %>%
    set_table_properties(layout = "autofit") %>%
    theme_vanilla()
}
add_footer_std <- function(ft, source_line, notes_line = NULL){
  ft <- flextable::add_footer_lines(ft, values = source_line)
  ft <- flextable::add_footer_lines(ft, values = "Authors’ calculations.")
  if (!is.null(notes_line)) ft <- flextable::add_footer_lines(ft, values = notes_line)
  ft
}
make_wide_cells <- function(df, model_order = NULL){
  out <- df %>%
    mutate(cell = to_cell(estimate, std_error, p_value, d = 3)) %>%
    dplyr::select(model, term, cell) %>%
    apply_term_labels() %>%
    tidyr::pivot_wider(names_from = model, values_from = cell)
  if (!is.null(model_order)) {
    keep <- intersect(model_order, names(out))
    out <- out %>% dplyr::select(term, dplyr::all_of(keep))
  }
  out
}

# -----------------------------------------------------------------------------
# 2) Paths (según tu estructura)
# -----------------------------------------------------------------------------
# Si existen: dir_tables y dir_clean en tu sesión.
stopifnot(exists("dir_tables"))

ccemg_csv <- file.path(dir_tables, "Table_CCE_MG_Coefficients.csv")
pcce_csv  <- file.path(dir_tables, "Table_CCE_Pooled_Coefficients.csv")
fd_csv    <- file.path(dir_tables, "Table_FD_DK_Coefficients.csv")

cd_csv    <- file.path(dir_tables, "Table_PesaranCD_FE.csv")
fe_dk_csv <- file.path(dir_tables, "Table_FE_DK_Coefficients.csv")

tp_country_csv <- file.path(dir_tables, "Table_TurningPoints_byCountry_bootstrap.csv")

phi_nody_csv <- file.path(dir_tables, "Table_ECM_phi_diagnostics_noDY.csv")
phi_std_csv  <- file.path(dir_tables, "Table_ECM_phi_diagnostics.csv")

# -----------------------------------------------------------------------------
# 3) Load CSVs (robusto)
# -----------------------------------------------------------------------------
load_csv <- function(p){
  if (!file.exists(p)) stop("No encuentro el archivo: ", p)
  readr::read_csv(p, show_col_types = FALSE)
}

pcce_table <- load_csv(pcce_csv)
fd_table   <- load_csv(fd_csv)

cd_table   <- load_csv(cd_csv)
fe_dk_tbl  <- load_csv(fe_dk_csv)

tp_country <- load_csv(tp_country_csv)

phi_diag <- if (file.exists(phi_nody_csv)) load_csv(phi_nody_csv) else load_csv(phi_std_csv)

# -----------------------------------------------------------------------------
# 4) TABLE 4: Robustness (Panel A PCCE pooled; Panel B FD)
# =============================================================================
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(readr)
  library(flextable)
  library(officer)
})

stopifnot(exists("dir_tables"))

# -----------------------------------------------------------------------------
# 0) Helpers robustos
# -----------------------------------------------------------------------------
sig_stars <- function(p){
  dplyr::case_when(
    is.na(p) ~ "",
    p < 0.01 ~ "***",
    p < 0.05 ~ "**",
    p < 0.10 ~ "*",
    TRUE     ~ ""
  )
}
fmt_num <- function(x, d = 3){
  ifelse(is.na(x), "", formatC(x, format = "f", digits = d))
}

# Versión ESCALAR (para usar dentro de mapply)
to_cell1 <- function(est, se, p, d = 3){
  if (is.na(est) || !is.finite(est)) return("")
  stars <- sig_stars(p)
  est_s <- fmt_num(est, d)
  if (is.na(se) || !is.finite(se)) return(paste0(est_s, stars))
  se_s  <- fmt_num(se, d)
  paste0(est_s, stars, "\n(", se_s, ")")
}

# Vectorizado (esto evita el error "condition has length > 1")
to_cell_vec <- function(est, se, p, d = 3){
  mapply(to_cell1, est, se, p, MoreArgs = list(d = d), USE.NAMES = FALSE)
}

term_labels <- c(
  "(Intercept)" = "Constant",
  "ln_GDPpc_centered" = "ln GDPpc (centered)",
  "ln_GDPpc_sq" = "[ln GDPpc (centered)]\u00B2",
  "ln_energy_pc" = "ln Energy use (pc)",
  "estabilidad_politica_vdem" = "Political stability (V-Dem index)",
  "gobernanza_vdem" = "Governance (V-Dem index)",
  "logpib_estabilidad" = "ln GDPpc \u00D7 Political stability",
  "energia_gobierno" = "ln Energy \u00D7 Governance"
)

apply_term_labels <- function(df){
  df %>% mutate(term = dplyr::recode(term, !!!term_labels, .default = term))
}

style_journal_ft <- function(ft){
  ft %>%
    fontsize(size = 10, part = "all") %>%
    font(fontname = "Times New Roman", part = "all") %>%
    align(align = "center", part = "header") %>%
    align(align = "left", part = "body", j = 1) %>%
    valign(valign = "top", part = "all") %>%
    set_table_properties(layout = "autofit") %>%
    theme_vanilla()
}

add_footer_std <- function(ft, source_line, notes_line = NULL){
  ft <- flextable::add_footer_lines(ft, values = source_line)
  ft <- flextable::add_footer_lines(ft, values = "Authors’ calculations.")
  if (!is.null(notes_line)) ft <- flextable::add_footer_lines(ft, values = notes_line)
  ft
}

# Estandariza nombres típicos -> estimate / std_error / p_value
standardize_cols <- function(df){
  nm <- names(df)
  
  # estimate
  if (!("estimate" %in% nm)) {
    cand <- intersect(nm, c("Estimate","est","coef","Coefficient","estimate_mean"))
    if (length(cand) > 0) df <- df %>% rename(estimate = !!cand[1])
  }
  
  # std_error
  if (!("std_error" %in% nm)) {
    cand <- intersect(nm, c("Std. Error","Std.Error","Std_Error","std.error","se","SE","StdErr"))
    if (length(cand) > 0) df <- df %>% rename(std_error = !!cand[1])
  }
  
  # p_value
  if (!("p_value" %in% nm)) {
    cand <- intersect(nm, c("Pr(>|t|)","Pr(>|z|)","p.value","pvalue","P>|t|","P>|z|","p"))
    if (length(cand) > 0) df <- df %>% rename(p_value = !!cand[1])
  }
  
  # si no hay p_value, lo dejamos NA (pero NO rompemos)
  if (!("p_value" %in% names(df))) df$p_value <- NA_real_
  if (!("std_error" %in% names(df))) df$std_error <- NA_real_
  
  df
}

make_wide_cells_safe <- function(df, model_order = NULL, digits = 3){
  df2 <- df %>%
    standardize_cols() %>%
    mutate(
      estimate  = as.numeric(estimate),
      std_error = as.numeric(std_error),
      p_value   = as.numeric(p_value),
      cell      = to_cell_vec(estimate, std_error, p_value, d = digits)
    ) %>%
    dplyr::select(model, term, cell) %>%
    apply_term_labels() %>%
    tidyr::pivot_wider(names_from = model, values_from = cell)
  
  if (!is.null(model_order)) {
    keep <- intersect(model_order, names(df2))
    df2 <- df2 %>% dplyr::select(term, dplyr::all_of(keep))
  }
  df2
}

load_csv <- function(p){
  if (!file.exists(p)) stop("No encuentro el archivo: ", p)
  readr::read_csv(p, show_col_types = FALSE)
}

# =============================================================================
# TABLE 4: Robustness (Panel A PCCE pooled; Panel B FD)
# =============================================================================
pcce_csv <- file.path(dir_tables, "Table_CCE_Pooled_Coefficients.csv")
fd_csv   <- file.path(dir_tables, "Table_FD_DK_Coefficients.csv")

pcce_table <- load_csv(pcce_csv)
fd_table   <- load_csv(fd_csv)

# Panel A: PCCE pooled (M2, M3)
tab4A <- pcce_table %>%
  filter(model %in% c("M2_PCCE","M3_PCCE")) %>%
  make_wide_cells_safe(model_order = c("M2_PCCE","M3_PCCE"), digits = 3)

# Panel B: First Differences DK (M2, M3)
tab4B <- fd_table %>%
  filter(model %in% c("M2_FD_DK","M3_FD_DK")) %>%
  make_wide_cells_safe(model_order = c("M2_FD_DK","M3_FD_DK"), digits = 3)

ft4A <- flextable(tab4A) %>% style_journal_ft()
ft4B <- flextable(tab4B) %>% style_journal_ft()

doc4 <- read_docx() %>%
  body_add_par("Table 4. Robustness checks", style = "heading 1") %>%
  body_add_par("Panel A. PCCE pooled (robustness to estimator under cross-sectional dependence)", style = "heading 2") %>%
  body_add_flextable(
    add_footer_std(
      ft4A,
      source_line = "Source: World Bank (WDI) and V-Dem v15.",
      notes_line  = "Notes: Coefficients with standard errors in parentheses. Significance: * p<0.10, ** p<0.05, *** p<0.01. PCCE pooled imposes slope homogeneity relative to CCE-MG."
    )
  ) %>%
  body_add_par("Panel B. First differences with Driscoll–Kraay SE (short-run interpretation)", style = "heading 2") %>%
  body_add_flextable(
    add_footer_std(
      ft4B,
      source_line = "Source: World Bank (WDI) and V-Dem v15.",
      notes_line  = "Notes: First differences remove unit roots but eliminate long-run interpretation. Coefficients with DK SE. Significance: * p<0.10, ** p<0.05, *** p<0.01."
    )
  )

print(doc4, target = file.path(dir_tables, "Table_4_Robustness.docx"))

# -----------------------------------------------------------------------------
# 5) APPENDIX A1: Diagnostic Tests (Panel A CD; Panel B Unit Roots IPS )
# -----------------------------------------------------------------------------
# =============================================================================
# APPENDIX A1 — Diagnostic Tests (CD + IPS)
# =============================================================================

doc_a1 <- read_docx()

# ----------------------------
# Panel A: Pesaran CD
# ----------------------------
cd_path <- file.path(dir_tables, "Table_PesaranCD_FE.csv")
cd_tbl  <- readr::read_csv(cd_path, show_col_types = FALSE)

ft_cd <- flextable(cd_tbl) %>%
  autofit()

ft_cd <- add_footer_lines(
  ft_cd,
  values = c(
    "Source: Authors’ calculations.",
    "Notes: Panel A reports Pesaran (2004) CD tests computed on FE residuals. Driscoll–Kraay standard errors correct inference but do not eliminate bias under strong cross-sectional dependence."
  )
)

doc_a1 <- doc_a1 %>%
  body_add_par("Table A1. Diagnostic Tests", style = "heading 1") %>%
  body_add_par("Panel A. Pesaran CD Test", style = "heading 2") %>%
  body_add_flextable(ft_cd)

# ----------------------------
# Panel B: IPS Unit Root Tests
# ----------------------------
ft_ips <- flextable(ips_tbl) %>%
  autofit()

ft_ips <- add_footer_lines(
  ft_ips,
  values = c(
    "Source: Authors’ calculations.",
    "Notes: Panel B reports Im–Pesaran–Shin (IPS) panel unit root tests. Unit root tests are provided for diagnostic purposes only."
  )
)

doc_a1 <- doc_a1 %>%
  body_add_par("Panel B. IPS Unit Root Tests", style = "heading 2") %>%
  body_add_flextable(ft_ips)

print(doc_a1, target = file.path(dir_tables, "Table_A1_Diagnostic_Tests.docx"))

message("OK: Appendix A1 (CD + IPS) exportado.")

# -----------------------------------------------------------------------------
# 6) APPENDIX A2: Robustness – FE with Driscoll–Kraay SE
# =============================================================================
fe_dk_csv <- file.path(dir_tables, "Table_FE_DK_Coefficients.csv")
fe_dk_tbl <- load_csv(fe_dk_csv)

a2_wide <- fe_dk_tbl %>%
  make_wide_cells_safe(model_order = c("M1_FE_DK","M2_FE_DK","M3_FE_DK"), digits = 3)

ft_a2 <- flextable(a2_wide) %>% style_journal_ft()
ft_a2 <- add_footer_std(
  ft_a2,
  source_line = "Source: World Bank (WDI) and V-Dem v15.",
  notes_line  = "Notes: FE in levels can be biased under I(1) regressors and strong common factors; DK corrects SE but not necessarily coefficient bias. Reported as robustness only. Significance: * p<0.10, ** p<0.05, *** p<0.01."
)

docA2 <- read_docx() %>%
  body_add_par("Appendix Table A2. Fixed effects with Driscoll–Kraay standard errors (robustness)", style = "heading 1") %>%
  body_add_flextable(ft_a2)

print(docA2, target = file.path(dir_tables, "Table_A2_FE_DK_Robustness.docx"))

message("OK: corregidas y exportadas Table 4 y Table A2 en /tables.")
# -----------------------------------------------------------------------------
# 7) APPENDIX A3: Institutional Indices (Correlation + Validation)
# -----------------------------------------------------------------------------
# Necesita panel_brics2; si no está en memoria, intenta cargar el RDS más reciente.
if (!exists("panel_brics2")) {
  if (exists("dir_clean")) {
    cand <- list.files(dir_clean, pattern = "^panel_BRICS_macro_vdem_.*\\.rds$", full.names = TRUE)
    if (length(cand) > 0) {
      cand <- cand[order(file.info(cand)$mtime, decreasing = TRUE)]
      panel_brics2 <- readRDS(cand[1])
    } else {
      stop("No encuentro panel_brics2 en memoria y tampoco un RDS panel_BRICS_macro_vdem_*.rds en dir_clean.")
    }
  } else {
    stop("No existe panel_brics2 en memoria y tampoco dir_clean definido.")
  }
}

# Panel A: correlación
inst_df <- panel_brics2 %>%
  dplyr::select(iso3c, year, estabilidad_politica_vdem, gobernanza_vdem) %>%
  tidyr::drop_na()

rho <- suppressWarnings(cor(inst_df$estabilidad_politica_vdem, inst_df$gobernanza_vdem, use = "complete.obs"))

a3A <- tibble::tibble(
  Pair = "Political stability (V-Dem) vs Governance (V-Dem)",
  Correlation = fmt_num(rho, 3),
  N = nrow(inst_df)
)

# Panel B: validación mínima
a3B <- inst_df %>%
  summarise(
    `Political stability: N` = sum(!is.na(estabilidad_politica_vdem)),
    `Political stability: Mean` = mean(estabilidad_politica_vdem, na.rm=TRUE),
    `Political stability: SD` = sd(estabilidad_politica_vdem, na.rm=TRUE),
    `Political stability: Min` = min(estabilidad_politica_vdem, na.rm=TRUE),
    `Political stability: Max` = max(estabilidad_politica_vdem, na.rm=TRUE),
    `Governance: N` = sum(!is.na(gobernanza_vdem)),
    `Governance: Mean` = mean(gobernanza_vdem, na.rm=TRUE),
    `Governance: SD` = sd(gobernanza_vdem, na.rm=TRUE),
    `Governance: Min` = min(gobernanza_vdem, na.rm=TRUE),
    `Governance: Max` = max(gobernanza_vdem, na.rm=TRUE)
  ) %>%
  pivot_longer(everything(), names_to="Metric", values_to="Value") %>%
  mutate(Value = ifelse(str_detect(Metric, ": N$"), as.character(round(Value)), fmt_num(Value, 3)))

ft_a3A <- flextable(a3A) %>% style_journal_ft()
ft_a3B <- flextable(a3B) %>% style_journal_ft()

ft_a3A <- add_footer_std(
  ft_a3A,
  source_line = "Source: V-Dem v15.",
  notes_line  = "Notes: Correlation computed over pooled BRICS 1990–2024 sample after index construction (z-score pooled as defined in the code)."
)
ft_a3B <- add_footer_std(
  ft_a3B,
  source_line = "Source: V-Dem v15.",
  notes_line  = "Notes: Validation panel reports basic distributional checks for the constructed indices."
)

docA3 <- read_docx() %>%
  body_add_par("Appendix Table A3. Institutional indices", style = "heading 1") %>%
  body_add_par("Panel A. Correlation", style = "heading 2") %>%
  body_add_flextable(ft_a3A) %>%
  body_add_par("Panel B. Validation summary", style = "heading 2") %>%
  body_add_flextable(ft_a3B)

print(docA3, target = file.path(dir_tables, "Table_A3_Institutional_Indices.docx"))

# -----------------------------------------------------------------------------
# 8) APPENDIX B1: Country-specific results (TPs + ECM phi diagnostics)
# -----------------------------------------------------------------------------
# Panel A: TPs por país (bootstrap)
b1A <- tp_country %>%
  mutate(
    tp_hat  = fmt_num(tp_hat, 0),
    tp_p025 = fmt_num(tp_p025, 0),
    tp_p975 = fmt_num(tp_p975, 0),
    gdp_min = fmt_num(gdp_min, 0),
    gdp_max = fmt_num(gdp_max, 0),
    b1 = fmt_num(b1, 3),
    b2 = fmt_num(b2, 3)
  ) %>%
  dplyr::select(iso3c, tp_hat, tp_p025, tp_p975, gdp_min, gdp_max, b2, within_support_country, note)

# Panel B: ECM phi diagnostics (elige el que ya existe)
b1B <- phi_diag %>%
  mutate(
    phi = fmt_num(phi, 3),
    abs_phi = fmt_num(abs_phi, 3)
  )

ft_b1A <- flextable(b1A) %>% style_journal_ft()
ft_b1B <- flextable(b1B) %>% style_journal_ft()

ft_b1A <- add_footer_std(
  ft_b1A,
  source_line = "Source: World Bank (WDI) and V-Dem v15.",
  notes_line  = "Notes: Turning points are computed from country-specific quadratic fits and parametric bootstrap intervals, then evaluated against each country’s observed GDPpc support. Values reported in PPP constant international dollars."
)
ft_b1B <- add_footer_std(
  ft_b1B,
  source_line = "Source: Authors’ calculations.",
  notes_line  = "Notes: ECM speed of adjustment (phi). Overshooting flagged when phi < -1 (annual frequency)."
)

docB1 <- read_docx() %>%
  body_add_par("Appendix Table B1. Country-specific results", style = "heading 1") %>%
  body_add_par("Panel A. Turning points by country (bootstrap)", style = "heading 2") %>%
  body_add_flextable(ft_b1A) %>%
  body_add_par("Panel B. ECM phi diagnostics", style = "heading 2") %>%
  body_add_flextable(ft_b1B)

print(docB1, target = file.path(dir_tables, "Table_B1_Country_Specific_Results.docx"))

message("OK: Table 4 + Appendices exportados a Word en: ", dir_tables)

###############################################################################
# FIGURES — EKC BRICS (1990–2024) | WDI + V-Dem
# ---------------------------------------------------------------------------
# Creates final figures agreed (MAIN + optional ROBUSTEZ) with:
#   - Color, aesthetic, and still B/W-print-safe (shapes + linetypes)
#   - PNG 300 dpi + PDF vector export
#   - No re-estimation: reads CSV outputs from /tables and plots
#
###############################################################################

# =========================
# 0) Packages
# =========================
pkgs <- c("tidyverse", "readr", "fs", "stringr", "glue", "ggplot2")
to_install <- pkgs[!pkgs %in% installed.packages()[, "Package"]]
if (length(to_install) > 0) install.packages(to_install)
invisible(lapply(pkgs, library, character.only = TRUE))

# =========================
# 1) Paths + Output options
# =========================
CONFIG <- list(
  out_fig  = "figures",
  in_tabs  = "tables",
  dpi      = 300,
  w_main   = 6.8,
  h_main   = 4.2,
  w_wide   = 6.8,
  h_wide   = 4.8
)
if (!dir.exists(CONFIG$out_fig)) dir.create(CONFIG$out_fig, recursive = TRUE)

# =========================
# 2) Flag-consistent BRICS palette + theme
# =========================
pal_brics <- c(
  BRA = "#1B9E77",  # Brazil green
  CHN = "#D73027",  # China red
  RUS = "#4575B4",  # Russia blue
  IND = "#F28E2B",  # India saffron/orange
  ZAF = "#006837"   # South Africa green (darker)
)

theme_journal_color <- function(base_size = 11, base_family = "serif") {
  theme_minimal(base_size = base_size, base_family = base_family) +
    theme(
      plot.title    = element_text(face = "bold", size = base_size + 2, color = "black"),
      plot.subtitle = element_text(size = base_size, color = "grey25"),
      axis.title    = element_text(size = base_size, color = "black"),
      axis.text     = element_text(size = base_size - 1, color = "black"),
      legend.position = "bottom",
      legend.title  = element_text(face = "plain"),
      legend.text   = element_text(color = "black"),
      panel.grid.major = element_line(color = "grey88", linewidth = 0.35),
      panel.grid.minor = element_blank(),
      plot.margin = margin(10, 10, 10, 10)
    )
}

save_png_pdf <- function(p, filename_stub, w, h, dpi = 300) {
  ggsave(file.path(CONFIG$out_fig, paste0(filename_stub, ".png")),
         p, width = w, height = h, dpi = dpi, units = "in", bg = "white")
  ggsave(file.path(CONFIG$out_fig, paste0(filename_stub, ".pdf")),
         p, width = w, height = h, units = "in")
  message("Saved: ", filename_stub, " (.png + .pdf)")
}

stop_if_missing_cols <- function(df, cols, df_name = "data") {
  miss <- setdiff(cols, names(df))
  if (length(miss) > 0) stop(glue("'{df_name}' missing columns: {paste(miss, collapse=', ')}"))
}
first_existing <- function(df, candidates, default = NULL) {
  hit <- candidates[candidates %in% names(df)]
  if (length(hit) == 0) return(default)
  df[[hit[1]]]
}

# =========================
# 3) Deterministic inputs
# =========================
path_cce      <- file.path(CONFIG$in_tabs, "Table_CCE_MG_Coefficients.csv")
path_tp_cy    <- file.path(CONFIG$in_tabs, "Table_TurningPoints_byCountry_bootstrap.csv")
path_tp_lr    <- file.path(CONFIG$in_tabs, "Table_TurningPoint_LR_CS_ECM.csv")
path_tp_fe    <- file.path(CONFIG$in_tabs, "Table_TurningPoint_panel_FE_illustrative.csv")
path_phi_cy   <- file.path(CONFIG$in_tabs, "Table_ECM_phi_country.csv")
# ---- helper: check required columns ----
stop_if_missing_cols <- function(df, cols, df_name = "data") {
  miss <- setdiff(cols, names(df))
  if (length(miss) > 0) stop(paste0("'", df_name, "' missing columns: ", paste(miss, collapse = ", ")))
}

# ---- helper: standardize coefficient table WITHOUT failing if cols don't exist ----
standardize_coef_table <- function(df) {
  df <- df %>%
    dplyr::rename_with(~stringr::str_replace_all(.x, "\\s+", "_"))
  
  df <- df %>%
    dplyr::rename(
      term     = dplyr::any_of(c("term", "variable", "regressor", "Parameter", "param", "name")),
      estimate = dplyr::any_of(c("estimate", "coef", "Coefficient", "Estimate", "b", "beta")),
      se       = dplyr::any_of(c("std_error", "std.error", "se", "Std_Error", "Std..Error", "StdError"))
    )
  
  df
}

###############################################################################
# FIG1 (MAIN): EKC curve from CCE-type benchmark (shape-only, normalized)
###############################################################################
path_cce <- file.path(CONFIG$in_tabs, "Table_CCE_MG_Coefficients.csv")

if (file.exists(path_cce)) {
  
  cce <- readr::read_csv(path_cce, show_col_types = FALSE) %>%
    standardize_coef_table()
  
  stop_if_missing_cols(cce, c("term", "estimate"), "CCE coefficients")
  
  b1_row <- cce %>%
    dplyr::filter(stringr::str_detect(term, stringr::regex(
      "ln.*GDP.*center|ln_GDPpc_centered|GDPpc_centered|lnGDP|log.*pib|pib.*center",
      ignore_case = TRUE
    ))) %>%
    dplyr::slice(1)
  
  b2_row <- cce %>%
    dplyr::filter(stringr::str_detect(term, stringr::regex(
      "GDP.*sq|GDP\\^2|ln_GDPpc_sq|GDPpc_sq|centered.*sq|pib.*sq",
      ignore_case = TRUE
    ))) %>%
    dplyr::slice(1)
  
  if (nrow(b1_row) == 0 || nrow(b2_row) == 0) {
    message("Fig1 skipped: cannot find GDP linear/squared terms in Table_CCE_MG_Coefficients.csv")
  } else {
    
    b1 <- b1_row$estimate[1]
    b2 <- b2_row$estimate[1]
    tp_x <- -b1 / (2 * b2)
    
    x_grid <- seq(-1.5, 1.5, length.out = 300)
    
    y_raw <- b1 * x_grid + b2 * x_grid^2
    y     <- y_raw - (b1 * 0 + b2 * 0)
    
    df_curve <- tibble::tibble(x = x_grid, y = y)
    
    add_band <- ("se" %in% names(cce)) && (!is.na(b1_row$se[1])) && (!is.na(b2_row$se[1]))
    if (add_band) {
      se_b1 <- b1_row$se[1]
      se_b2 <- b2_row$se[1]
      var_y <- (x_grid^2) * (se_b1^2) + (x_grid^4) * (se_b2^2)
      se_y  <- sqrt(var_y)
      
      df_curve <- df_curve %>%
        dplyr::mutate(
          y_lo = y - 1.96 * se_y,
          y_hi = y + 1.96 * se_y
        )
    }
    
    p1 <- ggplot2::ggplot(df_curve, ggplot2::aes(x = x, y = y)) +
      {if (add_band) ggplot2::geom_ribbon(ggplot2::aes(ymin = y_lo, ymax = y_hi),
                                          fill = "#9ecae1", alpha = 0.25)} +
      ggplot2::geom_line(linewidth = 1.05, color = "#2171b5") +
      ggplot2::geom_vline(xintercept = tp_x, linetype = "dashed",
                          linewidth = 0.7, color = "grey35") +
      ggplot2::labs(
        title = "Estimated EKC shape (CCE-type benchmark)",
        subtitle = "Centered ln(GDP pc). Curve normalized to emphasize curvature; dashed line indicates implied turning point.",
        x = "Centered ln(GDP per capita)",
        y = "Relative ln(CO\u2082 pc) (normalized shape)"
      ) +
      theme_journal_color()
    
    save_png_pdf(p1, "Fig1_MAIN_EKC_PanelCurve_CCE",
                 CONFIG$w_main, CONFIG$h_main, CONFIG$dpi)
  }
  
} else {
  message("Fig1 skipped: missing ", path_cce)
}

###############################################################################
# FIG2 (MAIN): Turning points by country + CI + support (range) — ROBUST VERSION
###############################################################################
path_tp_cy <- file.path(CONFIG$in_tabs, "Table_TurningPoints_byCountry_bootstrap.csv")

if (file.exists(path_tp_cy)) {
  
  tp <- readr::read_csv(path_tp_cy, show_col_types = FALSE) %>%
    dplyr::rename_with(~stringr::str_replace_all(.x, "\\s+", "_"))
  
  # 1) Create standardized columns ONLY if they exist (no failures)
  tp <- tp %>%
    dplyr::rename(
      iso3c = dplyr::any_of(c("iso3c", "country", "iso")),
      tp_hat = dplyr::any_of(c("tp_hat", "turning_point", "tp", "tp_point", "tp_value")),
      gdp_min = dplyr::any_of(c("gdp_min", "GDP_min", "min_gdp", "xmin", "support_min")),
      gdp_max = dplyr::any_of(c("gdp_max", "GDP_max", "max_gdp", "xmax", "support_max")),
      # CI candidates (if present)
      tp_p025 = dplyr::any_of(c("tp_p025", "tp_q025", "tp_lo", "tp_low", "ci_lo", "p025")),
      tp_p975 = dplyr::any_of(c("tp_p975", "tp_q975", "tp_hi", "tp_high", "ci_hi", "p975")),
      within_support_country = dplyr::any_of(c("within_support_country", "within_support", "within_support_ctry"))
    )
  
  # 2) Required columns to draw Fig2
  stop_if_missing_cols(tp, c("iso3c", "tp_hat", "gdp_min", "gdp_max"), "Turning points by country")
  
  # 3) If support flag not available, compute it
  if (!("within_support_country" %in% names(tp))) {
    tp <- tp %>% dplyr::mutate(within_support_country = (tp_hat >= gdp_min & tp_hat <= gdp_max))
  } else {
    # if exists but has many NA, fill NA with computed rule
    tp <- tp %>%
      dplyr::mutate(
        within_support_country = dplyr::if_else(
          is.na(within_support_country),
          (tp_hat >= gdp_min & tp_hat <= gdp_max),
          as.logical(within_support_country)
        )
      )
  }
  
  # 4) Country ordering (by TP value) + labels
  tp <- tp %>%
    dplyr::mutate(
      iso3c = forcats::fct_reorder(iso3c, tp_hat),  # REORDER
      support_label = dplyr::if_else(within_support_country, "Within support", "Out of support")
    )
  
  # 5) Build improved plot
  has_ci <- ("tp_p025" %in% names(tp)) && ("tp_p975" %in% names(tp)) &&
    (any(!is.na(tp$tp_p025)) && any(!is.na(tp$tp_p975)))
  
  p2 <- ggplot2::ggplot(tp, ggplot2::aes(y = iso3c)) +
    ggplot2::geom_segment(
      ggplot2::aes(x = gdp_min, xend = gdp_max, yend = iso3c, color = iso3c),
      linewidth = 2.2, alpha = 0.35
    ) +
    { if (has_ci)
      ggplot2::geom_errorbarh(
        ggplot2::aes(xmin = tp_p025, xmax = tp_p975),
        height = 0.15, linewidth = 0.6, color = "black", na.rm = TRUE
      )
    } +
    ggplot2::geom_point(
      ggplot2::aes(x = tp_hat, fill = iso3c, shape = support_label),
      size = 3.1, color = "black", stroke = 0.6
    ) +
    ggplot2::geom_text(
      ggplot2::aes(x = tp_hat, label = scales::dollar(tp_hat, accuracy = 1)),
      nudge_x = 5000, size = 3, hjust = 0
    ) +
    ggplot2::scale_shape_manual(values = c("Within support" = 21, "Out of support" = 1)) +
    ggplot2::scale_x_continuous(
      limits = c(0, 45000),
      labels = scales::label_comma(prefix = "$", accuracy = 1)
    ) +
    ggplot2::labs(
      title = "Turning points of the EKC by country (with support)",
      subtitle = "Colored bars show observed GDP per capita ranges; points show estimated turning points (CI when available).",
      x = "GDP per capita at turning point (PPP constant int.$)",
      y = NULL, color = NULL, fill = NULL, shape = NULL
    ) +
    theme_journal_color() +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_blank(),
      plot.margin = margin(10, 30, 10, 10)
    )
  
  # Guardar
  save_png_pdf(p2, "Fig2_MAIN_TurningPoints_ByCountry_Support",
               CONFIG$w_wide, CONFIG$h_main, CONFIG$dpi)  

  ggplot2::ggsave(
    filename = "figures/Fig2_TurningPoints_ByCountry_Support.png",
    plot = p2,
    width = 6.8,       # en pulgadas (≈ 174 mm) para una columna
    height = 4.2,      # razonable para que no sea muy alta
    units = "in",
    dpi = 300,         # resolución mínima para impresión
    bg = "white"
  )
  
} else {
  message("Fig2 skipped: missing ", path_tp_cy)
}

###############################################################################
# FIG3 (MAIN): FE turning point vs Long-run CS-ECM turning point (robust extract)
###############################################################################
path_tp_lr <- file.path(CONFIG$in_tabs, "Table_TurningPoint_LR_CS_ECM.csv")
path_tp_fe <- file.path(CONFIG$in_tabs, "Table_TurningPoint_panel_FE_illustrative.csv")

if (file.exists(path_tp_lr) && file.exists(path_tp_fe)) {
  
  tp_lr <- readr::read_csv(path_tp_lr, show_col_types = FALSE) %>%
    dplyr::rename_with(~stringr::str_replace_all(.x, "\\s+", "_"))
  tp_fe <- readr::read_csv(path_tp_fe, show_col_types = FALSE) %>%
    dplyr::rename_with(~stringr::str_replace_all(.x, "\\s+", "_"))
  
  message("Fig3 input columns (LR CS-ECM): ", paste(names(tp_lr), collapse = ", "))
  message("Fig3 input columns (FE illustrative): ", paste(names(tp_fe), collapse = ", "))
  
  # robust scalar extractor:
  get_tp_scalar <- function(df) {
    # common candidate names (add here if your pipeline uses a different label)
    cand <- c(
      "tp_hat", "turning_point", "turningpoint", "tp", "tp_value",
      "tp_gdp", "tp_gdppc", "tp_gdp_pc", "tp_gdp_percap",
      "gdp_tp", "gdp_at_tp", "gdp_at_turning_point",
      "turning_point_gdp", "tp_gdp_hat",
      "tp_level", "tp_point",
      # sometimes stored as log-space:
      "tp_ln", "tp_lngdp", "tp_x", "x_tp"
    )
    
    hit <- cand[cand %in% names(df)]
    if (length(hit) > 0) {
      val <- suppressWarnings(as.numeric(df[[hit[1]]][1]))
      if (!is.na(val)) return(val)
    }
    
    # fallback: first numeric column, first row
    num_cols <- names(df)[vapply(df, is.numeric, logical(1))]
    if (length(num_cols) > 0) {
      val <- suppressWarnings(as.numeric(df[[num_cols[1]]][1]))
      if (!is.na(val)) return(val)
    }
    
    # last resort: scan first row, coerce to numeric where possible
    row1 <- df[1, , drop = FALSE]
    for (nm in names(row1)) {
      val <- suppressWarnings(as.numeric(row1[[nm]][1]))
      if (!is.na(val)) return(val)
    }
    
    NA_real_
  }
  
  tp_lr_val <- get_tp_scalar(tp_lr)
  tp_fe_val <- get_tp_scalar(tp_fe)
  
  if (is.na(tp_lr_val) || is.na(tp_fe_val)) {
    message("Fig3 skipped: could not extract turning point scalar from one or both inputs.")
    message("Extracted tp_lr_val = ", tp_lr_val, " | tp_fe_val = ", tp_fe_val)
  } else {
    
    df3 <- tibble::tibble(
      estimator = factor(
        c("FE (within, illustrative)", "CS-ECM (long-run)"),
        levels = c("FE (within, illustrative)", "CS-ECM (long-run)")
      ),
      turning_point = c(tp_fe_val, tp_lr_val)
    )
    
    p3 <- ggplot2::ggplot(df3, ggplot2::aes(x = estimator, y = turning_point)) +
      ggplot2::geom_segment(
        ggplot2::aes(x = estimator, xend = estimator, y = 0, yend = turning_point),
        linewidth = 1.0, color = "grey70"
      ) +
      ggplot2::geom_point(size = 3.4, color = "black") +
      ggplot2::labs(
        title = "Turning point comparison: within FE vs long-run CS-ECM",
        subtitle = "A shift in the implied turning point illustrates how common factors and long-run dynamics can alter EKC conclusions.",
        x = NULL,
        y = "Turning point (GDP per capita, PPP constant int.$)"
      ) +
      theme_journal_color() +
      ggplot2::theme(panel.grid.major.x = ggplot2::element_blank())
    
    save_png_pdf(p3, "Fig3_MAIN_TurningPoints_FE_vs_LR_CS_ECM",
                 CONFIG$w_main, CONFIG$h_main, CONFIG$dpi)
  }
  
} else {
  message("Fig3 skipped: missing ", path_tp_lr, " or ", path_tp_fe)
}


###############################################################################
# FIGR1 (ROBUST): phi by country (convergence + overshooting) — ROBUST VERSION
###############################################################################
path_phi_cy <- file.path(CONFIG$in_tabs, "Table_ECM_phi_country.csv")

if (file.exists(path_phi_cy)) {
  
  phi <- readr::read_csv(path_phi_cy, show_col_types = FALSE) %>%
    dplyr::rename_with(~stringr::str_replace_all(.x, "\\s+", "_")) %>%
    # rename safely (no errors if a candidate doesn't exist)
    dplyr::rename(
      iso3c = dplyr::any_of(c("iso3c", "country_code", "iso", "code")),
      phi   = dplyr::any_of(c("phi", "ECM_phi", "ecm_phi", "phi_hat", "phi_value", "ecm_speed"))
    )
  
  stop_if_missing_cols(phi, c("iso3c", "phi"), "Phi by country")
  
  phi <- phi %>%
    dplyr::filter(!is.na(iso3c), !is.na(phi)) %>%
    dplyr::mutate(
      iso3c = factor(iso3c, levels = c("BRA","CHN","RUS","IND","ZAF")),
      class = dplyr::case_when(
        phi < -1 ~ "Overshooting (phi < -1)",
        phi < 0  ~ "Convergence (phi < 0)",
        TRUE     ~ "No convergence (phi >= 0)"
      ),
      class = factor(
        class,
        levels = c("Convergence (phi < 0)", "Overshooting (phi < -1)", "No convergence (phi >= 0)")
      )
    )
  
  pR1 <- ggplot2::ggplot(phi, ggplot2::aes(x = iso3c, y = phi)) +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.6, color = "grey45") +
    ggplot2::geom_point(ggplot2::aes(color = iso3c, shape = class),
                        size = 3.1, stroke = 0.8) +
    ggplot2::scale_color_manual(values = pal_brics) +
    ggplot2::scale_shape_manual(values = c(
      "Convergence (phi < 0)"       = 16,
      "Overshooting (phi < -1)"     = 17,
      "No convergence (phi >= 0)"   = 4
    )) +
    ggplot2::labs(
      title = "ECM convergence diagnostics by country",
      subtitle = "Phi measures speed of adjustment; overshooting flagged when phi < -1.",
      x = NULL,
      y = "Phi (ECM)",
      color = NULL,
      shape = NULL
    ) +
    theme_journal_color()
  
  save_png_pdf(pR1, "FigR1_ROBUST_phi_by_country",
               CONFIG$w_main, CONFIG$h_main, CONFIG$dpi)
  
} else {
  message("FigR1 skipped: missing ", path_phi_cy)
}

message("Done. Figures saved in: ", CONFIG$out_fig)
###############################################################################
# =========================
# FIG10: Panel EKC + Turning point (journal-ready)
# Uses: panel_brics columns:
#   GDP_pc_filled, ln_CO2pc, ln_GDPpc, ln_GDPpc_centered, ln_GDPpc_sq
library(tidyverse)
library(scales)

# -------------------------------------------------
# 0) Checks mínimos
# -------------------------------------------------
stopifnot(all(c(
  "GDP_pc_filled",
  "ln_CO2pc",
  "ln_GDPpc"
) %in% names(panel_brics)))

# -------------------------------------------------
# 1) Coeficientes escalares (desde listas)
# -------------------------------------------------
b1_est <- as.numeric(b1$est)
b2_est <- as.numeric(b2$est)

stopifnot(is.finite(b1_est), is.finite(b2_est))

# -------------------------------------------------
# 2) Media de ln(GDPpc) usada para centrar
# -------------------------------------------------
mu_ln_gdp <- mean(panel_brics$ln_GDPpc, na.rm = TRUE)

# -------------------------------------------------
# 3) Soporte observado en niveles (GDP pc)
# -------------------------------------------------
xmin <- min(panel_brics$GDP_pc_filled, na.rm = TRUE)
xmax <- max(panel_brics$GDP_pc_filled, na.rm = TRUE)

# -------------------------------------------------
# 4) Grid en niveles → métrica centrada del modelo
# -------------------------------------------------
x_grid_lvl <- seq(xmin, xmax, length.out = 400)
x_grid_c   <- log(x_grid_lvl) - mu_ln_gdp

# -------------------------------------------------
# 5) Curva EKC (FORMA) — común a A y B
# -------------------------------------------------
y_hat <- b1_est * x_grid_c + b2_est * (x_grid_c^2)

# normalización en x_c = 0 (forma pura)
y_hat <- y_hat - (b1_est * 0 + b2_est * 0)

# -------------------------------------------------
# 6) Opción A — EKC shape NORMALIZADA
# -------------------------------------------------
df_curve_A <- tibble(
  GDP_pc_filled = x_grid_lvl,
  yhat = y_hat
)

# -------------------------------------------------
# 7) Opción B — EKC shape RE-ANCLADA al nivel medio
# -------------------------------------------------
y_shift <- mean(panel_brics$ln_CO2pc, na.rm = TRUE)

df_curve_B <- tibble(
  GDP_pc_filled = x_grid_lvl,
  yhat = y_hat + y_shift
)

# chequeo rápido (útil en consola)
range(df_curve_A$yhat, na.rm = TRUE)
range(df_curve_B$yhat, na.rm = TRUE)
y_shift

# -------------------------------------------------
# 8) Turning point (centrado → niveles)
# -------------------------------------------------
tp_c   <- -b1_est / (2*b2_est)
tp_lvl <- exp(tp_c + mu_ln_gdp)

# -------------------------------------------------
# 9) FIG10A — EKC normalizada (conceptualmente estricta)
# -------------------------------------------------
p10_A <- ggplot() +
  annotate("rect", xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, alpha = 0.06) +
  geom_point(
    data = panel_brics,
    aes(x = GDP_pc_filled, y = ln_CO2pc),
    alpha = 0.18, size = 1.0
  ) +
  geom_line(
    data = df_curve_A,
    aes(x = GDP_pc_filled, y = yhat),
    linewidth = 1.15, color = "#2171b5"
  ) +
  geom_vline(
    xintercept = tp_lvl,
    linetype = "dashed",
    linewidth = 0.7,
    color = "grey35"
  ) +
  annotate(
    "label",
    x = tp_lvl,
    y = max(df_curve_A$yhat, na.rm = TRUE),
    label = paste0("TP panel = ", dollar(round(tp_lvl, 0))),
    hjust = -0.05,
    vjust = 1.2,
    label.size = 0.2
  ) +
  scale_x_continuous(labels = dollar) +
  labs(
    title = "Environmental Kuznets Curve (panel) and implied turning point",
    subtitle = "Dots: observed ln(CO\u2082pc). Line: normalized EKC shape from quadratic in centered ln(GDPpc).",
    x = "GDP per capita (PPP, constant international $)",
    y = "Observed ln(CO\u2082pc (points) / EKC shape (normalized line)"
  ) +
  theme_journal_color() +
  theme(
    plot.background  = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  )

# ===============================
# FIG10B — FINAL (colored points by country + TP)
# EKC panel + implied turning point
# ===============================

library(tidyverse)
library(scales)

# -------------------------------------------------
# 0) Checks mínimos
# -------------------------------------------------
stopifnot(all(c("GDP_pc_filled","ln_CO2pc","ln_GDPpc","iso3c") %in% names(panel_brics)))

# -------------------------------------------------
# 1) Paleta BRICS (banderas)
# -------------------------------------------------
pal_brics <- c(
  BRA = "#1B9E77",  # Brazil green
  CHN = "#D73027",  # China red
  RUS = "#4575B4",  # Russia blue
  IND = "#F28E2B",  # India saffron/orange
  ZAF = "#006837"   # South Africa green (darker)
)

# -------------------------------------------------
# 2) Coeficientes escalares (desde listas b1/b2)
# -------------------------------------------------
b1_est <- as.numeric(b1$est)
b2_est <- as.numeric(b2$est)
stopifnot(is.finite(b1_est), is.finite(b2_est))

# -------------------------------------------------
# 3) Media ln(GDPpc) para centrar y soporte observado
# -------------------------------------------------
mu_ln_gdp <- mean(panel_brics$ln_GDPpc, na.rm = TRUE)

xmin <- min(panel_brics$GDP_pc_filled, na.rm = TRUE)
xmax <- max(panel_brics$GDP_pc_filled, na.rm = TRUE)

# -------------------------------------------------
# 4) Grid en niveles → escala centrada del modelo
# -------------------------------------------------
x_grid_lvl <- seq(xmin, xmax, length.out = 400)
x_grid_c   <- log(x_grid_lvl) - mu_ln_gdp

# -------------------------------------------------
# 5) Curva EKC (forma) normalizada y re-anclada al nivel medio
# -------------------------------------------------
y_hat <- b1_est * x_grid_c + b2_est * (x_grid_c^2)
y_hat <- y_hat - (b1_est * 0 + b2_est * 0)   # normaliza en x_c=0

y_shift <- mean(panel_brics$ln_CO2pc, na.rm = TRUE)

df_curve_B <- tibble(
  GDP_pc_filled = x_grid_lvl,
  yhat = y_hat + y_shift
)

# -------------------------------------------------
# 6) Turning point (centrado → niveles)
# -------------------------------------------------
tp_c   <- -b1_est / (2*b2_est)
tp_lvl <- exp(tp_c + mu_ln_gdp)

# -------------------------------------------------
# 7) Plot (colored points)
# -------------------------------------------------
p10_B_col <- ggplot() +
  annotate("rect", xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, alpha = 0.06) +
  geom_point(
    data = panel_brics,
    aes(x = GDP_pc_filled, y = ln_CO2pc, color = iso3c),
    alpha = 0.35, size = 1.2
  ) +
  geom_line(
    data = df_curve_B,
    aes(x = GDP_pc_filled, y = yhat),
    linewidth = 1.15, color = "#2171b5"
  ) +
  geom_vline(
    xintercept = tp_lvl,
    linetype = "dashed",
    linewidth = 0.7,
    color = "grey35"
  ) +
  annotate(
    "label",
    x = tp_lvl,
    y = max(df_curve_B$yhat, na.rm = TRUE),
    label = paste0("TP panel = ", dollar(round(tp_lvl, 0))),
    hjust = -0.05,
    vjust = 1.2,
    label.size = 0.2
  ) +
  scale_x_continuous(labels = dollar) +
  scale_color_manual(values = pal_brics, breaks = c("BRA","CHN","IND","RUS","ZAF")) +
  labs(
    title = "Environmental Kuznets Curve (panel) and implied turning point",
    subtitle = "Dots: BRICS country-year observations. Line: EKC shape shifted to the sample mean.",
    x = "GDP per capita (PPP, constant international $)",
    y = "ln(CO\u2082 per capita)",
    color = NULL
  ) +
  theme_journal_color() +
  theme(
    plot.background  = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    legend.position = "bottom"
  )

# -------------------------------------------------
# 8) Guardado
# -------------------------------------------------
if (!dir.exists("figures")) dir.create("figures")

ggsave(
  "figures/Fig10B_EKC_Panel_TP_ShiftedMean_COLORS.png",
  p10_B_col,
  width = 6.8,
  height = 4.2,
  dpi = 300,
  units = "in",
  bg = "white"
)

ggsave(
  "figures/Fig10B_EKC_Panel_TP_ShiftedMean_COLORS.pdf",
  p10_B_col,
  width = 6.8,
  height = 4.2,
  units = "in"
)

# EPS (opcional, journal requirement)
ggsave(
  filename = "figures/Fig10B_EKC_Panel_TP_ShiftedMean_COLORS.eps",
  plot = p10_B_col,
  width = 6.8,
  height = 4.2,
  units = "in",
  device = cairo_ps
)
# =========================================================
# FIGURE 11 — Environmental Kuznets Curves by country (MG)
# =========================================================

library(tidyverse)
library(scales)

# -----------------------------
# 1) Cargar turning points por país
# -----------------------------
tp_country <- readr::read_csv(
  "tables/Table_TurningPoints_byCountry_bootstrap.csv",
  show_col_types = FALSE
)

stopifnot(all(c("iso3c","b1","b2","gdp_min","gdp_max",
                "tp_hat","within_support_country") %in% names(tp_country)))

# -----------------------------
# 2) Anclas por país (datos reales)
# -----------------------------
anchors <- panel_brics %>%
  group_by(iso3c) %>%
  summarise(
    mu_ln_gdp  = mean(ln_GDPpc,  na.rm = TRUE),
    mu_ln_co2  = mean(ln_CO2pc,  na.rm = TRUE),
    .groups = "drop"
  )

# -----------------------------
# 3) Construir curvas EKC por país
# -----------------------------
df_curve_country <- tp_country %>%
  left_join(anchors, by = "iso3c") %>%
  mutate(grid = map2(gdp_min, gdp_max, ~ seq(.x, .y, length.out = 250))) %>%
  unnest(grid) %>%
  group_by(iso3c) %>%
  mutate(
    x_c   = log(grid) - mu_ln_gdp,
    shape = b1 * x_c + b2 * (x_c^2),
    yhat  = shape - (b1*0 + b2*0) + mu_ln_co2
  ) %>%
  ungroup() %>%
  rename(GDP_pc = grid)

# -----------------------------
# 4) Turning points válidos (solo dentro de soporte)
# -----------------------------
tp_valid <- tp_country %>%
  left_join(anchors, by = "iso3c") %>%
  filter(within_support_country) %>%
  mutate(
    tp_label = paste0("TP = ", dollar(round(tp_hat, 0)))
  )

# Orden consistente
df_curve_country$iso3c <- factor(df_curve_country$iso3c,
                                 levels = c("BRA","CHN","IND","RUS","ZAF"))
tp_valid$iso3c <- factor(tp_valid$iso3c,
                         levels = c("BRA","CHN","IND","RUS","ZAF"))

# -----------------------------
# 5) Figura
# -----------------------------
p11 <- ggplot(df_curve_country, aes(x = GDP_pc, y = yhat)) +
  geom_line(color = "black", linewidth = 1) +
  geom_vline(
    data = tp_valid,
    aes(xintercept = tp_hat),
    linetype = "dashed", linewidth = 0.7
  ) +
  geom_text(
    data = tp_valid,
    aes(x = tp_hat, y = Inf, label = tp_label),
    vjust = 1.2, hjust = 0, size = 3.1
  ) +
  facet_wrap(~ iso3c, scales = "free_x", ncol = 3) +
  scale_x_continuous(labels = dollar) +
  labs(
    title = "Environmental Kuznets Curves by country and implied turning points",
    subtitle = "Mean Group estimates. Turning points shown only when within country-specific support.",
    x = "GDP per capita (PPP, constant international $)",
    y = "ln(CO\u2082 per capita)"
  ) +
  theme_minimal(base_size = 11, base_family = "serif") +
  theme(
    plot.title = element_text(face = "bold"),
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    plot.margin = margin(10,10,10,10)
  )

# -----------------------------
# 6) Exportar
# -----------------------------
if (!dir.exists("figures")) dir.create("figures")

ggsave(
  "figures/Fig11_EKC_MG_byCountry_TurningPoints.png",
  p11, width = 9.2, height = 4.8, dpi = 300, bg = "white"
)

ggsave(
  "figures/Fig11_EKC_MG_byCountry_TurningPoints.pdf",
  p11, width = 9.2, height = 4.8
)

# -----------------------------
# 4) Turning points: dentro vs fuera del soporte
# -----------------------------
tp_all <- tp_country %>%
  left_join(anchors, by = "iso3c") %>%
  mutate(
    tp_label = if_else(
      within_support_country,
      paste0("TP = ", dollar(round(tp_hat, 0))),
      paste0("TP (out) = ", dollar(round(tp_hat, 0)))
    ),
    tp_class = if_else(within_support_country, "Within support", "Out of support"),
    tp_class = factor(tp_class, levels = c("Within support", "Out of support"))
  )

# Orden consistente
df_curve_country$iso3c <- factor(df_curve_country$iso3c,
                                 levels = c("BRA","CHN","IND","RUS","ZAF"))
tp_all$iso3c <- factor(tp_all$iso3c,
                       levels = c("BRA","CHN","IND","RUS","ZAF"))

# -----------------------------
# 5) Figura (con TP para todos)
# -----------------------------
p11 <- ggplot(df_curve_country, aes(x = GDP_pc, y = yhat)) +
  geom_line(color = "black", linewidth = 1) +
  
  # TP dentro y fuera del soporte con estilos distintos
  geom_vline(
    data = tp_all,
    aes(xintercept = tp_hat, linetype = tp_class),
    linewidth = 0.7,
    color = "black"
  ) +
  scale_linetype_manual(values = c("Within support" = "dashed",
                                   "Out of support" = "dotted")) +
  
  # Etiquetas: arriba si está dentro; abajo si está fuera (para no estorbar)
  geom_text(
    data = tp_all %>% filter(within_support_country),
    aes(x = tp_hat, y = Inf, label = tp_label),
    vjust = 1.2, hjust = 0, size = 3.1
  ) +
  geom_text(
    data = tp_all %>% filter(!within_support_country),
    aes(x = tp_hat, y = -Inf, label = tp_label),
    vjust = -1.0, hjust = 0, size = 3.0
  ) +
  
  facet_wrap(~ iso3c, scales = "free_x", ncol = 3) +
  scale_x_continuous(labels = dollar) +
  labs(
    title = "Environmental Kuznets Curves by country and implied turning points",
    subtitle = "Mean Group estimates. Turning points inside support (dashed) vs outside support (dotted).",
    x = "GDP per capita (PPP, constant international $)",
    y = "ln(CO\u2082 per capita)",
    linetype = NULL
  ) +
  theme_minimal(base_size = 11, base_family = "serif") +
  theme(
    plot.title = element_text(face = "bold"),
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    plot.margin = margin(10,10,10,10)
  )

# -----------------------------
# 6) Exportar
# -----------------------------
if (!dir.exists("figures")) dir.create("figures")

ggsave("figures/Fig11B_EKC_MG_byCountry_TurningPoints.png",
       p11, width = 9.2, height = 4.8, dpi = 300, bg = "white")

ggsave("figures/Fig11B_EKC_MG_byCountry_TurningPoints.pdf",
       p11, width = 9.2, height = 4.8)

# =========================================================
# FIGURE A — Income–CO2 relationship by country (trajectories, facets)
# BRICS 1990–2024 | panel_brics already in memory
# Adds key-year labels: 1990, 2000, 2010, 2020, 2024
# =========================================================

library(tidyverse)
library(scales)

# -----------------------------
# 0) Paleta BRICS (bandera-ish)
# -----------------------------
pal_brics <- c(
  BRA = "#1B9E77",  # Brazil green
  CHN = "#D73027",  # China red
  RUS = "#4575B4",  # Russia blue
  IND = "#F28E2B",  # India orange/saffron
  ZAF = "#006837"   # South Africa darker green
)

# -----------------------------
# 1) Datos
# -----------------------------
stopifnot(all(c("iso3c","year","GDP_pc","ln_CO2pc","ln_GDPpc") %in% names(panel_brics)))

df_ic <- panel_brics %>%
  filter(iso3c %in% names(pal_brics)) %>%
  filter(!is.na(GDP_pc), !is.na(ln_CO2pc), !is.na(year)) %>%
  mutate(
    iso3c = factor(iso3c, levels = c("BRA","CHN","IND","RUS","ZAF"))
  ) %>%
  arrange(iso3c, year)

# Años clave a etiquetar (se ajusta automáticamente si no existe alguno)
key_years <- c(1990, 2000, 2010, 2020, 2024)

df_key <- df_ic %>%
  filter(year %in% key_years) %>%
  group_by(iso3c, year) %>%
  slice(1) %>%  # por si hubiera duplicados
  ungroup()

# Inicio/fin por país
df_start <- df_ic %>% group_by(iso3c) %>% slice_min(year, n = 1) %>% ungroup()
df_end   <- df_ic %>% group_by(iso3c) %>% slice_max(year, n = 1) %>% ungroup()

# -----------------------------
# 2) Tema “journal”
# -----------------------------
theme_journal <- function(base_size = 11, base_family = "serif") {
  theme_minimal(base_size = base_size, base_family = base_family) +
    theme(
      plot.title = element_text(face = "bold"),
      strip.text = element_text(face = "bold"),
      panel.grid.minor = element_blank(),
      legend.position = "none",
      plot.margin = margin(10,10,10,10)
    )
}

# -----------------------------
# 3) Plot (Figura A)
# -----------------------------
pA <- ggplot(df_ic, aes(x = GDP_pc, y = ln_CO2pc)) +
  
  # Trayectoria temporal (por país)
  geom_path(aes(color = iso3c, group = iso3c), linewidth = 0.95, alpha = 0.95) +
  
  # Puntos anuales (suaves)
  geom_point(aes(color = iso3c), size = 1.55, alpha = 0.55) +
  
  # Marcar inicio (círculo sólido) y fin (triángulo)
  geom_point(data = df_start, aes(color = iso3c), size = 2.8, shape = 16) +
  geom_point(data = df_end,   aes(color = iso3c), size = 2.8, shape = 17) +
  
  # Puntos para años clave (más visibles)
  geom_point(data = df_key, aes(color = iso3c), size = 2.2, shape = 21,
             fill = "white", stroke = 0.8) +
  
  # Etiquetas de años clave (pequeñas, sin fondo pesado)
  geom_text(
    data = df_key,
    aes(label = year),
    size = 3.0,
    nudge_y = 0.04,
    color = "black"
  ) +
  
  facet_wrap(~ iso3c, scales = "free_x", ncol = 3) +
  scale_color_manual(values = pal_brics) +
  scale_x_continuous(labels = dollar) +
  labs(
    title = "Income–CO₂ relationship by country (trajectories)",
    subtitle = "Dots are yearly observations; lines connect years. Key years labeled (1990, 2000, 2010, 2020, 2024).",
    x = "GDP per capita (PPP, constant international $)",
    y = "ln(CO\u2082 per capita)"
  ) +
  theme_journal()

# -----------------------------
# 4) Exportar
# -----------------------------
if (!dir.exists("figures")) dir.create("figures")

ggsave("figures/Fig9A_Income_CO2_Trajectories_Facets_KeyYears.png",
       pA, width = 9.2, height = 4.8, dpi = 300, bg = "white")

ggsave("figures/Fig9A_Income_CO2_Trajectories_Facets_KeyYears.pdf",
       pA, width = 9.2, height = 4.8)
