
rm(list = ls())

# ---------------------------------------------------------------------------
# NHANES August 2021 - August 2023  ("_L" files)
# Build an analysis frame for an additive model with Student-t errors:
#   log(triglycerides) ~ s(age) + s(BMI) + s(waist) + s(HDL) + s(glucose)
#                        + s(systolic BP) + s(income-poverty ratio)
#
# Requires: haven, dplyr  (optionally mgcv for the model at the end)
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# NHANES August 2021 - August 2023  ("_L" files)
# Build an analysis frame for an additive model with Student-t errors:
#   log(triglycerides) ~ s(age) + s(BMI) + s(waist) + s(HDL) + s(glucose)
#                        + s(systolic BP) + s(income-poverty ratio)
#
# Requires: haven, dplyr  (optionally mgcv for the model at the end)
# ---------------------------------------------------------------------------

library(haven)
library(dplyr)

base <- "https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2021/DataFiles/"

get_xpt <- function(name, ext = ".xpt") {
  url <- paste0(base, name, ext)
  message("downloading ", name)
  read_xpt(url)
}

# --- 1. Download components ------------------------------------------------
demo   <- get_xpt("DEMO_L")     # demographics + weights + masked design vars
bmx    <- get_xpt("BMX_L")      # body measures
trigly <- get_xpt("TRIGLY_L")   # LDL + triglycerides (fasting subsample)
hdl    <- get_xpt("HDL_L")      # HDL cholesterol
glu    <- get_xpt("GLU_L")      # fasting plasma glucose
bpxo   <- get_xpt("BPXO_L")     # oscillometric blood pressure
fastqx <- get_xpt("FASTQX_L")   # fasting questionnaire (fasting time)

# --- 1b. ALWAYS inspect names before selecting -----------------------------
# Variable names are NOT stable across cycles. In August 2021-August 2023 the
# triglyceride variable was renamed LBXTR -> LBXTLG because the glycerol-blanked
# assay was replaced. Run this before trusting any name below.
for (nm in c("demo","bmx","trigly","hdl","glu","bpxo","fastqx")) {
  cat("\n==", nm, "==\n"); print(names(get(nm)))
}

# --- 2. Merge on SEQN ------------------------------------------------------
# NOTE: verify each variable name against its own codebook page before relying
# on it. Names are stable across cycles for the core variables below, but the
# fasting-subsample weight in particular changes name between cycles.
dat <- demo %>%
  select(SEQN, RIDSTATR, RIAGENDR, RIDAGEYR, RIDRETH3, RIDEXMON,
         DMDBORN4, RIDEXPRG, INDFMPIR, WTMEC2YR, SDMVSTRA, SDMVPSU) %>%
  left_join(select(bmx,    SEQN, BMXBMI, BMXWAIST),        by = "SEQN") %>%
  left_join(select(trigly, SEQN, LBXTLG),                   by = "SEQN") %>%
  left_join(select(hdl,    SEQN, LBDHDD),                  by = "SEQN") %>%
  left_join(select(glu,    SEQN, LBXGLU),                  by = "SEQN") %>%
  left_join(select(bpxo,   SEQN, BPXOSY1, BPXOSY2, BPXOSY3), by = "SEQN") %>%
  left_join(select(fastqx, SEQN, PHAFSTHR),                by = "SEQN")

# Average the three oscillometric systolic readings
dat <- dat %>%
  mutate(SBP = rowMeans(cbind(BPXOSY1, BPXOSY2, BPXOSY3), na.rm = TRUE),
         SBP = ifelse(is.nan(SBP), NA, SBP))


# --- 3. Filters, withdat# --- 3. Filters, with a running count -------------------------------------
report <- function(d, label) { cat(sprintf("%-42s n = %5d\n", label, nrow(d))); d }

analysis <- dat %>%
  report("all participants") %>%
  filter(RIDSTATR == 2)                    %>% report("MEC examined") %>%
  filter(!is.na(LBXTLG))                    %>% report("has triglycerides (fasting sub)") %>%
  filter(RIDAGEYR >= 40)   %>% report("adults 40+") %>%
  filter(PHAFSTHR >= 12, PHAFSTHR <= 24)    %>% report("fasted 12-24 h") %>%
  filter(RIAGENDR == 1)                    %>% report("+ male only") %>%
  filter(DMDBORN4 == 1)                  %>% report("+ US-born") %>%
  # -------------------------------------------------------------------------
select(SEQN, LBXTLG, RIDAGEYR, BMXBMI, BMXWAIST, LBDHDD, LBXGLU, SBP, INDFMPIR,
       RIDRETH3, RIAGENDR) %>%
  na.omit()                                %>% report("complete cases (final)")

analysis <- analysis %>%
  mutate(log_tg = log(LBXTLG)) %>%
  rename(age = RIDAGEYR, bmi = BMXBMI, waist = BMXWAIST,
         hdl = LBDHDD, glucose = LBXGLU, sbp = SBP, pir = INDFMPIR,
         ethn = RIDRETH3, gender = RIAGENDR)

cat("\nfinal n =", nrow(analysis), " p = 7  ->  n/p =",
    round(nrow(analysis) / 7, 1), "\n")

df <- as.data.frame(analysis)
df$ethn <- as.factor(df$ethn)
contrasts(df$ethn) <- contr.treatment(levels(df$ethn), base = 3) 

save(df, file = "appl_data.RData")
