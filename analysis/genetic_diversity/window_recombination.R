library(data.table); library(dplyr)


ld <- fread(cmd = "zcat chr1.fastld.ld.gz",
            select = c("SNP_A", "BP_A", "R2")) %>%
  mutate(across(c(BP_A, R2), as.numeric))

# Per-SNP local LD score (mean r² of kept neighbors)
ld_per_snp <- ld %>%
  group_by(SNP_A, BP_A) %>%
  summarize(ld_score = mean(R2, na.rm = TRUE), .groups = "drop")


# Bin along chromosome
bin_size <- 100000
ld_binned <- ld_per_snp %>%
  mutate(bin = floor(BP_A / bin_size) * bin_size) %>%
  group_by(bin) %>%
  summarize(mean_ld = mean(ld_score, na.rm = TRUE),
            n_snps  = n(), .groups = "drop")

# Quick look
library(here)
png(here("TEMP.png"))
plot(ld_binned$bin/1e6, ld_binned$mean_ld, type="l",
     xlab="Chr1 position (Mb)", ylab="Mean local LD (r²)")
dev.off()
