
mdp <- read_table("mean_depth.txt", col_names = FALSE)
scafflen <- read_table("scaffold_lengths.txt", col_names = FALSE) %>% rename(length = X2)

gs <- read_table("genome_scaffolds.txt", col_names = FALSE)
mdp %>% filter(X4 == 0)
diff <- setdiff(gs$X1, mdp$X1)

nsite_per_scaff <- mdp %>% group_by(X1) %>% summarize(n = n(), meanDP = mean(X4))
df <- 
  left_join(nsite_per_scaff, scafflen) %>% 
  mutate(perc = n/length) %>% 
  arrange(perc)
write_tsv(df, "scaff_coverage.txt")

df %>% filter(perc < 0.001) %>% rename(scaffold = X1) %>% select(scaffold) %>% write_csv("low_genotype_scaffolds.csv")
write_csv(data.frame(scaffold = diff), "zero_genotype_scaffolds.csv")
