# Load required libraries

if(!require(dplyr)){
  install.packages(pkgs = 'dplyr', repos = 'https://stat.ethz.ch/CRAN/')
  library(dplyr)
}

if(!require(tidyr)){
  install.packages(pkgs = 'tidyr', repos = 'https://stat.ethz.ch/CRAN/')
  library(tidyr)
}

if(!require(readr)){
  install.packages(pkgs = 'tidyr', repos = 'https://stat.ethz.ch/CRAN/')
  library(readr)
}

if(!require(stringr)){
  install.packages(pkgs = 'stringr', repos = 'https://stat.ethz.ch/CRAN/')
  library(stringr)
}

if(!require(ggplot2)){
  install.packages(pkgs = 'ggplot2', repos = 'https://stat.ethz.ch/CRAN/')
  library(ggplot2)
}

if(!require(GUniFrac)){
  install.packages(pkgs = 'GUniFrac', repos = 'https://stat.ethz.ch/CRAN/')
  library(GUniFrac)
}

root <- "../.."
input.strains <- file.path(root, "results", "quantify_strains", "strain_quantification_matrix_count.tsv")
input.readcounts <- file.path(root, "results", "assign_taxonomy", "read_count_wide.tsv")
input.pal <- file.path(root, "results", "quantify_strains", "color_palette.tsv")
out.comp <- file.path(root, "analysis_results")
out.plots <- file.path(root, "analysis_plots")

### Create outdirs ###

cat("\nCreating directories\n")

dir.create(out.comp, recursive = TRUE, showWarnings = FALSE)
dir.create(out.plots, recursive = TRUE, showWarnings = FALSE)

### Inputs ###

cat("Reading inputs\n")

strains <- read.table(input.strains, sep = "\t", row.names = 1, check.names = F)
readcounts <- read_tsv(input.readcounts, show_col_types = F)
pal <- read_tsv(input.pal, show_col_types = F)

### Compare read count to strain count ###

strain_counts <- data.frame(
  SampleID = colnames(strains),
  genome_depth = colSums(strains)
) %>% 
  left_join(readcounts, by = "SampleID")

plot_counts <- ggplot(
  strain_counts,
  aes(
    x = AssignTaxonomy/1e3,
    y = genome_depth/1e3
  )
) +
  geom_smooth(method = 'lm', linetype = 2, linewidth = 0.2, color = "grey20") +
  geom_point(alpha=0.7) +
  scale_x_log10(breaks = c(30,seq(50,400,50))) +
  scale_y_log10(breaks = seq(10,110,10)) +
  labs(
    x = expression(paste("# reads after taxonomy (", 10^3, ")")),
    y = expression(paste("# genome equivalents (", 10^3, ")")),
  ) +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank()
  ) +
  coord_equal()
plot_counts

ggsave(file.path(out.plots, "read_genomes_counts.pdf"), plot_counts, device="pdf", width = 6, height = 6)

lm <- lm(genome_depth ~ AssignTaxonomy, strain_counts)
summary(lm)
qqnorm(lm$residuals)
qqline(lm$residuals) # not too skewed
# genome_depth = 0.276*reads

# each variable is a good proxy for the other


### Rarefy at different depths in genome equivalents ###

depths <- c(seq(1,9,2),seq(10,70,5))*1e3
set.seed(42)

rarefy_long <- function(mat=strains, rarefy_to, n_boot=50){
  avail_depth <- colSums(mat)
  if (length(avail_depth[avail_depth > rarefy_to]) < 2){
    cat(rarefy_to, ": requested depth is too high\n")
    return(NULL)
  } else {
    cat(rarefy_to, "\n")
    boot_list <- lapply(seq_len(n_boot), function(b) {
      
      raref <- t(Rarefy(t(mat), depth = rarefy_to)$otu.tab.rff)
      
      as.data.frame(raref) %>%
        mutate(
          Strain = rownames(raref),
          genome_depth_rarefaction = rarefy_to,
          bootstrap = b,
          .before = 1
        ) %>%
        pivot_longer(
          colnames(raref),
          names_to = "SampleID",
          values_to = "strain_genome_equivalents"
        )
    })
    
    do.call(rbind, boot_list)
  }
}

list_subsamples <- lapply(
  depths,
  function(x) rarefy_long(mat = strains, rarefy_to = x, n_boot = 100)
)

subsamples <- do.call(rbind, list_subsamples)

# order strains
subsamples$Strain <- factor(subsamples$Strain, levels = unique(pal$Strain), ordered = T)
# order samples by type, then increasing depth
readcounts <- readcounts %>% arrange(SampleType, AssignTaxonomy)
subsamples$SampleID <- factor(subsamples$SampleID, levels = readcounts$SampleID, ordered = T)

subsamples_summary <- subsamples %>%
  group_by(
    SampleID,
    Strain,
    genome_depth_rarefaction
  ) %>%
  summarise(
    median_ab = median(strain_genome_equivalents),
    lo = quantile(strain_genome_equivalents, 0.025),
    hi = quantile(strain_genome_equivalents, 0.975),
    sd = sd(strain_genome_equivalents),
    .groups = "drop"
  )

## plot

strain_ab_subsamples1 <- ggplot(
  subsamples_summary,
  aes(
    x = genome_depth_rarefaction/1e3,
    y = median_ab,
    color = Strain,
    fill = Strain
  )
) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.25, color = NA) +
  geom_line(linewidth = 0.5) +
  scale_x_continuous(breaks=seq(0,60,10)) +
  scale_y_log10(breaks=10^seq(0,5,1)) +
  scale_fill_manual(values = setNames(pal$color,pal$Strain), labels = setNames(pal$label,pal$Strain)) +
  scale_color_manual(values = setNames(pal$color,pal$Strain), labels = setNames(pal$label,pal$Strain)) +
  geom_vline(xintercept = 0.276*20*0.73, linetype = 2, linewidth = 0.4, color = "grey20") +
  geom_vline(xintercept = 0.276*40*0.73, linetype = 2, linewidth = 0.4, color = "grey20") +
  geom_vline(xintercept = 0.276*60*0.73, linetype = 2, linewidth = 0.4, color = "grey20") +
  labs(
    x = expression(paste("Subsample size (", 10^3, ")")),
    y = "Strain abundance (median ± 95% CI)"
  ) +
  theme_bw() +
  theme(
    legend.position = "top",
    legend.key.size = unit(0.6, 'cm'),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8),
    panel.grid.minor = element_blank()
  ) +
  facet_wrap(~ SampleID) +
  guides(colour = guide_legend(nrow = 5)) 
strain_ab_subsamples1

strain_ab_subsamples2 <- ggplot(
  subsamples_summary,
  aes(
    x = genome_depth_rarefaction/1e3,
    y = median_ab,
    color = Strain,
    fill = Strain
  )
) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.25, color = NA) +
  geom_line(linewidth = 0.5) +
  scale_x_log10() +
  scale_y_log10(breaks=10^seq(0,5,1)) +
  scale_fill_manual(values = setNames(pal$color,pal$Strain), labels = setNames(pal$label,pal$Strain)) +
  scale_color_manual(values = setNames(pal$color,pal$Strain), labels = setNames(pal$label,pal$Strain)) +
  geom_vline(xintercept = 0.276*20*0.73, linetype = 2, linewidth = 0.4, color = "grey20") +
  geom_vline(xintercept = 0.276*40*0.73, linetype = 2, linewidth = 0.4, color = "grey20") +
  geom_vline(xintercept = 0.276*60*0.73, linetype = 2, linewidth = 0.4, color = "grey20") +
  labs(
    x = expression(paste("Subsample size (", 10^3, ")")),
    y = "Strain abundance (median ± 95% CI)"
  ) +
  theme_bw() +
  theme(
    legend.position = "top",
    legend.key.size = unit(0.6, 'cm'),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8),
    panel.grid.minor = element_blank()
  ) +
  facet_wrap(~ SampleID) +
  guides(colour = guide_legend(nrow = 5)) 
strain_ab_subsamples2

ggsave(file.path(out.plots, "strain_ab_subsamples.pdf"), strain_ab_subsamples1, device="pdf", width = 12, height = 12)
ggsave(file.path(out.plots, "strain_ab_subsamples_loglog.pdf"), strain_ab_subsamples2, device="pdf", width = 12, height = 12)

## compute lm for each strain in each sample

lm_strain <- function(data, str_smpl){
  # extract data
  strain = str_smpl[[1]]
  sample = str_smpl[[2]]
  df <- data %>% 
    filter(Strain == strain & SampleID == sample)
  
  if (sum(df$strain_genome_equivalents) == 0){
    return(NULL)
  } 
  
  # run lm
  my_lm <- lm(strain_genome_equivalents ~ genome_depth_rarefaction, df)
  f <- summary(my_lm)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  # put everything in a df
  coefs <- data.frame(
    Strain = strain,
    SampleID = sample,
    intercept = unname(coef(my_lm)[1]),
    slope = unname(coef(my_lm)[2]),
    adj.r2 = summary(my_lm)$r.squared,
    pval = p
    )
  return(coefs)
}

combinations <- unique(subsamples[ ,c("Strain", "SampleID")])
combinations.list <- unname(split(combinations, seq(nrow(combinations))))

list_lms <- lapply(
  combinations.list,
  function(x) lm_strain(data = subsamples, str_smpl = x)
)

lms <- do.call(rbind, list_lms) %>% 
  mutate(log10_slope = log10(slope)) 

# order strains
lms$Strain <- factor(lms$Strain, levels = unique(pal$Strain), ordered = T)
lms$SampleID <- factor(lms$SampleID, levels = readcounts$SampleID, ordered = T)

# add final abundance at highest depth

final_ab <- subsamples_summary %>% 
  group_by(Strain, SampleID) %>% 
  summarize(
    max_median_ab = max(median_ab)
  ) %>% 
  mutate(log10_max_median_ab = log10(max_median_ab))

subsamples_summary <- subsamples_summary %>% 
  left_join(lms, by = c("SampleID", "Strain"))

lms <- lms %>% 
  left_join(final_ab, by=c("Strain", "SampleID"))

plot_lms <- ggplot(
  lms,
  aes(
    x = slope,
    y = adj.r2,
    color = max_median_ab
  )
) +
  geom_point(alpha=0.8) +
  scale_color_gradient(low="#09202E", high="#8ED1FB", trans = "log", breaks=10^c(1:4)) +
  scale_x_log10() +
  scale_y_continuous(breaks=seq(0,1,0.1)) +
  labs(
    x = "Slope",
    y = expression(R^2),
    color = "Median abundance at highest depth"
  ) +
  theme_bw() + 
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "top",
    legend.text = element_text(angle=40, hjust=0.7, vjust=0.8)
  )
plot_lms

ggsave(file.path(out.plots, "linear_models.pdf"), plot_lms, device="pdf", width = 8, height = 6)

# seems like both slope and r2 on their own is enough to distinguish contaminants

threshold = 0.8

subsamples_summary <- subsamples_summary %>%
  mutate(fit = if_else(adj.r2 >= threshold,"Linear fit", "Highly variable"))

strain_ab_subsamples3 <- ggplot(
  subsamples_summary %>% filter(!is.na(fit)),
  aes(
    x = genome_depth_rarefaction/1e3,
    y = median_ab,
    group = Strain
  )
) +
  geom_ribbon(aes(ymin = lo, ymax = hi, fill = fit), alpha = 0.25, color = NA) +
  geom_line(aes(color = fit), linewidth = 0.5) +
  scale_x_continuous(breaks=seq(0,60,10)) +
  scale_y_log10(breaks=10^seq(0,5,1)) +
  scale_color_manual(values=setNames(c("#00BFC4","#F8766D"),c("Linear fit", "Highly variable"))) +
  scale_fill_manual(values=setNames(c("#00BFC4","#F8766D"),c("Linear fit", "Highly variable"))) +
  geom_vline(xintercept = 0.276*20*0.73, linetype = 2, linewidth = 0.4, color = "grey20") +
  geom_vline(xintercept = 0.276*40*0.73, linetype = 2, linewidth = 0.4, color = "grey20") +
  geom_vline(xintercept = 0.276*60*0.73, linetype = 2, linewidth = 0.4, color = "grey20") +
  labs(
    x = expression(paste("Subsample size (", 10^3, ")")),
    y = "Strain abundance (median ± 95% CI)",
    color = "",
    fill = ""
  ) +
  theme_bw() +
  theme(
    legend.position = "top",
    legend.key.size = unit(0.6, 'cm'),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8),
    panel.grid.minor = element_blank()
  ) +
  facet_wrap(~ SampleID)
strain_ab_subsamples3

strain_ab_subsamples4 <- ggplot(
  subsamples_summary %>% filter(!is.na(fit)),
  aes(
    x = genome_depth_rarefaction/1e3,
    y = median_ab,
    group = Strain
  )
) +
  geom_ribbon(aes(ymin = lo, ymax = hi, fill = fit), alpha = 0.25, color = NA) +
  geom_line(aes(color = fit), linewidth = 0.5) +
  scale_x_continuous(breaks=seq(0,20,2), limits=c(0,20)) +
  scale_y_log10(breaks=10^seq(0,5,1)) +
  scale_color_manual(values=setNames(c("#00BFC4","#F8766D"),c("Linear fit", "Highly variable"))) +
  scale_fill_manual(values=setNames(c("#00BFC4","#F8766D"),c("Linear fit", "Highly variable"))) +
  geom_vline(xintercept = 0.276*20*0.73, linetype = 2, linewidth = 0.4, color = "grey20") +
  geom_vline(xintercept = 0.276*40*0.73, linetype = 2, linewidth = 0.4, color = "grey20") +
  geom_vline(xintercept = 0.276*60*0.73, linetype = 2, linewidth = 0.4, color = "grey20") +
  labs(
    x = expression(paste("Subsample size (", 10^3, ")")),
    y = "Strain abundance (median ± 95% CI)",
    color = "",
    fill = ""
  ) +
  theme_bw() +
  theme(
    legend.position = "top",
    legend.key.size = unit(0.6, 'cm'),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8),
    panel.grid.minor = element_blank()
  ) +
  facet_wrap(~ SampleID)
strain_ab_subsamples4

ggsave(file.path(out.plots, "strain_ab_subsamples_classified.pdf"), strain_ab_subsamples3, device="pdf", width = 12, height = 10)
ggsave(file.path(out.plots, "strain_ab_subsamples_classified_zoom.pdf"), strain_ab_subsamples4, device="pdf", width = 12, height = 10)

write.table(
  strain_counts,
  file.path(out.comp, "strain_read_counts.tsv"),
  sep = "\t",
  col.names = T,
  row.names = F,
  quote = F
)

write.table(
  lms,
  file.path(out.comp, "linear_models.tsv"),
  sep = "\t",
  col.names = T,
  row.names = F,
  quote = F
)

write.table(
  subsamples,
  file.path(out.comp, "subsamples_table.tsv"),
  sep = "\t",
  col.names = T,
  row.names = F,
  quote = F
)

write.table(
  subsamples_summary,
  file.path(out.comp, "subsamples_summary.tsv"),
  sep = "\t",
  col.names = T,
  row.names = F,
  quote = F
)

write.table(
  lms,
  file.path(out.comp, "linear_models.tsv"),
  sep = "\t",
  col.names = T,
  row.names = F,
  quote = F
)
