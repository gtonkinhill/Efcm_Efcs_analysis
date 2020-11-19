library(data.table)
library(tidyverse)
library(micropan)

df <- fread("~/Downloads/gene_presence_absence.csv") %>% as_tibble()
pa <- data.matrix(df[,15:ncol(df)]!='')

heap_df <- map_dfr(1:1000, ~{
  ppa <- pa[sample(nrow(pa), replace = FALSE), sample(ncol(pa), replace = FALSE)]
  cumlative <- rowSums(apply(ppa, 1, cumsum)>0)
  cumlative <- cumlative-cumlative[[1]]
  df <- tibble(N = 1:length(cumlative),
               naccessory = cumlative)
  res <- lm(nunique ~ logN, tibble(logN = log(1:length(cumlative)),
                                   nunique = log(cumlative+0.001))) 
  return(df %>% 
           add_column(logK=res$coefficients[[1]]) %>% 
           add_column(beta=res$coefficients[[2]]) %>%
           add_column(permutation=.x, .before = 1))
})

coefs <- heap_df[!duplicated(heap_df$permutation),]
q <- quantile(coefs$beta, c(0.5, 0.025, 0.975))


plotdf <- heap_df %>% group_by(N) %>%
  summarise(
    `accessory size` = mean(naccessory),
    std = sd(naccessory)
  )


ggplot(plotdf, aes(N, `accessory size`)) + 
  geom_ribbon(aes(ymin = `accessory size` - std,
                  ymax = `accessory size` + std),
              fill='#67a9cf', alpha=0.5) + 
  geom_line(size = 1, col='#b2182b') +
  theme_bw(base_size = 14) +
  xlab("Number of genomes") +
  ylab("Accessory size") +
  geom_label(x = 5, y = 90, 
             label = sprintf("α = %.2f (%.2f, %.2f)", q[[1]], q[[2]], q[[3]]),
             size=5)
  
ggsave("~/Downloads/E_faecalis_pangenome_alpha.pdf", width = 12, height = 7, device = cairo_pdf)
ggsave("~/Downloads/E_faecalis_pangenome_alpha.png", width = 12, height = 7)

