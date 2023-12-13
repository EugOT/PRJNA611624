library(tidyverse)
metrics_summary_SRR11278029 <- read_csv("/data/PRJNA611624/cellranger/SRR11278029/outs/metrics_summary.csv") %>% mutate(Run = "SRR11278029")
metrics_summary_SRR11278030 <- read_csv("/data/PRJNA611624/cellranger/SRR11278030/outs/metrics_summary.csv") %>% mutate(Run = "SRR11278030")
metrics_summary_SRR11278031 <- read_csv("/data/PRJNA611624/cellranger/SRR11278031/outs/metrics_summary.csv") %>% mutate(Run = "SRR11278031")
metrics_summary_SRR11278032 <- read_csv("/data/PRJNA611624/cellranger/SRR11278032/outs/metrics_summary.csv") %>% mutate(Run = "SRR11278032")
metrics_summary <-
    bind_rows(
        metrics_summary_SRR11278029,
        metrics_summary_SRR11278030,
        metrics_summary_SRR11278031,
        metrics_summary_SRR11278032)

metrics_summary |>
    select("Estimated Number of Cells", "Run")

write_tsv(metrics_summary, "/data/PRJNA611624/metrics_summary.tsv")

