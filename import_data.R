## Import Data

# Import Data
taxa_og <- read.csv(file = 'data/taxatable-filtered-absolute.tsv', check.names = FALSE, stringsAsFactors = FALSE, sep = "\t")
meta_og <- read.csv(file = 'data/meta.txt', check.names = FALSE, stringsAsFactors = FALSE, sep = "\t")

# Colors
diversigen_colors_extended <-c ("#1db17e","#0065bc","#ffc75c","#18181a","#aa6da3",
                                "#345511","#b4d2e7","#922D50", "#596157","#e3655b")
