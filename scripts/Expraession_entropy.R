library(ggpubr)
library(dplyr)
library(purrr)
library(BioQC)
library(ggplot2)
library(data.table)

args <- commandArgs(TRUE)

entrop <- function(a) {
  ent <- entropySpecificity(a, norm=TRUE)
  return(ent)
}

trans <- function(a) {
  a$target_id <- NULL
  a$type <- NULL
  a$sum <- NULL
  entrop(a)
}

type <- function(a) {
  col_name <- colnames(a)
  col_name <- col_name[! col_name %in% c("type", "target_id")]
  a$sum <-  rowSums(a[,col_name], na.rm=TRUE)
  a <- filter(a, sum >0)
  ent <- trans(a)
  a$entropy <- ent
  return(a)
}

plot <- function(a,b,c) {
  ggplot() +
    geom_density(data = a, aes(x = entropy, fill="Non"), alpha=0.4)+
    geom_density(data = b, aes(x = entropy, fill="Cons"), alpha=0.4)+
    geom_density(data = c, aes(x = entropy, fill="mRNA"), alpha=0.4)+
    scale_colour_manual("", 
                        breaks = c("Non", "Cons", "mRNA")) +
    xlab('entropy') + ylab('density') + guides(fill = guide_legend(title = "Types")) + 
    theme(text = element_text(size = 20))
  ggsave("data/output/expression/entropy.png", width = 20, height = 20)
}

log <- function(a) {
       col_name <- colnames(a)
       col_name <- col_name[! col_name %in% c("type", "target_id")]  
       col <- c()
       for (x in col_name) {
            col <- append(col, paste("log_", noquote(x), sep = ""))
            a[, paste("log_", noquote(x), sep = "")] <- log2(a[, noquote(x)])
       } 
       new_col <<- col
       return(a)
}

# install dataframe
df <- read.table(args[1], header = TRUE)

# log of columns
df_log <- log(df) 

# exprassin analysis
df_log[df_log == -Inf] <- 0
my_comparisons = list(c("consv", "non_cons"), c("consv", "prot"),  c("prot", "non_cons"))

ggboxplot(df_log, x="type", y= new_col , merge = "flip", ylab = "log2(TPM)", color = "type", palette = "jco")+  
  font("subtitle", size = 20)+
  font("caption", size = 20)+
  font("xlab", size = 20)+
  font("ylab", size = 20)+
  font("xy.text", size = 20) +theme(legend.text=element_text(size=20))
ggsave("data/output/expression/box_plot.png", width = 20, height = 20)

ggdensity(df_log, x = new_col , y = "..density..", facet.by = "type", merge = TRUE, xlab = "TPM", rug = TRUE, add = "median", palette = "jco") +  font("subtitle", size = 20)+
  font("caption", size = 20)+
  font("xlab", size = 20)+
  font("ylab", size = 20)+
  font("xy.text", size = 20) + theme(legend.text=element_text(size=20))
ggsave("data/output/expression/hist_expr.png", width = 20, height = 20)

# entropy analysis
non <-df[grepl("noncons", df$type),]
cons <-df[grepl("consv", df$type),]
mRNA <-df[grepl("mRNA", df$type),]

non_ent <- type(non)
cons_ent <- type(cons)
mRNA_ent <- type(mRNA)

plot(non_ent, cons_ent, mRNA_ent)




