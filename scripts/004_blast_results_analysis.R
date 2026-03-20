library(tidyverse)
library(purrr)
library(ggplot2)
library(ggridges)
library(gggenes)
library(writexl)

## Before running the script, create a folder named "blast_results"
## and placed there all the blast searches results files


### Create a function to sort the blast output files
# 1 = congruent direction, 0 = non congruent
mitogenome_palindrome_analysis = function(blast_filename, 
                                          match_minsize = 100){
  calculate_direction_congruence = function(x){
    congruence_value = abs(sum(sign(x$sstart - x$send)))
    res = (congruence_value == nrow(x)) * 1
    return(res)
  }
  
  calculate_congruence_score = function(x){
    congruence_value = abs(sum(sign(x$sstart - x$send)))
    return(congruence_value)
  }
  
  # Load blast results
  blast_res = read.table(blast_filename,
                         sep = "\t",
                         header = FALSE)
  colnames(blast_res) = c("qseqid","sseqid","pident","length", "mismatch",
                          "gapopen", "qstart", "qend","sstart","send",
                          "evalue","bitscore","qlen", "slen")
  # Remove short matches
  blast_res = subset(blast_res, length >= match_minsize)
  
  # Create nested tibble 
  blast_res_nest = blast_res %>% 
    group_by(qseqid) %>% 
    nest() %>%
    mutate(n_matches = map_dbl(data, ~nrow(.x)),
           query_len = map_dbl(data, ~max(.x$qlen)),
           avg_id = map_dbl(data, ~mean(.x$pident)))
  
  blast_res_nest$individual = rep(gsub("_blastn.txt","",blast_filename, fixed = TRUE), nrow(blast_res_nest))
  
  
  # Subset blast results to keep reads with > 1 match
  blast_res_nest_sub = subset(blast_res_nest,
                              n_matches > 1)
  
  blast_res_nest_sub = blast_res_nest_sub %>%
    mutate(
           direction_con = map_dbl(data,
                                   ~calculate_direction_congruence(.x)),
           congr_score = map_dbl(data,
                                 ~calculate_congruence_score(.x)))
  
  return(list(blast_res_nest = blast_res_nest,
              blast_res_nest_sub = blast_res_nest_sub))
  
}


### Process the results for all the individuals
Acu_gio_RIT094 = mitogenome_palindrome_analysis("./blast_results/Acu_gio_RIT094_all_blastn.txt", 
                                                match_minsize = 100)
Acu_gio_RPL027 = mitogenome_palindrome_analysis("./blast_results/Acu_gio_RPL027_all_blastn.txt", 
                                                match_minsize = 100)
Acu_gio_RSE010 = mitogenome_palindrome_analysis("./blast_results/Acu_gio_RSE010_all_blastn.txt", 
                                                match_minsize = 100)
Acu_mec_RSE086 = mitogenome_palindrome_analysis("./blast_results/Acu_mes_RSE086_all_blastn.txt", 
                                                match_minsize = 100)
Acu_mec_RSE087 = mitogenome_palindrome_analysis("./blast_results/Acu_mes_RSE087_all_blastn.txt", 
                                                match_minsize = 100)

# Merge the blast results for all individuals in the same table
all_results = list(Acu_gio_RIT094, 
                   Acu_gio_RPL027,
                   Acu_gio_RSE010,
                   Acu_mec_RSE086,
                   Acu_mec_RSE087)

# Extract the table with all the reads matches for every individuals
table_allmatches = do.call(rbind, lapply(all_results, FUN = function(x) x[[1]]))

# Extract the table with only the reads with >1 match and with the direction congruence
table_multimatches = do.call(rbind, lapply(all_results, FUN = function(x) x[[2]]))
table_multimatches_unnested = unnest(table_multimatches, cols = data)

### Create a plot with some random palindromic reads

# Helper function
order_matches_table = function(x, query_id){
  start_tab = x[,c("qstart","qend")]
  start_mins_order = order(apply(start_tab, MARGIN = 1, FUN = min))
  x = x[start_mins_order,]
  x$queryID = rep(query_id, nrow(x))
  x = x %>%
    mutate(start2 = ifelse(sstart > send, qend, qstart),
           end2 = ifelse(sstart > send, qstart, qend))
  return(x)
}


# Format data
table_allmatches_red = subset(table_allmatches, n_matches >= 1)
set.seed = 1479034 # This is needed to be sure that every run will give the same exact figure
                   # As reads are randomly sampled
table_allmatches_red = table_allmatches_red[sample(1:nrow(table_allmatches_red), 
                                                   size = 20,
                                                   prob = table_allmatches_red$n_matches**4),]

data_plot_genes = table_allmatches_red %>%
  mutate(plotdata = map2(.x = data, 
                         .y = qseqid, 
                         ~order_matches_table(.x, .y))) %>%
  select(c("individual","plotdata")) %>%
  unnest(cols = "plotdata") %>%
  mutate(direction = ifelse(start2 > end2, "Rev", "For"))

# Make the plot and save it as pdf
ggplot(data_plot_genes) +
  theme_grey()+
  theme(panel.background = element_blank(), 
        panel.grid.major.y = element_blank(), 
        panel.grid.minor.y = element_blank(), 
        panel.grid.minor.x = element_blank(), 
        panel.grid.major.x = element_blank(), 
        axis.ticks.y = element_blank(), 
        axis.line.x = element_line(colour = "grey20", 
                                   linewidth = 0.5),
        axis.ticks.x = element_line(colour = "grey20",
                                    linewidth = 0.5),
        strip.text = element_blank(),
        strip.background = element_blank()) +
  geom_segment(aes(x = 0, xend = qlen, 
                   y = queryID, yend = queryID), 
               linewidth = 1,
               col = "grey20") +
  geom_gene_arrow(aes(y = queryID,
                      xmin = start2,
                      xmax = end2,
                      fill = direction),
                  show.legend = FALSE)+
  scale_fill_manual(values = c("Rev" = "#e9a3c9", 
                               "For" = "#a1d76a")) +
  facet_wrap(~ queryID, scales = "free_y", ncol = 1) +
  labs(x = "", y = "",
       title = "mitogenome matches on some random reads")
                                           
ggsave("reads_allmito.pdf",
       height = 6*1.2,
       width = 12*1.2)



### Calculate statistics on number of matches reads and number of palindromic reads

# Get number of matched reads
matched_reads_num = as.data.frame.matrix(table(table_allmatches[,c("individual","n_matches")]))
# Get number of reads that are palindromic
individual_directions = as.data.frame.matrix(table(table_multimatches[,c("individual","direction_con")]))
colnames(individual_directions)[1] = "palindrom_seqs"
individual_directions = individual_directions[,1]
# Get tot number of reads per individual
reads_num = as.data.frame(table(table_allmatches[,c("individual")]))
colnames(reads_num)[2] = "tot_reads"

# Merge the stats together
stats_data = cbind(reads_num,
                   individual_directions,
                   matched_reads_num)
colnames(stats_data)[3] = "num_palindrom_seqs"
write_xlsx(stats_data, "stats_data.xlsx")

