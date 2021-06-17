require("dplyr")
require("emoa")
require("mco")
require("ggplot2")
library("plotly")

file_name <- './data/predictions_smith_waterman_PAM30_2021-06-16_23_18_07.csv'
                
df <- read.csv(file=file_name)

df_complemented_negated <- df %>%
  mutate(probability_complement = 1 - Probability,
         distance_neg = -pair_align_score) %>%
  select(Protein, Sequence, probability_complement, distance_neg, Info_peptide, Info_PepID) #%>% head() 

df_objectives <- df_complemented_negated %>%
  select(probability_complement, distance_neg) #%>% head() 

m_objectives <- t(as.matrix(df_objectives))

nondominated <- nondominated_points(m_objectives)

rank <- nds_rank(m_objectives) #, partial)

file_name_split <- unlist(strsplit(file_name, '_'))

pl <- ggplot(data=df_complemented_negated, aes(x=probability_complement, y=distance_neg)) +
  ggtitle("PAM30") +
  geom_point(aes(text=sprintf("Protein ID: %s Sequence: %s<br>Info_PepID: %s Info_peptide: %s", 
                              Protein, Sequence, 
                              Info_PepID, Info_peptide)))

# Plot Pareto fronts
for (idx in 1:max(rank)) {
  color <- 'blue'
  if (!idx %% 2) {
    color <- 'green'
  }
  
  # Prepare Pareto points
  pareto_front_points_indicies <- which(rank %in% c(idx))
  df_pareto_front_points <- df_objectives[pareto_front_points_indicies, 
                                          c('probability_complement', 'distance_neg')]
    
  # http://www.sthda.com/english/wiki/ggplot2-line-plot-quick-start-guide-r-software-and-data-visualization
  # http://www.sthda.com/english/wiki/ggplot2-add-straight-lines-to-a-plot-horizontal-vertical-and-regression-lines
  pl <- pl + geom_line(data=df_pareto_front_points, 
                       aes(x=probability_complement, y=distance_neg), 
                       color=color)
}

ggplotly(pl) 
