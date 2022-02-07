library(vegan)

samples_per_pool <- read_csv("sample_per_pool.csv")
samples_per_pool$Interval <- as.character(samples_per_pool$Interval)
samples_per_pool <-
  mutate(samples_per_pool, id = paste(Species, Interval, sep=""))
  


get_abundance <- function(df) {
  
  df <- 
    mutate(df, id = paste(Species, Interval, sep="")) %>%
    mutate(df, freq = viral_read_count/Sequences) %>%
    mutate(df, present = ifelse(freq>0, 1, 0)) %>%
    mutate(Interval=as.character(Interval)) %>%
    left_join(samples_per_pool) %>%
    mutate(Interval=as.integer(Interval)) %>%
    left_join(total_virus_reads) %>%
    mutate(df, normalised_reads = viral_read_count/`Number of individuals`)
  
  return(df)
  
}


get_comm_mat <- function(df) {
  
  df <- get_abundance(df)
  
  comm_mat <-
    df %>%
    select(virus_family, id, normalised_reads) %>%
    spread(virus_family, normalised_reads)
  
  comm_mat[is.na(comm_mat)] <- 0
  
  return(comm_mat)
}

get_comm_mat_s <- function(df) {
  
  df <- get_abundance(df)
  
  comm_mat <-
    df %>%
    select(virus_species, id, normalised_reads) %>%
    spread(virus_species, normalised_reads)
  
  comm_mat[is.na(comm_mat)] <- 0
  
  return(comm_mat)
}

get_comm_mat_n <- function(df) {
  
  df <- get_abundance(df)
  
  comm_mat <-
    df %>%
    select(virus_species, id, viral_read_count) %>%
    spread(virus_species, viral_read_count)
  
  comm_mat[is.na(comm_mat)] <- 0
  
  return(comm_mat)
}

get_comm_mat_c <- function(df) {
  
  df <- get_abundance(df)
  
  comm_mat <-
    df %>%
    select(contig_ID, id, normalised_reads) %>%
    spread(contig_ID, normalised_reads)
  
  comm_mat[is.na(comm_mat)] <- 0
  
  return(comm_mat)
}


midpoints = c("31/01/2017", "31/03/2017", "15/06/2017","15/09/2017", "12/12/2017")
midpoints[3:5]

get_diversity_stats <- function(comm_mat) {
  
  length <- length(comm_mat[1,])
  matrix <- comm_mat[, 2:length]
  species <- c(rep("AF", 3), rep("AS", 5), rep("MG", 5))
  time <- c(midpoints[3:5], midpoints, midpoints)
  H <- diversity(matrix, index = "shannon")
  species_richness <- specnumber(matrix)
  
  diversity.df <-
    data.frame(
      Interval = (time),
      Species = species,
      H = H,
      species_richness = species_richness
    )
  }

get_NMDS_dist <- function(comm_mat, distance) {
  
  length <- length(comm_mat[1, ])
  matrix = as.matrix(comm_mat[, 2:length])
  rand <- as.integer(round(runif(1, 1, 10000)))
  set.seed(rand)
  nmds = metaMDS(matrix, distance = distance, try = 20)
  
  species <- c(rep("AF", 3), rep("AS", 5), rep("MG", 5))
  time <- c(3:5, 1:5, 1:5)
  
  data.scores = as.data.frame(scores(nmds))
  data.scores$Sample = comm_mat$id
  data.scores$Time = dmy(time)
  data.scores$Species = species
  
  return(data.scores)
  
}

plot_diversity_stats <- function(comm_mat) {

  diversity <- get_diversity_stats(comm_mat)
  diversity$Time <- c(midpoints[3:5], rep(midpoints, 2))
  diversity.m <- melt(diversity, id.vars = c("Interval", "Species"))
  
  
  p1 <- ggplot(diversity, aes(dmy(Interval), species_richness, group=Species, colour=Species))+
    geom_point()+
    geom_line()+
    facet_wrap(~Species)+
    scale_color_brewer(palette="Set2")+
    theme_bw()+
    ylab("Species richness")+
    scale_x_date(limits =c(dmy("01/01/2017"), dmy("01/01/2018")))+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))+
    xlab("Time")
  
  p2 <- ggplot(diversity, aes(dmy(Interval), H, group=Species, colour=Species))+
    geom_point()+
    geom_line()+
    facet_wrap(~Species)+
    scale_color_brewer(palette="Set2")+
    theme_bw()+
    ylab("Shannon Diversity")+
    scale_x_date(limits =c(dmy("01/01/2017"), dmy("01/01/2018")))+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))+
    xlab("Time")
  
  figure <- ggarrange(p1, p2, align="hv", common.legend = TRUE, ncol =1, legend="right")
  
  return(figure)
  
}
