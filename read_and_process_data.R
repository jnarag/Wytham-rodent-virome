library(tidyverse)
library(stringr)
library(reshape2)
library(ggsci)
library(cowplot)
library(RColorBrewer)
library(ggpubr)


# blastx virus results table - blast information for each viral contig
blast_results_table <-
  read_csv("blastx-viruses-results.csv")


# blast results table filtered by blast score
blast_results_table.s <-
  subset(blast_results_table, blast_results_table$`e-value` < 1e-10)

# relational table of virus families and host range. Change the include column according to which viruses you're interested in. Currently "Y" corresponds to verterbrate viruses.
virus_families <- read_csv("families-host-range.csv")

# contigency table containing how many reads match a specific contig in all 13 samples
contingency_table <-
  read_csv("results-cleaning-blastx-contingency.csv")

annotate_and_filter_contingency_table <-
  function(contingency_table, blast_table) {
    contingency_table.s <-
      contingency_table %>%
      filter(`contig ID` %in% blast_table$`contig ID`)
    
    print(length(contingency_table.s$`contig ID`))
      
    taxonomy_list <- list()
    species_list <- list()
    type_list <- list()
    x <- list()
    
    for (i in 1:length(blast_table$`contig ID`)) {
      taxonomy_list[[i]] <- "NA"
      species_list[[i]] <- "NA"
      type_list[[i]] <- "NA"
      
      for (j in 1:length(virus_families$family)) {
        var <-
          str_detect(blast_table$taxonomy[i],
                     virus_families$family[j])
        
        if (var == TRUE) {
          type_list[[i]] <- virus_families$type[j]
          taxonomy_list[[i]] <- virus_families$family[j]
          
        }
        
        # find index corresponding to virus_family
        taxonomy_parts <-
          str_split(blast_table$taxonomy[i], ",")
        
        family_index <- 0
        
        for(f in 1:length(taxonomy_parts[[1]])) {
          
          if(taxonomy_parts[[1]][f] == virus_families$family[j]) {
            family_index = f
            species_list[[i]] = taxonomy_parts[[1]][f+1]
          }
        }
        
        
        
      }
      
      taxonomy_parts <-
        str_split(blast_table$taxonomy[i], ",")
      
      var2 <-
        str_detect(taxonomy_parts[[1]][2], "Viruses")
      
      if (var2 == TRUE) {
        
        # species_list[[i]] <-
        #   (taxonomy_parts[[1]][length(taxonomy_parts[[1]]) - 2])
        # print(taxonomy_parts[[1]][4])
        # taxonomy_list[[i]] <- blast_results_table.s$taxonomy[i]
        
        if (taxonomy_list[[i]] == "Tobaniviridae") {
          species_list[[i]] <-
            (taxonomy_parts[[1]][8])
        }
        
        if (taxonomy_list[[i]] == "Picornaviridae" ||
            taxonomy_list[[i]] == "Retroviridae" ||
            taxonomy_list[[i]] == "Dicistroviridae") {
          species_list[[i]] <-
            (taxonomy_parts[[1]][6])
          
        }
        if (taxonomy_list[[i]] == "Astroviridae" ||
            taxonomy_list[[i]] == "Myoviridae" ||
            taxonomy_list[[i]] == "Picobirnaviridae" ||
            taxonomy_list[[i]] == "Leviviridae"||
            taxonomy_list[[i]] == "Podoviridae"||
            taxonomy_list[[i]] == "Herelleviridae"||
            taxonomy_list[[i]] == "Narnaviridae"||
            taxonomy_list[[i]] == "Parvoviridae"||
            taxonomy_list[[i]] == "Reoviridae") {
          species_list[[i]] <-
            (taxonomy_parts[[1]][5])
        }
        if (taxonomy_list[[i]] == "Microviridae" ||
            taxonomy_list[[i]] == "Circoviridae") {
        species_list[[i]] <-
          (taxonomy_parts[[1]][4])
        }
        
        if (taxonomy_list[[i]] == "Coronaviridae") {
          species_list[[i]] <-
            "Betacoronavirus"
        }
      
        if(species_list[[i]] == "Chocolate lily virus A") {
          species_list[[i]] <-
            "Sadwavirus"
        }
        
      
        v_species <- (taxonomy_parts[[1]][length(taxonomy_parts[[1]])])
        
        if(v_species == "Apodemus agrarius picornavirus") {
          species_list[[i]] <-
            "Unclassified picornavirus 1"
          taxonomy_list[[i]] <- "Picornaviridae"
          type_list[[i]] <- "ss+RNA"
        }
        
        if(v_species == "Washington bat picornavirus" || v_species == "Washington bat picornavirus 2") {
          species_list[[i]] <-
            "Unclassified picornavirus 2"
          type_list[[i]] <- "ss+RNA"
        }

      }
      
    }
    
    print(length(as.character(type_list)))
    print(length(as.character(taxonomy_list)))
    print(length(as.character(species_list)))
    
    contingency_table.s$type <- as.character(type_list)
    contingency_table.s$virus_family <- as.character(taxonomy_list)
    contingency_table.s$species_list <- as.character(species_list)
    
    
    return(contingency_table.s)
    
  }

data <- annotate_and_filter_contingency_table(contingency_table, blast_results_table) %>%
  filter(virus_family != "NA")



# relative abundance by virus group/family
get_relAbundance <- function(list_of_viruses, data) {
  
  # data_subset to get read counts for all virus per species_time
  # d to get read counts for specific virus per species_time
  
  data_subset <-
    data %>%
    filter(virus_family %in% list_of_viruses)
  
  
  virus_family <- c("x")
  total_reads <- c(0)
  viral_reads <- c(0)
  interval <- c(0)
  species <- ("x")
  type <- c("x")
  

  
  for (f in 1:length(list_of_viruses)) {
    d <- data %>%
      filter(virus_family %in% list_of_viruses[f])
    
    if(dim(d)[1]>0) {
    for (i in 2:14) { # 13 columns corresponding to read data
      
      virus_type <- (as.character(unique(d[,15]))[1])

      total <- 0
      if (!is.na(as.numeric(unlist(data_subset[, i])))) {
        total <- sum(data_subset[, i])
      }
      
      read_count <- 0
      if (!is.na(as.numeric(unlist(d[, i])))) {
        read_count <- sum(d[, i])
        
      }
      
      colnames <- colnames(data_subset)
      parts <- str_split(colnames[i], "_")
      
      virus_family <- c(virus_family, list_of_viruses[f])
      total_reads <- c(total_reads, total)
      viral_reads <- c(viral_reads, read_count)
      species <- c(species, parts[[1]][4])
      interval <- c(interval, parts[[1]][5])
      type <- c(type, virus_type)
      
      
    }
    }
  }
  
  df <-
    data.frame(
      virus_family = virus_family,
      viral_read_count = viral_reads,
      total_read_count = total_reads,
      Species = species,
      Interval = interval,
      type = type
    )
  
  df <- na.omit(df[2:length(df[, 1]), ])
  
  return(df)
}

# relative abundance by virus species
get_relAbundance_s <- function(list_of_viruses, data) {

  data_subset <-
    data %>%
    filter(species_list %in% list_of_viruses)
  
  colnames <- colnames(data_subset)

  virus_family <- c("x")
  virus_species <- c("x")
  total_reads <- c(0)
  viral_reads <- c(0)
  interval <- c(0)
  species <- ("x")
  
  for (f in 1:length(list_of_viruses)) {
    d <- data%>%
      filter(species_list %in% list_of_viruses[f])
        
    if (dim(d)[1] > 0) {
      for (i in 2:14) {
        total <- 0
        if (!is.na(as.numeric(unlist(data_subset[, i])))) {
          total <- sum(data_subset[, i])
        }
        
        ssRNApos <- 0
        if (!is.na(as.numeric(unlist(d[, i])))) {
          ssRNApos <- sum(d[, i])
          
        }
        
        parts <- str_split(colnames[i], "_")
        
        virus_family <- c(virus_family, d$virus_family[1])
        virus_species <- c(virus_species, list_of_viruses[f])
        total_reads <- c(total_reads, total)
        viral_reads <- c(viral_reads, ssRNApos)
        species <- c(species, parts[[1]][4])
        interval <- c(interval, parts[[1]][5])
        
      }
    }
  }
  
  df <-
    data.frame(
      virus_family = virus_family,
      virus_species = virus_species,
      viral_read_count = viral_reads,
      total_read_count = total_reads,
      Species = species,
      Interval = interval
    )
  
  total_virus_reads.temp <- total_virus_reads %>%
    select(Interval, Sequences, Species)
  
  df <- na.omit(df[2:length(df[, 1]),]) %>%
    mutate(Interval=as.integer(Interval)) %>% 
    left_join(total_virus_reads.temp)
  
  
  return(df)
}

get_relAbundance_contigs <- function(list_of_viruses, data) {
  
  data_subset <-
    data %>%
    filter(virus_family %in% list_of_viruses)
  
  contig_ID <- c("x")
  virus_family <- c("x")
  virus_species <- c("x")
  total_reads <- c(0)
  viral_reads <- c(0)
  interval <- c(0)
  species <- ("x")
  
  list_of_contigs <- unique(data_subset$`contig ID`)
  
  for (f in 1:length(list_of_contigs)) {
    d <- data%>%
      filter(virus_family %in% list_of_viruses[f])
    
    if (dim(d)[1] > 0) {
      for (f in 1:length(list_of_contigs)) {
        d <- data%>%
          filter(`contig ID` %in% list_of_contigs[f])
        
        if (dim(d)[1] > 0) {
          for (i in 2:14) {
            total <- 0
            if (!is.na(as.numeric(unlist(data_subset[, i])))) {
              total <- sum(data_subset[, i])
            }
            
            
            read_count <- 0
            if (!is.na(as.numeric(unlist(d[, i])))) {
              read_count <- sum(d[, i])
              
            }
            
            colnames <- colnames(data_subset)
            parts <- str_split(colnames[i], "_")
            viral_reads <- c(viral_reads, read_count)
            contig_ID <- c(contig_ID, list_of_contigs[f])
            total_reads <- c(total_reads, total)
            species <- c(species, parts[[1]][4])
            interval <- c(interval, parts[[1]][5])
            
          }
        }
      }
    }
  }
  
  df <-
    data.frame(
      contig_ID = contig_ID,
      viral_read_count = viral_reads,
      total_read_count = total_reads,
      Species = species,
      Interval = interval
    )
  
  df <- na.omit(df[2:length(df[, 1]),])
  
  
  return(df)
}

get_relAbundance_c <- function(list_of_contigs, data) {
  
  data_subset <-
    data %>%
    filter(`contig ID` %in% list_of_contigs)
  
  contig_ID <- c("x")
  total_reads <- c(0)
  viral_reads <- c(0)
  interval <- c(0)
  species <- ("x")
  
  for (f in 1:length(list_of_contigs)) {
    d <- data%>%
      filter(`contig ID` %in% list_of_contigs[f])
    
    if (dim(d)[1] > 0) {
      for (i in 2:14) {
        total <- 0
        if (!is.na(as.numeric(unlist(data_subset[, i])))) {
          total <- sum(data_subset[, i])
        }
        
        
        read_count <- 0
        if (!is.na(as.numeric(unlist(d[, i])))) {
          read_count <- sum(d[, i])
          
        }
        
        colnames <- colnames(data_subset)
        parts <- str_split(colnames[i], "_")
        viral_reads <- c(viral_reads, read_count)
        contig_ID <- c(contig_ID, list_of_contigs[f])
        total_reads <- c(total_reads, total)
        species <- c(species, parts[[1]][4])
        interval <- c(interval, parts[[1]][5])
        
      }
    }
  }
  
  df <-
    data.frame(
      contig_ID = contig_ID,
      viral_read_count = viral_reads,
      total_read_count = total_reads,
      Species = species,
      Interval = interval
    )
  
  df <- na.omit(df[2:length(df[, 1]),])
  
  
  return(df)
}

# plot bar plot of rel abundance or read counts by virus family
plot_barplot <- function(df, relabundance) {
  y <- NULL
  if (relabundance == TRUE) {
    y <- df$viral_read_count / df$total_read_count
    ylab_string <- "Relative Abundance"
  }
  else{
    y <- df$viral_read_count
    ylab_string <- "Number of viral reads"
    
  }
  p <- 
    ggplot(df, aes(interval, y, fill = as.factor(virus_family))) +
    geom_bar(stat = "identity") + facet_wrap( ~ species, nrow = 1) +
    xlab("Time period") +
    ylab(ylab_string) +
    theme_bw() +
    theme(legend.position = "right") +
    scale_fill_brewer(palette = "Spectral", name = "")
  return(p)
}

# plot bar plot of rel abundance or read counts by virus subfamily
plot_barplot_s <- function(df, relabundance) {
  y <- NULL
  if (relabundance == TRUE) {
    y <- df$viral_read_count / df$total_read_count
    ylab_string <- "Relative Abundance"
  }
  else{
    y <- df$viral_read_count
    ylab_string <- "Number of viral reads"
    
  }
  p <- 
    ggplot(df, aes(interval, y, fill = as.factor(virus_species))) +
    geom_bar(stat = "identity") + facet_wrap( ~ species, nrow = 1) +
    xlab("Time period") +
    ylab(ylab_string) +
    theme_bw() +
    theme(legend.position = "right") +
    scale_fill_brewer(palette = "Spectral", name = "")
  return(p)
}


# all viruses rel abundance
all_virus <-
  virus_families$family
all_virus.df <- get_relAbundance(all_virus, data)




# enveloped viruses rel abundance
env_viruses <-
  c("Coronaviridae", "Paramyxoviridae", "Tobaniviridae", "Retroviridae")
enveloped_virus.df <- get_relAbundance(env_viruses, data)


# nonenveloped viruses rel abundance
nonenv_viruses <- 
  c("Astroviridae", "Hepeviridae", "Picornaviridae", "Reoviridae", "Picobirnaviridae") 
nonenveloped_virus.df <- get_relAbundance(nonenv_viruses, data)
nonenveloped_virus.df$virus_family <- factor(nonenveloped_virus.df$virus_family, levels=c(nonenv_viruses))


# veterbrate viruses rel abundance
verteb_viruses <- 
  c("Retroviridae", "Astroviridae", "Hepeviridae", "Picornaviridae", "Reoviridae", "Picobirnaviridae", "Coronaviridae", "Paramyxoviridae", "Tobaniviridae") 
verteb_viruses.df <- get_relAbundance(verteb_viruses, data)


# bacteriophage rel abundance
bacteriophage <-
  c("Herelleviridae","Podoviridae","Siphoviridae","Inoviridae","Myoviridae","Leviviridae", "Microviridae")
bacteriophage.df <- get_relAbundance(bacteriophage, data)
bacteriophage.df$virus_family <- factor(bacteriophage.df$virus_family, levels=c(bacteriophage))

bacteriophage2 <-  c("Herelleviridae","Myoviridae","Podoviridae","Siphoviridae","Inoviridae")
bacteriophage2.df <- get_relAbundance(bacteriophage2, data)
bacteriophage2.df$virus_family <- factor(bacteriophage2.df$virus_family, levels=c(bacteriophage))


# plant virus rel abundance
plant_viruses <- 
  c("Amalgaviridae","Betaflexiviridae","Botourmiaviridae","Bromoviridae","Closteroviridae","Potyviridae","Secoviridae","Solemoviridae", "Tymoviridae", "Virgaviridae", "Partitiviridae")
plant_virus.df <- get_relAbundance(plant_viruses, data)

# insect virus rel abundance

invert_viruses <- c("Lispiviridae", "Phasmaviridae", "Alphatetraviridae", "Dicistroviridae",
                    "Iflaviridae", "Mesoniviridae", "Polycipiviridae", "Nodaviridae", "Smacoviridae")

invert_virus.df <- get_relAbundance(invert_viruses, data)


# total virus reads over time
total_virus_reads <- read_csv("clean_viral_reads_stats.csv")
total_virus_reads <- total_virus_reads %>% 
  separate(Name, c(NA, NA, NA, "Species", "Interval"), sep="_") %>%
  rename(Sequences = `# Sequences`)
total_virus_reads$Interval <- as.integer(total_virus_reads$Interval)


all_virus_species <- data %>%
  filter(type %in% c("ss+RNA", "ss-RNA", "dsRNA", "reverse transcribing", "ssDNA", "dsDNA"))


env.v <- data %>% 
  filter(virus_family %in% env_viruses)

nonenv.v <- data %>% 
  filter(virus_family %in% nonenv_viruses)


vert <- data %>%
  filter(virus_family %in% verteb_viruses)

picorna <- data %>%
  filter(virus_family %in% c("Picornaviridae"))

phage <- data %>%
  filter(virus_family %in% bacteriophage)

all_virus_species.df  <- get_relAbundance_s(unique(all_virus_species$species_list), data)
verterb_virus_species.df  <- get_relAbundance_s(unique(vert$species_list), data)
picorna_species.df <- get_relAbundance_s(unique(picorna$species_list), data)
phage_species.df <- get_relAbundance_s(unique(phage$species_list), data)
enveloped_virus_species.df <- get_relAbundance_s(unique(env.v$species_list), data)
nonenveloped_virus_species.df  <- get_relAbundance_s(unique(nonenv.v$species_list), data)


midpoints = c("Jan 31", "Mar 31", "Jun 15","Sep 15", "Dec 12")
time <- rep(c((rep(midpoints,2)), midpoints[3:5]), length(unique(all_virus_species.df$virus_species)))
all_virus_species.df$date <- time
all_virus_species.df$date <- factor(all_virus_species.df$date, levels=(midpoints))

time <- rep(c((rep(midpoints,2)), midpoints[3:5]), length(unique(all_virus.df$virus_family)))
all_virus.df$date <- time
all_virus.df$date <- factor(all_virus.df$date, levels=(midpoints))


