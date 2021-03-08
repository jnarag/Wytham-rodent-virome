library(ggplot2)

library(wesanderson)
??wesanderson
pal = c('#31a354',"#EBCC2A","#e17700")

pico = read.csv('data/nerc_picorna.predictions_with_seq_names.csv', header = TRUE)
head(pico)
pico$Seq_IDs = as.factor(pico$Seq_IDs)

pico$Seq_IDs = factor(pico$Seq_IDs, levels = c('Hunnivirus (b)','Unclassified 1','Mosavirus',
                                         'Hunnivirus (a)','Unclassified 2 (a)',
                                         'Unclassified 2 (b)','Cardiovirus'))


pico$priority_category = as.factor(pico$priority_category)
pico$priority_category =factor(pico$priority_category, levels = c('Low', 'Medium','High'))
ggplot(pico, aes(Seq_IDs, calibrated_score_mean, col = priority_category))+
  geom_hline(yintercept=0.3, linetype="dashed", color = "black")+
  geom_point()+
  coord_flip()+
  theme_bw()+
  labs(x = '', y = 'probability infects humans', color = 'priority')+
  scale_y_continuous(expand = c(0.0,0),limits = c(0,1.03), breaks = seq(0,1,by= 0.25))+
  geom_errorbar(aes(ymin=calibrated_score_lower, 
                    ymax=calibrated_score_upper), width=0,
                position=position_dodge(0.05), alpha = 0.5)+
  theme(legend.position = c(0.8, 0.3))+
  scale_colour_manual(values = pal)
  
  
