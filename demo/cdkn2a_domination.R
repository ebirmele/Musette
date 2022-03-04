library(musette)
library(ComplexHeatmap)
library(tidyverse)
library(ggplot2)



data("tcga_bladder",package="musette")
bound=20  #number of solutions to generate
stepmode="strict"





blind_domination_step=c("ampli","dele") # alteration types  to be considered for the "blind" domination step
blind_percent=90 #  required percentage for the "blind" domination step
blind_distance=5000000 # maximal distance at which an alteration can blind-dominate another one

color_domination_step=c("ampli","dele") # same parameters for the "colored" domination.
red_percent=80  # two percentages (red and blue) have to be defined
blue_percent=80
color_distance=5000000

# choice of the reds
#reds= names(groups) %in% groupe_cis_1
#reds= groups=="basal"
reds= (groups == 'basal')
#reds= (groups == 'luminale I')
#reds= (groups == 'luminale II')
#reds= (groups %in% c('luminale I', 'luminale II'))
#reds= matrices$dele['CDKN2A',] #(CDKN2A perdu)
names(reds)=names(groups)
blues=!reds
names(blues)=names(groups)

ll=do.musette(matrices=matrices, reds=reds, blind_domination_step = blind_domination_step,
              color_domination_step = color_domination_step, blind_distance=blind_distance, 
              color_distance=color_distance, blind_percent=blind_percent, red_percent=red_percent, 
              blue_percent=blue_percent, chromosome = chromosome, longname=longname,
              pathways=pathways,position=position, bound=bound)

sol=ll$solutions            # dataframe of the solutions of the musette algorithm
alterations=ll$alterations  # dataframe of the alterations considered to generate the solutions



## liste des altérations dominées (au sens large) par CDKN2A

dominated <- alterations[which(alterations$name=="dele_CDKN2A"),]$followers$dele_CDKN2A

## vérification qu'ils sont tous sur le même chromosome: réponse oui
#sum(alterations[which(alterations$name%in%dominated),]$chromosome==9)
#length(dominated)

## où sont-ils: entre 17 319 452   et  35 923 318
#dominated_positions <- unlist(alterations[which(alterations$name%in%dominated),]$position)
#min(dominated_positions)
#max(dominated_positions)


## création d'un tibble qui contient les infos sur les deletions du chromosome 9

# on ne garde que les deletions 
single_alterations <- alterations[grepl("dele",alterations$name),] 

#puis que celles sur le chromosome 9 entre 10 000 000 et 50 000 000,
# et on ajoute une variable disant si la deletion est dominée par cdkn2a
cdkdata <- tibble(name = unlist(single_alterations$name),chromosome=unlist(single_alterations$chromosome),position=unlist(single_alterations$position)) %>%
            filter(chromosome==9) %>% 
            #filter(position>10000000) %>%
            filter(position<50000000) %>%
            mutate(dominated= name %in% dominated)


## ajout du nombre d'évhantillons pour lesquels chaque deletion est présente,
## en séparant celles où cdkn2a est aussi délété

ref_deletions <- matrices$dele[which(rownames(matrices$dele)=="CDKN2A"),]
common_deletions <- c()
specific_deletions <- c()

for (i in 1:length(cdkdata$name)){
  
  altname <- strsplit(cdkdata$name[i],"_")[[1]][2] #suppression du dele_ en début de nom
  alt_deletions <- matrices$dele[which(rownames(matrices$dele)==altname),]
  common_deletions <- c(common_deletions,sum(ref_deletions & alt_deletions))
  specific_deletions <- c(specific_deletions,sum(!ref_deletions & alt_deletions))
}

#position de cdkn2a et des zones pour dles dominés
cdkn2aposition <- cdkdata$position[cdkdata$name=="dele_CDKN2A"]

rect1min <- mean(c(cdkdata$position[61],cdkdata$position[62]))
rect1max <- mean(c(cdkdata$position[62],cdkdata$position[63]))
rect2min <- mean(c(cdkdata$position[73],cdkdata$position[74]))
rect2max <- mean(c(cdkdata$position[117],cdkdata$position[118]))
rect3min <- mean(c(cdkdata$position[118],cdkdata$position[119]))
rect3max <- mean(c(cdkdata$position[120],cdkdata$position[121]))
rect4min <- mean(c(cdkdata$position[125],cdkdata$position[126]))
rect4max <- mean(c(cdkdata$position[145],cdkdata$position[146]))
rect5min <- mean(c(cdkdata$position[150],cdkdata$position[151]))
rect5max <- mean(c(cdkdata$position[208],cdkdata$position[209]))




#dedoublement pour faire difference entre deletions communes avec cdkn2a
cdkcommon <- cdkdata %>% add_column(deletions=common_deletions,type="common")
cdkspecific <- cdkdata %>% add_column(deletions=specific_deletions,type="specific")

cdkdata <- bind_rows(cdkcommon,cdkspecific)



pdf(file = "domination.pdf")
ggplot(data=cdkdata, aes(x=position, y=deletions, group=type)) +
  #geom_line(aes(linetype=type))+
  geom_rect(aes(xmin = rect1min, ymin = -Inf, xmax = rect1max, ymax = Inf),
          fill = "lightblue") +
  geom_rect(aes(xmin = rect2min, ymin = -Inf, xmax = rect2max, ymax = Inf),
            fill = "lightblue") +
  geom_rect(aes(xmin = rect3min, ymin = -Inf, xmax = rect3max, ymax = Inf),
            fill = "lightblue") +
  geom_rect(aes(xmin = rect4min, ymin = -Inf, xmax = rect4max, ymax = Inf),
            fill = "lightblue") +
  geom_rect(aes(xmin = rect5min, ymin = -Inf, xmax = rect5max, ymax = Inf),
            fill = "lightblue") +
  geom_point(aes(shape=type, color = type))+
  geom_vline(xintercept = cdkn2aposition, linetype="dotted", 
             color = "red", size=0.5)  

dev.off()






