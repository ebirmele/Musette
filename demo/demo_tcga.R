library(musette)
library(ComplexHeatmap)

data("tcga_bladder",package="musette")
bound=20  #number of solutions to generate
stepmode="strict"





blind_domination_step=c("ampli","dele") # alteration types  to be considered for the "blind" domination step
blind_percent=90  #  required percentage for the "blind" domination step
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

#influence graph of the alterations. 
g <- influenceGraph(ll,20,TRUE)
visIgraph(g)


#Oncoplot of the best solution
solution.oncomatrix(ll,1,reds)



# Data export to csv files
csv_sol=musette2csv(sol)
csv_alterations=musette2csv(alterations)





demerge=function(l,i){ # de-fusion d'une solution, pas encore mis en place, a discuter
  if(i>length(l)) return(list(l))
  if(alterations[[l[[i]],"type"]]!="pseudo") return (demerge(l,i+1))
  aliases=alterations[[l[[i]],"merged_from"]]
  do.call(c,lapply(aliases,function(a){ 
    l[[i]]=a
    demerge(l,i+1)
  } 
  )
  )
}


