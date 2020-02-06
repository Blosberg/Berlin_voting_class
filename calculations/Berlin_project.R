##-------------- BV_SVM.R - By Bren Osberg, MDC Berlin, started October 2016.
#-- last updated on never.

##--- Does support vector classification on the voting districts in Berlin
## ==================================================================================================

# install.packages("tmap")
# install.packages("leaflet") #--- had to install these to get it working.
# install.packages("tmaptools")
# install.packages("tsne")



library("tmap")
library("leaflet")
library("tmaptools")

# CLEAN UP EVERYTHING:

# == shape data can be obtained from: 
# == https://www.wahlen-berlin.de/wahlen/BE2016/Wahlkreiseinteil/wahlkreiseinteil.asp?sel1=1253&sel2=1045
# == actual voting data obtained from:
# == https://www.wahlen-berlin.de/Wahlen/BE2016/afspraes/download/download.html


rm(list=ls())
setwd("/Users/blosberg/Desktop/Science/postdoc_MDC/Berlin_vote/calculations/")
source('./funcs_Berlin.R') #---get what functions you need.
source('/Users/blosberg/Desktop/Science/Code_library/funcs_general.R')
# source('/Users/blosberg/Desktop/Science/2016_postdoc_MDC/RCode_deconv/func_def.R') #---get what functions you need.

#takesep="cl12"
takesep="OW"
num_Districts=12

# ========= read in shape file ==================
Berlinshapefile <- "./source_data/RBS_OD_UWB_AGH2016/UWB.shp"
Berlingeo_fine <- read_shape(file=Berlinshapefile)
Coarse_shapefile  = "./source_data/Wkr_25833/Wkr_25833.shp"
Berlingeo_coarse <- read_shape(file=Coarse_shapefile)

# ========= read in vote results ================

Data_init = read.csv("./source_data/DL_BE_AB2016_V2_formatted.txt",sep="\t", header=TRUE,stringsAsFactors=TRUE)
# Data_init = read.csv("./source_data/DL_BE_AB2016_V2_formatted.txt",sep="\t", header=TRUE,stringsAsFactors=TRUE)
covariate_data_init = read.csv("./source_data/DL_BE_AH2016_covariates_Strukturdaten.csv",sep="\t", header=TRUE,stringsAsFactors=FALSE)
# covariate_data_init = read.table("./source_data/DL_BE_AH2016_covariates_Strukturdaten.csv",sep="\t", header=TRUE,stringsAsFactors=FALSE,                                  colClasses=c('character'))
# covariate_data_init[1:5,1:5]

#--- now get rid of all parties below 1%
dim_Data_init = dim(Data_init)

should_keep =  matrix(TRUE, 1, dim_Data_init[2]); names(should_keep) = names(Data_init)
party_votes =  matrix(0,    1, dim_Data_init[2]); names(party_votes) = names(Data_init)

Total_votes = sum(Data_init[,17])

THRESH=0.05
for( i in c(18:dim_Data_init[2]) )
{
  party_votes[i] = sum(Data_init[,i])
  if( party_votes[i]/Total_votes < THRESH )
  {
  should_keep[,i]=FALSE  
  }
}

Data     = Data_init[,should_keep]
num_parties=sum(as.integer(should_keep[18:dim_Data_init[2]]))
dim_Data = dim(Data)
#============== THIS IS NOW THE DATA MATRIX WITH ONLY THE PARTIES OVER 1% INCLUDED.       =======
#============== "Urn" and "Brief" results are still mixed. Next step is to seperate those =======
Data_Urn   = Data[ Data[,6]=="Urnenwahlbezirk", ]
Data_Brief = Data[ Data[,6]=="Briefwahlbezirk", ]

#--- now reorganize the data so that the points are in the same order as in the geo-shape file.
# --- first grab the last two entries displaced to the end of the shape file:
  dummy=dim(Data_Urn); num_ridings=dummy[1]; num_fields =dummy[2]
  temp_data = rearrange_to_match_shape_file(Data_Urn, Berlingeo_fine, num_Districts)

  second_to_last_index = which(as.integer(temp_data[,5]==608) * as.integer(temp_data[,3]==07) ==1 ) 
  second_to_last_entry = temp_data[second_to_last_index,]
  last_index   = which(as.integer(temp_data[,5]==116) * as.integer(temp_data[,3]==05) ==1 ) 
  last_entry   = temp_data[last_index,]  

  temp1 = temp_data[-c(second_to_last_index, last_index),]
  temp2 = rbind.data.frame(temp1, rbind(second_to_last_entry,last_entry))

  Data_Urn        = temp2
  vote_share_Urn  = Data_Urn[,18:(18+num_parties-1)]/rowSums(Data_Urn[,18:(18+num_parties-1)])
  
  #--- now do the same reorganization for the covariate data (i.e. put in same order as in the geo-shape file)
    # --- Again, first grab the last two entries displaced to the end of the shape file:
  
  temp_data = rearrange_covariate_to_match_shape_file(covariate_data_init, Berlingeo_fine)
  TEST_covmap = identical( paste( district_string(temp_data[,2]), temp_data[,4], sep=""  ),  paste( Berlingeo_fine$BEZ, Berlingeo_fine$UWB, sep=""  ) )
  covariate_data_rearranged = temp_data
  
#=====================================================================================================
#==== PERFORMED THESE CHECKS A SUFFICIENT NUMBER OF TIMES TO BE CONFIDENT THAT THEY PASS. SKIPPING THEM NOW FOR EFFICIENCY.
# # 
#  TEST_BEZMATCH =  identical(as.integer(as.character(temp2[,3])) ,as.integer(as.character(Berlingeo_fine$BEZ)) )
#  TEST_UWBMATCH =  identical(as.integer(as.character(temp2[,5])) ,as.integer(as.character(Berlingeo_fine$UWB)) )
# # 
#   temp_rot = get_rot_mat(as.character(temp2[,1]) , as.character(Data_Urn[,1])  )
#  if(  (sum(rowSums(temp_rot))==num_ridings && sum(sum(colSums(temp_rot) ))==num_ridings ) && (min(rowSums(temp_rot)) ==1 && min(sum(colSums(temp_rot) ) ) ) )
#    {
#    ROTCHECK=TRUE
#    } else
#    {
#    ROTCHECK=FALSE
#    } #--- this establishes that we have a one-to-one mapping in our rearranged data matrix
#  i=15
#  max(abs(as.integer(temp2[,i]) - temp_rot %*%  Data_Urn[,i]))
#     #--- (which is now in the same order as the polygons from the geometrical shape file) 
#     #--- for all the entries from the original data matrix
# 
# if( (TEST_BEZMATCH && TEST_UWBMATCH) &&   ROTCHECK )
#   { # if all of the above checks, pass, then we start using this as our data matrix.
#   Data_Urn = temp2
#   }else
#   {
#   stop("There's been a problem ordering the data entries to match the shape file. Halting.")
#   }
#=====================================================================================================
 


# qtm(Berlingeo_fine) # --- will make the plot of the city show up as a figure.
# str(Berlingeo_fine) # --- will print out text data of the contents.

# ----------------------------------------------------

k=3
results_Urn   = get_rotation_and_clusters( Data_Urn,   k=3) 
results_Brief = get_rotation_and_clusters( Data_Brief, k=3) 

Cluster=factor(results_Urn$kmeans_out$cluster)
levels(Cluster)=c("West","Ost","Integriert")

#==== THIS IS THE ACTUAL MAP COMMAND. COMMENTING THIS OUT FOR EFFICIENCY (SO WE DONT GENERATE THE FULL MAP EVERY TIME)
# Berlinmap_Urn <- append_data(Berlingeo_fine, data.frame(Cluster), key.shp = NULL, fixed.order=TRUE, key.data="Cluster", ignore.duplicates=TRUE)
# qtm(Berlinmap_Urn, "Cluster", fill.palette=c("deepskyblue","firebrick2","darkgreen") )

Berlinmap_Urn <- append_data(Berlingeo_fine, data.frame(covariate_data_rearranged), key.shp = NULL, fixed.order=TRUE, key.data="Cluster", ignore.duplicates=TRUE)

 
qtm(Berlinmap_Urn, "eutsche.45.60.Prozen")

  #================== FINISHED MAKING MAP---

plot(results_Urn$Westpoints_rotated[,1]*100, results_Urn$Westpoints_rotated[,2]*100, col="blue", pch=16,xlim=c(-27,25),ylim=c(-22,32), lwd=10, main="Vote space; East and West", xlab = "[s]",ylab="PC1")
points(results_Urn$Ostpoints_rotated[,1]*100, results_Urn$Ostpoints_rotated[,2]*100, col="red", pch=18, lwd=10)
legend("topleft",inset=.05, c( "West", "Ost"), col=c("blue","red"),pch=c(16,18,1,2,3),    bty = "n") 

plot(results_Urn$Datapoints_rotated[results_Urn$kmeans_out$cluster==1,1]*100,  results_Urn$Datapoints_rotated[results_Urn$kmeans_out$cluster==1,2]*100,xlim=c(-27,25),ylim=c(-22,32), lwd=3, main="Clusters", col="Blue",pch=1, xlab = "[s]",ylab="PC1")
points(results_Urn$Datapoints_rotated[results_Urn$kmeans_out$cluster==2,1]*100,  results_Urn$Datapoints_rotated[results_Urn$kmeans_out$cluster==2,2]*100,lwd=3, col="red",pch=2)
points(results_Urn$Datapoints_rotated[results_Urn$kmeans_out$cluster==3,1]*100,  results_Urn$Datapoints_rotated[results_Urn$kmeans_out$cluster==3,2]*100,lwd=3, col="Darkgreen",pch=4)

# legend("topleft",inset=.05, c( "West", "Ost", "cl. 1", "cl. 2", "cl. 3"), col=c("blue","red", "cyan", "pink","yellow"),pch=c(16,18,1,2,3)) 

##===========  NOW FOR THE BRIEF RESULTS =======

plot(results_Brief$Westpoints_rotated[,1]*100, results_Brief$Westpoints_rotated[,2]*100, col="blue", pch=16,xlim=c(-30,30),ylim=c(-40,35), main = "vote clusters -Brief \n x-axis along cluster means", xlab = "[s] (c2-c1)",ylab="PC1")
points(results_Brief$Ostpoints_rotated[,1]*100, results_Brief$Ostpoints_rotated[,2]*100, col="red", pch=18)
  c1_rotated = results_Brief$Datapoints_rotated[results_Brief$kmeans_out$cluster==1,]*100
  c2_rotated = results_Brief$Datapoints_rotated[results_Brief$kmeans_out$cluster==2,]*100
  c3_rotated = results_Brief$Datapoints_rotated[results_Brief$kmeans_out$cluster==3,]*100
    points(c1_rotated[,1],  c1_rotated[,2], col="cyan",pch=1)
    points(c2_rotated[,1],  c2_rotated[,2], col="pink",pch=2)
    points(c3_rotated[,1],  c3_rotated[,2], col="yellow",pch=3)
legend("bottomleft",inset=.05, c( "West", "Ost", "cl. 1", "cl. 2", "cl. 3"), col=c("blue","red", "cyan", "pink","yellow"),pch=c(16,18,1,2,3)) 
# points(cluster3points_Urn[,1],cluster4points_Urn[,2], col="brown",pch=3)
# rotated_cluster_means = t( (results_Brief$rot_mat) %*% t(results_Brief$kmeans_out$centers)  ) *100
points(mean(c1_rotated[,1]) , mean(c1_rotated[,2]),col="black",pch=20,lwd="3")
points(mean(c2_rotated[,1]) , mean(c2_rotated[,2]),col="black",pch=20,lwd="3")
points(mean(c3_rotated[,1]) , mean(c3_rotated[,2]),col="black",pch=20,lwd="3")

# 
# ALTERNATIVE PLOTTING ANGLES:
# plot(Ostpoints_rotated[,1], Ostpoints_rotated[,3], col="blue", pch=16,xlim=c(-0.3,0.3),ylim=c(-0.2,0.25))
# points(Westpoints_rotated[,1],Westpoints_rotated[,3],col="red", pch=18)
# 
# 
legend("topleft",inset=.05, c("Ost", "West"), col=c("blue","red"),pch=c(16,18)) 
# 
# plot(Ostpoints[,4], Ostpoints[,1], col="blue", pch=16, xlim=c(0,0.45),ylim=c(0,0.45),xlab = "Linke", ylab = "SPD")
# points(Westpoints[,4],Westpoints[,1],col="red", pch=18)
# legend("topright",inset=.05, c("Ost", "West"), col=c("blue","red"),pch=c(16,18)) 

barplot(100*results_Urn$kmeans_out$centers, main="Mean vote-share by cluster", cex.names=0.6, xlab="Party", srt=45, axis.lty=1, ylab="Vote share", col=c("darkblue","red","darkgreen"), beside=TRUE,    bty = "n")
legend(c(17),c(27), legend = c("West","Ost","Integriert"), fill = c("darkblue", "red","darkgreen"),    bty = "n")

# axis(1, at=1:8, labels=names()
#== ,xlim=c(-0.3,0.3),ylim=c(-0.15,0.03)
axis(side = 1, at = c(1:8))
axis(side=1, at=seq(1,8,2), labels=labs[seq(0,8,2)], cex.axis=0.35)
axis(side=1, at=seq(2,8,2), labels=labs[seq(1,7,2)], cex.axis=0.35)



tsne_districts = tsne::tsne(vote_share_Urn , initial_dims = num_parties, perplexity=50)

#===================== OW + CLUSTER PLOTS =============================
tsne_points_W=tsne_districts[Data_Urn[,9]=="W",]
tsne_points_O=tsne_districts[Data_Urn[,9]=="O",]
tsne_points_C1=tsne_districts[results_Urn$kmeans_out$cluster==1,]
tsne_points_C2=tsne_districts[results_Urn$kmeans_out$cluster==2,]
tsne_points_C3=tsne_districts[results_Urn$kmeans_out$cluster==3,]

plot(tsne_districts,col="black",lwd=2, ylab = "separation", xlab = "separation")
points(tsne_districts[Data_Urn[,9]=="W",],col="blue",lwd="3")
points(tsne_districts[Data_Urn[,9]=="O",],col="red",lwd="3")
points(tsne_points_C3,col="green",lwd="3")
legend("topleft",inset=.05, c( "West", "Ost","cluster3"), col=c("blue","red","green"),pch=c(1,1,1), bty = "n",lwd=3) 

#==================== WINNER PLOTS ==================================

# tsne_winners_plot
plot(tsne_districts,col="black",lwd=2, main="District winners", ylab = "separation", xlab = "separation")
#  
points(tsne_districts[which(max.col(Data_Urn[,18:(18+num_parties-1)])==1),],col="red", lwd=3)
points(tsne_districts[which(max.col(Data_Urn[,18:(18+num_parties-1)])==2),],col="blue", lwd=3)
points(tsne_districts[which(max.col(Data_Urn[,18:(18+num_parties-1)])==3),],col="green", lwd=3)
points(tsne_districts[which(max.col(Data_Urn[,18:(18+num_parties-1)])==4),],col="orange", lwd=3)
points(tsne_districts[which(max.col(Data_Urn[,18:(18+num_parties-1)])==6),],col="yellow", lwd=6)
points(tsne_districts[which(max.col(Data_Urn[,18:(18+num_parties-1)])==9),],col="purple", lwd=3)
legend("topleft",inset=.05, c( "SPD", "CDU","Green","Linke","FDP","AfD"), col=c("red","blue","green","orange","yellow","purple"),pch=c(1,1,1,1,1,1), bty = "n", lwd=3) 


which(max.col(Data_Urn[,18:(18+num_parties-1)])==10)

#  points(tsne_points_all[which(max.col(Data_Urn[,18:25])==2),],col="blue", lwd=2)
#  points(tsne_points_all[which(max.col(Data_Urn[,18:25])==3),],col="green", lwd=2)
#  Data_Urn[1,18:25]
#  points(tsne_points_all[which(max.col(Data_Urn[,18:25])==4),],col="orange", lwd=2)
#  points(tsne_points_all[which(max.col(Data_Urn[,18:25])==8),],col="purple", lwd=2)
#  legend("topleft",inset=.05, c( "SPD", "CDU","Green","Linke","AfD"), col=c("red","blue","green","orange","purple"),pch=c(1,1,1,1,1), bty = "n") 
