##############################################################################
# The heritability of promiscuity in the Lundy House Sparrows
##############################################################################

# Isabel Winney
# 20160115

rm(list=ls())

##############################################################################
# loading necessary packages
##############################################################################

library(MCMCglmm)
library(pedantics)
library(RODBC)
library(knitr)

##############################################################################
# what system, R, and package configurations am I using?
##############################################################################

#R version 3.0.2 (2013-09-25) -- "Frisbee Sailing"
#Copyright (C) 2013 The R Foundation for Statistical Computing
print(sessionInfo())

# output pasted on 20160118
#R version 3.0.2 (2013-09-25)
#Platform: i386-w64-mingw32/i386 (32-bit)
#
#locale:
#  [1] LC_COLLATE=English_United Kingdom.1252  LC_CTYPE=English_United Kingdom.1252   
#[3] LC_MONETARY=English_United Kingdom.1252 LC_NUMERIC=C                           
#[5] LC_TIME=English_United Kingdom.1252    
#
#attached base packages:
#  [1] grid      stats     graphics  grDevices utils     datasets  methods   base     
#
#other attached packages:
#  [1] pedantics_1.5    MasterBayes_2.52 kinship2_1.6.0   quadprog_1.5-5   genetics_1.3.8.1
#[6] mvtnorm_1.0-2    MASS_7.3-29      gtools_3.4.2     gdata_2.13.3     combinat_0.0-8  
#[11] MCMCglmm_2.22    ape_3.2          coda_0.17-1      Matrix_1.2-0    
#
#loaded via a namespace (and not attached):
#  [1] corpcor_1.6.7   cubature_1.1-2  lattice_0.20-31 nlme_3.1-111    tensorA_0.36    tools_3.0.2    



##############################################################################
# loading the pedigree
##############################################################################

# Load the most recent working house sparrow pedigree in to the file.
# This can be achieved with the script importSparrowPedigree

source("importSparrowPedigree.R")

# metadata is in importSparrowPedigree

head(sparrowpedigree)

# cut to the three first columns:
sparrowped <- sparrowpedigree[,1:3]
head(sparrowped)
tail(sparrowped)

# fix the pedigree (missing sires and dams are inserted at the top of the file)
{fixedsparrowped <- fixPedigree(sparrowped)

head(fixedsparrowped)
tail(fixedsparrowped)

summary(fixedsparrowped)
str(fixedsparrowped)}


# in previous work with the pedigree, I have imputed dams.
# This involves taking every individual with no assigned dam and
# assigning an unique dam to each of those individuals.
# at the moment, I will not be doing this in this script.

# The same can be done with the ASReml-R ready version of the pedigree:
{
  asrped <- sparrowpedigree[,5:7]
  head(asrped)
  tail(asrped)

  # fix the pedigree (missing sires and dams are inserted at the top of the file)
  fixedasrped <- fixPedigree(asrped)

  head(fixedasrped)
  tail(fixedasrped)

  summary(fixedasrped)
  str(fixedasrped)
}

##############################################################################
# loading the Lundy House Sparrow Database
##############################################################################

sparrowDB <- odbcConnectAccess('C:/Users/Issie/SkyDrive/PhD/SparrowDatabases/Database0.74_12Oct2015GenotypesUpdatedBMSBNS/SparrowDatabase0.74_Issie.mdb')

# view tables within the database
sqlTables(sparrowDB)


##############################################################################
# offspring data
##############################################################################

# Assumptions:
# 1   males are all dead --> the calculated phenotypes
#     are complete. This is not true but can be reconsidered later
#     later.
# 2   the 2013 pedigree is ok.
# 3   a brood has just one genetic mother. By looking
#     up the genetic mother of a brood, the mother that is
#     returned first is the one used in the analysis below. 
# 4   I can only make phenotypes using individuals that have been born
#     in a brood, otherwise I do not know their social and genetic parents
#     i.e. for a mother, I can only calculate her EPO and WPO from EPO/WPO
#     born in broods
# 5   All cohorts in the database are accurate (not true)


# I have made a list using the database from March 2015
# database called DataBaseV0.74-updated20150320-AlfredoCheckedMarch2015
# of all the birds associated with a brood, their cohorts,
# their broods, their parents.
# This was called birdbroodcohortparent-to2014-20150910.txt
# and now I am directly calling the database to extract the
# same data instead.

# This means that my measure of EPO and WPO totals is only
# based on these offspring and all the offspring that are
# caught in the winter etc I ignore.

# This is the Git version of an existing script in which I explored
# how to make the data and use the sparrow data set for an analysis
# of this.

offspring <- sqlQuery (sparrowDB,
                      "SELECT tblBirdID.BirdID, 
                              tblBirdID.Cohort,
                              tblBirdID.BroodRef, 
                              tblBroods.BroodName, 
                              tblBroods.NestboxRef, 
                              tblBroods.SocialDadID, 
                              tblBroods.SocialDadCertain, 
                              tblBroods.SocialMumID, 
                              tblBroods.SocialMumCertain
                      FROM tblBroods 
                      INNER JOIN tblBirdID 
                      ON tblBroods.BroodRef = tblBirdID.BroodRef;",
                      na.strings="NA")

head(offspring)
tail(offspring)
str(offspring)
summary(offspring)

# In how many cases do I know no parents versus both?
table(offspring$SocialDadCertain,
      offspring$SocialMumCertain)


# In 6433 cases out of 7628 I know both. That presumably includes
# cases where the parent was not known, rather than not certain

which(is.na(offspring$SocialDadID))
which(offspring$SocialDadCertain==0)

# Quick visual inspection shows a lot of shared rows.



#-------------------------------------------
# Add genetic parentage to offspring data
#-------------------------------------------


# This ensures that offspring with genetic AND social fathers can
# be assigned a status as extra-pair offspring or within-pair (EPO or WPO):

offspring$GeneticDadID <- sparrowped$sire[match(offspring$BirdID, 
                                                sparrowped$birdid)]

offspring$GeneticMumID <- sparrowped$dam[match(offspring$BirdID,
                                               sparrowped$birdid)]

head(offspring)
tail(offspring)

# Any missing years of genetic parents indicate the end of the pedigree.
# i.e. at the time of writing (20160118) the pedigree goes to 2013, and
# all offspring from 2014 and 2015 do not have genetic parents.

summary(offspring)


# where parentage is missing for the genetic male or the social male,
# we do not know whether an offspring is EPO or WPO. Therefore we 
# must remove these cases from the data set:

offspring2 <- offspring[-which(is.na(offspring$SocialDadID)),]
offspring3 <- offspring2[-which(is.na(offspring2$GeneticDadID)),]

# how much data was lost at each step:
length(offspring[,1])
length(offspring2[,1])
length(offspring3[,1])

# how many social fathers and genetic fathers are there?
length(unique(offspring3$SocialDadID))
length(unique(offspring3$GeneticDadID))


#-------------------------------------------
# Designate offspring as EPO or WPO
#-------------------------------------------

# WPO have the same social and genetic sire,
# EPO have different genetic and social sires.

offspring3$WPO <- ifelse(offspring3$SocialDadID==offspring3$GeneticDadID,
                         1, 0)
head(offspring3)
tail(offspring3)

offspring3$EPO <- ifelse(offspring3$SocialDadID!=offspring3$GeneticDadID,
                         1, 0)

table(offspring3$EPO, offspring3$WPO)

# 856 EPO. 4011 WPO. None identified as both. Good!

##############################################################################
# Additional consideration 1: exclude birds that are still alive
##############################################################################

# The consideration is that birds that are still alive are still 
# having EPO and WPO, so our phenotype is not completely accurate.
# What happens to the analysis without these birds?

##############################################################################
# Additional consideration 2: exclude birds that are not certain social parents
##############################################################################

# When we have imperfect sightings information in the field and/or infer
# social parentage later using the sightings and the genetic information,
# we tick 0 for 'SocialDadCertain' or 'SocialMumCertain' so that the parent
# can be marked as less reliable.

##############################################################################
# Close the database
##############################################################################

close(sparrowDB)
