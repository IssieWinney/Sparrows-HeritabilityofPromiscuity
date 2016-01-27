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
{
  fixedsparrowped <- fixPedigree(sparrowped)

  head(fixedsparrowped)
  tail(fixedsparrowped)

  summary(fixedsparrowped)
  str(fixedsparrowped)
  
  # re-name the pedigree to be animal, dam, sire
  names(fixedsparrowped) <- c("animal", "dam", "sire")
  head(fixedsparrowped)
}


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
  
  # re-name the pedigree to be animal, dam, sire
  names(fixedasrped) <- c("animal", "dam", "sire")
  head(fixedasrped)
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

#-------------------------------------------
# Load offspring data from database
#-------------------------------------------
{
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
}

#-------------------------------------------
# Add genetic parentage to offspring data
#-------------------------------------------
{

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
}

#-------------------------------------------
# Designate offspring as EPO or WPO
#-------------------------------------------
{
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
}

##############################################################################
# cohort data
##############################################################################

# extract the list of cohorts and Bird IDs

{
  sqlFetch(sparrowDB, "tblBirdID", max=10)

  birdcohort <- sqlQuery (sparrowDB,
                       "SELECT tblBirdID.BirdID, 
                       tblBirdID.Cohort
                       FROM tblBirdID;",
                      na.strings="NA")
  
  head(birdcohort)
  summary(birdcohort)
  # a missing cohort???
  birdcohort[which(is.na(birdcohort$Cohort)),]
  # this is a blank record from my last year :s
  # rest of the data seems fine
}


##############################################################################
# male phenotypes 
##############################################################################

# Now aggregate the data in to a phenotype for each male.
# Create a list of all GENETIC FATHERS in the data set and 
# the number of EPO and WPO that each father sires.
# This means that males that sire no offspring are NOT in 
# this data set

#-------------------------------------------
# Male dataset of genetic offspring per lifetime
#-------------------------------------------

{
  # make a list of unique males and count the number of WPO per male
  malephenotypes <- aggregate(offspring3$WPO, 
                              list(offspring3$GeneticDadID),
                              FUN=sum)
  head(malephenotypes)
  
  # now the number of EPO per male
  maleEPO <- aggregate(offspring3$EPO, 
                       list(offspring3$GeneticDadID),
                       FUN=sum)
  head(maleEPO)
  
  # add sensible names to the main data frame
  names(malephenotypes) <- c("animal", "WPO")
  head(malephenotypes)
  
  str(malephenotypes)
  str(maleEPO)
  
  # and put the two data frames together
  malephenotypes <- cbind(malephenotypes, maleEPO)
  
  head(malephenotypes)
  
  # check that all data points have the same male from the WPO
  # as EPO data set
  malephenotypes$animal==malephenotypes$Group.1
    
  # and that no males are different
  malephenotypes$animal!=malephenotypes$Group.1
  
    
  # good!
  
  # rename the data set with names that are useful later:
  names(malephenotypes) <- c("animal", "WPO", "maleID", "EPO")
  head(malephenotypes)
    
  
  # are any males in my file not in the pedigree?
  check1 <- fixedsparrowped$animal[match(malephenotypes$animal,
                                  fixedsparrowped$animal)]
  
  summary(check1)
  check1
  # wooo! So, no NA's means that all males are in the pedigree.
}

{
  # add maternal and paternal ID of each male:
  
  # convert factors to numeric:
  tail(fixedsparrowped)
  fixedsparrowped$animal <- as.numeric(levels(fixedsparrowped$animal))[fixedsparrowped$animal]
  fixedsparrowped$dam <- as.numeric(levels(fixedsparrowped$dam))[fixedsparrowped$dam]
  fixedsparrowped$sire <- as.numeric(levels(fixedsparrowped$sire))[fixedsparrowped$sire]
  tail(fixedsparrowped)
  
  # maternal ID
  head(malephenotypes)
  malephenotypes$maternalID <- fixedsparrowped$dam[match(malephenotypes$animal,
                                                         fixedsparrowped$animal)]
  
  # and paternal ID
  malephenotypes$paternalID <- fixedsparrowped$sire[match(malephenotypes$animal,
                                                         fixedsparrowped$animal)]
  
  head(malephenotypes)
  
  # do a check
  fixedsparrowped[which(fixedsparrowped$animal==500),]
  malephenotypes[which(malephenotypes$animal==500),]
  
  # that is correct!
  
  fixedsparrowped[which(fixedsparrowped$animal==4722),]
  malephenotypes[which(malephenotypes$animal==4722),]

  # also correct! Looking good!
  
  str(malephenotypes)
}

{
  # Add the cohort of each male:
  malephenotypes$Cohort <- birdcohort$Cohort[match(malephenotypes$animal,
                                                   birdcohort$BirdID)]
  
  head(malephenotypes)
  summary(malephenotypes)
  
  # check:
  sparrowpedigree[which(sparrowpedigree$birdid==4722),]
  malephenotypes[which(malephenotypes$animal==4722),]
  # good good.
}

{
  # set appropriate factors for analysis
  str(malephenotypes)
  malephenotypes$factoranimal <- as.factor(malephenotypes$animal)
  malephenotypes$factormaternalID <- as.factor(malephenotypes$maternalID)
  malephenotypes$factorpaternalID <- as.factor(malephenotypes$paternalID)
  malephenotypes$factorcohort <- as.factor(malephenotypes$Cohort)
  str(malephenotypes)
}

{
  # male lifespan. Lifespan can be defined in many ways:
  # last time as a social parent, last time seen in the genetic pedigree,
  # last time captured, last time sighted. The pedigree has a 'last seen 
  # alive' query that takes all of this information and integrates it to
  # give the most recent and most likely date seen alive.
  
  # the current database does not contain the most recent pedigree.
  # for now, for the purposes of this analysis, the last date seen
  # alive will be defined as the last time seen in the GENETIC PEDIGREE.
  
  # so, find the maximum cohort for a given sire in the pedigree:
  malelifespan <- aggregate(sparrowpedigree$cohort, 
                           list(sparrowpedigree$sire),
                           FUN=max)
  
  head(malelifespan)
  str(malelifespan)
  
  # match to the data frame:
  malephenotypes$lastyearPED <- malelifespan$x[match(malephenotypes$animal,
                                                     malelifespan$Group.1)]
  
  head(malephenotypes)
  
  # make lifespan variable by subtracting cohort from year last seen:
  malephenotypes$lifespanPED <- malephenotypes$lastyearPED - malephenotypes$Cohort
  summary(malephenotypes$lifespanPED)
  table(malephenotypes$lifespanPED)
  # one male reproduced in year zero...
  
  # make spare age variable:
  malephenotypes$lifespanPED6 <- malephenotypes$lifespanPED
  
  # assign >6 years old to 6 years old for sample size by finding the
  # >6 data points and replacing them with a 6:
  malephenotypes$lifespanPED6[which(malephenotypes$lifespanPED6>6)] <- 6
  
  # and assign the zero year old to be one:
  malephenotypes$lifespanPED6[which(malephenotypes$lifespanPED6==0)] <- 1
  
  # check
  summary(malephenotypes$lifespanPED6)
  table(malephenotypes$lifespanPED6)
  table(malephenotypes$lifespanPED)
  # numbers correct.
}

#-------------------------------------------
# Male dataset of genetic offspring per year
#-------------------------------------------

{
  # per year per GENETIC male identifier:
  offspring3$maleyear <- paste(offspring3$GeneticDadID, offspring3$Cohort, sep=".")
  
  head(offspring3)
  table(offspring3$maleyear)
  table(table(offspring3$maleyear))
  # up to 25 offspring!
  
  
  # aggregate number of WPO per male
  maleyearWPO <- aggregate(offspring3$WPO, 
                           list(offspring3$maleyear),
                           FUN=sum)
  head(maleyearWPO)
  
  # aggregate number of EPO per male
  maleyearEPO <- aggregate(offspring3$EPO, 
                           list(offspring3$maleyear),
                           FUN=sum)
  head(maleyearEPO)
  
  
  names(maleyearWPO) <- c("animal", "WPO")
  head(maleyearWPO)
  
  str(maleyearWPO)
  str(maleyearEPO)
  
  
  maleyear <- cbind(maleyearWPO, maleyearEPO)
  
  head(maleyear)
  
  # check that all data rows are matched to the right male:
  which(maleyear$animal==maleyear$Group.1)
  # no missing numbers --> good
  which(maleyear$animal!=maleyear$Group.1)
  # no missmatch --> correct
  
  # rename data set
  names(maleyear) <- c("maleyear1", "WPO", "maleyear", "EPO")
  head(maleyear)
  
  
  # add the male's ID:
  maleyear$animal <- offspring3$GeneticDadID[match(maleyear$maleyear, offspring3$maleyear)]
  head(maleyear)
  
  # add male ID for repeated measures:
  maleyear$maleid <- maleyear$animal
  
}
  
{
  # add maternal ID for maternal effects:
  maleyear$maternalID <- fixedsparrowped$dam[match(maleyear$animal,
                                              fixedsparrowped$animal)]
  head(maleyear)
  
  # and paternal ID
  maleyear$paternalID <- fixedsparrowped$sire[match(maleyear$animal,
                                                        fixedsparrowped$animal)]
  
  summary(maleyear)
}

{
  # Add the cohort of each male:
  maleyear$Cohort <- birdcohort$Cohort[match(maleyear$animal,
                                                   birdcohort$BirdID)]
  
  head(maleyear)
  summary(maleyear)
  
  # check:
  sparrowpedigree[which(sparrowpedigree$birdid==4722),]
  maleyear[which(maleyear$animal==4722),]
  # good good.
}

{
  # male age. Unlike with lifespan, age is easier to define as the
  # year breeding occurred minus the cohort of the male.
  
  # add the year breeding occurred to the data set:
  maleyear$year <- offspring3$Cohort[match(maleyear$maleyear, 
                                           offspring3$maleyear)]
  
  head(maleyear)
  tail(maleyear)
  str(maleyear)
  
  # and calculate age by subtracting cohort from year breeding occurred:
  maleyear$age <- maleyear$year - maleyear$Cohort
  
  head(maleyear)
  table(maleyear$age)
  
  # as with lifespan, there will be problems analysing the males >6
  # years old because of sample size:
  
  maleyear$age6 <- maleyear$age
  maleyear$age6[which(maleyear$age6>6)] <- 6
  
  # and the male reproducing in year zero is a database problem:
  maleyear$age6[which(maleyear$age6==0)] <- 1
  
  table(maleyear$age)
  table(maleyear$age6)
  # looks good for analysis
  
  str(maleyear)
}

{
  # set appropriate factors for analysis
  str(maleyear)
  maleyear$factoranimal <- as.factor(maleyear$animal)
  maleyear$factormaleID <- as.factor(maleyear$animal)
  maleyear$factormaternalID <- as.factor(maleyear$maternalID)
  maleyear$factorpaternalID <- as.factor(maleyear$paternalID)
  maleyear$factorcohort <- as.factor(maleyear$Cohort)
  maleyear$factoryear <- as.factor(maleyear$year)
  str(maleyear)
}

{
  # One analysis is to find out whether the heritability
  # of male behaviour changes as males become older.
  # for this, two sets of models will be needed. One where
  # heritability is calculated within each age of male through
  # an interaction of the age factor with the animal term, and
  # one where the heritability of promiscuity is calculated 
  # within a specific age class, with one model per age class.
  # Therefore, for this second set produce six sets of data
  # frames, one for each age class:
  
  maleyear.age1 <- maleyear[which(maleyear$age6==1),]
  maleyear.age2 <- maleyear[which(maleyear$age6==2),]
  maleyear.age3 <- maleyear[which(maleyear$age6==3),]
  maleyear.age4 <- maleyear[which(maleyear$age6==4),]
  maleyear.age5 <- maleyear[which(maleyear$age6==5),]
  maleyear.age6 <- maleyear[which(maleyear$age6==6),]
  
  summary(maleyear.age1)
  summary(maleyear.age2)
  summary(maleyear.age3)
  summary(maleyear.age4)
  summary(maleyear.age5)
  summary(maleyear.age6)
  
  length(maleyear.age1$age6)
  length(maleyear.age2$age6)
  length(maleyear.age3$age6)
  length(maleyear.age4$age6)
  length(maleyear.age5$age6)
  length(maleyear.age6$age6)
  
  table(maleyear$age6)
  # numbers of males pretty small when it gets to four and
  # above. But see how the model copes. It might be that such
  # a model is only possible with ages 1-3 or 1-4.
}

##############################################################################
# female phenotypes 
##############################################################################

# Aggregate WPO and EPO by GeneticMumID

{
  # for offspring3, are there any cases where the genetic mother is not known?
  which(is.na(offspring3$GeneticMumID))
  # yes, which means these individuals need to be removed from making the 
  # dataset for females:
  offspring4 <- offspring3[-which(is.na(offspring3$GeneticMumID)),]
  
  which(is.na(offspring4$GeneticMumID))
  summary(offspring4)
  
  # aggregate by GeneticMumID
  femalephenotypes <- aggregate(offspring4$WPO, 
                                list(offspring4$GeneticMumID),
                                FUN=sum)
  head(femalephenotypes)
  
  
  femaleEPO <- aggregate(offspring4$EPO, 
                         list(offspring4$GeneticMumID),
                         FUN=sum)
  head(femaleEPO)
  
  # aggregate WPO by GeneticMumID
  names(femalephenotypes) <- c("animal", "WPO")
  head(femalephenotypes)
  
  str(femalephenotypes)
  str(femaleEPO)
  
  # bind together data frames:
  femalephenotypes <- cbind(femalephenotypes, femaleEPO)
  
  head(femalephenotypes)
  
  which(femalephenotypes$animal==femalephenotypes$Group.1)
  which(femalephenotypes$animal!=femalephenotypes$Group.1)
  # good!
  
  
  names(femalephenotypes) <- c("animal", "WPO", "femaleID", "EPO")
  head(femalephenotypes)
  
  # are any females in my file not in the pedigree?
  femalephenotypes$check <- fixedsparrowped$animal[match(femalephenotypes$animal,
                                                  fixedsparrowped$animal)]
  
  summary(femalephenotypes$check)
  # no NA values implies all females are in the pedigree
}

{
  # add maternal and paternal ID of each female:
  
  # maternal ID
  head(femalephenotypes)
  femalephenotypes$maternalID <- fixedsparrowped$dam[match(femalephenotypes$animal,
                                                         fixedsparrowped$animal)]
  
  # and paternal ID
  femalephenotypes$paternalID <- fixedsparrowped$sire[match(femalephenotypes$animal,
                                                          fixedsparrowped$animal)]
  
  head(femalephenotypes)
  
  # do a check
  fixedsparrowped[which(fixedsparrowped$animal==980),]
  femalephenotypes[which(femalephenotypes$animal==980),]
  
  # correct!
  
  fixedsparrowped[which(fixedsparrowped$animal==6500),]
  femalephenotypes[which(femalephenotypes$animal==6500),]
  
  # also correct!
  
  str(femalephenotypes)
}

{
  # Add the cohort of each female:
  femalephenotypes$Cohort <- birdcohort$Cohort[match(femalephenotypes$animal,
                                                   birdcohort$BirdID)]
  
  head(femalephenotypes)
  summary(femalephenotypes)
  
  # check:
  sparrowpedigree[which(sparrowpedigree$birdid==7022),]
  femalephenotypes[which(femalephenotypes$animal==7022),]
  # that is interesting. The cohort of this one is not the same between
  # the birdcohort file and the pedigree.
}

{
  # set appropriate factors for analysis
  str(femalephenotypes)
  femalephenotypes$factoranimal <- as.factor(femalephenotypes$animal)
  femalephenotypes$factormaternalID <- as.factor(femalephenotypes$maternalID)
  femalephenotypes$factorpaternalID <- as.factor(femalephenotypes$paternalID)
  femalephenotypes$factorcohort <- as.factor(femalephenotypes$Cohort)
  str(femalephenotypes)
}

{
  # female lifespan. See discussion in the same section of making the
  # male data set to see the issue with defining lifespan.
  # for now, for the purposes of this analysis, the last date seen
  # alive will be defined as the last time seen in the GENETIC PEDIGREE.
  
  # so, find the maximum cohort for a given dam in the pedigree:
  femalelifespan <- aggregate(sparrowpedigree$cohort, 
                            list(sparrowpedigree$dam),
                            FUN=max)
  
  head(femalelifespan)
  str(femalelifespan)
  
  # match to the data frame:
  femalephenotypes$lastyearPED <- femalelifespan$x[match(femalephenotypes$animal,
                                                     femalelifespan$Group.1)]
  
  head(femalephenotypes)
  
  # make lifespan variable by subtracting cohort from year last seen:
  femalephenotypes$lifespanPED <- femalephenotypes$lastyearPED - femalephenotypes$Cohort
  summary(femalephenotypes$lifespanPED)
  table(femalephenotypes$lifespanPED)
  # very few females older than 5.
  
  # make spare age variable:
  femalephenotypes$lifespanPED5 <- femalephenotypes$lifespanPED
  
  # assign >5 years old to 5 years old:
  femalephenotypes$lifespanPED5[which(femalephenotypes$lifespanPED5>5)] <- 5
  
  
  # check
  summary(femalephenotypes$lifespanPED5)
  table(femalephenotypes$lifespanPED5)
  table(femalephenotypes$lifespanPED5)
  # numbers correct.
}

# do broods have a single genetic mother?
length(offspring4$GeneticMumID)
table(offspring4$SocialMumID==offspring4$GeneticMumID)
# so there are 15 cases where the genetic mother is not the social mother.
# this could be a problem with either the pedigree or the database but I
# am aware that the database being used for this analysis needs updating,
# and that some social parents will be corrected as a result.

offspring4[which((offspring4$SocialMumID==offspring4$GeneticMumID)==FALSE),]
# The K, L, and M broods are from my years, where I am aware that changes are
# needed.

#-------------------------------------------
# Female dataset of genetic offspring per year
#-------------------------------------------


{
  # per year per GENETIC female identifier:
  offspring4$femaleyear <- paste(offspring4$GeneticMumID, offspring4$Cohort, sep=".")
  
  head(offspring4)
  table(offspring4$femaleyear)
  table(table(offspring4$femaleyear))
  # up to 22 offspring this time.
  
  
  # aggregate number of WPO per female
  femaleyearWPO <- aggregate(offspring4$WPO, 
                           list(offspring4$femaleyear),
                           FUN=sum)
  head(femaleyearWPO)
  
  # aggregate number of EPO per female
  femaleyearEPO <- aggregate(offspring4$EPO, 
                           list(offspring4$femaleyear),
                           FUN=sum)
  head(femaleyearEPO)
  
  
  names(femaleyearWPO) <- c("animal", "WPO")
  head(femaleyearWPO)
  
  str(femaleyearWPO)
  str(femaleyearEPO)
  
  
  femaleyear <- cbind(femaleyearWPO, femaleyearEPO)
  
  head(femaleyear)
  
  # check that all data rows are matched to the right female:
  which(femaleyear$animal==femaleyear$Group.1)
  # no missing numbers --> good
  which(femaleyear$animal!=femaleyear$Group.1)
  # no missmatch --> correct
  
  # rename data set
  names(femaleyear) <- c("femaleyear1", "WPO", "femaleyear", "EPO")
  head(femaleyear)
  
  
  # add the female's ID:
  femaleyear$animal <- offspring4$GeneticMumID[match(femaleyear$femaleyear, offspring4$femaleyear)]
  head(femaleyear)
  
  # add female ID for repeated measures:
  femaleyear$femaleid <- femaleyear$animal
  
}

{
  # add maternal ID for maternal effects:
  femaleyear$maternalID <- fixedsparrowped$dam[match(femaleyear$animal,
                                                   fixedsparrowped$animal)]
  head(femaleyear)
  
  # and paternal ID
  femaleyear$paternalID <- fixedsparrowped$sire[match(femaleyear$animal,
                                                    fixedsparrowped$animal)]
  
  summary(femaleyear)
}

{
  # Add the cohort of each female:
  femaleyear$Cohort <- birdcohort$Cohort[match(femaleyear$animal,
                                             birdcohort$BirdID)]
  
  head(femaleyear)
  summary(femaleyear)
  
  # check:
  sparrowpedigree[which(sparrowpedigree$birdid==76),]
  femaleyear[which(femaleyear$animal==76),]
  # good good.
}

{
  # female age. Unlike with lifespan, age is easier to define as the
  # year breeding occurred minus the cohort of the female.
  
  # add the year breeding occurred to the data set:
  femaleyear$year <- offspring4$Cohort[match(femaleyear$femaleyear, 
                                             offspring4$femaleyear)]
  
  head(femaleyear)
  tail(femaleyear)
  str(femaleyear)
  
  # and calculate age by subtracting cohort from year breeding occurred:
  femaleyear$age <- femaleyear$year - femaleyear$Cohort
  
  head(femaleyear)
  table(femaleyear$age)
  
  # as with lifespan, merge ages >5
  
  femaleyear$age5 <- femaleyear$age
  femaleyear$age5[which(femaleyear$age5>5)] <- 5
  
  table(femaleyear$age)
  table(femaleyear$age5)
  # looks good for analysis
  
  str(femaleyear)
}

{
  # set appropriate factors for analysis
  str(femaleyear)
  femaleyear$factoranimal <- as.factor(femaleyear$animal)
  femaleyear$factorfemaleID <- as.factor(femaleyear$animal)
  femaleyear$factormaternalID <- as.factor(femaleyear$maternalID)
  femaleyear$factorpaternalID <- as.factor(femaleyear$paternalID)
  femaleyear$factorcohort <- as.factor(femaleyear$Cohort)
  femaleyear$factoryear <- as.factor(femaleyear$year)
  str(femaleyear)
}


#-------------------------------------------
# Female dataset of genetic offspring per brood
#-------------------------------------------

{
  # per brood per female:
  
  broodWPO <- aggregate(offspring4$WPO,
                        list(offspring4$BroodName),
                        FUN=sum)
  
  broodEPO <- aggregate(offspring4$EPO,
                        list(offspring4$BroodName),
                        FUN=sum)
  
  summary(broodWPO)
  summary(broodEPO)
  str(broodWPO)
  str(broodEPO)
  table(broodEPO$x)
  # some broods are completely or almost completely EPO
  # (ones with four or five EPO). I wonder how many of
  # these are misidentified parents.
  
  broodWE <- data.frame(BroodName=broodWPO$Group.1,
                        WPO=broodWPO$x,
                        EPO=broodEPO$x)
  
  summary(broodWE)
  
  hist(broodWE$EPO/(broodWE$EPO+broodWE$WPO))
  # most broods have no EPO, so this might be a hard
  # analysis to run.
  
  # how many are all EPO?
  table(broodWE$WPO)
  # 63
  
  # add female ID
  broodWE$MumID <- offspring4$GeneticMumID[match(broodWE$BroodName,
                                                        offspring4$BroodName)]
  
  length(unique(broodWE$MumID))
  length(unique(offspring4$GeneticMumID))
  # there are more mothers in the offspring data set than go in
  # the per brood data set. This is because there are a small
  # number of broods that have more than one genetic mother. For
  # example:
  offspring4[which(offspring4$BroodName=="A008"),]
  broodWE[which(broodWE$BroodName=="A008"),]
  
  # add the female's social partner.
  broodWE$DadID <- offspring4$SocialDadID[match(broodWE$BroodName,
                                                        offspring4$BroodName)]
  
  # make a pair ID from the male and female's ID
  broodWE$pairID <- paste(broodWE$MumID, broodWE$DadID, sep="plus")
  
  # year of reproduction
  broodWE$year <- offspring4$Cohort[match(broodWE$BroodName,
                                          offspring4$BroodName)]
  
  head(broodWE)
  
  
  # check that this is all matching the original data:
  broodWE[which(broodWE$BroodName=="J050"),]
  offspring4[which(offspring4$BroodName=="J050"),]
}

{
  # add maternal ID for maternal effects:
  broodWE$maternalID <- fixedsparrowped$dam[match(broodWE$MumID,
                                                     fixedsparrowped$animal)]
  head(broodWE)
  
  # and paternal ID
  broodWE$paternalID <- fixedsparrowped$sire[match(broodWE$MumID,
                                                      fixedsparrowped$animal)]
  
  # though maternal ID of the male will also be needed (to estimate
  # the effect from his mother):
  broodWE$DadmaternalID <- fixedsparrowped$dam[match(broodWE$DadID,
                                                  fixedsparrowped$animal)]
  
  summary(broodWE)
  
  # check against original data:
  broodWE[which(broodWE$BroodName=="J050"),]
  sparrowpedigree[which(sparrowpedigree$birdid==4037),]
  sparrowpedigree[which(sparrowpedigree$birdid==4753),]
}

{
  # Add the cohort of each female:
  broodWE$Cohort <- birdcohort$Cohort[match(broodWE$MumID,
                                               birdcohort$BirdID)]
  
  head(broodWE)
  summary(broodWE)
  
  # check:
  sparrowpedigree[which(sparrowpedigree$birdid==4739),]
  broodWE[which(broodWE$MumID==4739),]
  # good good.
}

{
  # female age. Age is defined as the
  # year breeding occurred minus the cohort of the female.
  
  str(broodWE)
  
  # calculate age by subtracting cohort from year breeding occurred:
  broodWE$age <- broodWE$year - broodWE$Cohort
  
  head(broodWE)
  table(broodWE$age)
  
  # merge ages >5
  
  broodWE$age5 <- broodWE$age
  broodWE$age5[which(broodWE$age5>5)] <- 5
  
  table(broodWE$age)
  table(broodWE$age5)
  # looks good for analysis
  
  str(broodWE)
}

{
  # set appropriate factors for analysis
  str(broodWE)
  broodWE$factoranimal <- as.factor(broodWE$MumID)
  broodWE$factorfemaleID <- as.factor(broodWE$MumID)
  broodWE$factormaternalID <- as.factor(broodWE$maternalID)
  broodWE$factorDadmaternalID <- as.factor(broodWE$DadmaternalID)
  broodWE$factorpaternalID <- as.factor(broodWE$paternalID)
  broodWE$factorcohort <- as.factor(broodWE$Cohort)
  broodWE$factoryear <- as.factor(broodWE$year)
  broodWE$factordadanimal <- factor(broodWE$DadID)
  broodWE$factorDadID <- factor(broodWE$DadID)
  broodWE$factorpairID <- factor(broodWE$pairID)
  str(broodWE)
}

{  
  # add the age of the social father:
  
  head(birdcohort)
  
  broodWE$socdadcohort <- birdcohort$Cohort[match(broodWE$DadID,
                                                  birdcohort$BirdID)]
  
  broodWE$socdadage <- broodWE$year - broodWE$socdadcohort
  
  summary(broodWE)
  table(broodWE$socdadage)
  
  # amalgamate >6
  broodWE$dadage6 <- broodWE$socdadage
  broodWE$dadage6[which(broodWE$dadage6>6)] <- 6
  
  # make the zero a one:
  broodWE$dadage6[which(broodWE$dadage6==0)] <- 1
  
  table(broodWE$dadage6)
}

##############################################################################
# Pruning the pedigree
##############################################################################

# each pedigree can be pruned to speed up the fitting of models.
# Pedigrees are pruned to only contain phenotyped individuals and their
# relatives. Uninformative individuals are removed.
{
  # prunePed requires a numeric bird ID so check that the fixed pedigree is numeric:
  # (https://stat.ethz.ch/pipermail/r-sig-mixed-models/2014q1/021604.html)
  
  str(fixedsparrowped)
  
  # now prune the pedigree
  prunedped.malelifetime <- prunePed(pedigree=fixedsparrowped, 
                                     keep=malephenotypes$animal)
  
  head(prunedped.malelifetime)
  
  # take a look at the pruned info:
  drawPedigree(prunedped.malelifetime)
  pedigreeinformation <- pedigreeStats(prunedped.malelifetime)
  pedStatSummary(pedigreeinformation)
  
  # for females:
  prunedped.femalelifetime <- prunePed(pedigree=fixedsparrowped, 
                                     keep=femalephenotypes$animal)
  
  head(prunedped.femalelifetime)
  
  
  # for both males and females ***for the broodWE data set
  # i.e. parents differ because 1 - in a couple of cases there
  # is more than one female per brood 2 - the maleID is the social
  # male of the brood***
  keepers <- c(broodWE$MumID, broodWE$DadID)
  prunedped.bothsexes <- prunePed(pedigree=fixedsparrowped, 
                                  keep=keepers)
  head(prunedped.bothsexes)
}

# and then the pedigrees can be converted to relationship matrices
# for the current version of MCMCglmm:
{
  invped.malelifetime <- inverseA(prunedped.malelifetime)$Ainv
  summary(invped.malelifetime)
  
  invped.femalelifetime <- inverseA(prunedped.femalelifetime)$Ainv
  summary(invped.femalelifetime)
  
  invped.bothsexes <- inverseA(prunedped.bothsexes)$Ainv
  summary(invped.bothsexes)
}

##############################################################################
# Plots of variables
##############################################################################
{
# how many EPO offspring do males have?
hist(malephenotypes$EPO, breaks=seq(-1,26,1))
# lots of zeroes.

# how many are EPO?
hist(malephenotypes$EPO/(malephenotypes$EPO+malephenotypes$WPO))
# lots of males with no EPO, some males with only EPO

# in the per year data?
hist(maleyear$EPO, breaks=seq(-1,16,1))
hist(maleyear$EPO/(maleyear$EPO+maleyear$WPO))

# what about in the female data set?
# here the main measure is proportion of EPO
hist(femalephenotypes$EPO/(femalephenotypes$EPO+femalephenotypes$WPO))
# much more variation compared to the males data set.


}
##############################################################################
# Priors
##############################################################################

# parameter expanded priors. Number before G is number of random effects. 
# .p for parameter expansion
{
  prior1G.p = list(R=list(V=1,nu=0.002),
                   G=list(G1=list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000)))
  
  prior2G.p = list(R=list(V=1,nu=0.002),
                   G=list(G1=list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000),
                          G2=list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000)))
  
  prior3G.p = list(R=list(V=1,nu=0.002),
                  G=list(G1=list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000),
                         G2=list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000),
                         G3=list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000)))
  
  prior4G.p = list(R=list(V=1,nu=0.002),
                   G=list(G1=list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000),
                          G2=list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000),
                          G3=list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000),
                          G4=list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000)))
  
  prior5G.p = list(R=list(V=1,nu=0.002),
                   G=list(G1=list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000),
                          G2=list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000),
                          G3=list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000),
                          G4=list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000),
                          G5=list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000)))
  
  prior6G.p = list(R=list(V=1,nu=0.002),
                   G=list(G1=list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000),
                          G2=list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000),
                          G3=list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000),
                          G4=list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000),
                          G5=list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000),
                          G6=list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000)))
}

# parameter expanded priors with a covariance structure on the first
# random effect
# .p for parameter expanded, special for covariance structure
{
  priorspecial4G.p = list(R=list(V=1,nu=0.002),
                          G=list(G1=list(V=diag(2), nu=2, alpha.mu=c(0,0), alpha.V=diag(2)*1000),
                                 G2=list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000),
                                 G3=list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000),
                                 G4=list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000)))
  
  priorspecial5G.p = list(R=list(V=1,nu=0.002),
                   G=list(G1=list(V=diag(2), nu=2, alpha.mu=c(0,0), alpha.V=diag(2)*1000),
                          G2=list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000),
                          G3=list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000),
                          G4=list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000),
                          G5=list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000)))
  
  priorspecial8G.p = list(R=list(V=1,nu=0.002),
                          G=list(G1=list(V=diag(2), nu=2, alpha.mu=c(0,0), alpha.V=diag(2)*1000),
                                 G2=list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000),
                                 G3=list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000),
                                 G4=list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000),
                                 G5=list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000),
                                 G6=list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000),
                                 G7=list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000),
                                 G8=list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000)))
}

# parameter expanded, covariance on the first TWO random effects
# .p for parameter expanded, 2special for two covariances
{
  prior2special6G.p = list(R=list(V=1,nu=0.002),
                          G=list(G1=list(V=diag(2), nu=2, alpha.mu=c(0,0), alpha.V=diag(2)*1000),
                                 G2=list(V=diag(2), nu=2, alpha.mu=c(0,0), alpha.V=diag(2)*1000),
                                 G3=list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000),
                                 G4=list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000),
                                 G5=list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000),
                                 G6=list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000)))
}

# for penetrance with male age models, a prior with a six-level
# covariance on the first random effect:
{
  prior4G.p.penetrance = list(R=list(V=1,nu=0.002),
                 G=list(G1=list(V=diag(6), nu=6, alpha.mu=c(0,0,0,0,0,0), alpha.V=diag(6)*1000),
                        G2=list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000),
                        G3=list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000),
                        G4=list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000)))
}


##############################################################################
# Repeatability
##############################################################################

{
  maler2EPOmulti <- MCMCglmm(cbind(EPO, WPO)~1,
                      ~factormaleID,
                      prior=prior1G.p,
                      data=maleyear,
                      family="multinomial2",
                      nitt=1000000,
                      thin=800,
                      burnin=200000)

  plot(maler2EPOmulti)
  # low to no repeatability for the male multinomial behaviour.
  # but this behaviour isn't well described by a multinomial
  # distribution in males anyway.
  autocorr(maler2EPOmulti$Sol)
  autocorr(maler2EPOmulti$VCV)
  # Good.
}

{
  maler2EPOpois <- MCMCglmm(EPO~1,
                        ~factormaleID,
                        prior=prior1G.p,
                        data=maleyear,
                        family="poisson",
                        nitt=1000000,
                        thin=800,
                        burnin=200000)
  plot(maler2EPOpois)
  # This is repeatable. Not strongly, and there is a bit 
  # of overlap with zero, but it is. Might improve if age
  # is added since EP behaviour changes with age. Could also
  # be contingent on the year of measurement of behaviour,
  # so a model with maleID, year, and a fixed factor of male
  # age would be good.
  autocorr(maler2EPOpois$Sol)
  autocorr(maler2EPOpois$VCV)
  # Below 0.1.
}

# male poisson with male age as fixed factor, and year of 
# reproduction as additional random effects:
{
  maler2EPOpois.age.yr <- MCMCglmm(EPO~factor(age6),
                            ~factormaleID + factoryear,
                            prior=prior2G.p,
                            data=maleyear,
                            family="poisson",
                            nitt=1000000,
                            thin=800,
                            burnin=200000)
  plot(maler2EPOpois.age.yr)
  # the repeatability per male is more clearly defined in this
  # case, now that the age of the male and the between-year
  # variation in the number of EPO is accounted for.
  autocorr(maler2EPOpois.age.yr$Sol)
  autocorr(maler2EPOpois.age.yr$VCV)
  # 0.8 a bit high but they are all low ish.
}

{
  femaler2EPOmulti <- MCMCglmm(cbind(EPO,WPO)~1,
                            ~factorfemaleID,
                            prior=prior1G.p,
                            data=femaleyear,
                            family="multinomial2",
                            nitt=1000000,
                            thin=800,
                            burnin=200000)
  plot(femaler2EPOmulti)
  # Nicely repeatable.
  autocorr(femaler2EPOmulti$Sol)
  autocorr(femaler2EPOmulti$VCV)
  # Good mixing.
}

# as with males, include extra factors, but this time to see
# whether repeatability remains or is a result of these 
# effects:
# female multinomial with female age as fixed factor, and year:
{
  femaler2EPOmulti.age.yr <- MCMCglmm(cbind(EPO,WPO)~factor(age5),
                              ~factorfemaleID + factoryear,
                              prior=prior2G.p,
                              data=femaleyear,
                              family="multinomial2",
                              nitt=1000000,
                              thin=800,
                              burnin=200000)
  plot(femaler2EPOmulti.age.yr)
  # Not much purpose to including age. Little effect of 
  # year of reproduction on this trait, suggesting the
  # trait is more likely to be a set promiscuity level
  # within each female (though not proof).
  autocorr(femaler2EPOmulti.age.yr$Sol)
  autocorr(femaler2EPOmulti.age.yr$VCV)
  # age2 to intercept has a value of -0.11. Too high.
  # Though negative so probably doesn't make a difference.
}

# but female age has virtually no effect --> just year:
{
  femaler2EPOmulti.yr <- MCMCglmm(cbind(EPO,WPO)~1,
                                      ~factorfemaleID + factoryear,
                                      prior=prior2G.p,
                                      data=femaleyear,
                                      family="multinomial2",
                                      nitt=1000000,
                                      thin=800,
                                      burnin=200000)
  plot(femaler2EPOmulti.yr)
  # female repeatable here, no effect of year of reproduction.
  autocorr(femaler2EPOmulti.yr$Sol)
  autocorr(femaler2EPOmulti.yr$VCV)
  # fine.
}

# as a per brood analysis for females:

{
  femaler2EPObrood <- MCMCglmm(cbind(EPO,WPO)~1,
                               ~factorfemaleID,
                               prior=prior1G.p,
                               data=broodWE,
                               family="multinomial2",
                               nitt=1000000,
                               thin=800,
                               burnin=200000)
  plot(femaler2EPObrood)
  # Variance values are double from before, but very similar
  # proportions. Females are repeatable between broods according
  # to this.
  autocorr(femaler2EPObrood$Sol)
  autocorr(femaler2EPObrood$VCV)
  # good.
}

# as with males, include extra factors, but this time to see
# whether repeatability remains or is a result of these 
# effects:
# female with year:
{
  femaler2EPObrood.yr <- MCMCglmm(cbind(EPO,WPO)~1,
                                      ~factorfemaleID + factoryear,
                                      prior=prior2G.p,
                                      data=broodWE,
                                      family="multinomial2",
                                      nitt=1000000,
                                      thin=800,
                                      burnin=200000)
  plot(femaler2EPObrood.yr)
  # Small amount of the femaleID variance has been absorbed
  # in to year. Year is not a significant value but is not
  # zero. Female is clearly repeatable.
  autocorr(femaler2EPObrood.yr$Sol)
  autocorr(femaler2EPObrood.yr$VCV)
}



##############################################################################
# Heritability of male behaviour
##############################################################################

# a model of heritability should have as a minimum the ID of the bird (if there
# are repeated measures), a pedigree linked bird ID, and a maternal ID.
# Strongly recommended: cohort

#-------------------------------------------
# Lifetime male behaviour as a multinomial variable
#-------------------------------------------

# with pedigree term, maternal ID, and male cohort
{
  maleh2EPO.multinomial.lifetime <- MCMCglmm(cbind(EPO, WPO)~1,
                               ~factoranimal + factormaternalID + factorcohort,
                               ginverse=list(factoranimal=invped.malelifetime),
                               prior=prior3G.p,
                               data=malephenotypes,
                               family="multinomial2",
                               nitt=1000000,
                               thin=800,
                               burnin=200000)
  
  plot(maleh2EPO.multinomial.lifetime)
  # strong maternal ID effect, little animal and cohort effects
  autocorr(maleh2EPO.multinomial.lifetime$Sol)
  autocorr(maleh2EPO.multinomial.lifetime$VCV)
  # animal:animal is a little high - this needs to be run for longer
}

# also accounting for male lifespan as a fixed factor between 1 and 6:
{
  maleh2EPO.multinomial.lifetime.lifespan <- MCMCglmm(cbind(EPO, WPO)~factor(lifespanPED6),
                                             ~factoranimal + factormaternalID + factorcohort,
                                             ginverse=list(factoranimal=invped.malelifetime),
                                             prior=prior3G.p,
                                             data=malephenotypes,
                                             family="multinomial2",
                                             nitt=1000000,
                                             thin=800,
                                             burnin=200000)
  
  plot(maleh2EPO.multinomial.lifetime.lifespan)
  # clearer cohort effects. Suggestion that two and three year olds
  # have fewer EPO than one year olds and one years are similar to 
  # four and more year olds in proportion of EPO. Perhaps this is because
  # Perhaps this is because we miss the WPO of some males and only detect
  # their existence through EPO, so this multinomial approach is not
  # ideal for males where WPO are not guaranteed.
  
  autocorr(maleh2EPO.multinomial.lifetime.lifespan$Sol)
  autocorr(maleh2EPO.multinomial.lifetime.lifespan$VCV)
  # autocorrelation good!
}

# additional consideration: whether incidental assortative mating by promiscuity
# leads to an artificially strong maternal effect. In this case, we would expect
# a paternal effect of comparable magnitude to the maternal effect.

{
  maleh2EPO.multinomial.lifetime.patID <- MCMCglmm(cbind(EPO, WPO)~1,
                                             ~factoranimal + factormaternalID + 
                                               factorpaternalID + factorcohort,
                                             ginverse=list(factoranimal=invped.malelifetime),
                                             prior=prior4G.p,
                                             data=malephenotypes,
                                             family="multinomial2",
                                             nitt=1000000,
                                             thin=800,
                                             burnin=200000)
  
  plot(maleh2EPO.multinomial.lifetime.patID)
  # with the addition of paternal ID, the maternal effect is much harder to 
  # estimate and the paternal effect is of comparable size, so maybe it is
  # an assortative mating thing. Cohort effects are clearer now.
  autocorr(maleh2EPO.multinomial.lifetime.patID$Sol)
  autocorr(maleh2EPO.multinomial.lifetime.patID$VCV)
  # fine
}

# paternal ID and lifespan...

#-------------------------------------------
# Lifetime male behaviour as a poisson variable
#-------------------------------------------

{
  maleh2EPO.poisson.lifetime <- MCMCglmm(EPO~1,
                                             ~factoranimal + factormaternalID + factorcohort,
                                             ginverse=list(factoranimal=invped.malelifetime),
                                             prior=prior3G.p,
                                             data=malephenotypes,
                                             family="poisson",
                                             nitt=1000000,
                                             thin=800,
                                             burnin=200000)
  
  plot(maleh2EPO.poisson.lifetime)
  # much more consistently estimated maternal ID and cohort in this
  # compared to the multinomial model
  autocorr(maleh2EPO.poisson.lifetime$Sol)
  autocorr(maleh2EPO.poisson.lifetime$VCV)
  # all fine enough.
}

# also accounting for male lifespan as a fixed factor between 1 and 6:
{
  maleh2EPO.poisson.lifetime.lifespan <- MCMCglmm(EPO~factor(lifespanPED6),
                                         ~factoranimal + factormaternalID + factorcohort,
                                         ginverse=list(factoranimal=invped.malelifetime),
                                         prior=prior3G.p,
                                         data=malephenotypes,
                                         family="poisson",
                                         nitt=1000000,
                                         thin=800,
                                         burnin=200000)
  
  plot(maleh2EPO.poisson.lifetime.lifespan)
  # YES! See how now the distribution is a poisson and therefore unbounded by
  # a male's WPO, the number of EPO increases in older males.
  # As before, maternal effects are strongest in this model.
  
  autocorr(maleh2EPO.poisson.lifetime.lifespan$Sol)
  autocorr(maleh2EPO.poisson.lifetime.lifespan$VCV)
  # looks good.
}

# with paternal ID included:

{
  maleh2EPO.poisson.lifetime.patID <- MCMCglmm(EPO~1,
                                         ~factoranimal + factormaternalID +
                                           factorpaternalID + factorcohort,
                                         ginverse=list(factoranimal=invped.malelifetime),
                                         prior=prior4G.p,
                                         data=malephenotypes,
                                         family="poisson",
                                         nitt=1000000,
                                         thin=800,
                                         burnin=200000)
  
  plot(maleh2EPO.poisson.lifetime.patID)
  # now, no paternal variance estimated. Maternal and cohort stay strong.
  autocorr(maleh2EPO.poisson.lifetime.patID$Sol)
  autocorr(maleh2EPO.poisson.lifetime.patID$VCV)
  # again, good sampling.
}

# and paternal ID with lifespan:
{
  maleh2EPO.poisson.lifetime.lifespan.patID <- MCMCglmm(EPO~factor(lifespanPED6),
                                                  ~factoranimal + factormaternalID +
                                                    factorpaternalID + factorcohort,
                                                  ginverse=list(factoranimal=invped.malelifetime),
                                                  prior=prior4G.p,
                                                  data=malephenotypes,
                                                  family="poisson",
                                                  nitt=1000000,
                                                  thin=800,
                                                  burnin=200000)
  
  plot(maleh2EPO.poisson.lifetime.lifespan.patID)
  # again, maternal ID matters out of all these random effects.
  # EPO goes up with lifespan.
  autocorr(maleh2EPO.poisson.lifetime.lifespan.patID$Sol)
  autocorr(maleh2EPO.poisson.lifetime.lifespan.patID$VCV)
  # fine
}

# problems with these models: males that are older can contribute more
# EPO. A per-year approach is preferable.

#-------------------------------------------
# Per year male behaviour as a multinomial variable
#-------------------------------------------


# with bird ID term, pedigree term, maternal ID, and male cohort
{
  maleh2EPO.multi.yr <- MCMCglmm(cbind(EPO, WPO)~1,
                                 ~factoranimal + factormaleID +
                                   factormaternalID + factorcohort,
                                 ginverse=list(factoranimal=invped.malelifetime),
                                 prior=prior4G.p,
                                 data=maleyear,
                                 family="multinomial2",
                                 nitt=1000000,
                                 thin=800,
                                 burnin=200000)
  plot(maleh2EPO.multi.yr)
  # This one struggles to estimate anything.
  autocorr(maleh2EPO.multi.yr$Sol)
  autocorr(maleh2EPO.multi.yr$VCV)
  # and not because there is autocorrelation
}

# with age as a fixed factor between 1 and 6
{
  maleh2EPO.multi.yr.age <- MCMCglmm(cbind(EPO, WPO)~factor(age6),
                                 ~factoranimal + factormaleID +
                                   factormaternalID + factorcohort,
                                 ginverse=list(factoranimal=invped.malelifetime),
                                 prior=prior4G.p,
                                 data=maleyear,
                                 family="multinomial2",
                                 nitt=1000000,
                                 thin=800,
                                 burnin=200000)
  plot(maleh2EPO.multi.yr.age)
  # Once age is accounted for, there is a clear cohort effect but not much more.
  autocorr(maleh2EPO.multi.yr.age$Sol)
  autocorr(maleh2EPO.multi.yr.age$VCV)
  # That will do
}

# with paternalID
{
  maleh2EPO.multi.yr.pat <- MCMCglmm(cbind(EPO, WPO)~1,
                                         ~factoranimal + factormaleID +
                                           factormaternalID + factorpaternalID +
                                           factorcohort,
                                         ginverse=list(factoranimal=invped.malelifetime),
                                         prior=prior5G.p,
                                         data=maleyear,
                                         family="multinomial2",
                                         nitt=1000000,
                                         thin=800,
                                         burnin=200000)
  plot(maleh2EPO.multi.yr.pat)
  # even less going on than before (though no male age in this model)
  autocorr(maleh2EPO.multi.yr.pat$Sol)
  autocorr(maleh2EPO.multi.yr.pat$VCV)
  # This looks fine but I think the chain behaviour is a little wobbly
}

# with paternalID and age
{
  maleh2EPO.multi.yr.pat.age <- MCMCglmm(cbind(EPO, WPO)~factor(age6),
                                     ~factoranimal + factormaleID +
                                       factormaternalID + factorpaternalID +
                                       factorcohort,
                                     ginverse=list(factoranimal=invped.malelifetime),
                                     prior=prior5G.p,
                                     data=maleyear,
                                     family="multinomial2",
                                     nitt=1000000,
                                     thin=800,
                                     burnin=200000)
  plot(maleh2EPO.multi.yr.pat.age)
  # Hum! A small amount of variance to the paternal ID here where before there
  # was none! Cohort is stronger. No effect is significant. So this means the
  # maternal effect on this trait breaks down / is harder to quantify in this data.
  autocorr(maleh2EPO.multi.yr.pat.age$Sol)
  autocorr(maleh2EPO.multi.yr.pat.age$VCV)
  # good.
}

# with year of breeding as a random effect, paternal ID, and age
{
  maleh2EPO.multi.yr.pat.age.year <- MCMCglmm(cbind(EPO, WPO)~factor(age6),
                                         ~factoranimal + factormaleID +
                                           factormaternalID + factorpaternalID +
                                           factorcohort + factoryear,
                                         ginverse=list(factoranimal=invped.malelifetime),
                                         prior=prior6G.p,
                                         data=maleyear,
                                         family="multinomial2",
                                         nitt=1000000,
                                         thin=800,
                                         burnin=200000)
  plot(maleh2EPO.multi.yr.pat.age.year)
  # Most things show nothing. Something loaded on paternal ID but not enough
  # to separate. Cohort easier to separate, still not easily distinguished
  # effect
  autocorr(maleh2EPO.multi.yr.pat.age.year$Sol)
  autocorr(maleh2EPO.multi.yr.pat.age.year$VCV)
  # factoryear:factoryear, -0.985 not ideal. I think the model is too
  # much for the data.
}

#-------------------------------------------
# Per year male behaviour as a poisson variable
#-------------------------------------------

# with bird ID term, pedigree term, maternal ID, and male cohort
{
  maleh2EPO.pois.yr <- MCMCglmm(EPO~1,
                                 ~factoranimal + factormaleID +
                                   factormaternalID + factorcohort,
                                 ginverse=list(factoranimal=invped.malelifetime),
                                 prior=prior4G.p,
                                 data=maleyear,
                                 family="poisson",
                                 nitt=1000000,
                                 thin=800,
                                 burnin=200000)
  
  plot(maleh2EPO.pois.yr)
  # as before, maternal ID is estimated where the other factors, including
  # cohort, are not.
  autocorr(maleh2EPO.pois.yr$Sol)
  autocorr(maleh2EPO.pois.yr$VCV)
  # fine.
}

# with age as a fixed factor between 1 and 6
{
  maleh2EPO.pois.yr.age <- MCMCglmm(EPO~factor(age6),
                                     ~factoranimal + factormaleID +
                                       factormaternalID + factorcohort,
                                     ginverse=list(factoranimal=invped.malelifetime),
                                     prior=prior4G.p,
                                     data=maleyear,
                                     family="poisson",
                                     nitt=1000000,
                                     thin=800,
                                     burnin=200000)
  plot(maleh2EPO.pois.yr.age)
  # again, age important but out of the random effects only mother is 
  # strong.
  autocorr(maleh2EPO.pois.yr.age$Sol)
  autocorr(maleh2EPO.pois.yr.age$VCV)
  # units:animal has a strong one. Should run for longer. Might be random.
}

# with paternalID
{
  maleh2EPO.pois.yr.pat <- MCMCglmm(EPO~1,
                                    ~factoranimal + factormaleID +
                                      factormaternalID + factorpaternalID +
                                      factorcohort,
                                    ginverse=list(factoranimal=invped.malelifetime),
                                    prior=prior5G.p,
                                    data=maleyear,
                                    family="poisson",
                                    nitt=1000000,
                                    thin=800,
                                    burnin=200000)
  plot(maleh2EPO.pois.yr.pat)
  # maternal ID the best estimated random effect, no other random 
  # effects estimated.
  autocorr(maleh2EPO.pois.yr.pat$Sol)
  autocorr(maleh2EPO.pois.yr.pat$VCV)
  # all good.
}

# with paternalID and age
{
  maleh2EPO.pois.yr.age.pat <- MCMCglmm(EPO~factor(age6),
                                        ~factoranimal + factormaleID +
                                          factormaternalID + factorpaternalID +
                                          factorcohort,
                                        ginverse=list(factoranimal=invped.malelifetime),
                                        prior=prior5G.p,
                                        data=maleyear,
                                        family="poisson",
                                        nitt=1000000,
                                        thin=800,
                                        burnin=200000)
  plot(maleh2EPO.pois.yr.age.pat)
  # the expected effects of increased EPO with age, as before the
  # only random effect is maternal. Not great estimation, plenty of
  # overlap with zero.
  autocorr(maleh2EPO.pois.yr.age.pat$Sol)
  autocorr(maleh2EPO.pois.yr.age.pat$VCV)
  # fine.
}

# with year of breeding as a random effect, paternal ID, and age
{
  maleh2EPO.pois.yr.age.pat.year <- MCMCglmm(EPO~factor(age6),
                                        ~factoranimal + factormaleID +
                                          factormaternalID + factorpaternalID +
                                          factorcohort + factoryear,
                                        ginverse=list(factoranimal=invped.malelifetime),
                                        prior=prior6G.p,
                                        data=maleyear,
                                        family="poisson",
                                        nitt=1000000,
                                        thin=800,
                                        burnin=200000)
  
  plot(maleh2EPO.pois.yr.age.pat.year)
  # maybe it is just going too far with the available data. Only
  # a maternal effect here as well. It might be that it is only
  # possible to present a full model in a supplement since the 
  # random effect estimations are weak.
  autocorr(maleh2EPO.pois.yr.age.pat.year$Sol)
  autocorr(maleh2EPO.pois.yr.age.pat.year$VCV)
  # as previously, all pass the test of being below 0.1.
}

# since very low signal, suggestive of overspecified model.
# try model with only male ID, animal, and maternal ID

{
  maleh2EPO.pois.yr.mincohort <- MCMCglmm(EPO~1,
                                ~factoranimal + factormaleID +
                                  factormaternalID,
                                ginverse=list(factoranimal=invped.malelifetime),
                                prior=prior3G.p,
                                data=maleyear,
                                family="poisson",
                                nitt=1000000,
                                thin=800,
                                burnin=200000)
  plot(maleh2EPO.pois.yr.mincohort)
  # maternal ID strongest as before.
  autocorr(maleh2EPO.pois.yr.mincohort$Sol)
  autocorr(maleh2EPO.pois.yr.mincohort$VCV)
  # fine
}

# and account for male age and year that EPO was measured to eliminate
# important sources of noise
{
  maleh2EPO.pois.yr.mincohort.year.age <- MCMCglmm(EPO~factor(age6),
                                          ~factoranimal + factormaleID +
                                            factormaternalID + factoryear,
                                          ginverse=list(factoranimal=invped.malelifetime),
                                          prior=prior4G.p,
                                          data=maleyear,
                                          family="poisson",
                                          nitt=1000000,
                                          thin=800,
                                          burnin=200000)
  plot(maleh2EPO.pois.yr.mincohort.year.age)
  # still strongest maternal ID.
  autocorr(maleh2EPO.pois.yr.mincohort.year.age$Sol)
  autocorr(maleh2EPO.pois.yr.mincohort.year.age$VCV)
  # some 0.8 autocorrelations.
}


#-------------------------------------------
# Per year male poisson - does promiscuity have penetrance?
#-------------------------------------------

# one possibility is that male promiscuity is heritable but 
# males cannot express their true levels of promiscuity until
# they are older. (this is confounded with the fact that older
# males could be better males and have a consistent EPO rate 
# for this or other reasons).

# so, calculate heritability for each age within males:
{
  maleh2EPO.age1 <- MCMCglmm(EPO~1,
                             ~factoranimal +  factormaternalID + factoryear,
                             ginverse=list(factoranimal=invped.malelifetime),
                             prior=prior3G.p,
                             data=maleyear.age1,
                             family="poisson",
                             nitt=1000000,
                             thin=800,
                             burnin=200000)
  
  plot(maleh2EPO.age1)
  # maternal ID comes out strongly.
}

{
  maleh2EPO.age2 <- MCMCglmm(EPO~1,
                             ~factoranimal +  factormaternalID + factoryear,
                             ginverse=list(factoranimal=invped.malelifetime),
                             prior=prior3G.p,
                             data=maleyear.age2,
                             family="poisson",
                             nitt=1000000,
                             thin=800,
                             burnin=200000)
  plot(maleh2EPO.age2)
}

{
  maleh2EPO.age3 <- MCMCglmm(EPO~1,
                             ~factoranimal +  factormaternalID + factoryear,
                             ginverse=list(factoranimal=invped.malelifetime),
                             prior=prior3G.p,
                             data=maleyear.age3,
                             family="poisson",
                             nitt=1000000,
                             thin=800,
                             burnin=200000)
  plot(maleh2EPO.age3)
  # more evidence of an animal term here. Still not strong.
}

{
  maleh2EPO.age4 <- MCMCglmm(EPO~1,
                             ~factoranimal +  factormaternalID + factoryear,
                             ginverse=list(factoranimal=invped.malelifetime),
                             prior=prior3G.p,
                             data=maleyear.age4,
                             family="poisson",
                             nitt=1000000,
                             thin=800,
                             burnin=200000)
  plot(maleh2EPO.age4)
}

{
  maleh2EPO.age5 <- MCMCglmm(EPO~1,
                             ~factoranimal +  factormaternalID + factoryear,
                             ginverse=list(factoranimal=invped.malelifetime),
                             prior=prior3G.p,
                             data=maleyear.age5,
                             family="poisson",
                             nitt=1000000,
                             thin=800,
                             burnin=200000)
  plot(maleh2EPO.age5)
}

{
  maleh2EPO.age6 <- MCMCglmm(EPO~1,
                             ~factoranimal +  factormaternalID + factoryear,
                             ginverse=list(factoranimal=invped.malelifetime),
                             prior=prior3G.p,
                             data=maleyear.age6,
                             family="poisson",
                             nitt=1000000,
                             thin=800,
                             burnin=200000)
  plot(maleh2EPO.age6)
  # nowt. I expected that though - the model has a sample size of about 30!
  autocorr(maleh2EPO.age6$Sol)
  autocorr(maleh2EPO.age6$VCV)
  # fine and dandy.
}

# and what should be a similar model: interact Va with age
# in the random effects (male ID and age cannot interact
# because each male has one data point per age):
# age has an us() structure because males can be measured at
# multiple ages, so there can be covariance:

{
  maleh2EPO.pois.byage <- MCMCglmm(EPO~factor(age6),
                                   ~us(factor(age6)):factoranimal + factormaleID +
                                     factormaternalID + factoryear,
                                   ginverse=list(factoranimal=invped.malelifetime),
                                   prior=prior4G.p.penetrance,
                                   data=maleyear,
                                   family="poisson",
                                   nitt=100000,
                                   thin=80,
                                   burnin=20000)
  plot(maleh2EPO.pois.byage)
  # so, hints of some Va at age 2-3, but is it right to have
  # all these covariances? Why have a genetic covariance between
  # age 1 and age 6, for example?
  # maternal ID is still there as the best estimated random effect.
}

{
  # with idh() variance structure:
  maleh2EPO.pois.byage.idh <- MCMCglmm(EPO~factor(age6),
                                   ~idh(factor(age6)):factoranimal + factormaleID +
                                     factormaternalID + factoryear,
                                   ginverse=list(factoranimal=invped.malelifetime),
                                   prior=prior4G.p.penetrance,
                                   data=maleyear,
                                   family="poisson",
                                   nitt=1000000,
                                   thin=800,
                                   burnin=200000)
}

#-------------------------------------------
# Lifetime female behaviour as a multinomial variable
#-------------------------------------------

# with pedigree term, maternal ID, and female cohort
{
  femaleh2EPO.multinomial.lifetime <- MCMCglmm(cbind(EPO, WPO)~1,
                                             ~factoranimal + factormaternalID + factorcohort,
                                             ginverse=list(factoranimal=invped.femalelifetime),
                                             prior=prior3G.p,
                                             data=femalephenotypes,
                                             family="multinomial2",
                                             nitt=1000000,
                                             thin=800,
                                             burnin=200000)
  plot(femaleh2EPO.multinomial.lifetime)
  # not much is apparent in this mode. It looks like the model
  # wants to calculate some variance for animal or maternal ID.
  # no variance is given to cohort.
  autocorr(femaleh2EPO.multinomial.lifetime$Sol)
  autocorr(femaleh2EPO.multinomial.lifetime$VCV)
  # ok.
}

# also accounting for female lifespan as a fixed factor between 1 and 5:
{
  femaleh2EPO.multinomial.lifetime.lifespan <- MCMCglmm(cbind(EPO, WPO)~factor(lifespanPED5),
                                                      ~factoranimal + factormaternalID + factorcohort,
                                                      ginverse=list(factoranimal=invped.femalelifetime),
                                                      prior=prior3G.p,
                                                      data=femalephenotypes,
                                                      family="multinomial2",
                                                      nitt=1000000,
                                                      thin=800,
                                                      burnin=200000)
  plot(femaleh2EPO.multinomial.lifetime.lifespan)
  # suggesting fewer EPO in longer lived females... likely due to a 
  # more accurate information on longer lived females.
  # slightly more ability to distinguish maternal ID. Suggests maternal
  # ID could matter.
  autocorr(femaleh2EPO.multinomial.lifetime.lifespan$Sol)
  # lifespan2:lifespan4 a randomly high autocorrelation value
  autocorr(femaleh2EPO.multinomial.lifetime.lifespan$VCV)
}

# additional consideration: whether incidental assortative mating by promiscuity
# leads to an artificially strong maternal effect. In this case, we would expect
# a paternal effect of comparable magnitude to the maternal effect.

{
  femaleh2EPO.multinomial.lifetime.patID <- MCMCglmm(cbind(EPO, WPO)~1,
                                                   ~factoranimal + factormaternalID + 
                                                     factorpaternalID + factorcohort,
                                                   ginverse=list(factoranimal=invped.femalelifetime),
                                                   prior=prior4G.p,
                                                   data=femalephenotypes,
                                                   family="multinomial2",
                                                   nitt=1000000,
                                                   thin=800,
                                                   burnin=200000)
  plot(femaleh2EPO.multinomial.lifetime.patID)
  # nothing distinguished at all any more.
  autocorr(femaleh2EPO.multinomial.lifetime.patID$Sol)
  # some high autocorrelation.
  autocorr(femaleh2EPO.multinomial.lifetime.patID$VCV)
}

# paternal ID and lifespan

{
  femaleh2EPO.multinomial.lifetime.patID.lifespan <- MCMCglmm(cbind(EPO, WPO)~factor(lifespanPED5),
                                                     ~factoranimal + factormaternalID + 
                                                       factorpaternalID + factorcohort,
                                                     ginverse=list(factoranimal=invped.femalelifetime),
                                                     prior=prior4G.p,
                                                     data=femalephenotypes,
                                                     family="multinomial2",
                                                     nitt=1000000,
                                                     thin=800,
                                                     burnin=200000)
  plot(femaleh2EPO.multinomial.lifetime.patID.lifespan)
  # very little. Something maternal but can't be estimated.
  autocorr(femaleh2EPO.multinomial.lifetime.patID.lifespan$Sol)
  autocorr(femaleh2EPO.multinomial.lifetime.patID.lifespan$VCV)
  # ok
}


#-------------------------------------------
# Per year female behaviour as a multinomial variable
#-------------------------------------------

# with pedigree term, female ID, maternal ID, and female cohort
{
  femaleh2EPO.multi.yr <- MCMCglmm(cbind(EPO, WPO)~1,
                                               ~factoranimal + factorfemaleID +
                                                 factormaternalID + factorcohort,
                                               ginverse=list(factoranimal=invped.femalelifetime),
                                               prior=prior4G.p,
                                               data=femaleyear,
                                               family="multinomial2",
                                               nitt=1000000,
                                               thin=800,
                                               burnin=200000)
  plot(femaleh2EPO.multi.yr)
  # looks like the model wants to partition a little to animal, female, and mother
  # but it is lacking the ability to partition to any.
  autocorr(femaleh2EPO.multi.yr$Sol)
  autocorr(femaleh2EPO.multi.yr$VCV)
  # ok
}

# with age as a fixed factor between 1 and 6
{
  femaleh2EPO.multi.yr.age <- MCMCglmm(cbind(EPO, WPO)~factor(age5),
                                   ~factoranimal + factorfemaleID +
                                     factormaternalID + factorcohort,
                                   ginverse=list(factoranimal=invped.femalelifetime),
                                   prior=prior4G.p,
                                   data=femaleyear,
                                   family="multinomial2",
                                   nitt=1000000,
                                   thin=800,
                                   burnin=200000)
  plot(femaleh2EPO.multi.yr.age)
  # no real effect of female age
  # still no clear separation with the random effects.
  autocorr(femaleh2EPO.multi.yr.age$Sol)
  # a bit high
  autocorr(femaleh2EPO.multi.yr.age$VCV)
  # one also a little high
}

# with paternalID
{
  femaleh2EPO.multi.yr.pat <- MCMCglmm(cbind(EPO, WPO)~1,
                                   ~factoranimal + factorfemaleID +
                                     factormaternalID + factorpaternalID +
                                     factorcohort,
                                   ginverse=list(factoranimal=invped.femalelifetime),
                                   prior=prior5G.p,
                                   data=femaleyear,
                                   family="multinomial2",
                                   nitt=1000000,
                                   thin=800,
                                   burnin=200000)
  plot(femaleh2EPO.multi.yr.pat)
  # not really many effects here either. Maybe all these models are
  # doing too much for the available data.
  autocorr(femaleh2EPO.multi.yr.pat$Sol)
  autocorr(femaleh2EPO.multi.yr.pat$VCV)
  # fine
}

# with paternalID and age
{
  femaleh2EPO.multi.yr.pat.age <- MCMCglmm(cbind(EPO, WPO)~factor(age5),
                                       ~factoranimal + factorfemaleID +
                                         factormaternalID + factorpaternalID +
                                         factorcohort,
                                       ginverse=list(factoranimal=invped.femalelifetime),
                                       prior=prior5G.p,
                                       data=femaleyear,
                                       family="multinomial2",
                                       nitt=1000000,
                                       thin=800,
                                       burnin=200000)
  plot(femaleh2EPO.multi.yr.pat.age)
  # nothing of importance.
  autocorr(femaleh2EPO.multi.yr.pat.age$Sol)
  autocorr(femaleh2EPO.multi.yr.pat.age$VCV)
  # ok
}

# with year of breeding as a random effect, paternal ID, and age
{
  femaleh2EPO.multi.yr.pat.age.year <- MCMCglmm(cbind(EPO, WPO)~factor(age5),
                                       ~factoranimal + factorfemaleID +
                                         factormaternalID + factorpaternalID +
                                         factorcohort + factoryear,
                                       ginverse=list(factoranimal=invped.femalelifetime),
                                       prior=prior6G.p,
                                       data=femaleyear,
                                       family="multinomial2",
                                       nitt=1000000,
                                       thin=800,
                                       burnin=200000)
  plot(femaleh2EPO.multi.yr.pat.age.year)
  # nothing, long tails suggestive of other values but no clear values.
  autocorr(femaleh2EPO.multi.yr.pat.age.year$Sol)
  # 0.8 autocorrelations.
  autocorr(femaleh2EPO.multi.yr.pat.age.year$VCV)
  # ok.
}

# given that nothing has come from these models, do a very
# simple version of female ID, animal, and mother --> maximum
# possible power.

{
  femaleh2EPO.multi.yr.mincohort <- MCMCglmm(cbind(EPO, WPO)~1,
                                   ~factoranimal + factorfemaleID +
                                     factormaternalID,
                                   ginverse=list(factoranimal=invped.femalelifetime),
                                   prior=prior3G.p,
                                   data=femaleyear,
                                   family="multinomial2",
                                   nitt=1000000,
                                   thin=800,
                                   burnin=200000)
  plot(femaleh2EPO.multi.yr.mincohort)
  # hum. Maternal ID almost expressed...
  autocorr(femaleh2EPO.multi.yr.mincohort$Sol)
  autocorr(femaleh2EPO.multi.yr.mincohort$VCV)
  # fine
}

# add year of breeding to model without cohort, as a possible confound
# for the measure of EPO:WPO:
{
  femaleh2EPO.multi.yr.mincohort.year <- MCMCglmm(cbind(EPO, WPO)~1,
                                             ~factoranimal + factorfemaleID +
                                               factormaternalID + factoryear,
                                             ginverse=list(factoranimal=invped.femalelifetime),
                                             prior=prior4G.p,
                                             data=femaleyear,
                                             family="multinomial2",
                                             nitt=1000000,
                                             thin=800,
                                             burnin=200000)
  plot(femaleh2EPO.multi.yr.mincohort.year)
  # hum. More evidence of a femaleID term than before, 
  # and similar for maternalID (still nothing to write
  # about but there you go. There are no imputed dams
  # in this).
  autocorr(femaleh2EPO.multi.yr.mincohort.year$Sol)
  autocorr(femaleh2EPO.multi.yr.mincohort.year$VCV)
  # ok
}


#-------------------------------------------
# Per brood female behaviour as a poisson variable
#-------------------------------------------

# with this being per brood, can and should account for the social
# male, his indirect genetic effect on the female, his phenotypic
# effect on the female, the male's year, and the year of reproduction.
# See Reid et al 2015 for the approach.
{
  broodEPO <- MCMCglmm(cbind(EPO, WPO)~factor(dadage6),
                        random=~str(factoranimal + factordadanimal) + 
                          factorfemaleID + factorDadID + factoryear,
                        family="multinomial2",
                        prior=priorspecial4G.p,
                        ginverse=list(factoranimal=invped.bothsexes,
                                      factordadanimal=invped.bothsexes),
                        data=broodWE,
                        nitt=1000000,
                        thin=800,
                        burnin=200000)
  plot(broodEPO)
  # results intriguing: some evidence for female additive genetics
  # for the trait, though not strong. No covariance. No male additive
  # genetics. Mum and Dad ID matter. Oh, but I forgot both mother's
  # mother and father's mother, to calculate maternal effects.
  # It also seems the social male's age doesn't affect this trait.
  autocorr(broodEPO$Sol)
  autocorr(broodEPO$VCV)
  # autocorrelation good.
}

# adding the maternal effects for the male and female.
# adding pairID to account for phenotypic similarity.
# adding a maternal shared effect for males and females
# in case maternal effects have the same outcome in each sex.

broodWE$pairMaternal <- paste(broodWE$maternalID, broodWE$DadmaternalID, sep="plus")
broodWE$factorpairMaternal <- as.factor(broodWE$pairMaternal)
{
  broodEPO.maternal.pair <- MCMCglmm(cbind(EPO, WPO)~factor(dadage6),
                       random=~str(factoranimal + factordadanimal) + 
                         factormaternalID + factorDadmaternalID + 
                         factorfemaleID + factorDadID + factoryear +
                         factorpairMaternal + factorpairID,
                       family="multinomial2",
                       prior=priorspecial8G.p,
                       ginverse=list(factoranimal=invped.bothsexes,
                                     factordadanimal=invped.bothsexes),
                       data=broodWE,
                       nitt=1000000,
                       thin=800,
                       burnin=200000)
  
  plot(broodEPO.maternal.pair)
  # again no effect of age of the social father
  # no real genetic effects or maternal, though the
  # stronger effect here is the animal and mother's mother.
  # suggestion of a shared maternal effect is strange. Better
  # estimated than individual maternal effects.
}

# but a problem with this is the number of NA values for maternal
# IDs.
# one solution is to remove the NA values and try the model again:
broodWE.fullmums1 <- broodWE[-which(is.na(broodWE$maternalID)),]
broodWE.fullmums <- broodWE.fullmums1[-which(is.na(broodWE.fullmums1$DadmaternalID)),]

{
  broodEPO.maternal.pair.mums <- MCMCglmm(cbind(EPO, WPO)~factor(dadage6),
                                     random=~str(factoranimal + factordadanimal) + 
                                       factormaternalID + factorDadmaternalID + 
                                       factorfemaleID + factorDadID + factoryear +
                                       factorpairMaternal + factorpairID,
                                     family="multinomial2",
                                     prior=priorspecial8G.p,
                                     ginverse=list(factoranimal=invped.bothsexes,
                                                   factordadanimal=invped.bothsexes),
                                     data=broodWE.fullmums,
                                     nitt=1000000,
                                     thin=800,
                                     burnin=200000)
  plot(broodEPO.maternal.pair.mums)
  # not a bit of a negative effect of older father, so females more
  # faithful to oldest partners but not a lot.
  # again, quite wishy-washy for effects. How did Jane Reid get any
  # results out of this?
}

##############################################################################
# Bivariate model: genetic or phenotypic covariance between male and female
# behaviour
##############################################################################


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
# Issues
##############################################################################

# To sort in the new pedigree:
# 1   broods that have more than one genetic mother
# 2   male of age zero
# 3   birds with a different cohort between the database and the pedigree.

##############################################################################
# Close the database
##############################################################################

close(sparrowDB)
