##############################################################################
# The heritability of promiscuity in the Lundy House Sparrows
##############################################################################

# Isabel Winney
# 20160115

rm(list=ls())

##############################################################################
# loading necessary packages and functions
##############################################################################

library(MCMCglmm)
library(pedantics)
library(RODBC)
library(knitr)

# function to check for autocorrelation in MCMCglmm models
source("C:/Users/Issie/SkyDrive/RFunctions/MCMCglmm-checkAutocorr-20160317.R")

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

sparrowDB <- odbcConnectAccess('C:/Users/Issie/SkyDrive/PhD/SparrowDatabases/Database0.74_20160306-ISW/SparrowDatabase0.74.mdb')

# view tables within the database
sqlTables(sparrowDB)

# now extract data for checking the pedigree: cohort data,
# sex data.

##############################################################################
# cohort data and sex data and last stage data
##############################################################################

# extract the list of cohorts and Bird IDs

{
  sqlFetch(sparrowDB, "tblBirdID", max=10)
  sqlFetch(sparrowDB, "sys_SexEstimates", max=10)
  
  birdcohort <- sqlQuery (sparrowDB,
                          "SELECT tblBirdID.BirdID, 
                                  tblBirdID.Cohort, 
                                  tblBirdID.LastStage, 
                                  sys_SexEstimates.SexEstimate, 
                                  tblBirdID.DeathDate,
                                  tblBirdID.BroodRef
                          FROM tblBirdID
                          INNER JOIN sys_SexEstimates 
                          ON tblBirdID.BirdID = sys_SexEstimates.BirdID;",
                          na.strings="NA")
  
  head(birdcohort)
  summary(birdcohort)
  # no missing cohorts
  # no missing last stage
  
  length(birdcohort[,1])
  
}

{
  # I cannot use the LastStage in the tblBirdID because this is
  # updated with pedigree information and sightings information.
  # Therefore everything observed as a parent in the pedigree
  # will have LastStage == 3.
  # I will instead try the last stage as capture. No parent should
  # have been an egg or chigg that was never captured alive.
  laststage <- sqlQuery (sparrowDB,
                          "SELECT tblCaptures.BirdID, 
                                  Max(tblCaptures.Stage) AS MaxOfStage,
                                  Min(tblCaptures.Stage) AS MinOfStage
                          FROM tblCaptures
                          GROUP BY tblCaptures.BirdID;",
                          na.strings="NA")
  
  head(laststage)
}

{
  # and lastly to pull out the first egg date of each brood
  # so that I can check that each genetic parent was alive at
  # the time of siring the offspring.
  # Doesn't have to be biologically true, it is a guideline
  # for which individuals to check. For example, females could
  # store sperm for a couple of weeks.
  
  firstegg <- sqlQuery (sparrowDB,
                         "SELECT  tblBroods.BroodRef, 
                                  tblBroodEvents.EventNumber, 
                                  tblBroodEvents.EventDate
                          FROM tblBroods 
                          INNER JOIN tblBroodEvents 
                          ON tblBroods.BroodRef = tblBroodEvents.BroodRef
                          WHERE (((tblBroodEvents.EventNumber)=0))
                          ORDER BY tblBroods.BroodRef;",
                         na.strings="NA")
  
  head(firstegg)
}

##############################################################################
# Checking the pedigree
##############################################################################

# how many parents have a LastStage ==1 or 0?
{
  dams <- data.frame(BirdID = unique(sparrowpedigree$dam),
                     dam=1)
  head(dams)
  
  sires <- data.frame(BirdID = unique(sparrowpedigree$sire),
                      sire=1)
  head(sires)
  
  damcheck <- merge(dams, laststage, by="BirdID", 
                    all.x=TRUE, incomparables = NA)
  
  sirecheck <- merge(sires, laststage, by="BirdID", 
                     all.x=TRUE, incomparables = NA)
  
  head(damcheck)
  tail(damcheck)
  head(sirecheck)
  tail(sirecheck)
  
  # now which are =<1
  which(damcheck$LastStage<2)
  which(sirecheck$LastStage<2)
  # well that didn't work. Probably because the last stage is
  # also from the pedigree so use my other data set of maximum
  # last stage from capture:
  
  
  damcheck <- merge(dams, laststage, by="BirdID", 
                    all.x=TRUE, incomparables = NA)
  
  sirecheck <- merge(sires, laststage, by="BirdID", 
                     all.x=TRUE, incomparables = NA)
  
  head(damcheck)
  tail(damcheck)
  head(sirecheck)
  tail(sirecheck)
  
  # now which are =<1
  which(damcheck$LastStage<2)
  which(sirecheck$LastStage<2)
  # ok. All last stages are >=2
} # nothing to change


# how many pedigree cohorts and database cohorts are different?
{
  # re-name the pedigree to have the same birdID column and
  # distinct cohort column from the birdcohort data set:
  names(sparrowpedigree) <- c("BirdID", "dam", "sire", 
                              "CohortPedigree", "Immigrant")
  head(sparrowpedigree)
  
  pedcheck <- merge(sparrowpedigree, birdcohort, by="BirdID", all.x=TRUE)
  
  head(pedcheck)
  tail(pedcheck)
  
  summary(pedcheck$CohortPedigree - pedcheck$Cohort)
  table(pedcheck$CohortPedigree - pedcheck$Cohort)
  # two NAs, 58 different cohorts.
  pedcheck[which(pedcheck$CohortPedigree - pedcheck$Cohort !=0),1:6]
}
# this needs looking at:
pedcheck[which(pedcheck$CohortPedigree - pedcheck$Cohort !=0),1:6]


# are any dams male? Are any sires female?
{
  pedcheck$damsex <- birdcohort$SexEstimate[match(pedcheck$dam,
                                                  birdcohort$BirdID)]
  
  head(pedcheck)
  summary(pedcheck$damsex)
  
  pedcheck$siresex <- birdcohort$SexEstimate[match(pedcheck$sire,
                                                   birdcohort$BirdID)]
  
  summary(pedcheck$siresex)
  # oops.
  which(pedcheck$siresex==0)
  # only one!
  pedcheck[which(pedcheck$siresex==0),]
  # so this one male is actually female
}


# lastly, are there any zombie parents?
{
  # fill in the hatch date for the individuals:
  pedcheck.hatched <- merge(pedcheck, firstegg, by="BroodRef", 
                            all.x=TRUE, incomparables = NA)
  
  head(pedcheck.hatched)
  tail(pedcheck.hatched)
  table(pedcheck.hatched$BroodRef)
  
  # and add the death date of the parents:
  pedcheck.hatched$damdeath <- birdcohort$DeathDate[match(pedcheck.hatched$dam,
                                                          birdcohort$BirdID)]
  
  pedcheck.hatched$siredeath <- birdcohort$DeathDate[match(pedcheck.hatched$sire,
                                                           birdcohort$BirdID)]
  
  head(pedcheck.hatched)
  pedcheck.hatched$damdeath
  pedcheck.hatched$siredeath
  
  # all dates of death should be after broods were started.
  # That means all these numbers should be negative:
  table(pedcheck.hatched$EventDate - pedcheck.hatched$damdeath)
  table(pedcheck.hatched$EventDate - pedcheck.hatched$siredeath)
  
  # but there are some positive dates for males, i.e. males that
  # died before they became fathers:
  pedcheck.hatched$siredeathdif <- pedcheck.hatched$EventDate - pedcheck.hatched$siredeath
  which(pedcheck.hatched$siredeathdif >0)
  pedcheck.hatched[which(pedcheck.hatched$siredeathdif >0),]
}
# most of these are 1192, which is reassuring.
# The others, I can believe 3394 was a father 2 days later.
# 35 days is a stretch.
# so double check the fathers for 3797, 3926, 5426, 5429, 5634.


{
  # how many nestlings in the pedigree do not have assigned parents:
  sparrowpedigree$laststage <- laststage$MaxOfStage[match(sparrowpedigree$BirdID,
                                                          laststage$BirdID)]
  
  sparrowpedigree$firststage <- laststage$MinOfStage[match(sparrowpedigree$BirdID,
                                                          laststage$BirdID)]
  
  
  nestlings <- subset(sparrowpedigree, sparrowpedigree$firststage<3)
  
  summary(nestlings)
  table(nestlings$laststage)
  
  # how many eggs have one or both parents missing:
  length(which(paste(nestlings$laststage, nestlings$dam)=="0 NA"))
  length(which(paste(nestlings$laststage, nestlings$sire)=="0 NA"))
  length(which(paste(nestlings$laststage, nestlings$dam, nestlings$sire)=="0 NA NA"))
  length(nestlings[,1])
  length(which(is.na(nestlings$dam)))
  length(which(is.na(nestlings$sire)))
}

##############################################################################
# offspring data
##############################################################################

# Assumptions:
# 1   males are all dead --> the calculated phenotypes
#     are complete. This is not true but can be reconsidered later.
# 2   the 2015 pedigree is good.
# 3   a brood has just one genetic mother. By looking
#     up the genetic mother of a brood, the mother that is
#     returned first is the one used in the analysis below. 
# 4   I can only make phenotypes using individuals that have been born
#     in a brood, otherwise I do not know their social and genetic parents
#     i.e. for a mother, I can only calculate her EPO and WPO from EPO/WPO
#     born in broods
# 5   All cohorts in the database are accurate (not true)


# In the past I used text files of data downloaded from the database,
# now I am calling the database directly.

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


# In 6468 cases out of 7628 I know both. That presumably includes
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
                                                sparrowped$id)]

offspring$GeneticMumID <- sparrowped$dam[match(offspring$BirdID,
                                               sparrowped$id)]

head(offspring)
tail(offspring)


summary(offspring)
# fantastic

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

  table(table(offspring3$GeneticDadID))
  # there are 48 sires that appear only once.
}

# remove cases where genetic and social mother do not match
{
  which(offspring3$SocialMumID!=offspring3$GeneticMumID)
  offspring3[which(offspring3$SocialMumID!=offspring3$GeneticMumID),]
  
  # the majority are single chicks in random broods, there is only one
  # case of a mother where the whole brood is assigned to a sister.
  
  # remove
  offspring4 <- offspring3[-which(offspring3$SocialMumID!=offspring3$GeneticMumID),]
  
  # check it makes sense:
  which(offspring4$SocialMumID!=offspring4$GeneticMumID)
  summary(offspring4)
}
#-------------------------------------------
# Designate offspring as EPO or WPO
#-------------------------------------------
{
# WPO have the same social and genetic sire,
# EPO have different genetic and social sires.

offspring4$WPO <- ifelse(offspring4$SocialDadID==offspring4$GeneticDadID,
                         1, 0)
head(offspring4)
tail(offspring4)

offspring4$EPO <- ifelse(offspring4$SocialDadID!=offspring4$GeneticDadID,
                         1, 0)

table(offspring4$EPO, offspring4$WPO)

# 1041 EPO. 4837 WPO. None identified as both. Good!
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
  malephenotypes <- aggregate(offspring4$WPO, 
                              list(offspring4$GeneticDadID),
                              FUN=sum)
  head(malephenotypes)
  
  # now the number of EPO per male
  maleEPO <- aggregate(offspring4$EPO, 
                       list(offspring4$GeneticDadID),
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
  fixedsparrowped[which(fixedsparrowped$animal==971),]
  malephenotypes[which(malephenotypes$animal==971),]
  
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
  sparrowpedigree[which(sparrowpedigree$id==4722),]
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
  malelifespan <- aggregate(sparrowpedigree$Cohort, 
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
  # great!
  
  # make spare age variable:
  malephenotypes$lifespanPED6 <- malephenotypes$lifespanPED
  
  # assign >6 years old to 6 years old for sample size by finding the
  # >6 data points and replacing them with a 6:
  malephenotypes$lifespanPED6[which(malephenotypes$lifespanPED6>6)] <- 6
  
  
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
  offspring4$maleyear <- paste(offspring4$GeneticDadID, offspring4$Cohort, sep=".")
  
  head(offspring4)
  table(offspring4$maleyear)
  table(table(offspring4$maleyear))
  # up to 25 offspring!
  
  
  # aggregate number of WPO per male
  maleyearWPO <- aggregate(offspring4$WPO, 
                           list(offspring4$maleyear),
                           FUN=sum)
  head(maleyearWPO)
  
  # aggregate number of EPO per male
  maleyearEPO <- aggregate(offspring4$EPO, 
                           list(offspring4$maleyear),
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
  maleyear$animal <- offspring4$GeneticDadID[match(maleyear$maleyear, offspring4$maleyear)]
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
  sparrowpedigree[which(sparrowpedigree$id==4722),]
  maleyear[which(maleyear$animal==4722),]
  # good good.
}

{
  # male age. Unlike with lifespan, age is easier to define as the
  # year breeding occurred minus the cohort of the male.
  
  # add the year breeding occurred to the data set:
  maleyear$year <- offspring4$Cohort[match(maleyear$maleyear, 
                                           offspring4$maleyear)]
  
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


{
  # Import foster brood information:
  summary(sqlFetch(sparrowDB, "tblFosterBroods"))
  
  # nine missing foster broods. Four missing bird IDs.
  # These would need to be removed, but this will be adequate as cross-foster data
  
  BirdFosterBrood <- sqlQuery(sparrowDB,
                              "SELECT   tblFosterBroods.BirdID, 
                              tblFosterBroods.FosterBrood, 
                              tblBroods.BroodRef, 
                              tblBroods.SocialDadID, 
                              tblBroods.SocialMumID
                              FROM tblFosterBroods 
                              INNER JOIN tblBroods ON tblFosterBroods.FosterBrood = tblBroods.BroodRef;",
                              na.strings="NA")
  
  summary(BirdFosterBrood)
  # four entries with no bird ID
  
  BirdFosterBrood <- BirdFosterBrood[-which(is.na(BirdFosterBrood$BirdID)),]
  
  summary(BirdFosterBrood)
  
  # There we go.
}

{
  # add the foster mother ID to the male phenotypes. Call
  # it social mother because later I will add the genetic 
  # mother to the ones that were not cross-fostered:
  
  maleyear$SocialMumID <- BirdFosterBrood$SocialMumID[match(maleyear$maleid,
                                                            BirdFosterBrood$BirdID)]
  
  summary(maleyear$SocialMumID)
  # so 466 out of 703 were not cross-fostered. Might make disentangling these
  # two mothers challenging.
  
  # replace the 570 NA social mothers with the genetic mother:
  
  for(i in 1:length(maleyear$SocialMumID)){
    if(is.na(maleyear$SocialMumID[i])){
      maleyear$SocialMumID[i] <- maleyear$maternalID[i]
    } else {
      maleyear$SocialMumID[i] <- maleyear$SocialMumID[i]
    }
  }
  
  summary(maleyear$SocialMumID)
  
  
  length(which(maleyear$SocialMumID==maleyear$maternalID))
  
  summary(maleyear$maternalID)
  
  # so the missing mothers are missing from the maternalID. I guess
  # they are founder males or so.
  
  # how do the maternal sibships differ with the cross-foster data?
  table(table(maleyear$maternalID))
  table(table(maleyear$SocialMumID))
  # some changes.
}

{
  # set as factor
  maleyear$factorSocialMumID <- as.factor(maleyear$SocialMumID)
}


##############################################################################
# female phenotypes 
##############################################################################

# Aggregate WPO and EPO by GeneticMumID

{
  # for offspring4, are there any cases where the genetic mother is not known?
  which(is.na(offspring4$GeneticMumID))
  # yes, which means these individuals need to be removed from making the 
  # dataset for females:
  offspring5 <- offspring4[-which(is.na(offspring4$GeneticMumID)),]
  
  which(is.na(offspring5$GeneticMumID))
  summary(offspring5)
  
  # aggregate by GeneticMumID
  femalephenotypes <- aggregate(offspring5$WPO, 
                                list(offspring5$GeneticMumID),
                                FUN=sum)
  head(femalephenotypes)
  
  
  femaleEPO <- aggregate(offspring5$EPO, 
                         list(offspring5$GeneticMumID),
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
  # I am consciously choosing the database cohort and not
  # the pedigree cohort, since I know there are discrepancies.
  femalephenotypes$Cohort <- birdcohort$Cohort[match(femalephenotypes$animal,
                                                   birdcohort$BirdID)]
  
  head(femalephenotypes)
  summary(femalephenotypes)
  
  # check:
  sparrowpedigree[which(sparrowpedigree$id==7022),]
  fixedsparrowped[which(fixedsparrowped$animal==7022),]
  femalephenotypes[which(femalephenotypes$animal==7022),]
  # that is interesting. This bird is in the fixed pedigree
  # only and not the original pedigree. And, this bird used
  # to have parents but does not any more.
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
  femalelifespan <- aggregate(sparrowpedigree$Cohort, 
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
length(offspring5$GeneticMumID)
table(offspring5$SocialMumID==offspring5$GeneticMumID)
# there should not be since I removed them earlier now.

#-------------------------------------------
# Female dataset of genetic offspring per year
#-------------------------------------------


{
  # per year per GENETIC female identifier:
  offspring5$femaleyear <- paste(offspring5$GeneticMumID, offspring5$Cohort, sep=".")
  
  head(offspring5)
  table(offspring5$femaleyear)
  table(table(offspring5$femaleyear))
  # up to 22 offspring this time.
  
  
  # aggregate number of WPO per female
  femaleyearWPO <- aggregate(offspring5$WPO, 
                           list(offspring5$femaleyear),
                           FUN=sum)
  head(femaleyearWPO)
  
  # aggregate number of EPO per female
  femaleyearEPO <- aggregate(offspring5$EPO, 
                           list(offspring5$femaleyear),
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
  femaleyear$animal <- offspring5$GeneticMumID[match(femaleyear$femaleyear, offspring5$femaleyear)]
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
  sparrowpedigree[which(sparrowpedigree$id==76),]
  femaleyear[which(femaleyear$animal==76),]
  # good good.
}

{
  # female age. Unlike with lifespan, age is easier to define as the
  # year breeding occurred minus the cohort of the female.
  
  # add the year breeding occurred to the data set:
  femaleyear$year <- offspring5$Cohort[match(femaleyear$femaleyear, 
                                             offspring5$femaleyear)]
  
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
  
  broodWPO <- aggregate(offspring5$WPO,
                        list(offspring5$BroodName),
                        FUN=sum)
  
  broodEPO <- aggregate(offspring5$EPO,
                        list(offspring5$BroodName),
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
  # 74
  
  # add female ID
  broodWE$MumID <- offspring5$GeneticMumID[match(broodWE$BroodName,
                                                        offspring5$BroodName)]
  
  length(unique(broodWE$MumID))
  length(unique(offspring5$GeneticMumID))
  # there are more mothers in the offspring data set than go in
  # the per brood data set. This is because there are a small
  # number of broods that have more than one genetic mother. For
  # example:
  offspring5[which(offspring5$BroodName=="A008"),]
  broodWE[which(broodWE$BroodName=="A008"),]
  
  # add the female's social partner.
  broodWE$DadID <- offspring5$SocialDadID[match(broodWE$BroodName,
                                                        offspring5$BroodName)]
  
  # make a pair ID from the male and female's ID
  broodWE$pairID <- paste(broodWE$MumID, broodWE$DadID, sep="plus")
  
  # year of reproduction
  broodWE$year <- offspring5$Cohort[match(broodWE$BroodName,
                                          offspring5$BroodName)]
  
  head(broodWE)
  
  
  # check that this is all matching the original data:
  broodWE[which(broodWE$BroodName=="J050"),]
  offspring5[which(offspring5$BroodName=="J050"),]
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
  sparrowpedigree[which(sparrowpedigree$id==4037),]
  sparrowpedigree[which(sparrowpedigree$id==4753),]
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
  
  
  table(broodWE$dadage6)
}

{
  # adding the social pair's shared maternal effect
  broodWE$pairMaternal <- paste(broodWE$maternalID, broodWE$DadmaternalID, sep="plus")
  broodWE$factorpairMaternal <- as.factor(broodWE$pairMaternal)
}

{
  # exlcuding the missing maternal IDs might be important in this data set:
  broodWE.fullmums1 <- broodWE[-which(is.na(broodWE$maternalID)),]
  broodWE.fullmums <- broodWE.fullmums1[-which(is.na(broodWE.fullmums1$DadmaternalID)),]
  
  summary(broodWE.fullmums)
  
}

##############################################################################
# combined phenotypes 
##############################################################################

# For hypotheses of the evolution of promiscuity to work, there needs to be
# shared genetic covariance between males and females. Therefore, I need a
# bivariate data frame to account for the essential random and fixed effects
# that impact each behaviour.
# I will use male EPO per year and female proportion EPO per year so that the
# model is not biased by the lifetime of the birds and so that there is more
# variation in the female data to help the model fit.

# from both data sets, I need bird ID, animal, maternal ID, and year effects.
# from the male data set, I need male age.

{
  # Load female data:
  
  head(femaleyear)

  f1 <- data.frame(animal=femaleyear$femaleid,
                       birdID=femaleyear$femaleid,
                       EPOf=femaleyear$EPO,
                       EPOm=NA,
                       WPOf=femaleyear$WPO,
                       WPOm=NA,
                       year=femaleyear$year,
                       maternalID=femaleyear$maternalID,
                       maleage=NA,
                       sex="F")

  head(f1)
  summary(f1)
  # some missing mothers. One case of 9 EPO
  f1[which(f1$EPOf==9),]
  # how strange. I suspect the father is misidentified.
  offspring[which(offspring$SocialMumID==4972),]
  # This is weird. There is only one social father. But he only
  # sires one offspring. I would have guessed an assignment error but
  # since he did sire one offspring I guess it is a very promiscuous
  # female.
  
  str(f1)
}

{
  # Load male data:
  
  head(maleyear)
  
  m1 <- data.frame(animal=maleyear$maleid,
                         birdID=maleyear$maleid,
                         EPOf=NA,
                         EPOm=maleyear$EPO,
                         WPOf=NA,
                         WPOm=maleyear$WPO,
                         year=maleyear$year,
                         maternalID=maleyear$maternalID,
                         maleage=maleyear$age6,
                         sex="M")
  
  head(m1)
  summary(m1)
  # 100 missing dams.
  str(m1)
}

{
  # combine the data frames:
  bothsexesyear <- rbind(f1, m1)
  
  head(bothsexesyear)
  tail(bothsexesyear)
  str(bothsexesyear)
  summary(bothsexesyear)
}

{
  # Was either of the male or female extra-pair? i.e. are extra-pair
  # offspring behaviourally different?
  
  head(bothsexesyear)
  head(offspring4)
  
  bothsexesyear$EPOstatus <- offspring4$EPO[match(bothsexesyear$birdID,
                                                  offspring4$BirdID)]
  
  summary(bothsexesyear$EPOstatus)
  # quite a few missing...
  
  # data set without these individuals...
  
  bothsexes.EPOstatus <- bothsexesyear[-which(is.na(bothsexesyear$EPOstatus)),]
  summary(bothsexes.EPOstatus)
}

##############################################################################
# Phenotype per offspring
##############################################################################

{
  # Add maternal ID of genetic mother and father
  
  offspring5$GenDadMaternalID <- fixedsparrowped$dam[match(offspring5$GeneticDadID,
                                                           fixedsparrowped$animal)]
  
  offspring5$GenMumMaternalID <- fixedsparrowped$dam[match(offspring5$GeneticMumID,
                                                           fixedsparrowped$animal)]
}

{
  # add male age. First, male cohort:
  
  offspring5$GenDadCohort <- birdcohort$Cohort[match(offspring5$GeneticDadID,
                                                     birdcohort$BirdID)]
  
  head(offspring5)
  
  # age by subtracting Dad's cohort from offspring age:
  offspring5$GenDadAge <- offspring5$Cohort - offspring5$GenDadCohort
  
  table(offspring5$GenDadAge)
  
  # correct father with age 0 to be age 1, and amalgamate age >=6
  
  offspring5$GenDadAge6 <- offspring5$GenDadAge
  offspring5$GenDadAge6[which(offspring5$GenDadAge6==0)] <- 1
  offspring5$GenDadAge6[which(offspring5$GenDadAge6>6)] <- 6
  
  table(offspring5$GenDadAge6)
}

{
  # set factors

  offspring5$factoryear <- as.factor(offspring5$Cohort)
  offspring5$factorGenDadID <- as.factor(offspring5$GeneticDadID)
  offspring5$factorGenDadanimal <- as.factor(offspring5$GeneticDadID)
  offspring5$factorGenMumID <- as.factor(offspring5$GeneticMumID)
  offspring5$factorGenMumanimal <- as.factor(offspring5$GeneticMumID)
  offspring5$factorGenDadMaternalID <- as.factor(offspring5$GenDadMaternalID)
  offspring5$factorGenMumMaternalID <- as.factor(offspring5$GenMumMaternalID)
  offspring5$factorGenDadAge6 <- as.factor(offspring5$GenDadAge6)
  
  summary(offspring5)
  str(offspring5)
}


##############################################################################
# Pruning the pedigree
##############################################################################

{
  # pedigree info
  drawPedigree(fixedsparrowped)
  pedigreeinformation <- pedigreeStats(fixedsparrowped)
  #
  #
  #
  #
  #
  #
  #
  #
  pedStatSummary(pedigreeinformation)
  
}

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
  #
  #
  #
  #
  #
  #
  #
  #
  pedStatSummary(pedigreeinformation)
  
  # for females:
  prunedped.femalelifetime <- prunePed(pedigree=fixedsparrowped, 
                                     keep=femalephenotypes$animal)
  
  head(prunedped.femalelifetime)
  
  
  # for both males and females ***for the broodWE data set
  # i.e. parents differ because 1 - in a couple of cases there
  # is more than one female per brood 2 - the maleID is the social
  # male of the brood***
  keepers <- c(unique(broodWE$MumID), unique(broodWE$DadID))
  prunedped.bothsexes <- prunePed(pedigree=fixedsparrowped, 
                                  keep=keepers)
  head(prunedped.bothsexes)
  
  
  # and for both males and females in the broodWE data set with
  # missing maternal IDs for males and females pruned out:
  keepers <- c(unique(broodWE.fullmums$MumID), unique(broodWE.fullmums$DadID))
  prunedped.bothsexes.fullmums <- prunePed(pedigree=fixedsparrowped, 
                                  keep=keepers)
  head(prunedped.bothsexes.fullmums)
  
  
  
  # for both males and females in the combined data set of per year
  # behaviour:
  keepers <- unique(bothsexesyear$birdID)
  prunedped.bothsexesyear <- prunePed(pedigree=fixedsparrowped, 
                                      keep=keepers)
  
  # pedigree info for both sexes:
  pedigreeinformation <- pedigreeStats(prunedped.bothsexesyear)
  #
  #
  #
  #
  #
  #
  #
  #
  pedStatSummary(pedigreeinformation)
}

#
#
#
#
#

# and then the pedigrees can be converted to relationship matrices
# for the current version of MCMCglmm:
{
  invped.malelifetime <- inverseA(prunedped.malelifetime)$Ainv
  summary(invped.malelifetime)
  
  invped.femalelifetime <- inverseA(prunedped.femalelifetime)$Ainv
  summary(invped.femalelifetime)
  
  invped.bothsexes <- inverseA(prunedped.bothsexes)$Ainv
  summary(invped.bothsexes)
  
  invped.bothsexes.fullmums <- inverseA(prunedped.bothsexes.fullmums)$Ainv
  summary(invped.bothsexes.fullmums)
  
  invped.bothsexesyear <- inverseA(prunedped.bothsexesyear)$Ainv
  summary(invped.bothsexesyear)
}

##############################################################################
# Plots of variables
##############################################################################
{
# how many EPO offspring do males have?
hist(malephenotypes$EPO, breaks=seq(-1,35,1))
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


# And when done per brood, are there many more zeros and ones?
hist(broodWE$EPO/(broodWE$EPO+broodWE$WPO))
# many more zeros with the rest of the data more evenly spread.

table(broodWE$EPO/(broodWE$EPO+broodWE$WPO))
# 999 zero values. How ironic.
}

{
  # in the analysis of propensity for an individual offspring to
  # be WPO or EPO, I need to know whether there would be any bias
  # inherent to correlating the male and female phenotypes, for
  # example since social pairs might pair for multiple broods and
  # have numerous WPO together, this could create a correlation
  # between phenotypes due to reproductive constraints rather than
  # assortative mating.
  
  # So, take a EPO/(total):
  # THIS IS DONE OVER THE LIFETIME, this is because although the 
  # phenotype becomes more accurate with more samples, this is the
  # level at which the intended model would characterise an individual's
  # phenotype and thus is what I need to compare between genetic pairs.
  femalephenotypes$EPO.WPO <- femalephenotypes$EPO/
    (femalephenotypes$EPO+femalephenotypes$WPO)
  femalephenotypes$EPO.WPO
  
  hist(femalephenotypes$EPO.WPO)
  
  malephenotypes$EPO.WPO <- malephenotypes$EPO/
    (malephenotypes$EPO+malephenotypes$WPO)
  malephenotypes$EPO.WPO
  
  hist(malephenotypes$EPO.WPO)
  # many suspicious cases of males with all EPO.
  
  # add to the offspring data frame:
  summary(offspring5)
  
  offspring5$femalephenotype <- femalephenotypes$EPO.WPO[match(offspring5$GeneticMumID,
                                                               femalephenotypes$femaleID)]
  
  offspring5$malephenotype <- malephenotypes$EPO.WPO[match(offspring5$GeneticDadID,
                                                               malephenotypes$maleID)]
  
  head(offspring5)
  summary(offspring5)
  
  # now plot phenotypes:
  
  plot(offspring5$femalephenotype, offspring5$malephenotype)
  cor(offspring5$femalephenotype, offspring5$malephenotype)
  
  # that actually looks good... This should work!
  table(offspring5$femalephenotype, offspring5$malephenotype)[1,]
  # there are 194 zero-zero values. Out of over 3k, that is
  # good.
}

# how many EPP per year?
{
  barplot(table(offspring5$EPO, offspring5$Cohort)[2,]/
    (table(offspring5$EPO, offspring5$Cohort)[2,]+table(offspring5$EPO, offspring5$Cohort)[1,]),
    ylim=c(0,0.3))
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
  
  
  prior7G.p = list(R=list(V=1,nu=0.002),
                   G=list(G1=list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000),
                          G2=list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000),
                          G3=list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000),
                          G4=list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000),
                          G5=list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000),
                          G6=list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000),
                          G7=list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000)))
  
  prior9G.p = list(R=list(V=1,nu=0.002),
                   G=list(G1=list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000),
                          G2=list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000),
                          G3=list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000),
                          G4=list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000),
                          G5=list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000),
                          G6=list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000),
                          G7=list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000),
                          G8=list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000),
                          G9=list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000)))
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

  # but potentially the maternal effects could also only manifest in
  # older males, so put a six-level factor on the second variance component
  # as well
  prior4G.p.penetranceon2 = list(R=list(V=1,nu=0.002),
                              G=list(G1=list(V=diag(6), nu=6, alpha.mu=c(0,0,0,0,0,0), alpha.V=diag(6)*1000),
                                     G2=list(V=diag(6), nu=6, alpha.mu=c(0,0,0,0,0,0), alpha.V=diag(6)*1000),
                                     G3=list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000),
                                     G4=list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000)))
  
}


# slightly informative inverse-Wishart priors, with nu=0.02:
# .strong used to denote the slightly informative prior
{
  prior4G.strong = list(R = list(V = 1,nu = 0.02),
                 G = list(G1 = list(V = 1, nu = 0.02),
                          G2 = list(V = 1, nu = 0.02),
                          G3 = list(V = 1, nu = 0.02),
                          G4 = list(V = 1, nu = 0.02)))
}

# binomial priors. These have residual variance fixed at 1. All
# other random effects are parameter expanded:
{
  priorbinomial4 = list(R = list(V = 1, fix = 1),
                           G = list(G1 = list(V = diag(2), nu = 2, alpha.mu = c(0,0), alpha.V = diag(2)*1000),
                                  G2 = list(V = diag(2), nu = 2, alpha.mu = c(0,0), alpha.V = diag(2)*1000),
                                  G3 = list(V = diag(2), nu = 2, alpha.mu = c(0,0), alpha.V = diag(2)*1000),
                                  G4 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000)))

  priorbinomial7 = list(R = list(V = 1, fix = 1),
                        G = list(G1 = list(V = diag(2), nu = 2, alpha.mu = c(0,0), alpha.V = diag(2)*1000),
                                 G2 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000),
                                 G3 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000),
                                 G4 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000),
                                 G5 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000),
                                 G6 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000),
                                 G7 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000)))
  
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
  # low repeatability for the male multinomial behaviour - poorly
  # estimated maleID term. But more strongly estimated than in 2013
  # data set.
  # This behaviour isn't well described by a multinomial
  # distribution in males anyway, since total reproduction varies between
  # males.
  autocorr(maler2EPOmulti$Sol)
  autocorr(maler2EPOmulti$VCV)
  # Good. But now I have written a function, I can check this:
  checkAutocorr(maler2EPOmulti, 0.1)
  checkAutocorr(maler2EPOmulti, 0.05)
  # Great!
  
  summary(maler2EPOmulti)
  # good sampling.
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
  # This is repeatable. Better estimate than with 2013 data
  # set (no overlap with zero). Might improve if age
  # is added since EP behaviour changes with age. Could also
  # be contingent on the year of measurement of behaviour,
  # so a model with maleID, year, and a fixed factor of male
  # age would be good.
  checkAutocorr(maler2EPOpois, 0.1) # excellent
  checkAutocorr(maler2EPOpois, 0.05) # one
  
  summary(maler2EPOpois)
  # Male ID term distinct from zero.
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
  # older males have more EPO, which we knew already.
  # the repeatability per male is more clearly defined in this
  # case, now that the age of the male and the between-year
  # variation in the number of EPO is accounted for.
  # year has a distinct peak, though overlaps zero
  checkAutocorr(maler2EPOpois.age.yr, 0.1) # excellent
  checkAutocorr(maler2EPOpois.age.yr, 0.05)
  # some moderate Sol autocorrelation due to factor(age)
  
  summary(maler2EPOpois.age.yr)
  # a dip in sampling for two year old males.
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
  checkAutocorr(femaler2EPOmulti, 0.1)
  checkAutocorr(femaler2EPOmulti, 0.05)
  # great stuff!
  
  summary(femaler2EPOmulti)
}

# as with males, include extra factors, but this time to see
# whether repeatability remains or is a result of these 
# effects:
# female multinomial with female age as fixed factor, and year as random:
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
  # Unlike in males, year is not contributing to this. It
  # suggests this is a female trait, whereas the male trait
  # requires facilitation.
  checkAutocorr(femaler2EPOmulti.age.yr, 0.1)
  checkAutocorr(femaler2EPOmulti.age.yr, 0.05)
  # fine
  
  summary(femaler2EPOmulti.age.yr)
  # Though actually it looks like there is a dip in EPO 
  # production in middle-aged females. Why? First thoughts
  # are sampling bias: do we have more unhatched eggs for 
  # some ages and more experienced females hatch more eggs
  # so we have good sampling?
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
  # female repeatable here, no clear effect of year of reproduction.
  # Really, age has as much reason to stay as year. i.e. none. Or maybe
  # less since there is a slight dip with age.
  checkAutocorr(femaler2EPOmulti.yr, 0.1)
  checkAutocorr(femaler2EPOmulti.yr, 0.05)
  # fine.
  
  summary(femaler2EPOmulti.yr)
  # very similar variance components to the previous model.
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
  checkAutocorr(femaler2EPObrood,0.1)
  # Oop. One high value
  
  autocorr(femaler2EPObrood$VCV)
  # Hum. All the other values are pretty low for that variable. 
  # I will see whether the autocorrelation is low in the following
  # models using the same data set.
  
  summary(femaler2EPObrood)
  # sampling is good.
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
  checkAutocorr(femaler2EPObrood.yr, 0.1)
  # No autocorrelation here in the more complex model, so I
  # guess the previous model has one high autocorrelation by
  # chance.
  
  checkAutocorr(femaler2EPObrood.yr, 0.05)
  
  summary(femaler2EPObrood.yr)
  # some loss of sampling for female ID. Potentially, this data
  # is too zero-inflated for the model to handle well, and this
  # inflates autocorrelation and reduces sampling.
}

# Since the year of breeding does not matter in females, but there
# is good female repeatability, is it female quality that matters 
# and would cohort be relevant?

{
  femaler2EPObrood.cohort <- MCMCglmm(cbind(EPO,WPO)~1,
                                  ~factorfemaleID + factorcohort,
                                  prior=prior2G.p,
                                  data=broodWE,
                                  family="multinomial2",
                                  nitt=1000000,
                                  thin=800,
                                  burnin=200000)
  
  plot(femaler2EPObrood.cohort)
  
  # very similar to year, cohort is not explaining variance in this
  # behaviour.
  
  checkAutocorr(femaler2EPObrood.cohort, 0.1)
  # this is fine
  
  summary(femaler2EPObrood.cohort)
  # and good sampling.
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
  # the units in the model are not properly estimated.
  # the model would have to be run longer.
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
  # here, units estimated but no maternal ID. Maybe either/or.
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
  # and again probems estimating units.
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
  # Here and onwards, the sample size isn't sufficient for this
  # type of analysis.
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
  autocorr(maleh2EPO.pois.byage$Sol)
  # age3:age2 0.11. Other high autocorrelations suggest this model
  # would need much longer.
  autocorr(maleh2EPO.pois.byage$VCV)
  # and the same with high autocorrelation here.
}

##########################################################################
##########################################################################
# At this point, I needed to install a new R version to do some
# other work.
# Here is the new session info that applies to models after this point:
print(sessionInfo())

R version 3.2.4 (2016-03-10)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows >= 8 x64 (build 9200)

locale:
  [1] LC_COLLATE=English_United Kingdom.1252  LC_CTYPE=English_United Kingdom.1252   
[3] LC_MONETARY=English_United Kingdom.1252 LC_NUMERIC=C                           
[5] LC_TIME=English_United Kingdom.1252    

attached base packages:
  [1] grid      stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
  [1] knitr_1.12.3     RODBC_1.3-12     pedantics_1.5    MasterBayes_2.52 kinship2_1.6.4  
[6] quadprog_1.5-5   genetics_1.3.8.1 mvtnorm_1.0-5    MASS_7.3-45      gtools_3.5.0    
[11] gdata_2.17.0     combinat_0.0-8   MCMCglmm_2.22.1  ape_3.4          coda_0.18-1     
[16] Matrix_1.2-4    

loaded via a namespace (and not attached):
  [1] lattice_0.20-33 corpcor_1.6.8   nlme_3.1-125    cubature_1.1-2  tools_3.2.4    
[6] tensorA_0.36 

##########################################################################
##########################################################################


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
  plot(maleh2EPO.pois.byage.idh)
  # again, suggestion age 2 has most power to detect animal term 
  # but not clearly detected. Male ID term does ok, though not
  # clearly estimated either. Maternal ID is the most clearly
  # estimated.
  autocorr(maleh2EPO.pois.byage.idh$Sol)
  autocorr(maleh2EPO.pois.byage.idh$VCV)
  # age1:animal under age4 a little high. Otherwise fine.
}

{
  # but there could also be penetrance of maternal effects.
  # And of BirdID, but there is no repeated measure within
  # years, so it is an across-year repeatability:
  
  maleh2EPO.pois.byageandmum.idh <- MCMCglmm(EPO~factor(age6),
                                       ~idh(factor(age6)):factoranimal + 
                                         idh(factor(age6)):factormaternalID + 
                                         factormaleID + factoryear,
                                       ginverse=list(factoranimal=invped.malelifetime),
                                       prior=prior4G.p.penetranceon2,
                                       data=maleyear,
                                       family="poisson",
                                       nitt=1000000,
                                       thin=800,
                                       burnin=200000)
  plot(maleh2EPO.pois.byageandmum.idh)
}


#-------------------------------------------
# Per year male poisson - is the maternal effect social or genetic?
#-------------------------------------------


{
  # analysis to partition maternal effects in to those from the social
  # versus the genetic mother:
  
  maleh2EPO.pois.twomums <- MCMCglmm(EPO~factor(age6),
                                     ~factoranimal + factormaleID +
                                       factormaternalID + factorSocialMumID +
                                       factoryear,
                                     ginverse=list(factoranimal=invped.malelifetime),
                                     prior=prior5G.p,
                                     data=maleyear,
                                     family="poisson",
                                     nitt=1000000,
                                     thin=800,
                                     burnin=200000)
  plot(maleh2EPO.pois.twomums)
}


{
  
  
  # Joel suggested tracing the matriline back to find the founding
  # female that represents a male or female and their phenotype.
  
  # one way to do this is with a hash table for the males and to
  # list the mothers of the mothers ad infinitum until they are
  # all NA. Then collapse the hash table down to the first non-
  # NA value:
  
  # create the mothers of mothers table:
  
  head(malephenotypes)
  motherhash <- data.frame(BirdID = malephenotypes$maleID,
                           Mum1 = malephenotypes$maternalID)
  
  summary(motherhash)
  # 47 NAs. Plenty to work on:
  
  motherhash$Mum2 <- fixedsparrowped$dam[match(motherhash$Mum1,
                                               fixedsparrowped$animal)]
  
  summary(motherhash)
  # 104
  
  
  motherhash$Mum3 <- fixedsparrowped$dam[match(motherhash$Mum2,
                                               fixedsparrowped$animal)]
  
  summary(motherhash)
  # 154
  
  
  
  motherhash$Mum4 <- fixedsparrowped$dam[match(motherhash$Mum3,
                                               fixedsparrowped$animal)]
  
  summary(motherhash)
  # 198
  
  
  
  motherhash$Mum5 <- fixedsparrowped$dam[match(motherhash$Mum4,
                                               fixedsparrowped$animal)]
  
  summary(motherhash)
  # 250
  
  
  
  motherhash$Mum6 <- fixedsparrowped$dam[match(motherhash$Mum5,
                                               fixedsparrowped$animal)]
  
  summary(motherhash)
  # 278
  
  
  motherhash$Mum7 <- fixedsparrowped$dam[match(motherhash$Mum6,
                                               fixedsparrowped$animal)]
  
  summary(motherhash)
  # 295
  
  
  motherhash$Mum8 <- fixedsparrowped$dam[match(motherhash$Mum7,
                                               fixedsparrowped$animal)]
  
  summary(motherhash)
  # 321
  
  
  
  motherhash$Mum9 <- fixedsparrowped$dam[match(motherhash$Mum8,
                                               fixedsparrowped$animal)]
  
  summary(motherhash)
  # 349
  
  
  motherhash$Mum10 <- fixedsparrowped$dam[match(motherhash$Mum9,
                                                fixedsparrowped$animal)]
  
  
  summary(motherhash)
  # 385
  
  motherhash$Mum11 <- fixedsparrowped$dam[match(motherhash$Mum10,
                                                fixedsparrowped$animal)]
  
  
  summary(motherhash)
  # just a couple more! 411 NAs!
  
  motherhash$Mum12 <- fixedsparrowped$dam[match(motherhash$Mum11,
                                                fixedsparrowped$animal)]
  
  
  summary(motherhash)
  # last one!
  
  motherhash$Mum13 <- fixedsparrowped$dam[match(motherhash$Mum12,
                                                fixedsparrowped$animal)]
  
  
  summary(motherhash)
  # done!
  
  
  for(i in 1:length(motherhash[,1])){
  # look up which column in motherhash is the most left missing mother:
  j <- min(which(is.na(motherhash[i,])))
  
  # if j == 2, this is the bird and there is no matriline. 
  # What I did before is return NA here, but actually I can return the
  # ID of the bird themselves.
  
  if(j == 2){
    motherhash$matriline[i] <- motherhash$BirdID[i]
    motherhash$j[i] <- j
  } else {
  # one before this mother is the desired matriline mum:
    motherhash$matriline[i] <- motherhash[i,j-1]
    motherhash$j[i] <- j
  }
  }
  
  head(motherhash)
  tail(motherhash)
  
  # that looks good! :)
  
  # how many matrilines are there?
  str(motherhash$matriline)
  motherhash$matrilinenumeric <- as.numeric(motherhash$matriline)
  length(unique(motherhash$matriline))
  head(motherhash)
  tail(motherhash)
  motherhash$matrilinenumeric
  
  
  # 71 matrilines, that is an extra 45 compared to when I replaced
  # birds with no matriline with 'NA'
  table(table(motherhash$matriline))
  
  # add to the maleyear data frame
  maleyear$matriline <- motherhash$matrilinenumeric[match(maleyear$maleid,
                                                   motherhash$BirdID)]
  summary(maleyear$matriline)
  
  maleyear$factorMatriline <- maleyear$matriline
}

{
  maleh2EPO.pois.matriline <- MCMCglmm(EPO~factor(age6),
                                     ~factoranimal + factormaleID +
                                       factormaternalID + factorMatriline +
                                       factoryear,
                                     ginverse=list(factoranimal=invped.malelifetime),
                                     prior=prior5G.p,
                                     data=maleyear,
                                     family="poisson",
                                     nitt=1000000,
                                     thin=800,
                                     burnin=200000)
  plot(maleh2EPO.pois.matriline)
}



#-------------------------------------------
# Lifetime female behaviour as a multinomial variable
#-------------------------------------------


########## Run in R 3.0 ###########


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
  posterior.mode(femaleh2EPO.multi.yr$VCV)
  mean(femaleh2EPO.multi.yr$VCV[,1])
  mean(femaleh2EPO.multi.yr$VCV[,2])
}


# how does this model change if I use a slightly informative 
# prior: an inverse-Wishart with nu=0.02
{
  femaleh2EPO.multi.strongprior <- MCMCglmm(cbind(EPO, WPO)~1,
                                                  ~factoranimal + factorfemaleID +
                                                    factormaternalID + factoryear,
                                                  ginverse=list(factoranimal=invped.femalelifetime),
                                                  prior=prior4G.strong,
                                                  data=femaleyear,
                                                  family="multinomial2",
                                                  nitt=1000000,
                                                  thin=800,
                                                  burnin=200000)
  plot(femaleh2EPO.multi.strongprior)
  # fairly similar, with year estimated though. So year benefits
  # from the stronger prior.
  posterior.mode(femaleh2EPO.multi.strongprior$VCV)
}

#-------------------------------------------
# Per brood female behaviour as a multinomial variable
#-------------------------------------------



# with this being per brood, can and should account for the social
# male, his indirect genetic effect on the female, his phenotypic
# effect on the female, the male's year, and the year of reproduction.
# See Reid et al 2015 for the approach.
{
  broodEPO.dadage <- MCMCglmm(cbind(EPO, WPO)~factor(dadage6),
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
  plot(broodEPO.dadage)
  # results intriguing: some evidence for female additive genetics
  # for the trait, though not strong. No covariance. No male additive
  # genetics. Mum and Dad ID matter. Oh, but I forgot both mother's
  # mother and father's mother, to calculate maternal effects.
  # It also seems the social male's age doesn't affect this trait.
  autocorr(broodEPO.dadage$Sol)
  # age3:age5 0.0999
  autocorr(broodEPO.dadage$VCV)
  # fine.
}

# adding the maternal effects for the male and female.
# adding pairID to account for phenotypic similarity.
# adding a maternal shared effect for males and females
# in case maternal effects have the same outcome in each sex.
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


########## From here onward, run in R 3.2.3 ##########

{
  broodEPO.maternal.paircov.mums <- MCMCglmm(cbind(EPO, WPO)~factor(dadage6),
                                     random=~str(factoranimal + factordadanimal) + 
                                       factormaternalID + factorDadmaternalID + 
                                       factorfemaleID + factorDadID + factoryear +
                                       factorpairMaternal + factorpairID,
                                     family="multinomial2",
                                     prior=priorspecial8G.p,
                                     ginverse=list(factoranimal=invped.bothsexes.fullmums,
                                                   factordadanimal=invped.bothsexes.fullmums),
                                     data=broodWE.fullmums,
                                     nitt=1000000,
                                     thin=800,
                                     burnin=200000)
  plot(broodEPO.maternal.paircov.mums)
  # a bit of a negative effect of older father, so females more
  # faithful to oldest partners but not a lot.
  # again, quite wishy-washy for effects.
  # Too much in a model? There isn't a lot that can be taken out
  # to help model fit except the genetic covariance could be knocked out.
}

{
  broodEPO.maternal.pair.mums <- MCMCglmm(cbind(EPO, WPO)~factor(dadage6),
                                          random=~factoranimal + factordadanimal + 
                                            factormaternalID + factorDadmaternalID + 
                                            factorfemaleID + factorDadID + factoryear +
                                            factorpairMaternal + factorpairID,
                                          family="multinomial2",
                                          prior=prior9G.p,
                                          ginverse=list(factoranimal=invped.bothsexes.fullmums,
                                                        factordadanimal=invped.bothsexes.fullmums),
                                          data=broodWE.fullmums,
                                          nitt=1000000,
                                          thin=800,
                                          burnin=200000)
  plot(broodEPO.maternal.pair.mums)
  # Nothing. Nothing at all. This is either over-stretching the data or 
  # something else is going on.
  
  # one problem could be my EPP rates.
  table(offspring4$EPO)
  # 17% EPO. 
  format(table(offspring4$EPO, offspring4$Cohort)[2,]/(table(offspring4$EPO, offspring4$Cohort)[1,]+
                                                  table(offspring4$EPO, offspring4$Cohort)[2,]),digits=2)
  
  # Per year rates vary from a low of 9% to a high of 24%.
  # Jane Reid has an average of 28% EPO (Reid et al 2014, Evolution)
  # So to have a comparable sample size of EPO I would need 1.5 times
  # her data set size.
  # Maybe there are too many zero broods to work on the per brood level
  # for females, and I have to work on the per year level for females to
  # have sufficient variation.
}


broodEPO.nogenetics <- MCMCglmm(cbind(EPO, WPO)~1,
                                           random=~factormaternalID + factorDadmaternalID + 
                                             factorfemaleID + factorDadID + factoryear +
                                             factorpairMaternal + factorpairID,
                                           family="multinomial2",
                                           prior=prior7G.p,
                                           ginverse=list(factoranimal=invped.bothsexes.fullmums,
                                                         factordadanimal=invped.bothsexes.fullmums),
                                           data=broodWE.fullmums,
                                           nitt=100000,
                                           thin=80,
                                           burnin=20000)

plot(broodEPO.nogenetics)


broodEPO.nogenetics.pair.dad <- MCMCglmm(cbind(EPO, WPO)~1,
                                random=~factorfemaleID + factoryear + 
                                  factorDadID + factorpairID,
                                family="multinomial2",
                                prior=prior4G.p,
                                ginverse=list(factoranimal=invped.bothsexes.fullmums,
                                              factordadanimal=invped.bothsexes.fullmums),
                                data=broodWE.fullmums,
                                nitt=1000000,
                                thin=800,
                                burnin=200000)

plot(broodEPO.nogenetics.pair.dad)

broodEPO.nogenetics.pair <- MCMCglmm(cbind(EPO, WPO)~1,
                                     random=~factorfemaleID + factoryear + 
                                       factorpairID,
                                     family="multinomial2",
                                     prior=prior3G.p,
                                     ginverse=list(factoranimal=invped.bothsexes.fullmums,
                                                   factordadanimal=invped.bothsexes.fullmums),
                                     data=broodWE.fullmums,
                                     nitt=100000,
                                     thin=80,
                                     burnin=20000)

##############################################################################
# Bivariate model: genetic or phenotypic covariance between male and female
# behaviour
##############################################################################

# For the bivariate model, I will use male behaviour within a year and
# female behaviour within a year. I will use female behaviour within
# a year because the alternative is per brood, but the problem with per
# brood is the high number of broods with zero offspring.

# model structure:
# two dependent variables: 1) proportion of a female's offspring
# within a year that are EPO, 2) total EPO of a male within a year.
# One fixed effect: male age on dependent variable 2 ONLY, 
# as a six-level factor.
# Trait specific intercepts (since the traits have fundamentally
#                           different distributions).
# and the following random effects:
# bird ID with a within-sex estimate i.e. idh() structure.
# animal with a within-sex AND between-sex estimate i.e. us() structure.
# maternal effect with a within-sex AND between-sex estimate i.e. us() structure.
# year of reproduction. Only one estimate because between-year fluctuations 
#       in the total number of EPO should affect both sexes equally since
#       both sexes data come from this one data set.
# Residual covariance within sex i.e. idh() structure.

# bird ID and residual are estimated without covariance, since no male is
# assayed for the female trait and vice versa.
# animal and maternal ID have covariance because males and females can be
# linked by the pedigree or by their mother's ID.
# year is explained previously.

{
  priorbiv4.p <- list(R=list(V=diag(2), nu=0.002),
                     G=list(G1=list(V=diag(2), nu=2, 
                                    alpha.mu=c(0,0), alpha.V=diag(2)*1000),
                            G2=list(V=diag(2), nu=2, 
                                    alpha.mu=c(0,0), alpha.V=diag(2)*1000),
                            G3=list(V=diag(2), nu=2, 
                                    alpha.mu=c(0,0), alpha.V=diag(2)*1000),
                            G4=list(V=1, nu=1, alpha.mu=0, alpha.V=1000)))
  
  priorbiv3.p <- list(R=list(V=diag(2), nu=0.002),
                      G=list(G1=list(V=diag(2), nu=2, 
                                     alpha.mu=c(0,0), alpha.V=diag(2)*1000),
                             G2=list(V=diag(2), nu=2, 
                                     alpha.mu=c(0,0), alpha.V=diag(2)*1000),
                             G3=list(V=1, nu=1, alpha.mu=0, alpha.V=1000)))
  
  bothsexesyear$age <- bothsexesyear$maleage
  bothsexesyear$age[is.na(bothsexesyear$age)] <- 1
  bothsexesyear$age
  table(bothsexesyear$age)
  table(bothsexesyear$maleage)
  table(bothsexesyear$sex)
  
  bivariateEPO <- MCMCglmm(cbind(EPOm, cbind(EPOf,WPOf))~trait +
                             factor(age),
                            random=~us(trait):animal + idh(trait):birdID +
                             us(trait):maternalID + year,
                           rcov=~idh(trait):units,
                            family=c("poisson", "multinomial2"),
                            prior=priorbiv4.p,
                            ginverse=list(animal=invped.bothsexesyear),
                            data=bothsexesyear,
                            nitt=1000000,
                            thin=800,
                            burnin=200000)
  # bear in mind, for interpretation the female intercept is 
  # at the reference level of age==1. But the male has a separate
  # intercept. Age is not interacted with trait, and there is
  # only variation in the male trait. So the age should fit to
  # just the male trait data.
  plot(bivariateEPO)
  
  autocorr(bivariateEPO$Sol)
}

# with the correct construct for specifying the level on which
# a factor should act: at.level():

# minus missing dams, to check whether the estimates would remain
# the same without MCMCglmm adding dams:

{
  bothsexes.minusmissingdams <- bothsexesyear[-which(is.na(bothsexesyear$maternalID)),]
  
  bivariateEPO.minusmissingdams <- MCMCglmm(cbind(EPOm, cbind(EPOf,WPOf))~trait +
                             factor(age),
                           random=~us(trait):animal + idh(trait):birdID +
                             us(trait):maternalID + year,
                           rcov=~idh(trait):units,
                           family=c("poisson", "multinomial2"),
                           prior=priorbiv4.p,
                           ginverse=list(animal=invped.bothsexesyear),
                           data=bothsexes.minusmissingdams,
                           nitt=1000000,
                           thin=800,
                           burnin=200000)
  plot(bivariateEPO.minusmissingdams)
}

##############################################################################
# Per offspring model
##############################################################################

# In this model, each offspring is a data point.
# The random effects include 
# 1   the identity of the genetic mother and father (no covariance)
# 2   the genetic contribution of the mother and father and their genetic covariance
# 3   the year the offspring was born to account for differences between years
#     in the number of EPO
# 4   the maternal effect of the genetic mother and father and their covariance
# 5   

offspring5$combinedGrandmothers <- paste(offspring5$GenDadMaternalID,
                                         offspring5$GenMumMaternalID,
                                         sep="plus")
# now it is important to remove NA grandmothers:

offspring5.minusGMums1 <- offspring5[-which(is.na(offspring5$GenDadMaternalID)),]
offspring5.minusGMums <- offspring5.minusGMums1[-which(is.na(offspring5.minusGMums1$GenMumMaternalID)),]

offspring5.minusGMums$factorGrandmothers <- as.factor(offspring5.minusGMums$combinedGrandmothers)

offspringEPO <- MCMCglmm(EPO~factor(GenDadAge6),
                           random=~str(factorGenDadanimal + factorGenMumanimal) + 
                           factorGenDadID + factorGenMumID + 
                           factorGenDadMaternalID + factorGenMumMaternalID + 
                           factorGrandmothers + factoryear,
                           family="categorical",
                           prior=priorbinomial7,
                           ginverse=list(factorGenDadanimal=invped.bothsexesyear,
                                         factorGenMumanimal=invped.bothsexesyear),
                           data=offspring5.minusGMums,
                           nitt=1000000,
                           thin=800,
                           burnin=200000)

plot(offspringEPO)
# this model is UNSTABLE. I wonder why.

# One possibility is to remove the father's age, which
# is known to affect EPO but not affect WPO production.
# That inter-dependency could be causing the instability.

# run a model without the father's age:

offspringEPO.minusdadage <- MCMCglmm(EPO~1,
                         random=~str(factorGenDadanimal + factorGenMumanimal) + 
                           factorGenDadID + factorGenMumID + 
                           factorGenDadMaternalID + factorGenMumMaternalID + 
                           factorGrandmothers + factoryear,
                         family="categorical",
                         prior=priorbinomial7,
                         ginverse=list(factorGenDadanimal=invped.bothsexesyear,
                                       factorGenMumanimal=invped.bothsexesyear),
                         data=offspring5.minusGMums,
                         nitt=1000000,
                         thin=800,
                         burnin=200000)

plot(offspringEPO.minusdadage)

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
# Additional consideration 3: random repeatability
##############################################################################

# Jon Slate asked what the repeatability of promiscuous behaviour would be
# in a completely randomised data set. This is to establish whether the
# repeatability of promiscuity that we detect is a real individual effect or
# an artefact.

# take the male phenotypes:

permute.data.males <- data.frame(EPO=maleyear$EPO)
permute.data.males$random <- rnorm(length(permute.data.males$EPO))

head(permute.data.males)

permute.data.males2 <- permute.data.males[order(permute.data.males$random),]
head(permute.data.males2)
tail(permute.data.males2)

permute.data.males2$maleid <- maleyear$maleid
head(permute.data.males2)
tail(permute.data.males2)

permute.data.males2$factormaleID <- permute.data.males2$maleid

# now run the model on the permuted data:
{
  maler2.permut <- MCMCglmm(EPO~1,
                            ~factormaleID,
                            prior=prior1G.p,
                            data=permute.data.males2,
                            family="poisson",
                            nitt=100000,
                            thin=80,
                            burnin=20000)
  plot(maler2.permut)
}

malerepeatability <- maler2.permut$VCV[,"factormaleID"]/
  (maler2.permut$VCV[,"factormaleID"] +
     maler2.permut$VCV[,"units"] +
     log(1/exp(maler2.permut$Sol[,"(Intercept)"])+1))

posterior.mode(malerepeatability)
HPDinterval(malerepeatability)


##############################################################################
# Additional consideration 4: how does male repeatability change when
# any 'floater' type males are removed?
##############################################################################

length(which(maleyear$WPO==0))

maleyear.minusfloaters <- maleyear[-which(maleyear$WPO==0),]

table(table(maleyear$maleid))

maler2EPOpois.minusfloaters <- MCMCglmm(EPO~1,
                          ~factormaleID,
                          prior=prior1G.p,
                          data=maleyear.minusfloaters,
                          family="poisson",
                          nitt=1000000,
                          thin=800,
                          burnin=200000)

plot(maler2EPOpois.minusfloaters)
plot(maler2EPOpois)



malerepeatability <- maler2EPOpois$VCV[,"factormaleID"]/
  (maler2EPOpois$VCV[,"factormaleID"] +
     maler2EPOpois$VCV[,"units"] +
     log(1/exp(maler2EPOpois$Sol[,"(Intercept)"])+1))
posterior.mode(malerepeatability)


malerepeatability2 <- maler2EPOpois.minusfloaters$VCV[,"factormaleID"]/
  (maler2EPOpois.minusfloaters$VCV[,"factormaleID"] +
     maler2EPOpois.minusfloaters$VCV[,"units"] +
     log(1/exp(maler2EPOpois.minusfloaters$Sol[,"(Intercept)"])+1))

posterior.mode(malerepeatability2)

# If we account for male age and year of reproduction, does repeatability go up or down?
malerepeatability3 <- maler2EPOpois.age.yr$VCV[,"factormaleID"]/
  (maler2EPOpois.age.yr$VCV[,"factormaleID"] +
     maler2EPOpois.age.yr$VCV[,"units"] + maler2EPOpois.age.yr$VCV[,"factoryear"] +
     log(1/exp(maler2EPOpois.age.yr$Sol[,"(Intercept)"])+1))
posterior.mode(malerepeatability3)

# But this has the interesting problem that repeatability is relative to the intercept
# so this repeatability is relative to one year old males (which doesn't make
# sense...). Our estimate of repeatability will change with age, e.g.:

malerepeatability4 <- maler2EPOpois.age.yr$VCV[,"factormaleID"]/
  (maler2EPOpois.age.yr$VCV[,"factormaleID"] +
     maler2EPOpois.age.yr$VCV[,"units"] + maler2EPOpois.age.yr$VCV[,"factoryear"] +
     log(1/exp(maler2EPOpois.age.yr$Sol[,"(Intercept)"] +
                 maler2EPOpois.age.yr$Sol[,"factor(age6)2"])+1))
posterior.mode(malerepeatability4)

# would be the estimate for two year old males, except males are only sampled in
# one year. Therefore, only the intercept-only repeatability makes sense to me...

##############################################################################
# Additional consideration 5: are EPO offspring more likely to be promiscuous?
##############################################################################

# Add the EPO status of the male or female as a fixed factor to the
# model:

{
  bivariateEPO.parentEPO <- MCMCglmm(cbind(EPOm, cbind(EPOf,WPOf))~trait +
                             factor(age) + trait:factor(EPOstatus),
                           random=~us(trait):animal + idh(trait):birdID +
                             us(trait):maternalID + year,
                           rcov=~idh(trait):units,
                           family=c("poisson", "multinomial2"),
                           prior=priorbiv4.p,
                           ginverse=list(animal=invped.bothsexesyear),
                           data=bothsexes.EPOstatus,
                           nitt=1000000,
                           thin=800,
                           burnin=200000)
  
  plot(bivariateEPO.parentEPO)
}

# We know there is no detectable genetic background to either trait,
# so it makes sense to try the same model without the animal term to
# see how the other variance components could change.

{
  bivariateEPO.parentEPO.nogenetics <- MCMCglmm(cbind(EPOm, cbind(EPOf,WPOf))~trait +
                                       factor(age) + trait:factor(EPOstatus),
                                     random=~idh(trait):birdID +
                                       us(trait):maternalID + year,
                                     rcov=~idh(trait):units,
                                     family=c("poisson", "multinomial2"),
                                     prior=priorbiv3.p,
                                     data=bothsexes.EPOstatus,
                                     nitt=1000000,
                                     thin=800,
                                     burnin=200000)
  
  plot(bivariateEPO.parentEPO.nogenetics)
  # that really doesn't make sense.
}

##############################################################################
# Issues
##############################################################################

# To find solutions and information for:
# 1   birds with a different cohort between the database and the pedigree.
# 2   whether eggs are more likely to be assigned as extra-pair offspring.
# 3   whether the results are stable when uncertain parents are removed.
# 4   whether results are stable once Aaron has finished the database changes.
# 5   genetic parents that are the wrong sex.
# 6   social parents that are the wrong sex.
# 7   dead genetic parents
# 8 

##############################################################################
# Close the database
##############################################################################

close(sparrowDB)
