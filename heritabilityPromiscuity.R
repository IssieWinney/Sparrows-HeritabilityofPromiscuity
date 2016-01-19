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
  # set appropriate factors for analysis
  str(maleyear)
  maleyear$factoranimal <- as.factor(maleyear$animal)
  maleyear$factormaternalID <- as.factor(maleyear$maternalID)
  maleyear$factorpaternalID <- as.factor(maleyear$paternalID)
  maleyear$factorcohort <- as.factor(maleyear$Cohort)
  str(maleyear)
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
  
  
}

# and then the pedigrees can be converted to relationship matrices
# for the current version of MCMCglmm:
{
  invped.malelifetime <- inverseA(prunedped.malelifetime)$Ainv
  summary(invped.malelifetime)
}

##############################################################################
# Plots of variables
##############################################################################

# how many EPO offspring do males have?
hist(malephenotypes$EPO, breaks=seq(-1,26,1))
# lots of zeroes.

# how many are EPO?
hist(malephenotypes$EPO/(malephenotypes$EPO+malephenotypes$WPO))
# lots of males with no EPO, some males with only EPO

# in the per year data?
hist(maleyear$EPO, breaks=seq(-1,16,1))
hist(maleyear$EPO/(maleyear$EPO+maleyear$WPO))

##############################################################################
# Priors
##############################################################################

# parameter expanded priors. Number before G is number of random effects. 
# .p for parameter expansion
{
  prior3G.p = list(R=list(V=1,nu=0.002),
                  G=list(G1=list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000),
                         G2=list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000),
                         G3=list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000)))
  
  prior4G.p = list(R=list(V=1,nu=0.002),
                   G=list(G1=list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000),
                          G2=list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000),
                          G3=list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000),
                          G4=list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000)))
}

##############################################################################
# Repeatability
##############################################################################

{
  maler2EPO <- MCMCglmm(cbind(EPO, WPO)~1,
                      ~factormale,
                      prior=prior1G,
                      data=maleyear,
                      family="multinomial2",
                      nitt=1000000,
                      thin=800,
                      burnin=200000)

  plot(maler2EPO)
}

##############################################################################
# Heritability
##############################################################################

# a model of heritability should have as a minimum the ID of the bird (if there
# are repeated measures), a pedigree linked bird ID, and a maternal ID.
# Strongly recommended: cohort

#-------------------------------------------
# Lifetime male behaviour as a multinomial variable
#-------------------------------------------

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
}

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

# problems with these models: males that are older can contribute more
# EPO. Male lifespan would need to go in a lifetime model.

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
