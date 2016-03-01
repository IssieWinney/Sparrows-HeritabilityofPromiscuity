##############################################################################
# loading the pedigree
##############################################################################

# Isabel Winney
# 20160115
# This is the file to import the Lundy house sparrow pedigree to a 
# project. The pedigree is updated every one to two years and the
# pedigree import procedure changes accordingly.

# This is for the pedigree up to and including 2015 data.

##############################################################################
# Manual changes to pedigree file
##############################################################################

# some individuals have been removed from the pedigree
# file for being their own mothers or being duplicated:
# 6250 is duplicated, and is its own mother
# 6385 is its own mother
# 6451 is its own mother
# 6526 is its own mother
# 6565 is its own mother

# These are observations 5073, 5205, 5256, 5326, 5361.
# I have manually removed these records from the
# original file.

# 4378, 5048, 6238, 6245
# These individuals were missing cohorts and were needed for 
# exploration and age analysis. Aded 2005 for 4378 (but this is
# unsure because capture date is a Shinichi 1st Jan date), 2009
# for 5048, and 2011 for 6238 and 6245


# One bird (bird ID 6833, asrid 5626) had a cohort of "xxx".
# this bird was tracked since it was a nestling, so the cohort
# was manually changed to 2013.


##############################################################################
# Loading the pedigree
##############################################################################

# read in the pedigree file:
sparrowpedigree <- read.table("C:/Users/Issie/SkyDrive/PhD/masterdatasheets/Pedigreeto2015-20160201.txt",
                   header=T, na.strings="NA")

# Manual additions:
# Find all blank cells and replace with NA.

# check the birds that were missing in the 2013 version are
# present in this version:
which(sparrowpedigree$id==5644)
which(sparrowpedigree$id==7641)
which(sparrowpedigree$id==7649)


head(sparrowpedigree)
tail(sparrowpedigree)
str(sparrowpedigree)

# are there any unusual cohorts:
unique(sparrowpedigree$Cohort)

# do any individuals turn up in both the sire and dam columns:
unique(sparrowpedigree$sire[match(sparrowpedigree$dam, sparrowpedigree$sire)])
# 4975 occurs as both a dam and a sire.

##############################################################################
# Metadata
##############################################################################

# birdid      the numerical identity of the bird, the same as the numerical ID
#             in the database (BirdID)
# dam         the genetic mother as assigned by microsatellite data (ID is her 
#             BirdID in the database)
# sire        the genetic father as assigned by microsatellite data (ID is his 
#             BirdID in the database)
# asrid       birdid but with birds numbered (nearly) consecutively. This is here for
#             analyses with ASReml-R which used to complain that BirdID had 
#             numbers missing from the sequence
# asrdam      genetic mother where her ID is her asrid
# asrsire     genetic father where his ID is his asrid
