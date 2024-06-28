#' Hybrid Blau Model Data Import
#'
#' @description This function runs the data import script that allows users to import data that is readable by the simulation
#'
#' @param dataDir Data object you would like the function to use. Needs to be a data frame
#' @param nDim Dimensions of your matrix. Currently the model is capable of handling two or three dimensional matrices (default = 3)
#' @param entities Number of organizations in the model (default = 16)
#' @param entlist list of the organization names/types
#' @param d1lab Label for dimension 1
#' @param d2lab Label for dimension 2
#' @param d3lab Label for dimension 3
#'
#' @return Object list of data and necessary objects for the simulation
#' @import
#' @export
#'
#' @examples
#' datalist <- HBSmDataImp(dataDir = HBSmexampledat,
#' nDim = 3,
#' entities=16,
#' entlist=orgnames,
#' d1lab = "educ",
#' d2lab = "age",
#' d3lab = "ic_income")

#Hybrid Model v6#
#Data Import Codes#
#Code updates and automation - 8/30/2020#
HBSmDataImp <- function(dataDir, nDim = 3, entities, entlist, d1lab, d2lab, d3lab = NULL) {
  if(is.null(nDim)){
    nDim <- 3}
  #assign('nDim', nDim, envir = .GlobalEnv)
  #assign('entlist', entlist, envir = .GlobalEnv)
  #assign('entities', entities, envir = .GlobalEnv)

#Part 1 - Import into dataframe#
Importdat <-data.frame()
#GSS data#
#Importdat <- read.csv("~/Dropbox/__Across Computer Doc Share/R Central/GSS Cleaned/GSS04_membershipdata-cleaned.csv", header = TRUE, sep = ",")
#Importdat <- read.csv(file = dataDir)
Importdat <- dataDir
#"C:/Users/calfr/OneDrive/Desktop/Backup/UrnModel/Data/GSS Cleaned/GSS 1974_new - Rready.csv"
#assign('Importdat', Importdat, envir = .GlobalEnv)
#Importdat <- read.csv("~/Dropbox/__Across Computer Doc Share/ONR Grant/Convocation 8 MP List/Conv 8 Blaunet data - session 7 - move.csv")
#Importdat <- read.csv("~/Dropbox/SPPA dataset/SPPA_clean.csv")

#Sperate into years - this is done for the SPPA dataset becasue I chose not the break up the data by year after cleaning it#
#Importdat_1982 = subset(Importdat, year == 1982)
#Importdat_1985 = subset(Importdat, year == 1985)
#Importdat_1992 = subset(Importdat, year == 1992)
#Importdat_1997 = subset(Importdat, year == 1997)
#Importdat_2002 = subset(Importdat, year == 2002)
#Importdat_2008 = subset(Importdat, year == 2008)
#Importdat_2012 = subset(Importdat, year == 2012)

#Ideal Point Data#
#Importdat <- read.csv("~/Dropbox/__Across Computer Doc Share/R Central/Ideal Point Data/Declaration and ideal point data clean/Conv 8 Session 3 ideal and dec - new.csv", header = TRUE, sep = ",")

#Changing the form of median_pt_rounded_v2 - #Line 16-34 for Ideal Point data only
#Importdat$median_pt.round.v2[Importdat$median_pt.round.v2 == -4.50] <- 1
#Importdat$median_pt.round.v2[Importdat$median_pt.round.v2 == -4.00] <- 2
#Importdat$median_pt.round.v2[Importdat$median_pt.round.v2 == -3.50] <- 3
#Importdat$median_pt.round.v2[Importdat$median_pt.round.v2 == -3.00] <- 4
#Importdat$median_pt.round.v2[Importdat$median_pt.round.v2 == -2.50] <- 5
#Importdat$median_pt.round.v2[Importdat$median_pt.round.v2 == -2.00] <- 6
#Importdat$median_pt.round.v2[Importdat$median_pt.round.v2 == -1.50] <- 7
#Importdat$median_pt.round.v2[Importdat$median_pt.round.v2 == -1.00] <- 8
#Importdat$median_pt.round.v2[Importdat$median_pt.round.v2 == -0.50] <- 9
#Importdat$median_pt.round.v2[Importdat$median_pt.round.v2 ==  0.00] <- 10
#Importdat$median_pt.round.v2[Importdat$median_pt.round.v2 ==  0.50] <- 11
#Importdat$median_pt.round.v2[Importdat$median_pt.round.v2 ==  1.00] <- 12
#Importdat$median_pt.round.v2[Importdat$median_pt.round.v2 ==  1.50] <- 13
#Importdat$median_pt.round.v2[Importdat$median_pt.round.v2 ==  2.00] <- 14
#Importdat$median_pt.round.v2[Importdat$median_pt.round.v2 ==  2.50] <- 15
#Importdat$median_pt.round.v2[Importdat$median_pt.round.v2 ==  3.00] <- 16
#Importdat$median_pt.round.v2[Importdat$median_pt.round.v2 ==  3.50] <- 17
#Importdat$median_pt.round.v2[Importdat$median_pt.round.v2 ==  4.00] <- 18
#Importdat$median_pt.round.v2[Importdat$median_pt.round.v2 ==  4.50] <- 19

#Importdat["loc_id"] #Just a random variable in the dataset

#Creation of demographic object
demographics <- list(d1lab, d2lab, d3lab)
#assign('demographics', demographics, envir = .GlobalEnv)

#Part 2 - Matrix Creation# - automation for number of dimensions - Dimension Import
#dim1var <- as.numeric(unlist(Importdat["edu.code"])) #For Ideal Point data
dim1var <- as.numeric(unlist(Importdat[d1lab]))
#dim1var <- as.numeric(unlist(Importdat["Educ"]))
#dim1var <- as.numeric(unlist(Importdat["Education"]))
dim1var <- dim1var + 1 # for if education has a zero in it. PROGRAM CAN NOT HANDLE ZERO VALUES IN DIMENSIONS!!!!!
#assign('dim1var', dim1var, envir = .GlobalEnv)

d1 <- range(dim1var, na.rm = TRUE)

d1
#assign('d1', d1, envir = .GlobalEnv)
dim2var <- as.numeric(unlist(Importdat[d2lab]))
#assign('dim2var', dim2var, envir = .GlobalEnv)
#dim2var <- as.numeric(unlist(Importdat["Age"]))

d2 <- range(dim2var, na.rm = TRUE)

d2
#assign('d2', d2, envir = .GlobalEnv)

#dim3var <- as.numeric(unlist(Importdat["dim3var"]))
#assistants <- as.numeric(unlist(Importdat["X..of.Assistants"]))
#assistants <- assistants + 1

#idealpoint <- as.numeric(unlist(Importdat["median_pt.round.v2"])) #for Ideal Point data


if (nDim == 3) {
  dim3var <- as.numeric(unlist(Importdat[d3lab]))
  #assign('dim3var', dim3var, envir = .GlobalEnv)
  d3 <- range(dim3var, na.rm = TRUE)
#d3 <- range(assistants, na.rm = TRUE)

d3
#assign('d3', d3, envir = .GlobalEnv)
}

#Manually edit ranges, not needed for all ranges, check data beforehand#
d1r <- d1[2]#-d1[1] #Education
#assign('d1r', d1r, envir = .GlobalEnv)
d2r <- d2[2]#-d2[1] #Age - needs adjusted to be out in the correct range (bottom of range is the youngest age in the space)
#assign('d2r', d2r, envir = .GlobalEnv)
if (nDim == 3) {d3r <- d3[2]
#assign('d3r', d3r, envir = .GlobalEnv)#-d3[1] #dim3var# #Ideal Point adds the -d3[1]
}

if (nDim == 2) {template_array <- array(0,dim=c(d1r,d2r), dimnames=NULL)} else
  {template_array <- array(0,dim=c(d1r,d2r,d3r), dimnames=NULL)}
#assign('template_array', template_array, envir = .GlobalEnv)
#Number of organizations/ entities in model#
#entities <- 13

#entities <- 15

#entities <- 9

#entities <- 20

OGorgs <- list()

for(i in 1:entities) {
  OGorgs[[i]] <- template_array
}
#assign('OGorgs', OGorgs, envir = .GlobalEnv)
#Because the range for the SPPA is slightly different for age (and dim3var in a prior model attempt), cordinate locations of individuals needs to be imported for

#Part 2 - Matrix Creation# - automation for number of dimensions - Dimension Import
#GSS uses lowercase for variable names, the Rada and SPPA use uppercase (if I am correct)#
#tempedulev_year <- as.numeric(unlist(Importdat_2012["educ"]))
tempedulev_year <- as.numeric(unlist(Importdat[d1lab]))
tempedulev_year <- tempedulev_year + 1 # for if education has a zero in it. PROGRAM CAN NOT HANDLE ZERO VALUES IN DIMENSIONS!!!!!
#assign('tempedulev_year', tempedulev_year, envir = .GlobalEnv)
#tempage_year <- as.numeric(unlist(Importdat_2012["age"]))
tempage_year <- as.numeric(unlist(Importdat[d2lab]))
#assign('tempage_year', tempage_year, envir = .GlobalEnv)
#income_year <- as.numeric(unlist(Importdat_2012["dim3var"]))
if (nDim == 3) {income_year <- as.numeric(unlist(Importdat[d3lab]))
#assign('income_year', income_year, envir = .GlobalEnv)
}

#Range check#
d1_check <- range(tempedulev_year, na.rm = TRUE)

d1_check
#assign('d1_check', d1_check, envir = .GlobalEnv)

d2_check <- range(tempage_year, na.rm = TRUE)

d2_check
#assign('d2_check', d2_check, envir = .GlobalEnv)

if (nDim == 3) {d3_check <- range(income_year, na.rm = TRUE)

d3_check
#assign('d3_check', d3_check, envir = .GlobalEnv)
}

coordims1 <- as.numeric(tempedulev_year)#-d1[1]) #y coord
#assign('coordims1', coordims1, envir = .GlobalEnv)
coordims2 <- as.numeric(tempage_year)#-d2[1]) #x coord
#assign('coordims2', coordims2, envir = .GlobalEnv)
if (nDim == 3) { coordims3 <- as.numeric(income_year)
#assign('coordims3', coordims3, envir = .GlobalEnv)#-d3[1]) #z coord #for ideal point data -d3[1] is needed
#coordims3 <- as.numeric(assistants)#-d3[1]) #z coord #for ideal point data -d3[1] is needed
#coordims3 <- as.numeric(idealpoint) #z
}


#Part 3#
#lengthterm = nrow(Importdat_2012) #number of observations in the data
lengthterm = nrow(Importdat) #number of observations in the data

p=1 #counter for the number of individuals in the dataset (lengthterm is the limit of this)

#Orgnames <- vector()
#Orgnames <- c("CLASSICAL","OPERA","TUNES_BAND","JAZZ","REGGAE_RAP","MOOD_DANCE","BLUES","COUNTRY","BLGRASS","ROCK_POP","FOLK","CHORAL_HYMNS","ETHNIC_LATIN")
#Orgnames <- c("CLASSICAL","OPERA","TUNES","JAZZ","REGGAE_RAP","BLUES","BAND_PARADE_BARBER","COUNTRY","BLGRASS","ROCK","MOOD_NEW","FOLK","CHORAL_HYMNS","ETHNIC_LATIN","OTHER")
Orgnames <- entlist
#assign('Orgnames', Orgnames, envir = .GlobalEnv)

#print(Orgnames[o])

#To switch between datasets, change the value here from o=22 to o=34. This is just updating the column index for the loop to start counting.
#o=20 #Indexing value to start at the correct row of the table(Importdat) for placing memberships
#o=34
#o=21
#o=8
o=1

#o=1

stopnum = o+entities #This should be the length of something...
#stopnum = 14 #This should be the length of something...
#stopnum = 14 #This should be the length of something...

i = 1



repeat{
  #temporg <- Importdat_2012[,Orgnames[o]]
  temporg <- Importdat[,Orgnames[o]]


  #assign('temporg', temporg, envir = .GlobalEnv)
  #temporg <- Importdat[,o]

  logInd = temporg > 0
  #assign('logInd', logInd, envir = .GlobalEnv)
  repeat{
    if(logInd[p] == TRUE){
      if (nDim == 2) {OGorgs[[i]][coordims1[p],coordims2[p]] = OGorgs[[i]][coordims1[p],coordims2[p]] + temporg[p]} else
      {OGorgs[[i]][coordims1[p],coordims2[p],coordims3[p]] = OGorgs[[i]][coordims1[p],coordims2[p],coordims3[p]] + temporg[p]}
      #assign('OGorgs', OGorgs, envir = .GlobalEnv)
    }
    p = p+1
    print(p) # for trouble shooting if import is incorrect
    if(p==lengthterm + 1)
      break
  }
  print(i)
  p=1
  i=i+1
  o=o+1

  if(o==stopnum)
    break
}

#Troubleshooting Orgs array#
#sum(OGorgs[[1]])
#sum(Importdat_1985$CLASSICAL)
#sum(OGorgs[[2]])
#sum(Importdat_1985$OPERA)

#loop to generate array of individual locations
p=1
totindi <- template_array

repeat{
  if (nDim == 2) {totindi[coordims1[p],coordims2[p]] = totindi[coordims1[p],coordims2[p]] + 1} else
  {totindi[coordims1[p],coordims2[p],coordims3[p]] = totindi[coordims1[p],coordims2[p],coordims3[p]] + 1}

  p = p+1
  print(p)
  if(p==lengthterm + 1)
    break
}

#Code below removes all but the entities, the orglist (may rename to be entities list), and coordinates for individual's
#locations. To keep additional objects, but them in the parenthesis.#
#ls()



#Check total resources - do every time you set up!
sum(Reduce("+",OGorgs))

sum(Reduce("+",totindi))
#assign('OGorgs', OGorgs, envir = .GlobalEnv)
#assign('totindi', totindi, envir = .GlobalEnv)
#Orgnames <- c("CLASSICAL","OPERA","TUNES_BAND","JAZZ","REGGAE_RAP","DANCE_MOOD","BLUES","COUNTRY","BLGRASS","ROCK_POP","FOLK","CHORAL_HYMNS","ETHNIC_LATIN")

#This setup code is for the SPPA dataset music consumption data#
#sum(Importdat_2012$CLASSICAL,Importdat_2012$OPERA,Importdat_2012$TUNES_BAND,Importdat_2012$JAZZ,Importdat_2012$REGGAE_RAP,Importdat_2012$MOOD_DANCE,Importdat_2012$BLUES,Importdat_2012$COUNTRY,Importdat_2012$BLGRASS,Importdat_2012$ROCK_POP,Importdat_2012$FOLK,Importdat_2012$CHORAL_HYMNS,Importdat_2012$ETHNIC_LATIN)


#sum(Importdat_2012$EVERY)
#sum(Importdat_2012$CLASSICAL)
#sum(Importdat_2012$OPERA)
#sum(Importdat_2012$TUNES_BAND)
#sum(Importdat_2012$JAZZ)
#sum(Importdat_2012$REGGAE_RAP)
#sum(Importdat_2012$MOOD_DANCE)
#sum(Importdat_2012$BLUES)
#sum(Importdat_2012$COUNTRY)
#sum(Importdat_2012$BLGRASS)
#sum(Importdat_2012$ROCK_POP)
#sum(Importdat_2012$FOLK)
#sum(Importdat_2012$CHORAL_HYMNS)
#sum(Importdat_2012$ETHNIC_LATIN)

#sum(OGorgs[[1]])
#sum(OGorgs[[2]])
#sum(OGorgs[[3]])
#sum(OGorgs[[4]])
#sum(OGorgs[[5]])
#sum(OGorgs[[6]])
#sum(OGorgs[[7]])
#sum(OGorgs[[8]])
#sum(OGorgs[[9]])
#sum(OGorgs[[10]])
#sum(OGorgs[[11]])
#sum(OGorgs[[12]])
#sum(OGorgs[[13]])

if (nDim == 2) {rm(list= ls()[!(ls() %in% c('OGorgs','entities','coordims1','coordims2','totindi','nDim'))])} else
{ rm(list= ls()[!(ls() %in% c('OGorgs','entities','coordims1','coordims2','coordims3','totindi','nDim'))]) }

if (nDim == 2) {
  result_list <- list(
    OGorgs = OGorgs,
    entities = entities,
    coordims1 = coordims1,
    coordims2 = coordims2,
    totindi = totindi,
    nDim = nDim
  )
} else {
  result_list <- list(
    OGorgs = OGorgs,
    entities = entities,
    coordims1 = coordims1,
    coordims2 = coordims2,
    coordims3 = coordims3,
    totindi = totindi,
    nDim = nDim
  )
}
return(result_list)
}
