#' Hybrid Blau Model Sim
#'
#' @description This function runs the main simulation model while allowing you to edit the simulations global parameters.
#'
#' @param glob Total number of runs that the model simulation will complete (default = 10)
#' @param chec Number of iterations before convergence check (default = 500). If no convergence is reached, the model will continue until it reaches the values set for maxr or until it finally does converge
#' @param maxr Number of runs until model is forced to stopped if no convergence is reached (default = 4000)
#' @param tol How little difference between org arrays there can be before the model is considered converged (default = 8)
#' @param nrange Range of the neighborhood cells. A neighborhood cell is considered to be the cells surrounding the focal cell that are within range if those value (default = 1). If r = 1 then there is one focal cell and 8 neighborhood cells.
#' @param convr Number of permutations for convergence check
#' @param imp Toggle for import data
#' @param bin Toggle for if any bin variables are on
#' @param HBSmdata Set this to equal the object you generated with the data import function
#' @param prefix prefix to your saveout data
#' @param year year of your saveout data
#' @param saveDir directory for where your data is saved to
#' @param notif Plays a sound when model is done. Set to TRUE or FALSE.
#'
#' @return The results of the model. The list outputlist has the relevant outputs of the model while objectdump contains objects used to run the model.
#' @import Matrix, beepr, cdfquantreg, stringr, truncnorm
#' @export
#'
#' @examples
#' HBSmSim(glob = 1,
#' chec = 200,
#' maxr = 500,
#' tol = 8,
#' nrange = 1,
#' convr = 15,
#' HBSmdata = datalist,
#' prefix="GSS",
#' year = "1974",
#' saveDir = saveDir)
HBSmSim <- function(glob = 10, chec = 500, maxr = 4000, tol = 8, nrange = 1, convr = 15, imp = 1, bin = 1, HBSmdata, prefix = "", year = "", saveDir, notif = TRUE) {

  if (!exists("nDim")) {
    nDim <- 3}
  if(is.null(glob)){
    glob <- 10}
  if(is.null(chec)){
    chec <- 500}
  if(is.null(maxr)){
    maxr <- 4000}
  if(is.null(tol)){
    tol <- 8}
  if(is.null(nrange)){
    nrange <- 1}
  if(is.null(convr)){
    convr <- 15}
  if(is.null(imp)){
    imp <- 1}
  if(is.null(bin)){
    bin <- 1}
  if(exists("entities")==F){
  entities <- 16}
  if(exists("OGorgs")==F){
    OGorgs <- list()}
  limit <- HBSmdata$entities + 1
  testlimit <- entities + 1
  mlc <- 0
  mem <- 0 #just for now until we fix mem. Structure still requires the object
  OGorgs <- HBSmdata$OGorgs
  entities <- HBSmdata$entities
  coordims1 <- HBSmdata$coordims1
  coordims3 <- HBSmdata$coordims2
  coordims3 <- HBSmdata$coordims3
  totindi <- HBSmdata$totindi
  nDim <- HBSmdata$nDim
  #assign('nDim', nDim, envir = .GlobalEnv)
  #assign('glob', glob, envir = .GlobalEnv) # total number of runs +1 - 26 standard #Playing with a different location for this#
  #assign('limit', entities + 1, envir = .GlobalEnv)  #limit for loops, set to number of orgs +1
  #assign('mlc', 0, envir = .GlobalEnv) #counter for the number of runs
  #assign('errorcount', 0, envir = .GlobalEnv) #Counter for the number of times a run is trapped in a loop or fails to converge (error runs in general)#
  #assign('itde', 1, envir = .GlobalEnv)
  #assign('chec', chec, envir = .GlobalEnv) #500 #750 #1000 #1500 #Number of iterations that pass before the model does a convergence check#
  #assign('maxr', maxr, envir = .GlobalEnv) #10000 #35000 #Number of iterations that pass until the model force stops and provides and error message#
  #assign('tol', tol, envir = .GlobalEnv) #Base is 1 #How little difference org arrays can have before convergence is considered reached#
  #assign('nrange', r, envir = .GlobalEnv) #range of cell neighborhood selection#
  #assign('convr', convr, envir = .GlobalEnv) #15 #10 #15 #20 #30 #The number of permutations that are considered when doing a convergence check#
  #mem effect weighting and setup#
  #Note, if the first observation point of a data set is being run, then both mem and import should be set to 0#
  #assign('mem', mem, envir = .GlobalEnv)
  #assign('import', imp, envir = .GlobalEnv)
  #assign('bin', bin, envir = .GlobalEnv) #set as 1 if there are bin variables, set to 0 is there is not bin variables #Change to a 1 for GSS data
  #assign('prefix', prefix, envir = .GlobalEnv)
  #assign('year', year, envir = .GlobalEnv)
  #assign('saveDir', saveDir, envir = .GlobalEnv)
  #assign('notif', notif, envir = .GlobalEnv)


  #urndataimp(dataDir = "C:/Users/calfr/OneDrive/Desktop/Backup/UrnModel/Data/GSS Cleaned/GSS 1974_new - Rready.csv")
  #Temp location for summary and parameter lists, don't know the best place to put these#
  exploitMhist <- list() #List for exploitation of members values#
  WCChist <- list() #List for weighted carrying capacity values#
  extenhist <- list() #List for extensivenes values#
  intencellhist <- list() #List for Cell Focused Intensiveness values#
  intenranghist <- list() #List for Range Focused Intensiveness values#
  sumorghist <- list() #List for sum of entity population#
  orgs_x_minhist <- list() #List of min values for x dimension#
  orgs_x_maxhist <- list() #List of max values for x dimension#
  orgs_y_minhist <- list() #List of min values for y dimension#
  orgs_y_maxhist <- list() #List of max values for y dimension#
  orgs_z_minhist <- list() #List of min values for z dimension#
  orgs_z_maxhist <- list() #List of max values for z dimension#
  orgs <- list()
  countlist <-list() #List created to keep count of the number of runs, more of a troubleshooting check. Will be part of the log and parameters
  locrun <- list()
  #assign('exploitMhist', exploitMhist, envir = .GlobalEnv)
  #assign('WCChist', WCChist, envir = .GlobalEnv)
  #assign('intencellhist', intencellhist, envir = .GlobalEnv)
  #assign('intenranghist', intenranghist, envir = .GlobalEnv)
  #assign('sumorghist', sumorghist, envir = .GlobalEnv)
  #assign('orgs_x_minhist', orgs_x_minhist, envir = .GlobalEnv)
  #assign('orgs_x_maxhist', orgs_x_maxhist, envir = .GlobalEnv)
  #assign('orgs_y_minhist', orgs_y_minhist, envir = .GlobalEnv)
  #assign('orgs_y_maxhist', orgs_y_maxhist, envir = .GlobalEnv)
  #assign('orgs_z_minhist', orgs_z_minhist, envir = .GlobalEnv)
  #assign('orgs_z_maxhist', orgs_z_maxhist, envir = .GlobalEnv)
  #assign('orgs', orgs, envir = .GlobalEnv)
  #assign('countlist', countlist, envir = .GlobalEnv)
  #assign('locrun', locrun, envir = .GlobalEnv)
  #assign('extenhist', extenhist, envir = .GlobalEnv)

  ###Needs commenting out for now mem effect###
  #if (mem == 1) {
  w <- .5 #.1 #.5 #.75 #1 #weight of the location of resource during the last observation period
  #prev_orgs is the intital observations#
  #if (import == 0){
  #prev_orgs <- readRDS("C:/Users/calfr/OneDrive/Desktop/Backup/UrnModel/Data/OGorgs_1974.rds") #saveout of r object will include year, make sure to change year on this to be correct. Will work on getting r to feed year into this labels later#
  #prev_orgs is the ecology space after the last permutation of the final iteration of the model#
  #if (import == 1){
  #prev_orgs <- readRDS("C:/Users/calfr/OneDrive/Desktop/Backup/UrnModel/Data/mem/orgs_1974.rds")
  #  }
  #}

  #assign('prev_orgs', prev_orgs, envir = .GlobalEnv)

  ###End of mem effect###

  ####Directory Setup####
  #The following code sets up the main directory for output, prefix setup,
  prefix <- prefix
  #prefix <- "DIPP"
  #prefix <- "Rada"
  #prefix <- "SPPA"

  year <- year

  saveDir <- saveDir #Mac
  dir.create(saveDir)
  #setwd(saveDir) #to make sure that the main directory works (trouble shooting)

  coordDir <- paste0("/coord")

  coordpath <- paste0(saveDir, coordDir)
  #assign('coordpath', coordpath, envir = .GlobalEnv)

  dir.create(file.path(coordpath))
  #assign('saveDir', saveDir, envir = .GlobalEnv)

  ####bin Variables####
  #Are there any bin variables? If so make this value a 1#
  #Bivar = 3 #Put the dimension that the bin variable is in here, need to uncomment for GSS data

  #Start time of program
  starttime <- Sys.time()

  ####Scripts loop####

  repeat{

    orgs <- OGorgs
    #assign('orgs', orgs, envir = .GlobalEnv)
    simerror <- 0 #Error value that, when increased to anything above 0, stops the current simulation
    #List Creation and array formation# - Merged with Simulation script
    #source('~/Dropbox/__Across Computer Doc Share/R Central/Code update - 472020/3 Dimension Codes/HBSm setup new m6 - 3D.R', echo=TRUE)

    #Simulation Script - Loop#
    #Save Objects List Creation#
    orglist <- list() #Parallel arrays for each entity are saved here and updated after each simulation run#
    sumorgs <- list() #Saves the sum of each entity for the run#
    totpops <- list() #The total population of the ecology for each run#
    simerror = 0
    errorcount = 0
    #adj<-0
    Coord = list() #List of coordinaes where any change happens#
    orgscoord <- list() #List of coordinates where changes happen for each entity in the ecology#
    Change = list() #List of changes made to each entity in the ecology#
    Changepos = list() #Running list of positive changes made to each entity for each run of the simulation#
    changpos = list() #List of positive changes made to each entities#
    Changeneg = list() #Running list of negative changes made to each entity for each run of the simulation#
    changneg = list() #List of negative changes made to entities#
    #assign('sumorgs', sumorgs, envir = .GlobalEnv)
    #assign('totpops', totpops, envir = .GlobalEnv)
    #assign('Coord', Coord, envir = .GlobalEnv)
    #assign('orgscoord', orgscoord, envir = .GlobalEnv)
    #assign('Change', Change, envir = .GlobalEnv)
    #assign('Changepos', Changepos, envir = .GlobalEnv)
    #assign('changpos', changpos, envir = .GlobalEnv)
    #assign('Changeneg', Changeneg, envir = .GlobalEnv)
    #assign('changneg', changneg, envir = .GlobalEnv)

    skiptrip = 0 #Count of times that skip logic is used#
    #assign('skiptrip', skiptrip, envir = .GlobalEnv)
    #totchange <-0
    legrun = 0 #The number of sets that have been run for the simulation run
    #assign('legrun', legrun, envir = .GlobalEnv)
    #Values for testing the update loops (mem and old)#
    oldupdate <- 0
    newupdate <- 0
    #assign('oldupdate', oldupdate, envir = .GlobalEnv)
    #assign('newupdate', newupdate, envir = .GlobalEnv)
    #Testing Resets#
    #orgs <- OGorgs
    #orglist <- list()

    #Simulatin mode code#
    #Setting of repeatition parameters for main loop(?)#
    Set <- 0 #Starting value for number of sets run#
    rundiff <- 0 #Creation of object for difference between runs, updated after each simulation set (each run of the simulation in sets of 100)#
    lowdiff <- 0 #Creation of object for counting number of times convergence criteria is reached#
    TR <- 0 #Object for counting total number of iterations the simulation carries out (times it selects a new focal cell)# - changed from 0 to 1
    CR <- 0
    convnum <- 0
    pass <- 0
    conversave <- list()

    IT <-0
    i <- as.numeric(1)
    a <-0
    p <-0

    #Load needed libraries#


    ####Main Loop####
    repeat{
      IT <-0
      #i <- as.numeric(1)
      #a <-0
      #p <-0

      ####Inner Loop####
      repeat{

        #Inner and total Interation counts, and printing count numbers #
        TR = TR+1
        IT = IT+1
        #assign('TR', TR, envir = .GlobalEnv)
        #assign('IT', IT, envir = .GlobalEnv)
        print(paste0("IT",IT))
        print(paste0("TR",TR))

        ###Selection Loop###
        repeat{
          #dims and focells#
          if (nDim == 3) {mod_dim <- dim(orgs[[1]][,,])
          #assign('mod_dim', mod_dim, envir = .GlobalEnv)
          }
          else {mod_dim <- dim(orgs[[1]][,])
          #assign('mod_dim', mod_dim, envir = .GlobalEnv)
          }

          cord_x <- sample(mod_dim[2],1)
          cord_y <- sample(mod_dim[1],1)
          #assign('cord_x', cord_x, envir = .GlobalEnv)
          #assign('cord_y', cord_y, envir = .GlobalEnv)
          if (nDim == 3) {cord_z <- sample(mod_dim[3],1)
          #assign('cord_z', cord_z, envir = .GlobalEnv)
          }


          #New neighborhood selection code#
          #x nrange designation#

          if (cord_x > nrange & cord_x <= mod_dim[2]) {nbhcord_x_neg <- (cord_x-nrange)
          }else {nbhcord_x_neg = 1}

          if (cord_x >= 1 + nrange & cord_x < mod_dim[2]-nrange) {nbhcord_x_pos <- (cord_x+nrange)
          }else {nbhcord_x_pos = mod_dim[2]}

          nbhcord_x_nrange <- (nbhcord_x_neg:nbhcord_x_pos)
          #assign('nbhcord_x_neg', nbhcord_x_neg, envir = .GlobalEnv)
          #assign('nbhcord_x_pos', nbhcord_x_pos, envir = .GlobalEnv)
          #assign('nbhcord_x_nrange', nbhcord_x_nrange, envir = .GlobalEnv)

          #y nrange designation#
          if (cord_y > nrange & cord_y <= mod_dim[1]) {nbhcord_y_neg <- (cord_y-nrange)
          }else {nbhcord_y_neg = 1}

          if (cord_y >= 1 + nrange & cord_y < mod_dim[1]-nrange) {nbhcord_y_pos <- (cord_y+nrange)
          }else {nbhcord_y_pos = mod_dim[1]}

          nbhcord_y_nrange <- (nbhcord_y_neg:nbhcord_y_pos)
          #assign('nbhcord_y_neg', nbhcord_y_neg, envir = .GlobalEnv)
          #assign('nbhcord_y_pos', nbhcord_y_pos, envir = .GlobalEnv)
          #assign('nbhcord_y_nrange', nbhcord_y_nrange, envir = .GlobalEnv)

          #z nrange designation#
          if (nDim == 3) {if (cord_z > nrange & cord_z <= mod_dim[3]) {nbhcord_z_neg <- (cord_z-nrange)
          }else {nbhcord_z_neg = 1}

            if (cord_z >= 1 + nrange & cord_z < mod_dim[3]-nrange) {nbhcord_z_pos <- (cord_z+nrange)
            }else {nbhcord_z_pos = mod_dim[3]}

            nbhcord_z_nrange <- (nbhcord_z_neg:nbhcord_z_pos)
            #assign('nbhcord_z_neg', nbhcord_z_neg, envir = .GlobalEnv)
            #assign('nbhcord_z_pos', nbhcord_z_pos, envir = .GlobalEnv)
            #assign('nbhcord_z_nrange', nbhcord_z_nrange, envir = .GlobalEnv)
            }

          #focal cell selection and neighborhood creation#
          #setup#
          focelist <- list()
          nbhlist <- list()
          Prob_org <- list()

          #Loop for subsetting array to create neighborhood#
          l <-1

          if (nDim == 3) {
            repeat{
              orgtemp <- array(Matrix::as.matrix(orgs[[l]]), dim=c(mod_dim[1],mod_dim[2],mod_dim[3]), dimnames = NULL)
              focelist[[l]] <- c(orgtemp[cord_y,cord_x,cord_z])
              nbhlist[[l]] <- c(sum(orgtemp[nbhcord_y_nrange,nbhcord_x_nrange,nbhcord_z_nrange]))
              l = l+1
              if (l == limit)
                break
            }
          }else
          {repeat{
            orgtemp <- array(as.matrix(orgs[[l]]), dim=c(mod_dim[1],mod_dim[2]), dimnames = NULL)
            focelist[[l]] <- c(orgtemp[cord_y,cord_x])
            nbhlist[[l]] <- c(sum(orgtemp[nbhcord_y_nrange,nbhcord_x_nrange]))
            l = l+1
            if (l == limit)
              break
          }
          }
          #assign('orgtemp', orgtemp, envir = .GlobalEnv)
          #assign('focelist', focelist, envir = .GlobalEnv)
          #assign('nbhlist', nbhlist, envir = .GlobalEnv)
          #Loop for creating a subset array of the same location in the ecology, but for a previous time period#
          #Note: will need the R-object for the previous array of the ecology imported#
          l <-1

          #  if(mem == 1){
          #   pre_focelist <- list()
          #   pre_nbhlist <- list()
          #   if (nDim == 3) {repeat{
          #     oldorgtemp <- array(as.matrix(prev_orgs[[l]]), dim=c(mod_dim[1],mod_dim[2],mod_dim[3]), dimnames = NULL)
          #     pre_focelist[[l]] <- as.numeric(oldorgtemp[cord_y,cord_x,cord_z])
          #     pre_nbhlist[[l]] <- as.numeric(sum(oldorgtemp[nbhcord_y_nrange,nbhcord_x_nrange,nbhcord_z_nrange]))
          #     l = l+1
          #     if (l == limit)
          #       break
          #   }
          #   } else {repeat{
          #     oldorgtemp <- array(as.matrix(prev_orgs[[l]]), dim=c(mod_dim[1],mod_dim[2]), dimnames = NULL)
          #     pre_focelist[[l]] <- as.numeric(oldorgtemp[cord_y,cord_x])
          #     pre_nbhlist[[l]] <- as.numeric(sum(oldorgtemp[nbhcord_y_nrange,nbhcord_x_nrange]))
          #     l = l+1
          #     if (l == limit)
          #       break
          #   }
          #   }
          #   #assign('oldorgtemp', oldorgtemp, envir = .GlobalEnv)
          #   #assign('pre_focelist', pre_focelist, envir = .GlobalEnv)
          #   #assign('pre_nbhlist', pre_nbhlist, envir = .GlobalEnv)
          # }

          #Reduces the focelist and nbhlist objects to a sum and then saves them as a different object#
          focsum <- Reduce("+",focelist)
          nbhsum <- Reduce("+",nbhlist)
          #assign('focsum', focsum, envir = .GlobalEnv)
          #assign('nbhsum', nbhsum, envir = .GlobalEnv)

          #Reduces the pre_focelist and pre_nbhlist objects to a sum and then saves them as a different object#
          #Contains new code#
          # if(mem == 1){
          #   pre_focsum <- Reduce("+",pre_focelist)
          #   if(length(pre_focsum) == 0){
          #     pre_focsum <- 0
          #   }
          #   pre_nbhsum <- Reduce("+",pre_nbhlist)
          #   if(length(pre_nbhsum) == 0){
          #     pre_nbhsum <- 0
          #   }
          #   #assign('pre_focsum', pre_focsum, envir = .GlobalEnv)
          #   #assign('pre_nbhsum', pre_nbhsum, envir = .GlobalEnv)
          # }


          ###Skip code for when nbhsum would sum to zero###
          #This code skips the summation process for the current cell and then triggers an alternate code for saving the cell cordinates
          #of locations the model made changes#
          if(nbhsum == 0){
            skiptrip <- skiptrip+1 #Counts the number of times skips are made, mostly for trouble shooting#

            l <- 1
            repeat{
              changpos[[l]] = 0
              changneg[[l]] = 0

              l=l+1

              if(l==limit)
                break
            }

            i <- 1

            tempsum <- list()

            repeat{
              tempsum[[i]] <- sum(orgs[[i]])
              i = i+1

              if(i == entities+1){
                break
              }
            }

            sumorgs[[TR]] <- c(tempsum)

            temptot <- Reduce("+", orgs)
            totpops[[TR]] <-Reduce("+",temptot)

            #Save-out of coordinates of each selected cell where changes are made#
            #In this case, becasue the model makes no changes, no coordinate
            Changepos[[TR]] <- changpos
            Changeneg[[TR]] <- changneg
            Change[[TR]] <- unlist(changpos) + unlist(changneg)
            #New code for changelist logic issue#
            #if(Change[[TR]] == 0){
            #  Change[[TR]] <- 0
            #}
            tchange = Change[[TR]]

            temporgscoord <- list()

            if(any(tchange > 0)){
              a = a+1

              if (nDim == 3) {Coord[[a]] <- c(cord_x,cord_y,cord_z)
              } else {Coord[[a]] <- c(cord_x,cord_y)}
            }

            e <- 1
            for(e in 1:entities){
              if(Change[[TR]][e] != 0){
                #c[[i]] = c[[i]]+1 #I honestly don't know what this code is for#
                if (nDim == 3) {temporgscoord[[e]] <- c(cord_y,cord_x,cord_z)}
                else {temporgscoord[[e]] <- c(cord_y,cord_x)} #For a 2D ecology#
              }
              else{
                if (nDim == 3) {temporgscoord[[e]] <- c(0,0,0)}
                else {temporgscoord[[e]] <- c(0,0)} #For a 2D ecology#
              }
            }
            orgscoord[[TR]] <- temporgscoord

            totchange = Reduce('+',Change)
            #}
            break
          }
          break
          #assign('totchange', totchange, envir = .GlobalEnv)
        }
        #assign('sumorgs', sumorgs, envir = .GlobalEnv)
        #assign('totpops', totpops, envir = .GlobalEnv)
        #assign('Changepos', Changepos, envir = .GlobalEnv)
        #assign('Changeneg', Changeneg, envir = .GlobalEnv)
        #assign('Change', Change, envir = .GlobalEnv)
        #assign('Coord', Coord, envir = .GlobalEnv)
        #assign('orgscoord', orgscoord, envir = .GlobalEnv)


        ###Simulation Loop###
        #This loop runs if nbhsum is greater then 0, otherwise the previous loop runs#
        if(nbhsum >0){

          #if(mem == 0 | pre_nbhsum == 0){
          if(mem == 0){
            l <-1
            oldupdate = oldupdate+1
            repeat{
              Prob_org[[l]] <-c(nbhlist[[l]]/nbhsum)
              Prob_org[[l]][is.na(Prob_org[[l]])] <- 0
              l=l+1
              if (l ==limit)
                break
            }
          }
          #Modified code below to include if statement "mem == 1" so that this is skipped if mem is set to 0#
          # if(mem == 1){
          #   if(pre_nbhsum == 0){
          #     l <-1
          #     oldupdate = oldupdate+1
          #     repeat{
          #       Prob_org[[l]] <-c(nbhlist[[l]]/nbhsum)
          #       Prob_org[[l]][is.na(Prob_org[[l]])] <- 0
          #       l=l+1
          #       if (l ==limit)
          #         break
          #     }
          #   }
          # }

          #Probability of being in an organization based off the prevelence of the organization in the neighborhood and the prevelence during the last model#
          #loops#
          # if(mem == 1){
          #   #if(mem == 1 & pre_nbhsum > 0){
          #   if(pre_nbhsum > 0){
          #     l <-1
          #     newupdate = newupdate+1
          #     repeat{
          #       Prob_org[[l]] <-c((nbhlist[[l]]+(pre_nbhlist[[l]]*w))/nbhsum) #Check and make sure this is correct
          #       #Prob_org[[l]] <-c((nbhlist[[l]]*((pre_nbhlist[[l]]/pre_nbhsum)*w))/(nbhsum))
          #       Prob_org[[l]][is.na(Prob_org[[l]])] <- 0
          #       l=l+1
          #       if (l ==limit)
          #         break
          #     }
          #   }
          # }

          Probmax <- max(unlist(Prob_org))
          Probmax[Probmax < 0] <- 0

          if(Probmax > 1){
            l <- 1
            repeat{
              Prob_org[[l]] = Prob_org[[l]]/Probmax
              #Prob_org[[l]] = Prob_org[[l]]/focsum
              l=l+1
              if (l ==limit)
                break
            }
          }
          #assign('Prob_org', Prob_org, envir = .GlobalEnv)
          #assign('Probmax', Probmax, envir = .GlobalEnv)

          probtest<- unlist(Prob_org)

          #List creation and paramter reset#
          stolist <- list()

          e <- 0
          l <- 1

          ####Stolcastic Rounding Nested Repeat Fuction####
          #Part 1 - Round#
          repeat{
            repeat{
              stol = (focsum*Prob_org[[l]])

              q <- abs(stol - trunc(stol))
              #Adjusted rounding code for scalable cell neighborhood#
              probrang <- c(1 - q, q)

              adj <- sample(0:1, size = 1, replace = FALSE, prob = probrang)

              if(stol <= 0) {adj = 0}

              stol = trunc(stol)+adj
              stolist[[l]] <-c(stol)
              l = l+1
              if(l==limit)
                break
              #}
            }
            #assign('probrang', probrang, envir = .GlobalEnv)
            #assign('adj', adj, envir = .GlobalEnv)

            #Part 2 - Check round and direct repeat(if criteria not met)#
            stolsum <- Reduce("+",stolist)
            if(stolsum == focsum) {e = e + 1}
            else{
              if(stolsum > focsum){
                repeat{
                  mintemp <- min(probtest[probtest > 0])
                  #assign('mintemp', mintemp, envir = .GlobalEnv)
                  loctemp <- match(mintemp,probtest)
                  #assign('loctemp', loctemp, envir = .GlobalEnv)
                  stolist[[loctemp]] <- stolist[[loctemp]] - 1
                  if(stolist[[loctemp]] == 0){
                    probtest[loctemp] <- 0
                  }
                  stolsum <- Reduce("+",stolist)
                  if(stolsum == focsum){
                    e = e + 1
                    break
                  }
                }
              }

              l <- 1
            }

            if(e == 1)
              break
          }
          #assign('stolsum', stolsum, envir = .GlobalEnv)
          #assign('probtest', probtest, envir = .GlobalEnv)
          #assign('stolist', stolist, envir = .GlobalEnv)


          #Inner and total Interation counts, and printing count numbers #
          #TR = TR+1
          #IT = IT+1
          #print(paste0("IT",IT))
          #print(paste0("TR",TR))

          #Changes Lists#
          ctempval <-list()

          l = 1

          repeat{
            if (nDim == 3) {ctempval[[l]] <- c(stolist[[l]]-orgs[[l]][cord_y,cord_x,cord_z])}
            else {ctempval[[l]] <- c(stolist[[l]]-orgs[[l]][cord_y,cord_x])} #For a 2D ecology#

            if(ctempval[[l]] > 0) {changpos[[l]] = ctempval[[l]]}
            else {changpos[[l]] = 0}

            if(ctempval[[l]] < 0) {changneg[[l]] = ctempval[[l]]}
            else {changneg[[l]] = 0}

            l=l+1
            if(l==limit)
              break
          }
          #assign('ctempval', ctempval, envir = .GlobalEnv)

          #Updating cells#
          l<-1
          repeat{
            if (nDim == 3) {orgs[[l]][cord_y,cord_x,cord_z] <- stolist[[l]]}
            else {orgs[[l]][cord_y,cord_x] <- stolist[[l]]} #For a 2D ecology#
            orgs[[l]][is.na(orgs[[l]])] <- 0
            l=l+1
            if(l ==limit)
              break
          }
          #assign('orgs', orgs, envir = .GlobalEnv)

          #Sum of orgs and save of coordinates#
          tempsum <- list()

          i <- 1
          repeat{
            tempsum[[i]] <- sum(orgs[[i]])
            i = i+1

            if(i == entities+1)
              break
          }
          #assign('tempsum', tempsum, envir = .GlobalEnv)

          sumorgs[[TR]] <- c(tempsum)

          temptot <- Reduce("+", orgs)
          totpops[[TR]] <-Reduce("+",temptot)
          #assign('sumorgs', sumorgs, envir = .GlobalEnv)
          #assign('temptot', temptot, envir = .GlobalEnv)
          #assign('totpops', totpops, envir = .GlobalEnv)

          #Save-out of coordinates#
          #This saves out coordinate for locations where changes were made in the model#
          Changepos[[TR]] <- changpos
          Changeneg[[TR]] <- changneg
          Change[[TR]] <- unlist(changpos) + unlist(changneg)

          tchange = Change[[TR]]

          temporgscoord <- list()

          if(any(tchange > 0)){
            a = a+1
            if (nDim == 2) {Coord[[a]] <- c(cord_x,cord_y)} #For a 2D ecology#
            else {Coord[[a]] <- c(cord_x,cord_y,cord_z)}
          }

          e <- 1

          for(e in 1:entities){
            if(Change[[TR]][e] != 0){
              #c[[i]] = c[[i]]+1 #I honestly don't know what this code is for#
              if (nDim == 3) {temporgscoord[[e]] <- c(cord_y,cord_x,cord_z)}
              else {temporgscoord[[e]] <- c(cord_y,cord_x)} #For a 2D ecology#
            }
            else{
              if (nDim == 3) {temporgscoord[[e]] <- c(0,0,0)}
              else {temporgscoord[[e]] <- c(0,0)} #For a 2D ecology#
            }
          }
          orgscoord[[TR]] <- temporgscoord

          totchange = Reduce('+',Change)
          #break
        }
        #assign('totchange', totchange, envir = .GlobalEnv)
        #assign('orgscoord', orgscoord, envir = .GlobalEnv)
        #assign('Change', Change, envir = .GlobalEnv)
        #assign('Coord', Coord, envir = .GlobalEnv)
        #assign('tchange', tchange, envir = .GlobalEnv)
        #assign('Changepos', Changepos, envir = .GlobalEnv)
        #assign('Changeneg', Changeneg, envir = .GlobalEnv)

        ####Convergence Loop####
        if(IT == chec){
          Set = Set + 1
          print(paste0("Set",Set))

          Runsum = unlist(Changepos[[TR]])
          thisrun = Reduce('+',Runsum)
          if(Set > 2){
            rundiff = abs(thisrun - lastrun)
          }
          legrun <- Set
          #legrun[Set] = thisrun
          lastrun = thisrun
        }

        ###Print out of differences between run sets (the chec object). Useful for checking if changes are being made###
        ###If no changes are made, there is a problem with the intital data or code###
        #print(paste0("rundiff"," ",rundiff))

        if(TR %% 10==0){ #was TR#
          CR = CR+1
          orglist[[CR]] <- orgs#[[i]]
        }
        #assign('orglist', orglist, envir = .GlobalEnv)
        #assign('CR', CR, envir = .GlobalEnv)
        #IT <- 1000

        ####Convergence Code####
        #This code is done in multiple parts, which consolidate the ecology by each organization and then compare it to the tol.#
        if(IT == chec){
          temptest <- list()
          aveorgtemp <- list()
          aveorgtemp2 <- list()
          l <- 1
          #leng <- 1

          repeat{
            #Saves the length of the orglist vector#
            lenglim <- length(orglist)
            #Creates a starting values that is an adjustment of the lenght of the orglist list so that only the last n permutations are considered#
            # N here is convr#
            leng <- (lenglim - convr)
            #if(lenglim <= convr){
            # leng <- lenglim}

            #if(lenglim >= convr){
            #leng <- lenglim - convr}

            #Create a list with all permuations of the ecology that are within the specified convergence nrange#
            repeat{
              aveorgtemp <- c(orglist[leng])
              aveorgtemp2[[lenglim - leng]] <- aveorgtemp[[1]][[l]]
              leng = leng+1
              if(leng==lenglim)#+1)
                break
            }
            #Averaging the ecology by each cell. Here an ecology is created where the number of resources drawn from a cell are the
            #average across the convergence nrange#
            #Note, each organziation gets its own ecology here, so the next section does comparisons cross all orgs#
            #What was this code doing - G#
            aveorgtemp3 <- Reduce(`+`, aveorgtemp2) / length(aveorgtemp2)
            temptest[[l]] <- aveorgtemp3

            l=l+1
            if(l==limit)
              break
          }
          #assign('temptest', temptest, envir = .GlobalEnv)
          #assign('aveorgtemp3', aveorgtemp3, envir = .GlobalEnv)
          #assign('aveorgtemp2', aveorgtemp2, envir = .GlobalEnv)
          #assign('aveorgtemp', aveorgtemp, envir = .GlobalEnv)
          #assign('leng', leng, envir = .GlobalEnv)
          #assign('lenglim', lenglim, envir = .GlobalEnv)
          #}

          #Var creation for the convergence logic#
          converg <- list()
          convergcount <- list()

          l <- 1
          convnum <- convnum + 1

          #This loop creates the object for comparison for convergence by taking the temptest list, which has the averages from the
          #Convergence nrange being considered and creates a vectors of times that the averages are not nearly the same as the values in
          #The orgs list (where permuted results are saved). A count of mismatches is created and then compared to the tol specified by the user#
          repeat{
            converg[[l]] <- all.equal.raw(temptest[[l]],orgs[[l]],tol = 0, check.attributes = TRUE)
            convergcount[[l]] <- stringr::str_remove(converg[[l]]," element mismatch")
            convergcount[[l]] <- stringr::str_remove(convergcount[[l]],"es")
            convergcount[[l]][convergcount[[l]] == "TRUE"] <- 0
            convergcount[[l]] <- as.numeric(convergcount[[l]])
            conversave[[convnum]] <- convergcount
            l=l+1
            if(l==limit)
              break
          }

          pass <- 0

          #Comparison of the mismatches to the tol. This loops does the comparison and also counts the number of passes
          # or times that the mismatch is lower or equal to the tol. If the count of passed entities equals entities, then the model is considered convergenced.#
          l <- 1
          repeat{
            if(convergcount[[l]] <= tol) {pass = pass+1}
            else {pass = pass}
            l=l+1
            if(l==limit)
              break
          }
          #Number of times
          print(paste0("Number of orgs converged ",pass))

        }
        if(IT == chec){
          break
          #IT <- 0
          #break
        }
      }
      #assign('converg', converg, envir = .GlobalEnv)
      #assign('conversave', conversave, envir = .GlobalEnv)
      #assign('convergcount', convergcount, envir = .GlobalEnv)
      #assign('convnum', convnum, envir = .GlobalEnv)
      #Logic code for passing or failing convergence check#
      if(pass==entities|TR>=maxr){
        if(pass==entities){
          print("The model has converged!")
        }
        if(TR>=maxr){
          print("Model is close to reaching a pre-defined limit for iterations and has not converged. Simulation had stopped")
          simerror = simerror + 1
          errorcount = errorcount + 1
        }
        break
      }
    }

    if(simerror >= 1 & simerror <= 3){
      #Save Objects List Creation#
      orglist <- list() #Parallel arrays for each entity are saved here and updated after each simulation run#
      sumorgs <- list() #Saves the sum of each entity for the run#
      totpops <- list() #The total population of the ecology for each run#
      simerror = 0
      errorcount = 0
      #adj<-0
      Coord = list() #List of coordinaes where any change happens#
      orgscoord <- list() #List of coordinates where changes happen for each entity in the ecology#
      Change = list() #List of changes made to each entity in the ecology#
      Changepos = list() #Running list of positive changes made to each entity for each run of the simulation#
      changpos = list() #List of positive changes made to each entities#
      Changeneg = list() #Running list of negative changes made to each entity for each run of the simulation#
      changneg = list() #List of negative changes made to entities#
      #assign('sumorgs', sumorgs, envir = .GlobalEnv)
      #assign('totpops', totpops, envir = .GlobalEnv)
      #assign('Coord', Coord, envir = .GlobalEnv)
      #assign('orgscoord', orgscoord, envir = .GlobalEnv)
      #assign('Change', Change, envir = .GlobalEnv)
      #assign('Changepos', Changepos, envir = .GlobalEnv)
      #assign('changpos', changpos, envir = .GlobalEnv)
      #assign('Changeneg', Changeneg, envir = .GlobalEnv)
      #assign('changneg', changneg, envir = .GlobalEnv)

      skiptrip = 0 #Count of times that skip logic is used#
      #assign('skiptrip', skiptrip, envir = .GlobalEnv)
      #totchange <-0
      legrun = 0 #The number of sets that have been run for the simulation run
      #assign('legrun', legrun, envir = .GlobalEnv)
      #Values for testing the update loops (mem and old)#
      oldupdate <- 0
      newupdate <- 0
      #assign('oldupdate', oldupdate, envir = .GlobalEnv)
      #assign('newupdate', newupdate, envir = .GlobalEnv)
      #Testing Resets#
      #orgs <- OGorgs
      #orglist <- list()

      #Simulatin mode code#
      #Setting of repeatition parameters for main loop(?)#
      Set <- 0 #Starting value for number of sets run#
      rundiff <- 0 #Creation of object for difference between runs, updated after each simulation set (each run of the simulation in sets of 100)#
      lowdiff <- 0 #Creation of object for counting number of times convergence criteria is reached#
      TR <- 0 #Object for counting total number of iterations the simulation carries out (times it selects a new focal cell)# - changed from 0 to 1
      CR <- 0
      convnum <- 0
      pass <- 0
      conversave <- list()

      IT <-0
      i <- as.numeric(1)
      a <-0
      p <-0


      ####Main Loop####
      repeat{
        IT <-0
        #i <- as.numeric(1)
        #a <-0
        #p <-0

        ####Inner Loop####
        repeat{

          #Inner and total Interation counts, and printing count numbers #
          TR = TR+1
          IT = IT+1
          #assign('TR', TR, envir = .GlobalEnv)
          #assign('IT', IT, envir = .GlobalEnv)
          print(paste0("IT",IT))
          print(paste0("TR",TR))

          ###Selection Loop###
          repeat{
            #dims and focells#
            if (nDim == 3) {mod_dim <- dim(orgs[[1]][,,])
            #assign('mod_dim', mod_dim, envir = .GlobalEnv)
            }
            else {mod_dim <- dim(orgs[[1]][,])
            #assign('mod_dim', mod_dim, envir = .GlobalEnv)
            }

            cord_x <- sample(mod_dim[2],1)
            cord_y <- sample(mod_dim[1],1)
            #assign('cord_x', cord_x, envir = .GlobalEnv)
            #assign('cord_y', cord_y, envir = .GlobalEnv)
            if (nDim == 3) {cord_z <- sample(mod_dim[3],1)
            #assign('cord_z', cord_z, envir = .GlobalEnv)
            }


            #New neighborhood selection code#
            #x nrange designation#

            if (cord_x > nrange & cord_x <= mod_dim[2]) {nbhcord_x_neg <- (cord_x-nrange)
            }else {nbhcord_x_neg = 1}

            if (cord_x >= 1 + nrange & cord_x < mod_dim[2]-nrange) {nbhcord_x_pos <- (cord_x+nrange)
            }else {nbhcord_x_pos = mod_dim[2]}

            nbhcord_x_nrange <- (nbhcord_x_neg:nbhcord_x_pos)
            #assign('nbhcord_x_neg', nbhcord_x_neg, envir = .GlobalEnv)
            #assign('nbhcord_x_pos', nbhcord_x_pos, envir = .GlobalEnv)
            #assign('nbhcord_x_nrange', nbhcord_x_nrange, envir = .GlobalEnv)

            #y nrange designation#
            if (cord_y > nrange & cord_y <= mod_dim[1]) {nbhcord_y_neg <- (cord_y-nrange)
            }else {nbhcord_y_neg = 1}

            if (cord_y >= 1 + nrange & cord_y < mod_dim[1]-nrange) {nbhcord_y_pos <- (cord_y+nrange)
            }else {nbhcord_y_pos = mod_dim[1]}

            nbhcord_y_nrange <- (nbhcord_y_neg:nbhcord_y_pos)
            #assign('nbhcord_y_neg', nbhcord_y_neg, envir = .GlobalEnv)
            #assign('nbhcord_y_pos', nbhcord_y_pos, envir = .GlobalEnv)
            #assign('nbhcord_y_nrange', nbhcord_y_nrange, envir = .GlobalEnv)

            #z nrange designation#
            if (nDim == 3) {if (cord_z > nrange & cord_z <= mod_dim[3]) {nbhcord_z_neg <- (cord_z-nrange)
            }else {nbhcord_z_neg = 1}

              if (cord_z >= 1 + nrange & cord_z < mod_dim[3]-nrange) {nbhcord_z_pos <- (cord_z+nrange)
              }else {nbhcord_z_pos = mod_dim[3]}

              nbhcord_z_nrange <- (nbhcord_z_neg:nbhcord_z_pos)
              #assign('nbhcord_z_neg', nbhcord_z_neg, envir = .GlobalEnv)
              #assign('nbhcord_z_pos', nbhcord_z_pos, envir = .GlobalEnv)
              #assign('nbhcord_z_nrange', nbhcord_z_nrange, envir = .GlobalEnv)
              }

            #focal cell selection and neighborhood creation#
            #setup#
            focelist <- list()
            nbhlist <- list()
            Prob_org <- list()

            #Loop for subsetting array to create neighborhood#
            l <-1

            if (nDim == 3) {
              repeat{
                orgtemp <- array(as.matrix(orgs[[l]]), dim=c(mod_dim[1],mod_dim[2],mod_dim[3]), dimnames = NULL)
                focelist[[l]] <- c(orgtemp[cord_y,cord_x,cord_z])
                nbhlist[[l]] <- c(sum(orgtemp[nbhcord_y_nrange,nbhcord_x_nrange,nbhcord_z_nrange]))
                l = l+1
                if (l == limit)
                  break
              }
            }else
            {repeat{
              orgtemp <- array(as.matrix(orgs[[l]]), dim=c(mod_dim[1],mod_dim[2]), dimnames = NULL)
              focelist[[l]] <- c(orgtemp[cord_y,cord_x])
              nbhlist[[l]] <- c(sum(orgtemp[nbhcord_y_nrange,nbhcord_x_nrange]))
              l = l+1
              if (l == limit)
                break
            }
            }
            #assign('orgtemp', orgtemp, envir = .GlobalEnv)
            #assign('focelist', focelist, envir = .GlobalEnv)
            #assign('nbhlist', nbhlist, envir = .GlobalEnv)
            #Loop for creating a subset array of the same location in the ecology, but for a previous time period#
            #Note: will need the R-object for the previous array of the ecology imported#
            l <-1

            # if(mem == 1){
            #   pre_focelist <- list()
            #   pre_nbhlist <- list()
            #   if (nDim == 3) {repeat{
            #     oldorgtemp <- array(as.matrix(prev_orgs[[l]]), dim=c(mod_dim[1],mod_dim[2],mod_dim[3]), dimnames = NULL)
            #     pre_focelist[[l]] <- as.numeric(oldorgtemp[cord_y,cord_x,cord_z])
            #     pre_nbhlist[[l]] <- as.numeric(sum(oldorgtemp[nbhcord_y_nrange,nbhcord_x_nrange,nbhcord_z_nrange]))
            #     l = l+1
            #     if (l == limit)
            #       break
            #   }
            #   } else {repeat{
            #     oldorgtemp <- array(as.matrix(prev_orgs[[l]]), dim=c(mod_dim[1],mod_dim[2]), dimnames = NULL)
            #     pre_focelist[[l]] <- as.numeric(oldorgtemp[cord_y,cord_x])
            #     pre_nbhlist[[l]] <- as.numeric(sum(oldorgtemp[nbhcord_y_nrange,nbhcord_x_nrange]))
            #     l = l+1
            #     if (l == limit)
            #       break
            #   }
            #   }
            #   #assign('oldorgtemp', oldorgtemp, envir = .GlobalEnv)
            #   #assign('pre_focelist', pre_focelist, envir = .GlobalEnv)
            #   #assign('pre_nbhlist', pre_nbhlist, envir = .GlobalEnv)
            # }

            #Reduces the focelist and nbhlist objects to a sum and then saves them as a different object#
            focsum <- Reduce("+",focelist)
            nbhsum <- Reduce("+",nbhlist)
            #assign('focsum', focsum, envir = .GlobalEnv)
            #assign('nbhsum', nbhsum, envir = .GlobalEnv)

            #Reduces the pre_focelist and pre_nbhlist objects to a sum and then saves them as a different object#
            #Contains new code#
            # if(mem == 1){
            #   pre_focsum <- Reduce("+",pre_focelist)
            #   if(length(pre_focsum) == 0){
            #     pre_focsum <- 0
            #   }
            #   pre_nbhsum <- Reduce("+",pre_nbhlist)
            #   if(length(pre_nbhsum) == 0){
            #     pre_nbhsum <- 0
            #   }
            #   #assign('pre_focsum', pre_focsum, envir = .GlobalEnv)
            #   #assign('pre_nbhsum', pre_nbhsum, envir = .GlobalEnv)
            # }


            ###Skip code for when nbhsum would sum to zero###
            #This code skips the summation process for the current cell and then triggers an alternate code for saving the cell cordinates
            #of locations the model made changes#
            if(nbhsum == 0){
              skiptrip <- skiptrip+1 #Counts the number of times skips are made, mostly for trouble shooting#

              l <- 1
              repeat{
                changpos[[l]] = 0
                changneg[[l]] = 0

                l=l+1

                if(l==limit)
                  break
              }

              i <- 1

              tempsum <- list()

              repeat{
                tempsum[[i]] <- sum(orgs[[i]])
                i = i+1

                if(i == entities+1){
                  break
                }
              }

              sumorgs[[TR]] <- c(tempsum)

              temptot <- Reduce("+", orgs)
              totpops[[TR]] <-Reduce("+",temptot)

              #Save-out of coordinates of each selected cell where changes are made#
              #In this case, becasue the model makes no changes, no coordinate
              Changepos[[TR]] <- changpos
              Changeneg[[TR]] <- changneg
              Change[[TR]] <- unlist(changpos) + unlist(changneg)
              #New code for changelist logic issue#
              #if(Change[[TR]] == 0){
              #  Change[[TR]] <- 0
              #}
              tchange = Change[[TR]]

              temporgscoord <- list()

              if(any(tchange > 0)){
                a = a+1

                if (nDim == 3) {Coord[[a]] <- c(cord_x,cord_y,cord_z)
                } else {Coord[[a]] <- c(cord_x,cord_y)}
              }

              e <- 1
              for(e in 1:entities){
                if(Change[[TR]][e] != 0){
                  #c[[i]] = c[[i]]+1 #I honestly don't know what this code is for#
                  if (nDim == 3) {temporgscoord[[e]] <- c(cord_y,cord_x,cord_z)}
                  else {temporgscoord[[e]] <- c(cord_y,cord_x)} #For a 2D ecology#
                }
                else{
                  if (nDim == 3) {temporgscoord[[e]] <- c(0,0,0)}
                  else {temporgscoord[[e]] <- c(0,0)} #For a 2D ecology#
                }
              }
              orgscoord[[TR]] <- temporgscoord

              totchange = Reduce('+',Change)
              #}
              break
            }
            break
            #assign('totchange', totchange, envir = .GlobalEnv)
          }
          #assign('sumorgs', sumorgs, envir = .GlobalEnv)
          #assign('totpops', totpops, envir = .GlobalEnv)
          #assign('Changepos', Changepos, envir = .GlobalEnv)
          #assign('Changeneg', Changeneg, envir = .GlobalEnv)
          #assign('Change', Change, envir = .GlobalEnv)
          #assign('Coord', Coord, envir = .GlobalEnv)
          #assign('orgscoord', orgscoord, envir = .GlobalEnv)


          ###Simulation Loop###
          #This loop runs if nbhsum is greater then 0, otherwise the previous loop runs#
          if(nbhsum >0){

            #if(mem == 0 | pre_nbhsum == 0){
            if(mem == 0){
              l <-1
              oldupdate = oldupdate+1
              repeat{
                Prob_org[[l]] <-c(nbhlist[[l]]/nbhsum)
                Prob_org[[l]][is.na(Prob_org[[l]])] <- 0
                l=l+1
                if (l ==limit)
                  break
              }
            }
            #Modified code below to include if statement "mem == 1" so that this is skipped if mem is set to 0#
            # if(mem == 1){
            #   if(pre_nbhsum == 0){
            #     l <-1
            #     oldupdate = oldupdate+1
            #     repeat{
            #       Prob_org[[l]] <-c(nbhlist[[l]]/nbhsum)
            #       Prob_org[[l]][is.na(Prob_org[[l]])] <- 0
            #       l=l+1
            #       if (l ==limit)
            #         break
            #     }
            #   }
            # }

            #Probability of being in an organization based off the prevelence of the organization in the neighborhood and the prevelence during the last model#
            #loops#
            # if(mem == 1){
            #   #if(mem == 1 & pre_nbhsum > 0){
            #   if(pre_nbhsum > 0){
            #     l <-1
            #     newupdate = newupdate+1
            #     repeat{
            #       Prob_org[[l]] <-c((nbhlist[[l]]+(pre_nbhlist[[l]]*w))/nbhsum) #Check and make sure this is correct
            #       #Prob_org[[l]] <-c((nbhlist[[l]]*((pre_nbhlist[[l]]/pre_nbhsum)*w))/(nbhsum))
            #       Prob_org[[l]][is.na(Prob_org[[l]])] <- 0
            #       l=l+1
            #       if (l ==limit)
            #         break
            #     }
            #   }
            # }

            Probmax <- max(unlist(Prob_org))
            Probmax[Probmax < 0] <- 0

            if(Probmax > 1){
              l <- 1
              repeat{
                Prob_org[[l]] = Prob_org[[l]]/Probmax
                #Prob_org[[l]] = Prob_org[[l]]/focsum
                l=l+1
                if (l ==limit)
                  break
              }
            }
            #assign('Prob_org', Prob_org, envir = .GlobalEnv)
            #assign('Probmax', Probmax, envir = .GlobalEnv)

            probtest<- unlist(Prob_org)

            #List creation and paramter reset#
            stolist <- list()

            e <- 0
            l <- 1

            ####Stolcastic Rounding Nested Repeat Fuction####
            #Part 1 - Round#
            repeat{
              repeat{
                stol = (focsum*Prob_org[[l]])

                q <- abs(stol - trunc(stol))
                #Adjusted rounding code for scalable cell neighborhood#
                probrang <- c(1 - q, q)

                adj <- sample(0:1, size = 1, replace = FALSE, prob = probrang)

                if(stol <= 0) {adj = 0}

                stol = trunc(stol)+adj
                stolist[[l]] <-c(stol)
                l = l+1
                if(l==limit)
                  break
                #}
              }
              #assign('probrang', probrang, envir = .GlobalEnv)
              #assign('adj', adj, envir = .GlobalEnv)

              #Part 2 - Check round and direct repeat(if criteria not met)#
              stolsum <- Reduce("+",stolist)
              if(stolsum == focsum) {e = e + 1}
              else{
                if(stolsum > focsum){
                  repeat{
                    mintemp <- min(probtest[probtest > 0])
                    #assign('mintemp', mintemp, envir = .GlobalEnv)
                    loctemp <- match(mintemp,probtest)
                    #assign('loctemp', loctemp, envir = .GlobalEnv)
                    stolist[[loctemp]] <- stolist[[loctemp]] - 1
                    if(stolist[[loctemp]] == 0){
                      probtest[loctemp] <- 0
                    }
                    stolsum <- Reduce("+",stolist)
                    if(stolsum == focsum){
                      e = e + 1
                      break
                    }
                  }
                }

                l <- 1
              }

              if(e == 1)
                break
            }
            #assign('stolsum', stolsum, envir = .GlobalEnv)
            #assign('probtest', probtest, envir = .GlobalEnv)
            #assign('stolist', stolist, envir = .GlobalEnv)


            #Inner and total Interation counts, and printing count numbers #
            #TR = TR+1
            #IT = IT+1
            #print(paste0("IT",IT))
            #print(paste0("TR",TR))

            #Changes Lists#
            ctempval <-list()

            l = 1

            repeat{
              if (nDim == 3) {ctempval[[l]] <- c(stolist[[l]]-orgs[[l]][cord_y,cord_x,cord_z])}
              else {ctempval[[l]] <- c(stolist[[l]]-orgs[[l]][cord_y,cord_x])} #For a 2D ecology#

              if(ctempval[[l]] > 0) {changpos[[l]] = ctempval[[l]]}
              else {changpos[[l]] = 0}

              if(ctempval[[l]] < 0) {changneg[[l]] = ctempval[[l]]}
              else {changneg[[l]] = 0}

              l=l+1
              if(l==limit)
                break
            }
            #assign('ctempval', ctempval, envir = .GlobalEnv)

            #Updating cells#
            l<-1
            repeat{
              if (nDim == 3) {orgs[[l]][cord_y,cord_x,cord_z] <- stolist[[l]]}
              else {orgs[[l]][cord_y,cord_x] <- stolist[[l]]} #For a 2D ecology#
              orgs[[l]][is.na(orgs[[l]])] <- 0
              l=l+1
              if(l ==limit)
                break
            }
            #assign('orgs', orgs, envir = .GlobalEnv)

            #Sum of orgs and save of coordinates#
            tempsum <- list()

            i <- 1
            repeat{
              tempsum[[i]] <- sum(orgs[[i]])
              i = i+1

              if(i == entities+1)
                break
            }
            #assign('tempsum', tempsum, envir = .GlobalEnv)

            sumorgs[[TR]] <- c(tempsum)

            temptot <- Reduce("+", orgs)
            totpops[[TR]] <-Reduce("+",temptot)
            #assign('sumorgs', sumorgs, envir = .GlobalEnv)
            #assign('temptot', temptot, envir = .GlobalEnv)
            #assign('totpops', totpops, envir = .GlobalEnv)

            #Save-out of coordinates#
            #This saves out coordinate for locations where changes were made in the model#
            Changepos[[TR]] <- changpos
            Changeneg[[TR]] <- changneg
            Change[[TR]] <- unlist(changpos) + unlist(changneg)

            tchange = Change[[TR]]

            temporgscoord <- list()

            if(any(tchange > 0)){
              a = a+1
              if (nDim == 2) {Coord[[a]] <- c(cord_x,cord_y)} #For a 2D ecology#
              else {Coord[[a]] <- c(cord_x,cord_y,cord_z)}
            }

            e <- 1

            for(e in 1:entities){
              if(Change[[TR]][e] != 0){
                #c[[i]] = c[[i]]+1 #I honestly don't know what this code is for#
                if (nDim == 3) {temporgscoord[[e]] <- c(cord_y,cord_x,cord_z)}
                else {temporgscoord[[e]] <- c(cord_y,cord_x)} #For a 2D ecology#
              }
              else{
                if (nDim == 3) {temporgscoord[[e]] <- c(0,0,0)}
                else {temporgscoord[[e]] <- c(0,0)} #For a 2D ecology#
              }
            }
            orgscoord[[TR]] <- temporgscoord

            totchange = Reduce('+',Change)
            #break
          }
          #assign('totchange', totchange, envir = .GlobalEnv)
          #assign('orgscoord', orgscoord, envir = .GlobalEnv)
          #assign('Change', Change, envir = .GlobalEnv)
          #assign('Coord', Coord, envir = .GlobalEnv)
          #assign('tchange', tchange, envir = .GlobalEnv)
          #assign('Changepos', Changepos, envir = .GlobalEnv)
          #assign('Changeneg', Changeneg, envir = .GlobalEnv)

          ####Convergence Loop####
          if(IT == chec){
            Set = Set + 1
            print(paste0("Set",Set))

            Runsum = unlist(Changepos[[TR]])
            thisrun = Reduce('+',Runsum)
            if(Set > 2){
              rundiff = abs(thisrun - lastrun)
            }
            legrun <- Set
            #legrun[Set] = thisrun
            lastrun = thisrun
          }

          ###Print out of differences between run sets (the chec object). Useful for checking if changes are being made###
          ###If no changes are made, there is a problem with the intital data or code###
          #print(paste0("rundiff"," ",rundiff))

          if(TR %% 10==0){ #was TR#
            CR = CR+1
            orglist[[CR]] <- orgs#[[i]]
          }
          #assign('orglist', orglist, envir = .GlobalEnv)
          #assign('CR', CR, envir = .GlobalEnv)
          #IT <- 1000

          ####Convergence Code####
          #This code is done in multiple parts, which consolidate the ecology by each organization and then compare it to the tol.#
          if(IT == chec){
            temptest <- list()
            aveorgtemp <- list()
            aveorgtemp2 <- list()
            l <- 1
            #leng <- 1

            repeat{
              #Saves the length of the orglist vector#
              lenglim <- length(orglist)
              #Creates a starting values that is an adjustment of the lenght of the orglist list so that only the last n permutations are considered#
              # N here is convr#
              leng <- (lenglim - convr)
              #if(lenglim <= convr){
              # leng <- lenglim}

              #if(lenglim >= convr){
              #leng <- lenglim - convr}

              #Create a list with all permuations of the ecology that are within the specified convergence nrange#
              repeat{
                aveorgtemp <- c(orglist[leng])
                aveorgtemp2[[lenglim - leng]] <- aveorgtemp[[1]][[l]]
                leng = leng+1
                if(leng==lenglim)#+1)
                  break
              }
              #Averaging the ecology by each cell. Here an ecology is created where the number of resources drawn from a cell are the
              #average across the convergence nrange#
              #Note, each organziation gets its own ecology here, so the next section does comparisons cross all orgs#
              #What was this code doing - G#
              aveorgtemp3 <- Reduce(`+`, aveorgtemp2) / length(aveorgtemp2)
              temptest[[l]] <- aveorgtemp3

              l=l+1
              if(l==limit)
                break
            }
            #assign('temptest', temptest, envir = .GlobalEnv)
            #assign('aveorgtemp3', aveorgtemp3, envir = .GlobalEnv)
            #assign('aveorgtemp2', aveorgtemp2, envir = .GlobalEnv)
            #assign('aveorgtemp', aveorgtemp, envir = .GlobalEnv)
            #assign('leng', leng, envir = .GlobalEnv)
            #assign('lenglim', lenglim, envir = .GlobalEnv)
            #}

            #Var creation for the convergence logic#
            converg <- list()
            convergcount <- list()

            l <- 1
            convnum <- convnum + 1

            #This loop creates the object for comparison for convergence by taking the temptest list, which has the averages from the
            #Convergence nrange being considered and creates a vectors of times that the averages are not nearly the same as the values in
            #The orgs list (where permuted results are saved). A count of mismatches is created and then compared to the tol specified by the user#
            repeat{
              converg[[l]] <- all.equal.raw(temptest[[l]],orgs[[l]],tol = 0, check.attributes = TRUE)
              convergcount[[l]] <- stringr::str_remove(converg[[l]]," element mismatch")
              convergcount[[l]] <- stringr::str_remove(convergcount[[l]],"es")
              convergcount[[l]][convergcount[[l]] == "TRUE"] <- 0
              convergcount[[l]] <- as.numeric(convergcount[[l]])
              conversave[[convnum]] <- convergcount
              l=l+1
              if(l==limit)
                break
            }

            pass <- 0

            #Comparison of the mismatches to the tol. This loops does the comparison and also counts the number of passes
            # or times that the mismatch is lower or equal to the tol. If the count of passed entities equals entities, then the model is considered convergenced.#
            l <- 1
            repeat{
              if(convergcount[[l]] <= tol) {pass = pass+1}
              else {pass = pass}
              l=l+1
              if(l==limit)
                break
            }
            #Number of times
            print(paste0("Number of orgs converged ",pass))

          }
          if(IT == chec){
            break
            #IT <- 0
            #break
          }
        }
        #assign('converg', converg, envir = .GlobalEnv)
        #assign('conversave', conversave, envir = .GlobalEnv)
        #assign('convergcount', convergcount, envir = .GlobalEnv)
        #assign('convnum', convnum, envir = .GlobalEnv)
        #Logic code for passing or failing convergence check#
        if(pass==entities|TR>=maxr){
          if(pass==entities){
            print("The model has converged!")
          }
          if(TR>=maxr){
            print("Model is close to reaching a pre-defined limit for iterations and has not converged. Simulation had stopped")
            simerror = simerror + 1
            errorcount = errorcount + 1
          }
          break
        }
      }
    }

    if(simerror >= 1 & simerror <= 3){
      #Save Objects List Creation#
      orglist <- list() #Parallel arrays for each entity are saved here and updated after each simulation run#
      sumorgs <- list() #Saves the sum of each entity for the run#
      totpops <- list() #The total population of the ecology for each run#
      simerror = 0
      errorcount = 0
      #adj<-0
      Coord = list() #List of coordinaes where any change happens#
      orgscoord <- list() #List of coordinates where changes happen for each entity in the ecology#
      Change = list() #List of changes made to each entity in the ecology#
      Changepos = list() #Running list of positive changes made to each entity for each run of the simulation#
      changpos = list() #List of positive changes made to each entities#
      Changeneg = list() #Running list of negative changes made to each entity for each run of the simulation#
      changneg = list() #List of negative changes made to entities#
      #assign('sumorgs', sumorgs, envir = .GlobalEnv)
      #assign('totpops', totpops, envir = .GlobalEnv)
      #assign('Coord', Coord, envir = .GlobalEnv)
      #assign('orgscoord', orgscoord, envir = .GlobalEnv)
      #assign('Change', Change, envir = .GlobalEnv)
      #assign('Changepos', Changepos, envir = .GlobalEnv)
      #assign('changpos', changpos, envir = .GlobalEnv)
      #assign('Changeneg', Changeneg, envir = .GlobalEnv)
      #assign('changneg', changneg, envir = .GlobalEnv)

      skiptrip = 0 #Count of times that skip logic is used#
      #assign('skiptrip', skiptrip, envir = .GlobalEnv)
      #totchange <-0
      legrun = 0 #The number of sets that have been run for the simulation run
      #assign('legrun', legrun, envir = .GlobalEnv)
      #Values for testing the update loops (mem and old)#
      oldupdate <- 0
      newupdate <- 0
      #assign('oldupdate', oldupdate, envir = .GlobalEnv)
      #assign('newupdate', newupdate, envir = .GlobalEnv)
      #Testing Resets#
      #orgs <- OGorgs
      #orglist <- list()

      #Simulatin mode code#
      #Setting of repeatition parameters for main loop(?)#
      Set <- 0 #Starting value for number of sets run#
      rundiff <- 0 #Creation of object for difference between runs, updated after each simulation set (each run of the simulation in sets of 100)#
      lowdiff <- 0 #Creation of object for counting number of times convergence criteria is reached#
      TR <- 0 #Object for counting total number of iterations the simulation carries out (times it selects a new focal cell)# - changed from 0 to 1
      CR <- 0
      convnum <- 0
      pass <- 0
      conversave <- list()

      IT <-0
      i <- as.numeric(1)
      a <-0
      p <-0


      ####Main Loop####
      repeat{
        IT <-0
        #i <- as.numeric(1)
        #a <-0
        #p <-0

        ####Inner Loop####
        repeat{

          #Inner and total Interation counts, and printing count numbers #
          TR = TR+1
          IT = IT+1
          #assign('TR', TR, envir = .GlobalEnv)
          #assign('IT', IT, envir = .GlobalEnv)
          print(paste0("IT",IT))
          print(paste0("TR",TR))

          ###Selection Loop###
          repeat{
            #dims and focells#
            if (nDim == 3) {mod_dim <- dim(orgs[[1]][,,])
            #assign('mod_dim', mod_dim, envir = .GlobalEnv)
            }
            else {mod_dim <- dim(orgs[[1]][,])
            #assign('mod_dim', mod_dim, envir = .GlobalEnv)
            }

            cord_x <- sample(mod_dim[2],1)
            cord_y <- sample(mod_dim[1],1)
            #assign('cord_x', cord_x, envir = .GlobalEnv)
            #assign('cord_y', cord_y, envir = .GlobalEnv)
            if (nDim == 3) {cord_z <- sample(mod_dim[3],1)
            #assign('cord_z', cord_z, envir = .GlobalEnv)
            }


            #New neighborhood selection code#
            #x nrange designation#

            if (cord_x > nrange & cord_x <= mod_dim[2]) {nbhcord_x_neg <- (cord_x-nrange)
            }else {nbhcord_x_neg = 1}

            if (cord_x >= 1 + nrange & cord_x < mod_dim[2]-nrange) {nbhcord_x_pos <- (cord_x+nrange)
            }else {nbhcord_x_pos = mod_dim[2]}

            nbhcord_x_nrange <- (nbhcord_x_neg:nbhcord_x_pos)
            #assign('nbhcord_x_neg', nbhcord_x_neg, envir = .GlobalEnv)
            #assign('nbhcord_x_pos', nbhcord_x_pos, envir = .GlobalEnv)
            #assign('nbhcord_x_nrange', nbhcord_x_nrange, envir = .GlobalEnv)

            #y nrange designation#
            if (cord_y > nrange & cord_y <= mod_dim[1]) {nbhcord_y_neg <- (cord_y-nrange)
            }else {nbhcord_y_neg = 1}

            if (cord_y >= 1 + nrange & cord_y < mod_dim[1]-nrange) {nbhcord_y_pos <- (cord_y+nrange)
            }else {nbhcord_y_pos = mod_dim[1]}

            nbhcord_y_nrange <- (nbhcord_y_neg:nbhcord_y_pos)
            #assign('nbhcord_y_neg', nbhcord_y_neg, envir = .GlobalEnv)
            #assign('nbhcord_y_pos', nbhcord_y_pos, envir = .GlobalEnv)
            #assign('nbhcord_y_nrange', nbhcord_y_nrange, envir = .GlobalEnv)

            #z nrange designation#
            if (nDim == 3) {if (cord_z > nrange & cord_z <= mod_dim[3]) {nbhcord_z_neg <- (cord_z-nrange)
            }else {nbhcord_z_neg = 1}

              if (cord_z >= 1 + nrange & cord_z < mod_dim[3]-nrange) {nbhcord_z_pos <- (cord_z+nrange)
              }else {nbhcord_z_pos = mod_dim[3]}

              nbhcord_z_nrange <- (nbhcord_z_neg:nbhcord_z_pos)
              #assign('nbhcord_z_neg', nbhcord_z_neg, envir = .GlobalEnv)
              #assign('nbhcord_z_pos', nbhcord_z_pos, envir = .GlobalEnv)
              #assign('nbhcord_z_nrange', nbhcord_z_nrange, envir = .GlobalEnv)
              }

            #focal cell selection and neighborhood creation#
            #setup#
            focelist <- list()
            nbhlist <- list()
            Prob_org <- list()

            #Loop for subsetting array to create neighborhood#
            l <-1

            if (nDim == 3) {
              repeat{
                orgtemp <- array(as.matrix(orgs[[l]]), dim=c(mod_dim[1],mod_dim[2],mod_dim[3]), dimnames = NULL)
                focelist[[l]] <- c(orgtemp[cord_y,cord_x,cord_z])
                nbhlist[[l]] <- c(sum(orgtemp[nbhcord_y_nrange,nbhcord_x_nrange,nbhcord_z_nrange]))
                l = l+1
                if (l == limit)
                  break
              }
            }else
            {repeat{
              orgtemp <- array(as.matrix(orgs[[l]]), dim=c(mod_dim[1],mod_dim[2]), dimnames = NULL)
              focelist[[l]] <- c(orgtemp[cord_y,cord_x])
              nbhlist[[l]] <- c(sum(orgtemp[nbhcord_y_nrange,nbhcord_x_nrange]))
              l = l+1
              if (l == limit)
                break
            }
            }
            #assign('orgtemp', orgtemp, envir = .GlobalEnv)
            #assign('focelist', focelist, envir = .GlobalEnv)
            #assign('nbhlist', nbhlist, envir = .GlobalEnv)
            #Loop for creating a subset array of the same location in the ecology, but for a previous time period#
            #Note: will need the R-object for the previous array of the ecology imported#
            l <-1

            # if(mem == 1){
            #   pre_focelist <- list()
            #   pre_nbhlist <- list()
            #   if (nDim == 3) {repeat{
            #     oldorgtemp <- array(as.matrix(prev_orgs[[l]]), dim=c(mod_dim[1],mod_dim[2],mod_dim[3]), dimnames = NULL)
            #     pre_focelist[[l]] <- as.numeric(oldorgtemp[cord_y,cord_x,cord_z])
            #     pre_nbhlist[[l]] <- as.numeric(sum(oldorgtemp[nbhcord_y_nrange,nbhcord_x_nrange,nbhcord_z_nrange]))
            #     l = l+1
            #     if (l == limit)
            #       break
            #   }
            #   } else {repeat{
            #     oldorgtemp <- array(as.matrix(prev_orgs[[l]]), dim=c(mod_dim[1],mod_dim[2]), dimnames = NULL)
            #     pre_focelist[[l]] <- as.numeric(oldorgtemp[cord_y,cord_x])
            #     pre_nbhlist[[l]] <- as.numeric(sum(oldorgtemp[nbhcord_y_nrange,nbhcord_x_nrange]))
            #     l = l+1
            #     if (l == limit)
            #       break
            #   }
            #   }
            #   #assign('oldorgtemp', oldorgtemp, envir = .GlobalEnv)
            #   #assign('pre_focelist', pre_focelist, envir = .GlobalEnv)
            #   #assign('pre_nbhlist', pre_nbhlist, envir = .GlobalEnv)
            # }

            #Reduces the focelist and nbhlist objects to a sum and then saves them as a different object#
            focsum <- Reduce("+",focelist)
            nbhsum <- Reduce("+",nbhlist)
            #assign('focsum', focsum, envir = .GlobalEnv)
            #assign('nbhsum', nbhsum, envir = .GlobalEnv)

            #Reduces the pre_focelist and pre_nbhlist objects to a sum and then saves them as a different object#
            #Contains new code#
            # if(mem == 1){
            #   pre_focsum <- Reduce("+",pre_focelist)
            #   if(length(pre_focsum) == 0){
            #     pre_focsum <- 0
            #   }
            #   pre_nbhsum <- Reduce("+",pre_nbhlist)
            #   if(length(pre_nbhsum) == 0){
            #     pre_nbhsum <- 0
            #   }
            #   #assign('pre_focsum', pre_focsum, envir = .GlobalEnv)
            #   #assign('pre_nbhsum', pre_nbhsum, envir = .GlobalEnv)
            # }


            ###Skip code for when nbhsum would sum to zero###
            #This code skips the summation process for the current cell and then triggers an alternate code for saving the cell cordinates
            #of locations the model made changes#
            if(nbhsum == 0){
              skiptrip <- skiptrip+1 #Counts the number of times skips are made, mostly for trouble shooting#

              l <- 1
              repeat{
                changpos[[l]] = 0
                changneg[[l]] = 0

                l=l+1

                if(l==limit)
                  break
              }

              i <- 1

              tempsum <- list()

              repeat{
                tempsum[[i]] <- sum(orgs[[i]])
                i = i+1

                if(i == entities+1){
                  break
                }
              }

              sumorgs[[TR]] <- c(tempsum)

              temptot <- Reduce("+", orgs)
              totpops[[TR]] <-Reduce("+",temptot)

              #Save-out of coordinates of each selected cell where changes are made#
              #In this case, becasue the model makes no changes, no coordinate
              Changepos[[TR]] <- changpos
              Changeneg[[TR]] <- changneg
              Change[[TR]] <- unlist(changpos) + unlist(changneg)
              #New code for changelist logic issue#
              #if(Change[[TR]] == 0){
              #  Change[[TR]] <- 0
              #}
              tchange = Change[[TR]]

              temporgscoord <- list()

              if(any(tchange > 0)){
                a = a+1

                if (nDim == 3) {Coord[[a]] <- c(cord_x,cord_y,cord_z)
                } else {Coord[[a]] <- c(cord_x,cord_y)}
              }

              e <- 1
              for(e in 1:entities){
                if(Change[[TR]][e] != 0){
                  #c[[i]] = c[[i]]+1 #I honestly don't know what this code is for#
                  if (nDim == 3) {temporgscoord[[e]] <- c(cord_y,cord_x,cord_z)}
                  else {temporgscoord[[e]] <- c(cord_y,cord_x)} #For a 2D ecology#
                }
                else{
                  if (nDim == 3) {temporgscoord[[e]] <- c(0,0,0)}
                  else {temporgscoord[[e]] <- c(0,0)} #For a 2D ecology#
                }
              }
              orgscoord[[TR]] <- temporgscoord

              totchange = Reduce('+',Change)
              #}
              break
            }
            break
            #assign('totchange', totchange, envir = .GlobalEnv)
          }
          #assign('sumorgs', sumorgs, envir = .GlobalEnv)
          #assign('totpops', totpops, envir = .GlobalEnv)
          #assign('Changepos', Changepos, envir = .GlobalEnv)
          #assign('Changeneg', Changeneg, envir = .GlobalEnv)
          #assign('Change', Change, envir = .GlobalEnv)
          #assign('Coord', Coord, envir = .GlobalEnv)
          #assign('orgscoord', orgscoord, envir = .GlobalEnv)


          ###Simulation Loop###
          #This loop runs if nbhsum is greater then 0, otherwise the previous loop runs#
          if(nbhsum >0){

            #if(mem == 0 | pre_nbhsum == 0){
            if(mem == 0){
              l <-1
              oldupdate = oldupdate+1
              repeat{
                Prob_org[[l]] <-c(nbhlist[[l]]/nbhsum)
                Prob_org[[l]][is.na(Prob_org[[l]])] <- 0
                l=l+1
                if (l ==limit)
                  break
              }
            }
            #Modified code below to include if statement "mem == 1" so that this is skipped if mem is set to 0#
            # if(mem == 1){
            #   if(pre_nbhsum == 0){
            #     l <-1
            #     oldupdate = oldupdate+1
            #     repeat{
            #       Prob_org[[l]] <-c(nbhlist[[l]]/nbhsum)
            #       Prob_org[[l]][is.na(Prob_org[[l]])] <- 0
            #       l=l+1
            #       if (l ==limit)
            #         break
            #     }
            #   }
            # }

            #Probability of being in an organization based off the prevelence of the organization in the neighborhood and the prevelence during the last model#
            #loops#
            # if(mem == 1){
            #   #if(mem == 1 & pre_nbhsum > 0){
            #   if(pre_nbhsum > 0){
            #     l <-1
            #     newupdate = newupdate+1
            #     repeat{
            #       Prob_org[[l]] <-c((nbhlist[[l]]+(pre_nbhlist[[l]]*w))/nbhsum) #Check and make sure this is correct
            #       #Prob_org[[l]] <-c((nbhlist[[l]]*((pre_nbhlist[[l]]/pre_nbhsum)*w))/(nbhsum))
            #       Prob_org[[l]][is.na(Prob_org[[l]])] <- 0
            #       l=l+1
            #       if (l ==limit)
            #         break
            #     }
            #   }
            # }

            Probmax <- max(unlist(Prob_org))
            Probmax[Probmax < 0] <- 0

            if(Probmax > 1){
              l <- 1
              repeat{
                Prob_org[[l]] = Prob_org[[l]]/Probmax
                #Prob_org[[l]] = Prob_org[[l]]/focsum
                l=l+1
                if (l ==limit)
                  break
              }
            }
            #assign('Prob_org', Prob_org, envir = .GlobalEnv)
            #assign('Probmax', Probmax, envir = .GlobalEnv)

            probtest<- unlist(Prob_org)

            #List creation and paramter reset#
            stolist <- list()

            e <- 0
            l <- 1

            ####Stolcastic Rounding Nested Repeat Fuction####
            #Part 1 - Round#
            repeat{
              repeat{
                stol = (focsum*Prob_org[[l]])

                q <- abs(stol - trunc(stol))
                #Adjusted rounding code for scalable cell neighborhood#
                probrang <- c(1 - q, q)

                adj <- sample(0:1, size = 1, replace = FALSE, prob = probrang)

                if(stol <= 0) {adj = 0}

                stol = trunc(stol)+adj
                stolist[[l]] <-c(stol)
                l = l+1
                if(l==limit)
                  break
                #}
              }
              #assign('probrang', probrang, envir = .GlobalEnv)
              #assign('adj', adj, envir = .GlobalEnv)

              #Part 2 - Check round and direct repeat(if criteria not met)#
              stolsum <- Reduce("+",stolist)
              if(stolsum == focsum) {e = e + 1}
              else{
                if(stolsum > focsum){
                  repeat{
                    mintemp <- min(probtest[probtest > 0])
                    #assign('mintemp', mintemp, envir = .GlobalEnv)
                    loctemp <- match(mintemp,probtest)
                    #assign('loctemp', loctemp, envir = .GlobalEnv)
                    stolist[[loctemp]] <- stolist[[loctemp]] - 1
                    if(stolist[[loctemp]] == 0){
                      probtest[loctemp] <- 0
                    }
                    stolsum <- Reduce("+",stolist)
                    if(stolsum == focsum){
                      e = e + 1
                      break
                    }
                  }
                }

                l <- 1
              }

              if(e == 1)
                break
            }
            #assign('stolsum', stolsum, envir = .GlobalEnv)
            #assign('probtest', probtest, envir = .GlobalEnv)
            #assign('stolist', stolist, envir = .GlobalEnv)


            #Inner and total Interation counts, and printing count numbers #
            #TR = TR+1
            #IT = IT+1
            #print(paste0("IT",IT))
            #print(paste0("TR",TR))

            #Changes Lists#
            ctempval <-list()

            l = 1

            repeat{
              if (nDim == 3) {ctempval[[l]] <- c(stolist[[l]]-orgs[[l]][cord_y,cord_x,cord_z])}
              else {ctempval[[l]] <- c(stolist[[l]]-orgs[[l]][cord_y,cord_x])} #For a 2D ecology#

              if(ctempval[[l]] > 0) {changpos[[l]] = ctempval[[l]]}
              else {changpos[[l]] = 0}

              if(ctempval[[l]] < 0) {changneg[[l]] = ctempval[[l]]}
              else {changneg[[l]] = 0}

              l=l+1
              if(l==limit)
                break
            }
            #assign('ctempval', ctempval, envir = .GlobalEnv)

            #Updating cells#
            l<-1
            repeat{
              if (nDim == 3) {orgs[[l]][cord_y,cord_x,cord_z] <- stolist[[l]]}
              else {orgs[[l]][cord_y,cord_x] <- stolist[[l]]} #For a 2D ecology#
              orgs[[l]][is.na(orgs[[l]])] <- 0
              l=l+1
              if(l ==limit)
                break
            }
            #assign('orgs', orgs, envir = .GlobalEnv)

            #Sum of orgs and save of coordinates#
            tempsum <- list()

            i <- 1
            repeat{
              tempsum[[i]] <- sum(orgs[[i]])
              i = i+1

              if(i == entities+1)
                break
            }
            #assign('tempsum', tempsum, envir = .GlobalEnv)

            sumorgs[[TR]] <- c(tempsum)

            temptot <- Reduce("+", orgs)
            totpops[[TR]] <-Reduce("+",temptot)
            #assign('sumorgs', sumorgs, envir = .GlobalEnv)
            #assign('temptot', temptot, envir = .GlobalEnv)
            #assign('totpops', totpops, envir = .GlobalEnv)

            #Save-out of coordinates#
            #This saves out coordinate for locations where changes were made in the model#
            Changepos[[TR]] <- changpos
            Changeneg[[TR]] <- changneg
            Change[[TR]] <- unlist(changpos) + unlist(changneg)

            tchange = Change[[TR]]

            temporgscoord <- list()

            if(any(tchange > 0)){
              a = a+1
              if (nDim == 2) {Coord[[a]] <- c(cord_x,cord_y)} #For a 2D ecology#
              else {Coord[[a]] <- c(cord_x,cord_y,cord_z)}
            }

            e <- 1

            for(e in 1:entities){
              if(Change[[TR]][e] != 0){
                #c[[i]] = c[[i]]+1 #I honestly don't know what this code is for#
                if (nDim == 3) {temporgscoord[[e]] <- c(cord_y,cord_x,cord_z)}
                else {temporgscoord[[e]] <- c(cord_y,cord_x)} #For a 2D ecology#
              }
              else{
                if (nDim == 3) {temporgscoord[[e]] <- c(0,0,0)}
                else {temporgscoord[[e]] <- c(0,0)} #For a 2D ecology#
              }
            }
            orgscoord[[TR]] <- temporgscoord

            totchange = Reduce('+',Change)
            #break
          }
          #assign('totchange', totchange, envir = .GlobalEnv)
          #assign('orgscoord', orgscoord, envir = .GlobalEnv)
          #assign('Change', Change, envir = .GlobalEnv)
          #assign('Coord', Coord, envir = .GlobalEnv)
          #assign('tchange', tchange, envir = .GlobalEnv)
          #assign('Changepos', Changepos, envir = .GlobalEnv)
          #assign('Changeneg', Changeneg, envir = .GlobalEnv)

          ####Convergence Loop####
          if(IT == chec){
            Set = Set + 1
            print(paste0("Set",Set))

            Runsum = unlist(Changepos[[TR]])
            thisrun = Reduce('+',Runsum)
            if(Set > 2){
              rundiff = abs(thisrun - lastrun)
            }
            legrun <- Set
            #legrun[Set] = thisrun
            lastrun = thisrun
          }

          ###Print out of differences between run sets (the chec object). Useful for checking if changes are being made###
          ###If no changes are made, there is a problem with the intital data or code###
          #print(paste0("rundiff"," ",rundiff))

          if(TR %% 10==0){ #was TR#
            CR = CR+1
            orglist[[CR]] <- orgs#[[i]]
          }
          #assign('orglist', orglist, envir = .GlobalEnv)
          #assign('CR', CR, envir = .GlobalEnv)
          #IT <- 1000

          ####Convergence Code####
          #This code is done in multiple parts, which consolidate the ecology by each organization and then compare it to the tol.#
          if(IT == chec){
            temptest <- list()
            aveorgtemp <- list()
            aveorgtemp2 <- list()
            l <- 1
            #leng <- 1

            repeat{
              #Saves the length of the orglist vector#
              lenglim <- length(orglist)
              #Creates a starting values that is an adjustment of the lenght of the orglist list so that only the last n permutations are considered#
              # N here is convr#
              leng <- (lenglim - convr)
              #if(lenglim <= convr){
              # leng <- lenglim}

              #if(lenglim >= convr){
              #leng <- lenglim - convr}

              #Create a list with all permuations of the ecology that are within the specified convergence nrange#
              repeat{
                aveorgtemp <- c(orglist[leng])
                aveorgtemp2[[lenglim - leng]] <- aveorgtemp[[1]][[l]]
                leng = leng+1
                if(leng==lenglim)#+1)
                  break
              }
              #Averaging the ecology by each cell. Here an ecology is created where the number of resources drawn from a cell are the
              #average across the convergence nrange#
              #Note, each organziation gets its own ecology here, so the next section does comparisons cross all orgs#
              #What was this code doing - G#
              aveorgtemp3 <- Reduce(`+`, aveorgtemp2) / length(aveorgtemp2)
              temptest[[l]] <- aveorgtemp3

              l=l+1
              if(l==limit)
                break
            }
            #assign('temptest', temptest, envir = .GlobalEnv)
            #assign('aveorgtemp3', aveorgtemp3, envir = .GlobalEnv)
            #assign('aveorgtemp2', aveorgtemp2, envir = .GlobalEnv)
            #assign('aveorgtemp', aveorgtemp, envir = .GlobalEnv)
            #assign('leng', leng, envir = .GlobalEnv)
            #assign('lenglim', lenglim, envir = .GlobalEnv)
            #}

            #Var creation for the convergence logic#
            converg <- list()
            convergcount <- list()

            l <- 1
            convnum <- convnum + 1

            #This loop creates the object for comparison for convergence by taking the temptest list, which has the averages from the
            #Convergence nrange being considered and creates a vectors of times that the averages are not nearly the same as the values in
            #The orgs list (where permuted results are saved). A count of mismatches is created and then compared to the tol specified by the user#
            repeat{
              converg[[l]] <- all.equal.raw(temptest[[l]],orgs[[l]],tol = 0, check.attributes = TRUE)
              convergcount[[l]] <- stringr::str_remove(converg[[l]]," element mismatch")
              convergcount[[l]] <- stringr::str_remove(convergcount[[l]],"es")
              convergcount[[l]][convergcount[[l]] == "TRUE"] <- 0
              convergcount[[l]] <- as.numeric(convergcount[[l]])
              conversave[[convnum]] <- convergcount
              l=l+1
              if(l==limit)
                break
            }

            pass <- 0

            #Comparison of the mismatches to the tol. This loops does the comparison and also counts the number of passes
            # or times that the mismatch is lower or equal to the tol. If the count of passed entities equals entities, then the model is considered convergenced.#
            l <- 1
            repeat{
              if(convergcount[[l]] <= tol) {pass = pass+1}
              else {pass = pass}
              l=l+1
              if(l==limit)
                break
            }
            #Number of times
            print(paste0("Number of orgs converged ",pass))

          }
          if(IT == chec){
            break
            #IT <- 0
            #break
          }
        }
        #assign('converg', converg, envir = .GlobalEnv)
        #assign('conversave', conversave, envir = .GlobalEnv)
        #assign('convergcount', convergcount, envir = .GlobalEnv)
        #assign('convnum', convnum, envir = .GlobalEnv)
        #Logic code for passing or failing convergence check#
        if(pass==entities|TR>=maxr){
          if(pass==entities){
            print("The model has converged!")
          }
          if(TR>=maxr){
            print("Model is close to reaching a pre-defined limit for iterations and has not converged. Simulation had stopped")
            simerror = simerror + 1
            errorcount = errorcount + 1
          }
          break
        }
      }
    }

    if(simerror >= 3){

      print("Simulation Script stuck in loop")
    }

    if(simerror < 3){

      mlc = mlc+1
      #assign('mlc', mlc, envir = .GlobalEnv)
      #Evaluation Script - Loop#
      indices <- list()

      orgs_y_min <- list()
      orgs_y_max <- list()
      orgs_lengy <- list()

      orgs_x_min <- list()
      orgs_x_max <- list()
      orgs_lengx <- list()
      if (nDim == 3) {
        orgs_z_min <- list()
        orgs_z_max <- list()
        orgs_lengz <- list()
      }
      rang_orgs <- list()

      Bivar_count_hi = list() #Blank if no bin var
      Bivar_count_lo = list() #Blank if no bin var

      #Code updates and automation - 4/30/2020#

      #Extensiveness Loops and Indices#

      #Indices Creation#
      #Indices are used as markers of where resources are in the space#
      l <-1

      repeat{
        if (nDim == 3) {
          orgtemp2 <- array(as.matrix(orgs[[l]]), dim=c(mod_dim[[1]],mod_dim[[2]],mod_dim[[3]]), dimnames = NULL)
        } else
        {
          orgtemp2 <- array(as.matrix(orgs[[l]]), dim=c(mod_dim[[1]],mod_dim[[2]]), dimnames = NULL)
        }
        indices[[l]] <- which(orgtemp2 !=0 ,arr.ind = TRUE)
        l=l+1
        if(l ==limit)
          break
      }
      #assign('indices', indices, envir = .GlobalEnv)
      #assign('orgtemp2', orgtemp2, envir = .GlobalEnv)
      #Orgs range list creation#
      #Creates list for the min and max of each dimension. If an impossible value is the result (such as division by 0), then replaces value with 0. #

      l <-1 #Steps for each org to work through indices list#
      a <-1 #Steps for each org in the below list#
      repeat{
        orgs_x_min[[a]] <- c(min(indices[[l]][,1]))
        orgs_x_min[[a]][orgs_x_min[[a]] == "Inf" | orgs_x_min[[a]] == "-Inf"] <- 0
        orgs_x_max[[a]] <- c(max(indices[[l]][,1]))
        orgs_x_max[[a]][orgs_x_max[[a]] == "Inf" | orgs_x_max[[a]] == "-Inf"] <- 0
        orgs_lengx[[a]] <- c((orgs_x_max[[a]]-orgs_x_min[[a]])+1)
        #assign('orgs_x_min', orgs_x_min, envir = .GlobalEnv)
        #assign('orgs_x_max', orgs_x_max, envir = .GlobalEnv)
        #assign('orgs_lengx', orgs_lengx, envir = .GlobalEnv)

        orgs_y_min[[a]] <- c(min(indices[[l]][,2]))
        orgs_y_min[[a]][orgs_y_min[[a]] == "Inf" | orgs_y_min[[a]] == "-Inf"] <- 0
        orgs_y_max[[a]] <- c(max(indices[[l]][,2]))
        orgs_x_max[[a]][orgs_y_max[[a]] == "Inf" | orgs_y_max[[a]] == "-Inf"] <- 0
        orgs_lengy[[a]] <- c((orgs_y_max[[a]]-orgs_y_min[[a]])+1)
        #assign('orgs_y_min', orgs_y_min, envir = .GlobalEnv)
        #assign('orgs_y_max', orgs_y_max, envir = .GlobalEnv)
        #assign('orgs_lengy', orgs_lengy, envir = .GlobalEnv)

        #rang_orgs[[a]] <- c(orgs_lengx[[a]]*orgs_lengy[[a]]) #For a 2D ecology#
        if (nDim == 3) {
          orgs_z_min[[a]] <- c(min(indices[[l]][,3]))
          orgs_z_min[[a]][orgs_z_min[[a]] == "Inf" | orgs_z_min[[a]] == "-Inf"] <- 0
          orgs_z_max[[a]] <- c(max(indices[[l]][,3]))
          orgs_z_max[[a]][orgs_z_max[[a]] == "Inf" | orgs_z_max[[a]] == "-Inf"] <- 0
          orgs_lengz[[a]] <- c((orgs_z_max[[a]]-orgs_z_min[[a]])+1)
          rang_orgs[[a]] <- c(orgs_lengx[[a]]*orgs_lengy[[a]]*orgs_lengz[[a]])
          #assign('orgs_z_min', orgs_z_min, envir = .GlobalEnv)
          #assign('orgs_z_max', orgs_z_max, envir = .GlobalEnv)
          #assign('orgs_lengz', orgs_lengz, envir = .GlobalEnv)
        } else {
          rang_orgs[[a]] <- c(orgs_lengx[[a]]*orgs_lengy[[a]])
        }
        #assign('rang_orgs', rang_orgs, envir = .GlobalEnv)

        a=a+1
        l=l+1
        if(l ==limit)
          break
      }

      #Total Eco list creation#
      totaeco <- Reduce('+',orgs)
      ecoindices<- which(totaeco!=0, arr.ind = TRUE)
      eco_y_min <- min(ecoindices[,1],na.rm = TRUE)
      eco_y_max <- max(ecoindices[,1],na.rm = TRUE)
      eco_y <- c(eco_y_min,eco_y_max)
      eco_y_leng <- (eco_y_max-eco_y_min)+1
      #assign('totaeco', totaeco, envir = .GlobalEnv)
      #assign('ecoindices', ecoindices, envir = .GlobalEnv)
      #assign('eco_y_min', eco_y_min, envir = .GlobalEnv)
      #assign('eco_y_max', eco_y_max, envir = .GlobalEnv)
      #assign('eco_y', eco_y, envir = .GlobalEnv)
      #assign('eco_y_leng', eco_y_leng, envir = .GlobalEnv)

      eco_x_min <- min(ecoindices[,2],na.rm = TRUE)
      eco_x_max <- max(ecoindices[,2],na.rm = TRUE)
      eco_x <- c(eco_x_min,eco_x_max)
      eco_x_leng <- (eco_x_max-eco_x_min)+1
      #assign('eco_x_min', eco_x_min, envir = .GlobalEnv)
      #assign('eco_x_max', eco_x_max, envir = .GlobalEnv)
      #assign('eco_x', eco_x, envir = .GlobalEnv)
      #assign('eco_x_leng', eco_x_leng, envir = .GlobalEnv)

      if (nDim == 2) { ecorange <- eco_y_leng*eco_x_leng
      } else {
        eco_z_min <- min(ecoindices[,3],na.rm = TRUE)
        eco_z_max <- max(ecoindices[,3],na.rm = TRUE)
        eco_z <- c(eco_z_min,eco_z_max)
        eco_z_leng <- (eco_z_max-eco_z_min)+1
        #assign('eco_z_min', eco_z_min, envir = .GlobalEnv)
        #assign('eco_z_max', eco_z_max, envir = .GlobalEnv)
        #assign('eco_z', eco_z, envir = .GlobalEnv)
        #assign('eco_z_leng', eco_z_leng, envir = .GlobalEnv)

        ecorange <- eco_y_leng*eco_x_leng*eco_z_leng
      }
      #assign('ecorange', ecorange, envir = .GlobalEnv)
      #Extensiveness Calculation#
      #Setup of ext list and reset of limit#
      l <-1

      ext_org <- list()

      repeat{
        ext_org[[l]] <- rang_orgs[[l]]/ecorange

        l=l+1
        if(l ==limit)
          break
      }
      #assign('ext_org', ext_org, envir = .GlobalEnv)

      #Extensiveness Loops#
      #Reset of limit and creation of a repetition limit object#
      l <- 1
      replim <- list()

      repeat{
        if (nDim == 3) {
          replim[[l]] <- c(length(indices[[l]])/3)
        }
        if (nDim == 2) {
          replim[[l]] <- c(length(indices[[l]])/2) #For a 2D ecology#
        }
        l = l+1
        if(l==limit)
          break
      }
      #assign('replim', replim, envir = .GlobalEnv)

      #Reset of Parameters and Intensiveness list creation#
      l <- 1
      r <- 1
      cjk <- list()
      xjk <- list()
      top <- list()
      caltop <- list()

      #Calculation of cjk, xjk, and compression of top of equation for later calculations#
      repeat{
        cords <- indices[[l]]
        if(length(cords) != 0){
          if (nDim == 3) {
            for(r in 1:replim[[l]]){
              cjk[[r]] <- c(totaeco[cords[r,1],cords[r,2],cords[r,3]])
              xjk[[r]] <- c(orgs[[l]][cords[r,1],cords[r,2],cords[r,3]])
              #print(r)
            }
          }
          if (nDim == 2) {
            #For a 2D ecology#
            for(r in 1:replim[[l]]){
              cjk[[r]] <- c(totaeco[cords[r,1],cords[r,2]])
              xjk[[r]] <- c(orgs[[l]][cords[r,1],cords[r,2]])
              #print(l)
            }
          }
          if(r==replim[[l]]){
            for(r in 1:replim[[l]]){
              caltop[[r]] <- c(xjk[[r]]/cjk[[r]])
            }

            top[[l]] <- Reduce('+',caltop)
            cjk <- list()
            xjk <- list()
            caltop <- list()
            l=l+1}
          else{l=l}

          if(l==limit)
            break
        }
        else {top[[l]] <- 1
        l=l+1
        if(l==limit)
          break
        }
      }
      #assign('cjk', cjk, envir = .GlobalEnv)
      #assign('xjk', xjk, envir = .GlobalEnv)
      #assign('top', top, envir = .GlobalEnv)
      #assign('caltop', caltop, envir = .GlobalEnv)

      #Intensiveness Calculation range#
      l <- 1
      inten_rang <- list()
      repeat{
        inten_rang[[l]] <- c(top[[l]]/rang_orgs[[l]])
        l = l+1
        if(l==limit)
          break
      }
      #assign('inten_rang', inten_rang, envir = .GlobalEnv)

      #Intensiveness Calculation Cell#
      l <- 1
      inten_cell <- list()
      repeat{
        inten_cell[[l]] <- c(top[[l]]/replim[[l]])
        l = l+1
        if(l==limit)
          break
      }
      #assign('inten_cell', inten_cell, envir = .GlobalEnv)

      #bin summary code#
      l <- 1
      if(bin == 1){
        repeat{
          if (nDim == 3) {
            Bivar_count_lo[l] <- c(sum(orgs[[l]][,,1]))
            Bivar_count_hi[l] <- c(sum(orgs[[l]][,,2]))
          }
          if (nDim == 2) {
            Bivar_count_lo[l] <- c(sum(orgs[[l]][,1]))
            Bivar_count_hi[l] <- c(sum(orgs[[l]][,2]))
          }
          l = l+1
          if(l==limit)
            break
        }
      }
      #assign('Bivar_count_lo', Bivar_count_lo, envir = .GlobalEnv)
      #assign('Bivar_count_hi', Bivar_count_hi, envir = .GlobalEnv)

      #Export Script#
      #setwd(coordpath)

      #totiter <- TR-1
      totiter <- TR

      e <- 1
      t <- 1

      repeat{
        tempcoordloc <- list()
        orgtitle <- paste0("orgcoords ",e)
        for(t in 1:totiter){
          tempcoordloc[t] <- c(orgscoord[[t]][e])
        }
        write.csv(tempcoordloc, file = paste0(prefix," run ",mlc," ",orgtitle, ".csv"))
        e <- e+1
        t <- 1
        if(e == entities+1){
          break
        }
      }
      #assign('tempcoordloc', tempcoordloc, envir = .GlobalEnv)
      #assign('e', e, envir = .GlobalEnv)
      #Creation of subfolder for each run#
      runDir <- paste0("/run",mlc)
      runpath <- paste0(saveDir, runDir)
      dir.create(file.path(runpath))

      #setwd(runpath)

      write.csv(inten_rang, file = paste0(prefix," range intensiveness.csv"))
      write.csv(inten_cell, file = paste0(prefix," cell intensiveness.csv"))
      write.csv(ext_org, file = paste0(prefix," exstensiveness.csv"))

      #sumorgs saveout#
      for(e in 1:entities){
        sumtitle <- paste0("sumorg ",e)
        write.csv(sumorgs[[TR]][e], file = paste0(prefix," ", sumtitle, ".csv"))
      }

      #orglist saveout#
      for(e in 1:entities){
        orgtitle <- paste0("orglist ",e)
        write.csv(orglist[[CR]][e], file = paste0(prefix," ", orgtitle, ".csv"))
      }

      #total population count saveout#
      write.csv(totpops,file = paste0(prefix," totpops.csv"))

      #total ecology image saveout#
      write.csv(totaeco,file = paste0(prefix," totaeco.csv"))

      #total change count saveout#
      write.csv(totchange,file = paste0(prefix," totchange.csv"))

      #General coordinate saveout#
      write.csv(Coord, file = paste0(prefix," coordinate.csv"))

      #Dimension max and min saveout#
      write.csv(orgs_x_min,file = paste0(prefix," orgs_x_min.csv"))
      write.csv(orgs_x_max,file = paste0(prefix," orgs_x_max.csv"))

      write.csv(orgs_y_min,file = paste0(prefix," orgs_y_min.csv"))
      write.csv(orgs_y_max,file = paste0(prefix," orgs_y_max.csv"))

      if (nDim == 3) {
        write.csv(orgs_z_min,file = paste0(prefix," orgs_z_min.csv"))
        write.csv(orgs_z_max,file = paste0(prefix," orgs_z_max.csv"))
      }
      #Export bin variable counts#
      if(bin == 1 ){
        write.csv(Bivar_count_hi, file = paste0(prefix," orgs_binary_max.csv"))
        write.csv(Bivar_count_lo, file = paste0(prefix," orgs_binary_min.csv"))
      }

      #Summary Script - saveout at final run#

      simDir <- paste0("/AS - After Simulation")

      simpath <- paste0(saveDir, simDir)

      dir.create(file.path(simpath))

      simDir2 <- paste0("/AS - After Simulation/means")

      simpath2 <- paste0(saveDir, simDir2)

      dir.create(file.path(simpath2))

      #setwd(simpath2)

      #Loops for averaging output values#
      #Loop for averaging extensiveness values and save out of R#
      extenhist[[mlc]] <- c(t(ext_org))

      if (mlc == glob){
        extenlist <- matrix(unlist(extenhist), ncol = entities, byrow = TRUE)
        extenlist[extenlist == "-Inf" | extenlist == "Inf"] <- 0
        exten_means <- colMeans(extenlist)
        #setwd(saveDir)
        write.csv(extenlist,file = paste0(prefix," ",year," exten.csv"))
        #setwd(simpath2)
        write.csv(exten_means,file = paste0(prefix," ",year," exten_means.csv"))
        #assign('extenlist', extenlist, envir = .GlobalEnv)
        #assign('exten_means', exten_means, envir = .GlobalEnv)
      }

      #Loop for averaging cell focused intensiveness values and save out of R#
      intencellhist[[mlc]] <- c(t(inten_cell))
      #assign('intencellhist', intencellhist, envir = .GlobalEnv)

      if (mlc == glob){
        intencelllist <- matrix(unlist(intencellhist), ncol = entities, byrow = TRUE)
        intencelllist[intencelllist == "-Inf" | intencelllist == "Inf"] <- 0
        intencell_means <- colMeans(intencelllist)
        #setwd(saveDir)
        write.csv(intencelllist,file = paste0(prefix," ",year," intencell.csv"))
        #setwd(simpath2)
        write.csv(intencell_means,file = paste0(prefix," ",year," intencell_means.csv"))
        #assign('intencelllist', intencelllist, envir = .GlobalEnv)
        #assign('intencell_means', intencell_means, envir = .GlobalEnv)
      }

      #Loop for averaging range focused intensiveness values and save out of R#
      intenranghist[[mlc]] <- c(t(inten_rang))
      #assign('intenranghist', intenranghist, envir = .GlobalEnv)

      if (mlc == glob){
        intenranglist <- matrix(unlist(intenranghist), ncol = entities, byrow = TRUE)
        intenranglist[intenranglist == "-Inf" | intenranglist == "Inf"] <- 0
        intenrang_means <- colMeans(intenranglist)
        #setwd(saveDir)
        write.csv(intenranglist,file = paste0(prefix," ",year," intenrang.csv"))
        #setwd(simpath2)
        write.csv(intenrang_means,file = paste0(prefix," ",year," intenrang_means.csv"))
        #assign('intenranglist', intenranglist, envir = .GlobalEnv)
        #assign('intenrang_means', intenrang_means, envir = .GlobalEnv)
      }


      #Loop for averaging sum of entity values and save out of R#
      sumorghist[[mlc]] <- c(sumorgs[[TR-1]])
      #assign('sumorghist', sumorghist, envir = .GlobalEnv)

      if (mlc == glob){
        sumorglist <- matrix(unlist(sumorghist), ncol = entities, byrow = TRUE)
        sumorglist[sumorglist == "-Inf" | sumorglist == "Inf"] <- 0
        sumorghist_means <- colMeans(sumorglist)
        #setwd(saveDir)
        write.csv(sumorghist,file = paste0(prefix," ",year," sumorghist.csv"))
        #setwd(simpath2)
        write.csv(sumorghist_means,file = paste0(prefix," ",year," sumorghist_means.csv"))
        #assign('sumorglist', sumorglist, envir = .GlobalEnv)
        #assign('sumorghist_means', sumorghist_means, envir = .GlobalEnv)
      }


      #Average minimun and maximun values for each dimension in the Blau space#
      #X dimension#
      orgs_x_minhist[[mlc]] <- c(t(orgs_x_min))
      #assign('orgs_x_minhist', orgs_x_minhist, envir = .GlobalEnv)

      if (mlc == glob){
        orgs_x_minlist <- matrix(unlist(orgs_x_minhist), ncol = entities, byrow = TRUE)
        orgs_x_minlist[orgs_x_minlist == "-Inf" | orgs_x_minlist == "Inf"] <- 0
        orgs_x_min_means <- colMeans(orgs_x_minlist)
        #setwd(saveDir)
        write.csv(orgs_x_minlist,file = paste0(prefix," ",year," orgs_x_min.csv"))
        #setwd(simpath2)
        write.csv(orgs_x_min_means,file = paste0(prefix," ",year," orgs_x_min_means.csv"))
        #assign('orgs_x_minlist', orgs_x_minlist, envir = .GlobalEnv)
        #assign('orgs_x_min_means', orgs_x_min_means, envir = .GlobalEnv)
      }


      orgs_x_maxhist[[mlc]] <- c(t(orgs_x_max))
      #assign('orgs_x_maxhist', orgs_x_maxhist, envir = .GlobalEnv)

      if (mlc == glob){
        orgs_x_maxlist <- matrix(unlist(orgs_x_maxhist), ncol = entities, byrow = TRUE)
        orgs_x_maxlist[orgs_x_maxlist == "-Inf" | orgs_x_maxlist == "Inf"] <- 0
        orgs_x_max_means <- colMeans(orgs_x_maxlist)
        #setwd(saveDir)
        write.csv(orgs_x_maxlist,file = paste0(prefix," ",year," orgs_x_max.csv"))
        #setwd(simpath2)
        write.csv(orgs_x_max_means,file = paste0(prefix," ",year," orgs_x_max_means.csv"))
        #assign('orgs_x_maxlist', orgs_x_maxlist, envir = .GlobalEnv)
        #assign('orgs_x_max_means', orgs_x_max_means, envir = .GlobalEnv)
      }


      #Y dimension#
      orgs_y_minhist[[mlc]] <- c(t(orgs_y_min))
      #assign('orgs_y_minhist', orgs_y_minhist, envir = .GlobalEnv)


      if (mlc == glob){
        orgs_y_minlist <- matrix(unlist(orgs_y_minhist), ncol = entities, byrow = TRUE)
        orgs_y_minlist[orgs_y_minlist == "-Inf" | orgs_y_minlist == "Inf"] <- 0
        orgs_y_min_means <- colMeans(orgs_y_minlist)
        #setwd(saveDir)
        write.csv(orgs_y_minlist,file = paste0(prefix," ",year," orgs_y_min.csv"))
        #setwd(simpath2)
        write.csv(orgs_y_min_means,file = paste0(prefix," ",year," orgs_y_min_means.csv"))
        #assign('orgs_y_minlist', orgs_y_minlist, envir = .GlobalEnv)
        #assign('orgs_y_min_means', orgs_y_min_means, envir = .GlobalEnv)
      }


      orgs_y_maxhist[[mlc]] <- c(t(orgs_y_max))
      #assign('orgs_y_maxhist', orgs_y_maxhist, envir = .GlobalEnv)

      if (mlc == glob){
        orgs_y_maxlist <- matrix(unlist(orgs_y_maxhist), ncol = entities, byrow = TRUE)
        orgs_y_maxlist[orgs_y_maxlist == "-Inf" | orgs_y_maxlist == "Inf"] <- 0
        orgs_y_max_means <- colMeans(orgs_y_maxlist)
        #setwd(saveDir)
        write.csv(orgs_y_maxlist,file = paste0(prefix," ",year," orgs_y_max.csv"))
        #setwd(simpath2)
        write.csv(orgs_y_max_means,file = paste0(prefix," ",year," orgs_y_max_means.csv"))
        #assign('orgs_y_maxlist', orgs_y_maxlist, envir = .GlobalEnv)
        #assign('orgs_y_max_means', orgs_y_max_means, envir = .GlobalEnv)
      }


      #Z dimension#
      if (nDim == 3) {
        orgs_z_minhist[[mlc]] <- c(t(orgs_z_min))
        #assign('orgs_z_minhist', orgs_z_minhist, envir = .GlobalEnv)

        if (mlc == glob){
          orgs_z_minlist <- matrix(unlist(orgs_z_minhist), ncol = entities, byrow = TRUE)
          orgs_z_minlist[orgs_z_minlist == "-Inf" | orgs_z_minlist == "Inf"] <- 0
          orgs_z_min_means <- colMeans(orgs_z_minlist)
          #setwd(saveDir)
          write.csv(orgs_z_minlist,file = paste0(prefix," ",year," orgs_z_min.csv"))
          #setwd(simpath2)
          write.csv(orgs_z_min_means,file = paste0(prefix," ",year," orgs_z_min_means.csv"))
          #assign('orgs_z_minlist', orgs_x_minlist, envir = .GlobalEnv)
          #assign('orgs_z_min_means', orgs_z_min_means, envir = .GlobalEnv)
        }


        orgs_z_maxhist[[mlc]] <- c(t(orgs_z_max))
        #assign('orgs_z_maxhist', orgs_z_maxhist, envir = .GlobalEnv)

        if (mlc == glob){
          orgs_z_maxlist <- matrix(unlist(orgs_z_maxhist), ncol = entities, byrow = TRUE)
          orgs_z_maxlist[orgs_z_maxlist == "-Inf" | orgs_z_maxlist == "Inf"] <- 0
          orgs_z_max_means <- colMeans(orgs_z_maxlist)
          #setwd(saveDir)
          write.csv(orgs_z_maxlist,file = paste0(prefix," ",year," orgs_z_max.csv"))
          #setwd(simpath2)
          write.csv(orgs_z_max_means,file = paste0(prefix," ",year," orgs_z_max_means.csv"))
          #assign('orgs_z_maxlist', orgs_z_maxlist, envir = .GlobalEnv)
          #assign('orgs_z_max_means', orgs_z_max_means, envir = .GlobalEnv)
        }

      }
      #Parameter save and output#
      locrun[[mlc]] <- c(TR-1)

      #assign('locrun', locrun, envir = .GlobalEnv)

      if (mlc == glob){
        total_runs <- mlc
        write.csv(total_runs,file = paste0(prefix," ",year," total_runs.csv"))
        write.csv(locrun,file = paste0(prefix," ",year," locrun_per_globrun.csv"))
      }

      #Carrying capacity and explotation script#



      simDir <- paste0("/AS - After Simulation")

      simpath <- paste0(saveDir, simDir)

      dir.create(file.path(simpath))

      #setwd(simpath)

      #write.csv(Bivar_count_hi, file = paste0(prefix," orgs_binary_max.csv"))

      #setwd("~/Desktop/Temp Folder")

      #year <- "1987"

      #aveorglist <- list()
      #aveorglist <- orglist
      #aveorglist <- mean(aveorglist[[1:9800]][[1]])

      #aveorglist <- lapply(orglist,mean(orglist[[]][[1]]))
      ####
      e <- 1
      i <- 1

      aveorglist_temp <- list()
      aveorglist_matrix <- list()

      repeat{
        for(e in 1:length(orglist)) {
          aveorglist_temp[[e]] <- orglist[[e]][[i]]
        }
        #assign('aveorglist_temp', aveorglist_temp, envir = .GlobalEnv)

        aveorglist_matrix[[i]] <- Reduce("+",aveorglist_temp)/length(orglist)
        #assign('aveorglist_matrix', aveorglist_matrix, envir = .GlobalEnv)

        #aveorglist_matrix[[i]] <- round(aveorglist_matrix[[i]])
        i <- i + 1
        if(i == entities + 1)
          break
        #print(e)
        #print(i)
        e <- 1
      }
      ####

      #Average simulated population means#
      sumorg_sim <- list()

      i <- 0

      repeat{
        i = i+1
        sumorg_sim[[i]] <- c(sum(aveorglist_matrix[[i]]))

        if(i == entities){
          break
        }
      }
      #assign('sumorg_sim', sumorg_sim, envir = .GlobalEnv)

      totalmatrix <- Reduce('+',orgs)
      #assign('totalmatrix', totalmatrix, envir = .GlobalEnv)

      write.csv(totalmatrix, file = paste0(prefix," ",year," ",mlc," totalmatrix - simulated.csv"))

      total <- sum(totalmatrix)
      #assign('total', total, envir = .GlobalEnv)

      resource_sim <- list()

      i <- 0

      repeat{
        i = i+1
        resource_sim[[i]] <- Reduce("+",totalmatrix[aveorglist_matrix[[i]] > 0])


        if(i == entities){
          break
        }
      }

      q <- 1

      repeat{
        if(sum(resource_sim[[q]]) == 0){
          resource_sim[[q]] <- 0
        }
        q = q+1
        if(q == entities)
          break
      }
      #assign('resource_sim', resource_sim, envir = .GlobalEnv)

      write.csv(resource_sim, file = paste0(prefix," ",year," ",mlc," resources - Simulated.csv"))

      population <- list()

      i <- 0

      repeat{
        i = i + 1
        population[[i]] <- Reduce("+",aveorglist_matrix[[i]])

        if (i == entities){
          break
        }
      }
      #assign('population', population, envir = .GlobalEnv)

      write.csv(population, file = paste0(prefix," ",year," ",mlc," population - Simulated.csv"))

      capacity <- list()
      wcapacity <- list()

      weights <- list()
      weights <- totalmatrix/ totindi
      weights[is.nan(weights)] <- 0
      weights[weights == Inf] <- 0


      i <- 0
      repeat{
        i <- i + 1
        capacity[[i]] <- totalmatrix
        capacity[[i]][aveorglist_matrix[[i]] == 0] <- 0
        wcapacity[[i]] <- capacity[[i]]*weights

        if(i == entities){
          break
        }
      }
      #assign('weights', weights, envir = .GlobalEnv)
      #assign('capacity', capacity, envir = .GlobalEnv)
      #assign('wcapacity', wcapacity, envir = .GlobalEnv)

      sim_wcarry_capacity <- list()

      i <- 0

      repeat{
        i <- i + 1
        sim_wcarry_capacity[[i]] <- sum(wcapacity[[i]])

        if(i == entities){
          break
        }
      }
      #assign('sim_wcarry_capacity', sim_wcarry_capacity, envir = .GlobalEnv)

      write.csv(sim_wcarry_capacity, file = paste0(prefix," ", year," ",mlc," wcarry capacity - simulated.csv"))

      exploitation <- list()
      exploitation_value_members <- list()
      exploitation_value_space <- list()
      exploitation_value_space_range <- list()

      i <- 0

      repeat{
        i <- i + 1
        exploitation[[i]] <- aveorglist_matrix[[i]] - wcapacity[[i]]
        exploitation_value_members[i] <- abs(Reduce("+",exploitation[[i]])/sum(totalmatrix*weights))
        exploitation_value_space[i] <- exploitation_value_members[[i]]/sum(wcapacity[[i]]>0)
        exploitation_value_space_range[i] <- abs(sum(exploitation[[i]])/sum(wcapacity[[i]]>0))

        if(i == entities){
          break
        }
      }
      #assign('exploitation', exploitation, envir = .GlobalEnv)
      #assign('exploitation_value_members', exploitation_value_members, envir = .GlobalEnv)
      #assign('exploitation_value_space', exploitation_value_space, envir = .GlobalEnv)
      #assign('exploitation_value_space_range', exploitation_value_space_range, envir = .GlobalEnv)

      write.csv(exploitation_value_members, file = paste0(prefix," ",year," ",mlc," sim Exploitation value members - Simulated.csv"))
      write.csv(exploitation_value_space, file = paste0(prefix," ",year," ",mlc," sim Exploitation value space - Simulated.csv"))
      write.csv(exploitation_value_space_range, file = paste0(prefix," ",year," ",mlc," sim Exploitation value space range - Simulated.csv"))

      exploitMhist[[mlc]] <- c(t(exploitation_value_members))
      WCChist[[mlc]] <- c(t(sim_wcarry_capacity))
      #assign('exploitMhist', exploitMhist, envir = .GlobalEnv)
      #assign('WCChist', WCChist, envir = .GlobalEnv)

      countlist <- c(mlc)
    }
    #assign('countlist', countlist, envir = .GlobalEnv)

    if(mlc == glob | simerror >= 3){

      break
      #endtime <- Sys.time
    }
  }

  ####Saveout of discribitve metrics of Simulated and Observed Ecologies####
  if(simerror < 3){
    #setwd("~/Desktop/Temp Folder")
    #year for saveout file names#
    #year <- "1987"

    OGDir<-paste0("/OG - before simulations")
    OGpath<-paste0(saveDir,OGDir)

    dir.create(OGpath)

    #setwd(OGpath)

    #Evaluation list creation#
    #This are lists that values for the new metrics get saved into#

    indices <- list()

    orgs_y_min <- list()
    orgs_y_max <- list()
    orgs_lengy <- list()

    orgs_x_min <- list()
    orgs_x_max <- list()
    orgs_lengx <- list()

    if (nDim == 3) {
      orgs_z_min <- list()
      orgs_z_max <- list()
      orgs_lengz <- list()
    }
    rang_orgs <- list()

    Bivar_count_hi = list() #Blank if no bin var
    Bivar_count_lo = list() #Blank if no bin var

    #Extensiveness Loops and Indices#

    #For instances which you want to generate discriptive metrics on intital observations#
    #dim <- dim(orgs[[1]][,])
    #orgs <- OGorgs

    #Indices Creation#
    #Indices are used as markers of where resources are in the space#
    l <-1

    repeat{
      if (nDim == 3) {
        orgtemp2 <- array(as.matrix(OGorgs[[l]]), dim=c(mod_dim[[1]],mod_dim[[2]],mod_dim[[3]]), dimnames = NULL)
      }
      if (nDim == 2) {
        orgtemp2 <- array(as.matrix(OGorgs[[l]]), dim=c(mod_dim[[1]],mod_dim[[2]]), dimnames = NULL) #For a 2D ecology#
      }
      indices[[l]] <- which(orgtemp2 !=0 ,arr.ind = TRUE)

      l=l+1
      if(l ==limit)
        break
    }
    #assign('indices', indices, envir = .GlobalEnv)
    #assign('orgtemp2', orgtemp2, envir = .GlobalEnv)
    #Orgs range list creation#
    #Creates list for the min and max of each dimension. If an impossible value is the result (such as division by 0), then replaces value with 0. #

    l <-1 #Steps for each org to work through indices list#
    a <-1 #Steps for each org in the below list#
    repeat{
      orgs_x_min[[a]] <- c(min(indices[[l]][,1]))
      orgs_x_min[[a]][orgs_x_min[[a]] == "Inf" | orgs_x_min[[a]] == "-Inf"] <- 0
      orgs_x_max[[a]] <- c(max(indices[[l]][,1]))
      orgs_x_max[[a]][orgs_x_max[[a]] == "Inf" | orgs_x_max[[a]] == "-Inf"] <- 0
      orgs_lengx[[a]] <- c((orgs_x_max[[a]]-orgs_x_min[[a]])+1)
      #assign('orgs_x_min', orgs_x_min, envir = .GlobalEnv)
      #assign('orgs_x_max', orgs_x_max, envir = .GlobalEnv)
      #assign('orgs_lengx', orgs_lengx, envir = .GlobalEnv)

      orgs_y_min[[a]] <- c(min(indices[[l]][,2]))
      orgs_y_min[[a]][orgs_y_min[[a]] == "Inf" | orgs_y_min[[a]] == "-Inf"] <- 0
      orgs_y_max[[a]] <- c(max(indices[[l]][,2]))
      orgs_x_max[[a]][orgs_y_max[[a]] == "Inf" | orgs_y_max[[a]] == "-Inf"] <- 0
      orgs_lengy[[a]] <- c((orgs_y_max[[a]]-orgs_y_min[[a]])+1)
      #assign('orgs_y_min', orgs_y_min, envir = .GlobalEnv)
      #assign('orgs_y_max', orgs_y_max, envir = .GlobalEnv)
      #assign('orgs_lengy', orgs_lengy, envir = .GlobalEnv)

      rang_orgs[[a]] <- c(orgs_lengx[[a]]*orgs_lengy[[a]])
      if (nDim == 3) {
        orgs_z_min[[a]] <- c(min(indices[[l]][,3]))
        orgs_z_min[[a]][orgs_z_min[[a]] == "Inf" | orgs_z_min[[a]] == "-Inf"] <- 0
        orgs_z_max[[a]] <- c(max(indices[[l]][,3]))
        orgs_z_max[[a]][orgs_z_max[[a]] == "Inf" | orgs_z_max[[a]] == "-Inf"] <- 0
        orgs_lengz[[a]] <- c((orgs_z_max[[a]]-orgs_z_min[[a]])+1)
        #assign('orgs_z_min', orgs_z_min, envir = .GlobalEnv)
        #assign('orgs_z_max', orgs_z_max, envir = .GlobalEnv)
        #assign('orgs_lengz', orgs_lengz, envir = .GlobalEnv)

        rang_orgs[[a]] <- c(orgs_lengx[[a]]*orgs_lengy[[a]]*orgs_lengz[[a]])
      } else {
        rang_orgs[[a]] <- c(orgs_lengx[[a]]*orgs_lengy[[a]])
      }
      #assign('rang_orgs', rang_orgs, envir = .GlobalEnv)
      a=a+1
      l=l+1
      if(l ==limit)
        break
    }

    #Total Eco list creation#
    totaeco <- Reduce('+',OGorgs)
    ecoindices<- which(totaeco!=0, arr.ind = TRUE)

    eco_y_min <- min(ecoindices[,1],na.rm = TRUE)
    eco_y_max <- max(ecoindices[,1],na.rm = TRUE)
    eco_y <- c(eco_y_min,eco_y_max)
    eco_y_leng <- (eco_y_max-eco_y_min)+1
    #assign('totaeco', totaeco, envir = .GlobalEnv)
    #assign('ecoindices', ecoindices, envir = .GlobalEnv)
    #assign('eco_y_min', eco_y_min, envir = .GlobalEnv)
    #assign('eco_y_max', eco_y_max, envir = .GlobalEnv)
    #assign('eco_y', eco_y, envir = .GlobalEnv)
    #assign('eco_y_leng', eco_y_leng, envir = .GlobalEnv)

    eco_x_min <- min(ecoindices[,2],na.rm = TRUE)
    eco_x_max <- max(ecoindices[,2],na.rm = TRUE)
    eco_x <- c(eco_x_min,eco_x_max)
    eco_x_leng <- (eco_x_max-eco_x_min)+1
    #assign('eco_x_min', eco_x_min, envir = .GlobalEnv)
    #assign('eco_x_max', eco_x_max, envir = .GlobalEnv)
    #assign('eco_x', eco_x, envir = .GlobalEnv)
    #assign('eco_x_leng', eco_x_leng, envir = .GlobalEnv)
    if (nDim == 3) {
      eco_z_min <- min(ecoindices[,3],na.rm = TRUE)
      eco_z_max <- max(ecoindices[,3],na.rm = TRUE)
      eco_z <- c(eco_z_min,eco_z_max)
      eco_z_leng <- (eco_z_max-eco_z_min)+1
      #assign('eco_z_min', eco_z_min, envir = .GlobalEnv)
      #assign('eco_z_max', eco_z_max, envir = .GlobalEnv)
      #assign('eco_z', eco_z, envir = .GlobalEnv)
      #assign('eco_z_leng', eco_z_leng, envir = .GlobalEnv)

      ecorange <- eco_y_leng*eco_x_leng*eco_z_leng
    }
    if (nDim == 2) {
      ecorange <- eco_y_leng*eco_x_leng #For a 2D ecology#
    }
    #assign('ecorange', ecorange, envir = .GlobalEnv)

    #Extensiveness Calculation#
    #Setup of ext list and reset of limit#
    l <-1

    ext_org <- list()

    repeat{
      ext_org[[l]] <- rang_orgs[[l]]/ecorange

      l=l+1
      if(l ==limit)
        break
    }
    #assign('ext_org', ext_org, envir = .GlobalEnv)
    #Extensiveness Loops#
    #Reset of limit and creation of a repetition limit object#
    l <- 1
    replim <- list()

    repeat{
      if (nDim == 3) {
        replim[[l]] <- c(length(indices[[l]])/3)
      }
      if (nDim == 2) {
        replim[[l]] <- c(length(indices[[l]])/2) #For a 2D ecology#
      }
      l = l+1
      if(l==limit)
        break
    }
    #assign('replim', replim, envir = .GlobalEnv)
    #Reset of Parameters and Intensiveness list creation#
    l <- 1
    r <- 1
    cjk <- list()
    xjk <- list()
    top <- list()
    caltop <- list()

    #Calculation of cjk, xjk, and compression of top of equation for later calculations#
    repeat{
      cords <- indices[[l]]
      if(length(cords) != 0){
        if (nDim == 3) {
          for(r in 1:replim[[l]]){
            cjk[[r]] <- c(totaeco[cords[r,1],cords[r,2],cords[r,3]])
            xjk[[r]] <- c(OGorgs[[l]][cords[r,1],cords[r,2],cords[r,3]]) #I think it's here!!! OGorgs
            #print(r)
          }
        }

        if (nDim == 2) {
          #For a 2D ecology#
          for(r in 1:replim[[l]]){
            cjk[[r]] <- c(totaeco[cords[r,1],cords[r,2]])
            xjk[[r]] <- c(OGorgs[[l]][cords[r,1],cords[r,2]]) #I think it's here!!!
            #print(r)
          }
        }

        if(r==replim[[l]]){
          for(r in 1:replim[[l]]){
            caltop[[r]] <- c(xjk[[r]]/cjk[[r]])
          }

          top[[l]] <- Reduce('+',caltop)
          cjk <- list()
          xjk <- list()
          caltop <- list()
          l=l+1}
        else{l=l}

        if(l==limit)
          break
      }
      else {top[[l]] <- 1
      l=l+1
      if(l==limit)
        break
      }
    }
    #assign('cjk', cjk, envir = .GlobalEnv)
    #assign('xjk', xjk, envir = .GlobalEnv)
    #assign('top', top, envir = .GlobalEnv)
    #assign('caltop', caltop, envir = .GlobalEnv)

    #Intensiveness Calculation range#
    inten_rang <- list()

    l <- 1

    repeat{
      inten_rang[[l]] <- c(top[[l]]/rang_orgs[[l]])
      l = l+1
      if(l==limit)
        break
    }
    #assign('inten_rang', inten_rang, envir = .GlobalEnv)

    #Intensiveness Calculation Cell#
    inten_cell <- list()

    l <- 1

    repeat{
      inten_cell[[l]] <- c(top[[l]]/replim[[l]])
      l = l+1
      if(l==limit)
        break
    }
    #assign('inten_cell', inten_cell, envir = .GlobalEnv)

    #bin summary code#
    l <- 1
    if(bin == 1){
      repeat{
        if (nDim == 3) {
          Bivar_count_lo[l] <- c(sum(OGorgs[[l]][,,1]))
          Bivar_count_hi[l] <- c(sum(OGorgs[[l]][,,2]))
        }
        if (nDim == 2) {
          Bivar_count_lo[l] <- c(sum(OGorgs[[l]][,1]))
          Bivar_count_hi[l] <- c(sum(OGorgs[[l]][,2]))
        }
        l = l+1
        if(l==limit)
          break
      }
    }
    #assign('Bivar_count_lo', Bivar_count_lo, envir = .GlobalEnv)
    #assign('Bivar_count_hi', Bivar_count_hi, envir = .GlobalEnv)

    #Sum of orgs population codes#
    OGorgsum <- list()

    i <- 1

    repeat{
      OGorgsum[[i]] <- sum(OGorgs[[i]])
      i = i+1

      if(i == entities+1){
        break
      }
    }
    #assign('OGorgsum', OGorgsum, envir = .GlobalEnv)

    #Export of CSV's with values for discriptive metrics and general ecology discriptives#
    write.csv(inten_rang, file = paste0(prefix," ",year," OG range intensiveness.csv"))
    write.csv(inten_cell, file = paste0(prefix," ",year," OG cell intensiveness.csv"))
    write.csv(ext_org, file = paste0(prefix," ",year," OG exstensiveness.csv"))
    write.csv(OGorgsum, file = paste0(prefix," ",year," OG orgsum.csv"))

    #Creation of dataframe to turn to csv
    output_df <- data.frame(
      inten_rang = unlist(inten_rang),
      inten_cell = unlist(inten_cell),
      ext_org = unlist(ext_org),
      OGorgsum = unlist(OGorgsum)
    )

    #assign('output_df', output_df, envir = .GlobalEnv)
    # Write the dataframe to a CSV file
    write.csv(output_df, file = paste0(prefix," ",year," descriptive output.csv"))

    #sumorgs saveout#
    for(e in 1:entities){
      sumtitle <- paste0("OG sumorg ",e)
      write.csv(sumorgs[[TR/10]][e], file = paste0(prefix," ",year," ", sumtitle, ".csv"))
    }

    #orglist saveout#
    for(e in 1:entities){
      orgtitle <- paste0("OG orglist ",e)
      write.csv(orglist[[TR/10]][e], file = paste0(prefix," ",year," ", orgtitle, ".csv"))
    }

    #total ecology image saveout#
    write.csv(totaeco,file = paste0(prefix," ",year," OG totaeco.csv"))

    #Dimension max and min saveout#
    write.csv(orgs_x_min,file = paste0(prefix," ",year," OG orgs_x_min.csv"))
    write.csv(orgs_x_max,file = paste0(prefix," ",year," OG orgs_x_max.csv"))

    write.csv(orgs_y_min,file = paste0(prefix," ",year," OG orgs_y_min.csv"))
    write.csv(orgs_y_max,file = paste0(prefix," ",year," OG orgs_y_max.csv"))
    if (nDim == 3) {
      write.csv(orgs_z_min,file = paste0(prefix," ",year," OG orgs_z_min.csv"))
      write.csv(orgs_z_max,file = paste0(prefix," ",year," OG orgs_z_max.csv"))
    }
    #Export bin variable counts#
    if(bin == 1 ){
      write.csv(Bivar_count_hi, file = paste0(prefix," ",year, " OG orgs_binary_max.csv"))
      write.csv(Bivar_count_lo, file = paste0(prefix," ",year, " OG orgs_binary_min.csv"))
    }



    #setwd(saveDir)


    #setwd(OGpath)

    #write.csv(Bivar_count_hi, file = paste0(prefix," orgs_binary_max.csv"))

    #setwd("~/Desktop/Temp Folder")

    #year <- "1974"
    #Observation Means#
    #This section, utilizing sumorg, is not used for anything in the seed or simulation versions of the
    #capacity and exploitation code. I'm going to quick turn it into a list which will iterate over the entities variable#
    sumorg_OG <- list()

    i <- 0

    repeat{
      i = i+1
      sumorg_OG[[i]] <- c(sum(OGorgs[[i]]))
      if(i == entities){
        break
      }
    }
    #assign('sumorg_OG', sumorg_OG, envir = .GlobalEnv)

    #Summation for the total ecology (all organizations included, but location of resources is saved)#
    totalmatrix <- Reduce('+',OGorgs)

    total <- sum(totalmatrix)

    #assign('totalmatrix', totalmatrix, envir = .GlobalEnv)
    #assign('total', total, envir = .GlobalEnv)

    resource_OG <- list()

    i <- 0

    repeat{
      i = i+1
      resource_OG[[i]] <- Reduce("+",totalmatrix[OGorgs[[i]]!= 0])
      if(i == entities){
        break
      }
    }

    q <- 1

    repeat{
      if(sum(resource_OG[[q]]) == 0){
        resource_OG[[q]] <- 0
      }
      q = q+1
      if(q == entities)
        break
    }
    #assign('resource_OG', resource_OG, envir = .GlobalEnv)

    write.csv(resource_OG, file = paste0(prefix," ",year," OG resources.csv"))

    capacity <- list()
    wcapacity <- list()

    weights <- list()
    weights <- totalmatrix/ totindi
    weights[is.nan(weights)] <- 0
    weights[weights == Inf] <- 0

    #assign('weights', weights, envir = .GlobalEnv)


    i <- 0
    repeat{
      i <- i + 1
      capacity[[i]] <- totalmatrix
      capacity[[i]][OGorgs[[i]] == 0] <- 0
      wcapacity[[i]] <- capacity[[i]]*weights
      if(i == entities){
        break
      }
    }
    #assign('capacity', capacity, envir = .GlobalEnv)
    #assign('wcapacity', wcapacity, envir = .GlobalEnv)

    OG_wcarry_capacity <- list()
    OG_capacity <- list()

    i <- 0

    repeat{
      i <- i + 1
      OG_wcarry_capacity[[i]] <- sum(wcapacity[[i]])
      OG_capacity[[i]] <- sum(capacity[[i]])
      if(i == entities){
        break
      }
    }
    #assign('OG_wcarry_capacity', OG_wcarry_capacity, envir = .GlobalEnv)
    #assign('OG_capacity', OG_capacity, envir = .GlobalEnv)

    #Save Exploitation values out#
    write.csv(OG_wcarry_capacity, file = paste0(prefix," ", year," OG wcarry capacity.csv"))
    write.csv(OG_capacity, file = paste0(prefix," ", year," OG carry capacity.csv"))

    exploitation <- list()
    exploitation_value <- list()
    exploitation_value_members_weighted <- list()
    exploitation_value_members_unweighted <- list()
    exploitation_value_space <- list()
    exploitation_value_space_range <- list()

    i <- 0

    repeat{
      i <- i + 1
      exploitation[[i]] <- OGorgs[[i]] - wcapacity[[i]]
      exploitation_value_members[i] <- abs(Reduce("+",exploitation[[i]])/sum(totalmatrix*weights))
      exploitation_value_space[i] <- exploitation_value_members[[i]]/sum(wcapacity[[i]]>0)
      exploitation_value_space_range[i] <- abs(sum(exploitation[[i]])/sum(wcapacity[[i]]>0))

      if(i == entities){
        break
      }
    }
    #assign('exploitation', exploitation, envir = .GlobalEnv)
    #assign('exploitation_value_members', exploitation_value_members, envir = .GlobalEnv)
    #assign('exploitation_value_space', exploitation_value_space, envir = .GlobalEnv)
    #assign('exploitation_value_space_range', exploitation_value_space_range, envir = .GlobalEnv)

    #Original values for Exploitation#
    OG_Exploitation_members <- exploitation_value_members
    OG_Exploitation_space <- exploitation_value_space
    OG_Exploitation_space_range <- exploitation_value_space_range
    #assign('OG_Exploitation_members', OG_Exploitation_members, envir = .GlobalEnv)
    #assign('OG_Exploitation_space', OG_Exploitation_space, envir = .GlobalEnv)
    #assign('OG_Exploitation_space_range', OG_Exploitation_space_range, envir = .GlobalEnv)

    #Save Exploitation values out#
    write.csv(exploitation_value_members, file = paste0(prefix," ",year," OG Exploitation value members.csv"))
    write.csv(exploitation_value_space, file = paste0(prefix," ",year," OG Exploitation value space.csv"))
    write.csv(exploitation_value_space_range, file = paste0(prefix," ",year," OG Exploitation value space range.csv"))

    #New Directory to make things cleaner, saveDir is still used for the tables that will be needed for GOF test#
    simDir <- paste0("/AS - After Simulation")

    simpath <- paste0(saveDir, simDir)

    dir.create(file.path(simpath))

    simDir2 <- paste0("/AS - After Simulation/means")

    simpath2 <- paste0(saveDir, simDir2)

    dir.create(file.path(simpath2))

    #setwd(simpath2)

    #Setup matrix objects#
    exploitation_mean <- matrix()
    wcarry_mean <- matrix()

    #Mean for simulated results of exploitation metric#
    exploithist_matrix <- matrix(unlist(exploitMhist), nrow = 25, ncol = 16)

    #assign('exploithist_matrix', exploithist_matrix, envir = .GlobalEnv)

    i <- 0

    repeat{
      i = i+1
      exploitation_mean[i] <- mean(exploithist_matrix[,i])
      if(i == entities){
        break
      }
    }
    #assign('exploitation_mean', exploitation_mean, envir = .GlobalEnv)


    #Mean for simulated results of weighted carrying capacity metric#
    wcarry_matrix <- matrix(unlist(WCChist), nrow = 25, ncol = 16)

    #assign('wcarry_matrix', wcarry_matrix, envir = .GlobalEnv)


    i <- 0

    repeat{
      i = i+1
      wcarry_mean[i] <- mean(wcarry_matrix[,i])
      if(i == entities){
        break
      }
    }
    #assign('wcarry_mean', wcarry_mean, envir = .GlobalEnv)

    write.csv(wcarry_mean, file = paste0(prefix," ", year," wcarry capacity - simulated.csv"))
    write.csv(exploitation_mean, file = paste0(prefix," ",year," Exploitation value members - Simulated.csv"))



    #Extensiveness GOF t-test#
    exten_GOF_means <- matrix()
    exten_GOF_SD <- matrix()
    exten_GOF_min <- matrix()
    exten_GOF_max <- matrix()
    exten_adj_mean <- matrix()
    exten_adj_var <- matrix()
    exten_adj_SD <- matrix()
    exten_t_score <- matrix()
    exten_p_value <- matrix()

    #Each for loop utilizes the truncnorm pkg#
    for (ln in 1:entities){
      exten_GOF_means[ln] <- mean(extenlist[,ln])
      exten_GOF_SD[ln] <- sd(extenlist[,ln])
      exten_GOF_min[ln] <- min(extenlist[,ln])
      exten_GOF_max[ln] <- max(extenlist[,ln])
      exten_adj_mean[ln] <- truncnorm::etruncnorm(a=exten_GOF_min[ln], b=exten_GOF_max[ln], mean=exten_GOF_means[ln], sd=exten_GOF_SD[ln])
      exten_adj_var[ln] <- truncnorm::vtruncnorm(a=exten_GOF_min[ln], b=exten_GOF_max[ln], mean=exten_GOF_means[ln], sd=exten_GOF_SD[ln])
      if(is.nan(exten_adj_mean[ln]) == "TRUE"){
        exten_adj_mean[ln] <- exten_GOF_means[ln]
      }
      if(is.nan(exten_adj_var[ln]) == "TRUE"){
        exten_adj_var[ln] <- var(extenlist[,ln])
      }
      exten_adj_SD[ln] <- sqrt(exten_adj_var[ln])
      exten_t_score[ln] <- (exten_adj_mean[ln]-as.numeric(ext_org[ln]))/(exten_adj_SD[ln]/sqrt(glob))
      if(is.nan(exten_t_score[ln]) == "TRUE"){
        exten_t_score[ln] <- 0
      }
      exten_p_value[ln] <- 2*pt(exten_t_score[ln],glob-1,lower.tail = FALSE)
    }
    #assign('exten_GOF_means', exten_GOF_means, envir = .GlobalEnv)
    #assign('exten_GOF_SD', exten_GOF_SD, envir = .GlobalEnv)
    #assign('exten_GOF_min', exten_GOF_min, envir = .GlobalEnv)
    #assign('exten_GOF_max', exten_GOF_max, envir = .GlobalEnv)
    #assign('exten_adj_mean', exten_adj_mean, envir = .GlobalEnv)
    #assign('exten_adj_var', exten_adj_var, envir = .GlobalEnv)
    #assign('exten_adj_SD', exten_adj_SD, envir = .GlobalEnv)
    #assign('exten_t_score', exten_t_score, envir = .GlobalEnv)
    #assign('exten_p_value', exten_p_value, envir = .GlobalEnv)

    #Cell Focused Intensiveness GOF#
    intencell_GOF_means <- matrix()
    intencell_GOF_SD <- matrix()
    intencell_GOF_min <- matrix()
    intencell_GOF_max <- matrix()
    intencell_adj_mean <- matrix()
    intencell_adj_var <- matrix()
    intencell_adj_SD <- matrix()
    intencell_t_score <- matrix()
    intencell_p_value <- matrix()

    for (ln in 1:entities){
      intencell_GOF_means[ln] <- mean(intencelllist[,ln])
      intencell_GOF_SD[ln] <- sd(intencelllist[,ln])
      intencell_GOF_min[ln] <- min(intencelllist[,ln])
      intencell_GOF_max[ln] <- max(intencelllist[,ln])
      intencell_adj_mean[ln] <- truncnorm::etruncnorm(a=intencell_GOF_min[ln], b=intencell_GOF_max[ln], mean=intencell_GOF_means[ln], sd=intencell_GOF_SD[ln])
      intencell_adj_var[ln] <- truncnorm::vtruncnorm(a=intencell_GOF_min[ln], b=intencell_GOF_max[ln], mean=intencell_GOF_means[ln], sd=intencell_GOF_SD[ln])
      if(is.nan(intencell_adj_mean[ln]) == "TRUE"){
        intencell_adj_mean[ln] <- intencell_GOF_means[ln]
      }
      if(is.nan(intencell_adj_var[ln]) == "TRUE"){
        intencell_adj_var[ln] <- var(intencelllist[,ln])
      }
      intencell_adj_SD[ln] <- sqrt(intencell_adj_var[ln])
      intencell_t_score[ln] <- (intencell_adj_mean[ln]-as.numeric(inten_cell[ln]))/(intencell_adj_SD[ln]/sqrt(glob))
      if(is.nan(intencell_t_score[ln]) == "TRUE"){
        intencell_t_score[ln] <- 0
      }
      intencell_p_value[ln] <- 2*pt(intencell_t_score[ln],glob-1,lower.tail = FALSE)
    }
    #assign('intencell_GOF_means', intencell_GOF_means, envir = .GlobalEnv)
    #assign('intencell_GOF_SD', intencell_GOF_SD, envir = .GlobalEnv)
    #assign('intencell_GOF_min', intencell_GOF_min, envir = .GlobalEnv)
    #assign('intencell_GOF_max', intencell_GOF_max, envir = .GlobalEnv)
    #assign('intencell_adj_mean', intencell_adj_mean, envir = .GlobalEnv)
    #assign('intencell_adj_var', intencell_adj_var, envir = .GlobalEnv)
    #assign('intencell_adj_SD', intencell_adj_SD, envir = .GlobalEnv)
    #assign('intencell_t_score', intencell_t_score, envir = .GlobalEnv)
    #assign('intencell_p_value', intencell_p_value, envir = .GlobalEnv)

    #Range Focused Intensiveness GOF#
    intenrang_GOF_means <- matrix()
    intenrang_GOF_SD <- matrix()
    intenrang_GOF_min <- matrix()
    intenrang_GOF_max <- matrix()
    intenrang_adj_mean <- matrix()
    intenrang_adj_var <- matrix()
    intenrang_adj_SD <- matrix()
    intenrang_t_score <- matrix()
    intenrang_p_value <- matrix()

    for (ln in 1:entities){
      intenrang_GOF_means[ln] <- mean(intenranglist[,ln])
      intenrang_GOF_SD[ln] <- sd(intenranglist[,ln])
      intenrang_GOF_min[ln] <- min(intenranglist[,ln])
      intenrang_GOF_max[ln] <- max(intenranglist[,ln])
      intenrang_adj_mean[ln] <- truncnorm::etruncnorm(a=intenrang_GOF_min[ln], b=intenrang_GOF_max[ln], mean=intenrang_GOF_means[ln], sd=intenrang_GOF_SD[ln])
      intenrang_adj_var[ln] <- truncnorm::vtruncnorm(a=intenrang_GOF_min[ln], b=intenrang_GOF_max[ln], mean=intenrang_GOF_means[ln], sd=intenrang_GOF_SD[ln])
      if(is.nan(intenrang_adj_mean[ln]) == "TRUE"){
        intenrang_adj_mean[ln] <- intenrang_GOF_means[ln]
      }
      if(is.nan(intenrang_adj_var[ln]) == "TRUE"){
        intenrang_adj_var[ln] <- var(intenranglist[,ln])
      }
      intenrang_adj_SD[ln] <- sqrt(intenrang_adj_var[ln])
      intenrang_t_score[ln] <- (intenrang_adj_mean[ln]-as.numeric(inten_rang[ln]))/(intenrang_adj_SD[ln]/sqrt(glob))
      if(is.nan(intenrang_t_score[ln]) == "TRUE"){
        intenrang_t_score[ln] <- 0
      }
      intenrang_p_value[ln] <- 2*pt(intenrang_t_score[ln],glob-1,lower.tail = FALSE)
    }
    #assign('intenrang_GOF_means', intenrang_GOF_means, envir = .GlobalEnv)
    #assign('intenrang_GOF_SD', intenrang_GOF_SD, envir = .GlobalEnv)
    #assign('intenrang_GOF_min', intenrang_GOF_min, envir = .GlobalEnv)
    #assign('intenrang_GOF_max', intenrang_GOF_max, envir = .GlobalEnv)
    #assign('intenrang_adj_mean', intenrang_adj_mean, envir = .GlobalEnv)
    #assign('intenrang_adj_var', intenrang_adj_var, envir = .GlobalEnv)
    #assign('intenrang_adj_SD', intenrang_adj_SD, envir = .GlobalEnv)
    #assign('intenrang_t_score', intenrang_t_score, envir = .GlobalEnv)
    #assign('intenrang_p_value', intenrang_p_value, envir = .GlobalEnv)

    #setwd(saveDir)

    write.csv(exten_p_value, file = paste0(prefix," ",year," GOF_extensiveness.csv"))
    write.csv(intencell_p_value, file = paste0(prefix," ",year," GOF_Intensiveness_cell.csv"))
    write.csv(intenrang_p_value, file = paste0(prefix," ",year," GOF_Intensiveness_range.csv"))

    #Exploitation and Carrying Capacity Import - only for old runs, should be commented out#
    #The commented code below is being kept in for legacy purposes#

    #ln <- 1
    #Exploitation_M_gof <- matrix()
    #for(rn in 1:glob){
    #  exploitation_M_temp <- read.csv(file = paste0("~/Dropbox/SPPA/Model Output/Sunbelt Runs/SPPA ",year," tol 25 crange 10/AS - After Simulation/SPPA ",year," ",rn," Exploitation value members - Simulated.csv"))
    #Exploitation_M[ln] <- read.csv(file = paste0("~/Dropbox/SPPA/Model Output/Sunbelt Runs/SPPA 1982 tol 25 crange 10/AS - After Simulation/SPPA 1982 ",ln," Exploitation value members - Simulated.csv"))
    #  if(rn == 1){Exploitation_M_gof <- exploitation_M_temp}
    #  if(rn > 1){Exploitation_M_gof[rn,] <- exploitation_M_temp}
    #if(ln > 1){Exploitation_M <- merge(Exploitation_M, exploitation_temp)}
    #Exploitation_M[ln] <- as.matrix(exploitation_temp)
    #  }

    #unlist(Exploitation_M_gof)
    #Exploitation_M_gof <- subset(Exploitation_M_gof, select = -X)
    #names(Exploitation_M_gof) <- NULL

    #Wcarry_capacity_gof <- matrix()
    #for(rn in 1:glob){
    #  wcarry_capacity_temp <- read.csv(file = paste0("~/Dropbox/SPPA/Model Output/Sunbelt Runs/SPPA ",year," tol 25 crange 10/AS - After Simulation/SPPA ",year," ",rn," wcarry capacity - Simulated.csv"))
    #Exploitation_M[ln] <- read.csv(file = paste0("~/Dropbox/SPPA/Model Output/Sunbelt Runs/SPPA 1982 tol 25 crange 10/AS - After Simulation/SPPA 1982 ",ln," Exploitation value members - Simulated.csv"))
    #  if(rn == 1){Wcarry_capacity_gof <- wcarry_capacity_temp}
    #  if(rn > 1){Wcarry_capacity_gof[rn,] <- wcarry_capacity_temp}
    #if(ln > 1){Exploitation_M <- merge(Exploitation_M, exploitation_temp)}
    #Exploitation_M[ln] <- as.matrix(exploitation_temp)
    #}

    #unlist(Wcarry_capacity_gof)
    #Wcarry_capacity_gof <- subset(Wcarry_capacity_gof, select = -X)
    #names(Wcarry_capacity_gof) <- NULL

    #Change matrix name so I don't need to redo the code below (important for legacy issues, but that's about all)#
    Exploitation_M_gof <- exploithist_matrix

    #Exploitation GOF#

    Exploitation_M_GOF_means <- matrix()
    Exploitation_M_GOF_SD <- matrix()
    Exploitation_M_GOF_min <- matrix()
    Exploitation_M_GOF_max <- matrix()
    Exploitation_M_adj_mean <- matrix()
    Exploitation_M_adj_var <- matrix()
    Exploitation_M_adj_SD <- matrix()
    Exploitation_M_t_score <- matrix()
    Exploitation_M_p_value <- matrix()

    for (ln in 1:entities){
      Exploitation_M_GOF_means[ln] <- mean(Exploitation_M_gof[,ln])
      Exploitation_M_GOF_SD[ln] <- sd(Exploitation_M_gof[,ln])
      Exploitation_M_GOF_min[ln] <- min(Exploitation_M_gof[,ln])
      Exploitation_M_GOF_max[ln] <- max(Exploitation_M_gof[,ln])
      Exploitation_M_adj_mean[ln] <- truncnorm::etruncnorm(a=Exploitation_M_GOF_min[ln], b=Exploitation_M_GOF_max[ln], mean=Exploitation_M_GOF_means[ln], sd=Exploitation_M_GOF_SD[ln])
      Exploitation_M_adj_var[ln] <- truncnorm::vtruncnorm(a=Exploitation_M_GOF_min[ln], b=Exploitation_M_GOF_max[ln], mean=Exploitation_M_GOF_means[ln], sd=Exploitation_M_GOF_SD[ln])
      if(is.nan(Exploitation_M_adj_mean[ln]) == "TRUE"){
        Exploitation_M_adj_mean[ln] <- Exploitation_M_GOF_means[ln]
      }
      if(is.nan(Exploitation_M_adj_var[ln]) == "TRUE"){
        Exploitation_M_adj_var[ln] <- var(Exploitation_M_gof[,ln])
      }
      Exploitation_M_adj_SD[ln] <- sqrt(Exploitation_M_adj_var[ln])
      Exploitation_M_t_score[ln] <- (Exploitation_M_adj_mean[ln]-as.numeric(OG_Exploitation_members[ln]))/(Exploitation_M_adj_SD[ln]/sqrt(glob))
      if(is.nan(Exploitation_M_t_score[ln]) == "TRUE"){
        Exploitation_M_t_score[ln] <- 0
      }
      Exploitation_M_p_value[ln] <- 2*pt(Exploitation_M_t_score[ln],glob-1,lower.tail = FALSE)
    }
    #assign('Exploitation_M_GOF_means', Exploitation_M_GOF_means, envir = .GlobalEnv)
    #assign('Exploitation_M_GOF_SD', Exploitation_M_GOF_SD, envir = .GlobalEnv)
    #assign('Exploitation_M_GOF_min', Exploitation_M_GOF_min, envir = .GlobalEnv)
    #assign('Exploitation_M_GOF_max', Exploitation_M_GOF_max, envir = .GlobalEnv)
    #assign('Exploitation_M_adj_mean', Exploitation_M_adj_mean, envir = .GlobalEnv)
    #assign('Exploitation_M_adj_var', Exploitation_M_adj_var, envir = .GlobalEnv)
    #assign('Exploitation_M_adj_SD', Exploitation_M_adj_SD, envir = .GlobalEnv)
    #assign('Exploitation_M_t_score', Exploitation_M_t_score, envir = .GlobalEnv)
    #assign('Exploitation_M_p_value', Exploitation_M_p_value, envir = .GlobalEnv)

    #Change matrix name so I don't need to redo the code below (important for legacy issues, but that's about all)#
    Wcarry_capacity_gof <- wcarry_matrix

    #Weighted Carrying Capacity GOF#

    Wcarry_capacity_GOF_means <- matrix()
    Wcarry_capacity_GOF_SD <- matrix()
    Wcarry_capacity_GOF_min <- matrix()
    Wcarry_capacity_GOF_max <- matrix()
    Wcarry_capacity_adj_mean <- matrix()
    Wcarry_capacity_adj_var <- matrix()
    Wcarry_capacity_adj_SD <- matrix()
    Wcarry_capacity_t_score <- matrix()
    Wcarry_capacity_p_value <- matrix()

    for (ln in 1:entities){
      Wcarry_capacity_GOF_means[ln] <- mean(Wcarry_capacity_gof[,ln])
      Wcarry_capacity_GOF_SD[ln] <- sd(Wcarry_capacity_gof[,ln])
      Wcarry_capacity_GOF_min[ln] <- min(Wcarry_capacity_gof[,ln])
      Wcarry_capacity_GOF_max[ln] <- max(Wcarry_capacity_gof[,ln])
      Wcarry_capacity_adj_mean[ln] <- truncnorm::etruncnorm(a=Wcarry_capacity_GOF_min[ln], b=Wcarry_capacity_GOF_max[ln], mean=Wcarry_capacity_GOF_means[ln], sd=Wcarry_capacity_GOF_SD[ln])
      Wcarry_capacity_adj_var[ln] <- truncnorm::vtruncnorm(a=Wcarry_capacity_GOF_min[ln], b=Wcarry_capacity_GOF_max[ln], mean=Wcarry_capacity_GOF_means[ln], sd=Wcarry_capacity_GOF_SD[ln])
      if(is.nan(Wcarry_capacity_adj_mean[ln]) == "TRUE"){
        Wcarry_capacity_adj_mean[ln] <- Wcarry_capacity_GOF_means[ln]
      }
      if(is.nan(Wcarry_capacity_adj_var[ln]) == "TRUE"){
        Wcarry_capacity_adj_var[ln] <- var(Wcarry_capacity_gof[,ln])
      }
      Wcarry_capacity_adj_SD[ln] <- sqrt(Wcarry_capacity_adj_var[ln])
      Wcarry_capacity_t_score[ln] <- (Wcarry_capacity_adj_mean[ln]-as.numeric(OG_wcarry_capacity[ln]))/(Wcarry_capacity_adj_SD[ln]/sqrt(glob))
      if(is.nan(Wcarry_capacity_t_score[ln]) == "TRUE"){
        Wcarry_capacity_t_score[ln] <- 0
      }
      Wcarry_capacity_p_value[ln] <- 2*pt(Wcarry_capacity_t_score[ln],glob-1,lower.tail = FALSE)
    }
    #assign('Wcarry_capacity_GOF_means', Wcarry_capacity_GOF_means, envir = .GlobalEnv)
    #assign('Wcarry_capacity_M_GOF_SD', Exploitation_GOF_SD, envir = .GlobalEnv)
    #assign('Wcarry_capacity_GOF_min', Wcarry_capacity_GOF_min, envir = .GlobalEnv)
    #assign('Wcarry_capacity_GOF_max', Wcarry_capacity_GOF_max, envir = .GlobalEnv)
    #assign('Wcarry_capacity_adj_mean', Wcarry_capacity_adj_mean, envir = .GlobalEnv)
    #assign('Wcarry_capacity_adj_var', Wcarry_capacity_adj_var, envir = .GlobalEnv)
    #assign('Wcarry_capacity_adj_SD', Wcarry_capacity_adj_SD, envir = .GlobalEnv)
    #assign('Wcarry_capacity_t_score', Wcarry_capacity_t_score, envir = .GlobalEnv)
    #assign('Wcarry_capacity_p_value', Wcarry_capacity_p_value, envir = .GlobalEnv)

    #setwd(saveDir)
    write.csv(Exploitation_M_p_value, file = paste0(prefix," ",year," GOF_Exploitation_M_p_value.csv"))
    write.csv(Wcarry_capacity_p_value, file = paste0(prefix," ",year," GOF_Wcarry_capacity_p_value.csv"))
  }

  ####Saveout of OGorgs array and orgs array as r objects (saved as r objects for easy future import)####
  #setwd(saveDir)
  saveRDS(OGorgs, file = "OGorgs_1974.rds") #Original array from observations#
  saveRDS(orgs, file = "orgs_1974.rds") #data after final iterations. This is final iteration because I'm getting it after the final run, but could more this code and rework it to be after a run of the model#

  ####total run time script####
  endtime <- Sys.time()
  totaltime <- difftime(endtime, starttime, units = c("mins"))
  print(starttime)
  print(endtime)
  print(totaltime)

  ###Move objects to environment###

  # if (nDim == 2) {
  #   objectdump <- list(c(aveorglist_temp, aveorgtemp, aveorgtemp2, changneg, changpos, cjk,
  #                        ctempval, orgs_lengx, orgs_lengy, orgs_x_max, orgs_x_min, orgs_x_maxhist,
  #                        orgs_x_maxlist, orgs_x_minlist, orgs_x_minhist, orgs_y_max, orgs_y_min,
  #                        orgs_y_maxhist, orgs_y_maxlist, orgs_y_minlist, orgs_y_minhist, tempcoordloc,
  #                        tempsum, temptest, Changeneg, Changepos, ctempval))
  #
  #   outputlist <- list(c(aveorglist_matrix, Bivar_count_hi, Bivar_count_lo, capacity, Change, converg,
  #                        convergcount, conversave, Coord, demographics, ecoindices, exploitation,
  #                        exploitation_value_members, exploitation_value_space, exploitation_value_space_range,
  #                        exploitMhist, ext_org, extenhist, extenlist, focelist, totpops, wcapacity, sumorgs,
  #                        sumorg_sim, sim_wcarry_capacity, output_df, resource_sim, OGorgs, OGorgsum,
  #                        orglist, orgs))
  # }
  # if (nDim == 3) {
  #   objectdump <- list(c(aveorglist_temp, aveorgtemp, aveorgtemp2, caltop, changneg, changpos, cjk,
  #                        ctempval, orgs_lengx, orgs_lengy, orgs_lengz, orgs_z_max, orgs_z_min, orgs_z_maxhist,
  #                        orgs_z_maxlist, orgs_z_minlist, orgs_z_minhist, orgs_x_max, orgs_x_min, orgs_x_maxhist,
  #                        orgs_x_maxlist, orgs_x_minlist, orgs_x_minhist, orgs_y_max, orgs_y_min,
  #                        orgs_y_maxhist, orgs_y_maxlist, orgs_y_minlist, orgs_y_minhist, tempcoordloc,
  #                        tempsum, temptest, Changeneg, Changepos, ctempval, Importdat))
  #
  #   outputlist <- list(c(aveorglist_matrix, Bivar_count_hi, Bivar_count_lo, capacity, Change, converg,
  #                        convergcount, conversave, Coord, demographics, ecoindices, exploitation,
  #                        exploitation_value_members, exploitation_value_space, exploitation_value_space_range,
  #                        exploitMhist, ext_org, extenhist, extenlist, focelist, totpops, wcapacity, sumorgs,
  #                        sumorg_sim, sim_wcarry_capacity, output_df, resource_sim, OGorgs, OGorgsum,
  #                        orglist, orgs))
  # }

##Objects/values I am unsure of## ~Gage

  if (bin==99) {
    stringr::str_replace()
    cdfquantreg::Ambdata
  }

  #bivar_count_hi, bivar_count_lo, converg, convergcount, extenhist, Importdat,

  if (notif == TRUE){
    beepr::beep()
  }



  ####Cleaning Scripts####
  #Running the code below clears all put the initial data from R#
  #Remove all but intital imports. Useful if you want to rerun the models, but do not want to reimport the data#
  #NOTE: THIS WILL REMOVE EVERYTHING AFTER THE FINAL MODEL RUN (mlc value)!
  #ls()
  #rm(list= ls()[!(ls() %in% c('OGorgs'))])
  testlist <- list(
    limit = limit,
    testlimit = testlimit
  )
  return(testlist)
}


