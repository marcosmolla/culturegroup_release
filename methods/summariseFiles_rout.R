suppressMessages(library(dplyr))
library(reshape2)

path <- getwd()

if(file.exists(paste(path,"output/",sep="/"))){
  ## Define functions
  loadData <- function(PATH, PNPR=FALSE, JULIA=FALSE){
    summaryFilePath <- paste(path, paste(tail(strsplit(path, split="/")[[1]],1),"summary",sep="_"), sep="/") # Set path to save results file
    pnpr2_data <- do.call(rbind, lapply(list.files(PATH), function(f){ # for all simulations do
      load(paste(PATH, f, sep="")) # Load file
      return(setup) # Return the setup file (containts summarised results for all files)
    }))
    save(pnpr2_data, file=summaryFilePath) # Save new file or append to previous
  }

  loadRepertoireData <- function(PATH){
    summaryFilePath3 <- paste(path, paste(tail(strsplit(path, split="/")[[1]],1),"traitHistogram",sep="_"), sep="/") # Set path to save results file
    traitHistogram_c <- do.call(rbind, lapply(list.files(PATH), function(f){ # for all simulations do
      load(paste(PATH, f, sep="")) # Load file
      return(
        data.frame(
        # State values
        # neibhood=setup[1,"neibhood"],
        payoffMethod=setup[1,"payoffMethod"],
        innovationSuccess=setup[1,"innovationSuccess"],
        indSuccessRate=setup[1,"indSuccessRate"],
        # ID values
        seed=sample(x=LETTERS, size=10, replace=T) %>% paste(collapse=""),
        # Results
        traits=sort(colSums(repertoir!=0), decreasing=TRUE),
        position=1:ncol(repertoir))
      ) # Return the setup file (containts summarised results for all files)
    }))

    traitHistogram_m <- do.call(rbind, lapply(list.files(PATH), function(f){ # for all simulations do
      load(paste(PATH, f, sep="")) # Load file
      return(
        data.frame(
        # State values
        # neibhood=setup[1,"neibhood"],
        payoffMethod=setup[1,"payoffMethod"],
        innovationSuccess=setup[1,"innovationSuccess"],
        indSuccessRate=setup[1,"indSuccessRate"],
        # ID values
        seed=sample(x=LETTERS, size=10, replace=T) %>% paste(collapse=""),
        # Results
        traits=sort(colMeans(repertoir), decreasing=TRUE),
        position=1:ncol(repertoir))
      ) # Return the setup file (containts summarised results for all files)
    }))

    save(traitHistogram_c, traitHistogram_m, file=summaryFilePath3) # Save new file or append to previous
  }

  loadTimeData <- function(PATH){
    summaryFilePath2 <- paste(path, paste(tail(strsplit(path, split="/")[[1]],1),"timeData",sep="_"), sep="/") # Set path to save results file
    timeData <- do.call(rbind, lapply(list.files(PATH), function(f){ # for all simulations do
      load(paste(PATH, f, sep="")) # Load file
      return(
        data.frame(
        # State values
        # neibhood=setup[1,"neibhood"],
        payoffMethod=setup[1,"payoffMethod"],
        innovationSuccess=setup[1,"innovationSuccess"],
        densityDependence=setup[1,"densityDependence"],
        indSuccessRate=setup[1,"indSuccessRate"],
        totalGenerations=setup[1,"generation"],
        nSeedTraitsTotal=setup[1,"nSeedTraitsTotal"],
        startSeedTraits=setup[1,"startSeedTraits"],
        # ID values
        generations=1:setup[1,"generation"],
        seed=lapply(1:setup[1,"repe"], function(y) sample(x=LETTERS, size=10, replace=T)) %>% lapply(paste,collapse="") %>% unlist %>% rep(each=setup[1,"generation"]),
        # Results
        recPN=melt(recPNl)[,"value"],
        recPR=melt(recPRl)[,"value"],
        recNTraits=melt(recNTraitsl)[,"value"],
        recClust=melt(recClustl)[,"value"],
        recDeg=melt(recDegl)[,"value"],
        recPath=melt(recPathl)[,"value"],
        recTraitL=melt(recTraitLl)[,"value"],
        recTraitN=melt(recTraitNl)[,"value"],
        recBetaDiv=melt(recBetaDivl)[,"value"],
        recTraitDiversity=melt(recTraitDiversityl)[,"value"]#,
        # recLearnedIno=melt(recLearnedInol)[,"value"],
        # recLearnedSoc=melt(recLearnedSocl)[,"value"]
        )
      ) # Return the setup file (containts summarised results for all files)
    }))
    save(timeData, file=summaryFilePath2) # Save new file or append to previous
  }


  ## Run summarising function
  loadData(PATH=paste(path,"output/",sep="/"), PNPR=TRUE, JULIA=TRUE)
  loadRepertoireData(PATH=paste(path,"output/",sep="/"))
  loadTimeData(PATH=paste(path,"output/",sep="/"))
  # Zip output folder
  #zip(zipfile=paste(path,"output",sep="/"),files=paste(path,"output/",sep="/"), flags="-q -r")
  zip(zipfile="output",files="output/", flags="-q -r")
  # Remove unziped folder
  unlink(paste(path,"output/",sep="/"), recursive=TRUE)
  # Remove Slurm files
  system(paste("find ", path,"/ -name \"slurm*\" -delete", sep=""))
  print("Done!")
} else {
  print("There is no output directory to summarise (probably already done that!) :D ")
}
