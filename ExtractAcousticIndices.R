#----------------------------------- About -------------------------------------
# Script for extracting the following 7 Acoustic Indices from field and VSE-
# based audio recordings of 6 sites: Acoustic Complexity Index (ACI), Acoustic
# Diversity Index (ADI), Acoustic Eveness (AEve), Bioacoustic Index (Bio),
# Acoustic Entropy (H), Median of the Amplitude Envelope (M), Normalised 
# Difference Soundscape Index (NDSI).


# V3.0, 04.12.2023

#-------------------------1. Import Packages and Audio -------------------------
# setwd("your_directory_path_here") Be sure to set the working directory to where the relevant data have been stored.

library(soundecology)
library(seewave)
library(tuneR)

# Names of files for each site:
site1_files <- c("SilwoodFieldRecToUse-Mic1-LP1.wav","SilwoodReRecVertToUse-Mic1-LP1.wav","SilwoodReRec45ToUse-Mic1-LP1.wav","SilwoodReRecHToUse-Mic1-LP1.wav")
site2_files <- c("SilwoodFieldRecToUse-Mic1-LP2.wav","SilwoodReRecVertToUse-Mic1-LP2.wav","SilwoodReRec45ToUse-Mic1-LP2.wav","SilwoodReRecHToUse-Mic1-LP2.wav")
site3_files <- c("SilwoodFieldRecToUse-Mic1-LP3.wav","SilwoodReRecVertToUse-Mic1-LP3.wav","SilwoodReRec45ToUse-Mic1-LP3.wav","SilwoodReRecHToUse-Mic1-LP3.wav")
site4_files <- c("SilwoodFieldRecToUse-Mic1-LP4.wav","SilwoodReRecVertToUse-Mic1-LP4.wav","SilwoodReRec45ToUse-Mic1-LP4.wav","SilwoodReRecHToUse-Mic1-LP4.wav")
site5_files <- c("SilwoodFieldRecToUse-Mic1-LP5.wav","SilwoodReRecVertToUse-Mic1-LP5.wav","SilwoodReRec45ToUse-Mic1-LP5.wav","SilwoodReRecHToUse-Mic1-LP5.wav")
site6_files <- c("SilwoodFieldRecToUse-Mic1-LP6.wav","SilwoodReRecVertToUse-Mic1-LP6.wav","SilwoodReRec45ToUse-Mic1-LP6.wav","SilwoodReRecHToUse-Mic1-LP6.wav")


allWaves <- vector(mode='list', length=6) #Create list in which to store lists of Wave files of (re-)recordings from each of the 6 sites.
mins <- vector(mode='character', length=6) #Vector to store shortest length in seconds of the (re-)recordings from each site.

fs <- 16000 #Sampling frequency

# Import the wave audio files:
allWaves[[1]] = lapply(site1_files, readWave)
allWaves[[2]] = lapply(site2_files, readWave)
allWaves[[3]] = lapply(site3_files, readWave)
allWaves[[4]] = lapply(site4_files, readWave)
allWaves[[5]] = lapply(site5_files, readWave)
allWaves[[6]] = lapply(site6_files, readWave)

# Get the length in s of the shortest of the re-recordings from each site:
mins[1] <- min(lengths(allWaves[[1]]))/fs
mins[2] <- min(lengths(allWaves[[2]]))/fs
mins[3] <- min(lengths(allWaves[[3]]))/fs
mins[4] <- min(lengths(allWaves[[4]]))/fs
mins[5] <- min(lengths(allWaves[[5]]))/fs
mins[6] <- min(lengths(allWaves[[6]]))/fs


#-----------------------2. Cut audio into 30s segments--------------------------

#First create lists for each site which will contain 4 lists each, one for 
# indices calculated for each of the 4 (re-)recordings of each site (Field,
# VSE Vertical, VSE 45°, VSE Horizontal):
site1_30 <- vector(mode='list', length=4)
site2_30 <- vector(mode='list', length=4)
site3_30 <- vector(mode='list', length=4)
site4_30 <- vector(mode='list', length=4)
site5_30 <- vector(mode='list', length=4)
site6_30 <- vector(mode='list', length=4)

# We use 20 windows as the windows are 30s long and each recording is
# approximately 10 minutes long:
for (i in 1:4){
  site1_30[[i]] <- vector(mode='list', length=20) 
  site2_30[[i]] <- vector(mode='list', length=20)
  site3_30[[i]] <- vector(mode='list', length=20)
  site4_30[[i]] <- vector(mode='list', length=20)
  site5_30[[i]] <- vector(mode='list', length=20)
  site6_30[[i]] <- vector(mode='list', length=20)
}

# Now cut the wave audio files into 30 s windows:
for (i in 1:4){
  for (j in 1:20){
    if (j < 20){ # Get the first 19 windows:
      site1_30[[i]][[j]] = cutw(allWaves[[1]][[i]], fs, channel=1, from = 30*(j-1), to = 30*j, choose = FALSE,plot = FALSE, marks = TRUE, output="Wave")
      site2_30[[i]][[j]] = cutw(allWaves[[2]][[i]], fs, channel=1, from = 30*(j-1), to = 30*j, choose = FALSE,plot = FALSE, marks = TRUE, output="Wave")
      site3_30[[i]][[j]] = cutw(allWaves[[3]][[i]], fs, channel=1, from = 30*(j-1), to = 30*j, choose = FALSE,plot = FALSE, marks = TRUE, output="Wave")
      site4_30[[i]][[j]] = cutw(allWaves[[4]][[i]], fs, channel=1, from = 30*(j-1), to = 30*j, choose = FALSE,plot = FALSE, marks = TRUE, output="Wave")
      site5_30[[i]][[j]] = cutw(allWaves[[5]][[i]], fs, channel=1, from = 30*(j-1), to = 30*j, choose = FALSE,plot = FALSE, marks = TRUE, output="Wave")
      site6_30[[i]][[j]] = cutw(allWaves[[6]][[i]], fs, channel=1, from = 30*(j-1), to = 30*j, choose = FALSE,plot = FALSE, marks = TRUE, output="Wave")
    } else { # For the final window, use however many seconds are remaining in the shortest recording of each site (usually less than 30 s):
      site1_30[[i]][[j]] = cutw(allWaves[[1]][[i]], fs, channel=1, from = 30*(j-1), to = as.numeric(mins[1]), choose = FALSE,plot = FALSE, marks = TRUE, output="Wave")
      site2_30[[i]][[j]] = cutw(allWaves[[2]][[i]], fs, channel=1, from = 30*(j-1), to = as.numeric(mins[2]), choose = FALSE,plot = FALSE, marks = TRUE, output="Wave")
      site3_30[[i]][[j]] = cutw(allWaves[[3]][[i]], fs, channel=1, from = 30*(j-1), to = as.numeric(mins[3]), choose = FALSE,plot = FALSE, marks = TRUE, output="Wave")
      site4_30[[i]][[j]] = cutw(allWaves[[4]][[i]], fs, channel=1, from = 30*(j-1), to = as.numeric(mins[4]), choose = FALSE,plot = FALSE, marks = TRUE, output="Wave")
      site5_30[[i]][[j]] = cutw(allWaves[[5]][[i]], fs, channel=1, from = 30*(j-1), to = as.numeric(mins[5]), choose = FALSE,plot = FALSE, marks = TRUE, output="Wave")
      site6_30[[i]][[j]] = cutw(allWaves[[6]][[i]], fs, channel=1, from = 30*(j-1), to = as.numeric(mins[6]), choose = FALSE,plot = FALSE, marks = TRUE, output="Wave")
    }
  }
}


#---------------3. Calculate Acoustic Indices and Save to CSVs------------------


#Create vectors to store the calculated indices 
aci_30 <- vector(mode='numeric', length=20)
adi_30 <- vector(mode='numeric', length=20)
aeev_30 <- vector(mode='numeric', length=20)
bio_30 <- vector(mode='numeric', length=20)
ndsi_30 <- vector(mode='numeric', length=20)
H_30 <- vector(mode='numeric', length=20)
M_30 <- vector(mode='numeric', length=20)


site <- 1
max.freq <- fs/2 #Max frequency is the Nyquist frequency (half of sampling frequency).

rerecNames <- c("FieldVert","LabVert","Lab45deg","LabH")

#Now loop through all sites, (re-)recordings, and window sizes to get the acoustic indices for all, and save the results to csv 
for (k in 1:6){
  
  if (k == 1){
    current_site30 <- site1_30
  }   else if (k == 2){
    current_site30 <- site2_30
  }  else if (k == 3){
    current_site30 <- site3_30
  }  else if (k == 4){
    current_site30 <- site4_30
  }  else if (k == 5){
    current_site30 <- site5_30
  }  else {
    current_site30 <- site6_30
  }
  
  # Extract the 7 acoustic indices:
  for (i in 1:4){
    for (j in 1:20){

      aci.list = acoustic_complexity(current_site30[[i]][[j]], min_freq = NA, max_freq = max.freq)
      aci_30[j] = as.numeric(aci.list$AciTotAll_left_bymin)
        
      adi.list = acoustic_diversity(current_site30[[i]][[j]], max_freq = max.freq)
      adi_30[j] = as.numeric(adi.list$adi_left)
        
      aeev.list = acoustic_evenness(current_site30[[i]][[j]], max_freq = max.freq)
      aeev_30[j] = as.numeric(aeev.list$aei_left)
        
      bio.list = bioacoustic_index(current_site30[[i]][[j]], max_freq = max.freq)
      bio_30[j] = as.numeric(bio.list$left_area)
        
      ndsi.list = ndsi(current_site30[[i]][[j]], bio_max = max.freq)
      ndsi_30[j] = as.numeric(ndsi.list$ndsi_left)
        
      H_30[j] = as.numeric(H(current_site30[[i]][[j]], fs))
      M_30[j] = as.numeric(M(current_site30[[i]][[j]], fs))


    }
    indices30 <- cbind(aci_30,adi_30,aeev_30,bio_30,ndsi_30,H_30,M_30) # Combine the results
    write.csv(indices30,paste("Site",k,rerecNames[i],"30sAI.csv",sep = ""),row.names = FALSE) # Write results to CSV files (1 for each recording)
  }
}

