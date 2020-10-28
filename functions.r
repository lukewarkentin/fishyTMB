# -----------------------------------
# Functions
# -----------------------------------

# Function to test if data is already present; if false, write data to csv
get_SR_dat <- function(x) {
  if(file.exists(paste0("data_in/", x))) {  # test if data is already present in folder
          warning("File already exists") # if file exists, give warning message
         # if file does not exist, download .csv file from SalmonLRP_RetroEval github repository
      } else { dat <- read.csv("https://raw.githubusercontent.com/Pacific-salmon-assess/SalmonLRP_RetroEval/withChum/SCChumStudy/DataOut/SRdatWild.csv", stringsAsFactors = FALSE)
         dat1 <- dat[, names(dat) %in% c("Year", "CU", "Escape", "Recruit")] # select columns to keep
         names(dat1) <- c("year", "CU", "spawners", "recruits") # rename columns
         write.csv(dat1, paste0("data_in/", x), row.names = FALSE)
      } # end of else statement
}

# Function to test whether TMB files are already compiled, if not, don't recompile. Then load model
loadTMB <- function(x) {
  modpath <- paste0("TMB/", x) # model path string
  if(file.exists(paste0(modpath, ".dll"))) { # test if TMB .cpp file is already present
        warning("Already compiled: .dll file already exists") # if file exists, give warning
      } else { compile(paste0(modpath, ".cpp")) # compile .cpp file of model
  } # end of else statement 
  dyn.load(dynlib(modpath)) # load .dll file of model
}



