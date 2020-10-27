# -----------------------------------
# Functions
# -----------------------------------

# Function to test if data is already present; if false, write data to csv
save_dat <- function(x) {
  if(file.exists(paste0("data_in/", x))) {  # test if data is already present in folder
          warning("File already exists") # if file exists, give warning message
         # if file does not exist, download .csv file from SalmonLRP_RetroEval github repository
      } else { dat <- read.csv("https://raw.githubusercontent.com/Pacific-salmon-assess/SalmonLRP_RetroEval/withChum/SCChumStudy/DataOut/SRdatWild.csv", stringsAsFactors = FALSE)
         dat1 <- dat[, names(dat) %in% c("Year", "CU", "Escape", "Recruit")] # select columns to keep
         names(dat1) <- c("year", "CU", "spawners", "recruits") # rename columns
         write.csv(dat1, paste0("data_in/", x), row.names = FALSE)
      } # end of else statement
}

# 
