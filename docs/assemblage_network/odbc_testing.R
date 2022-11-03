library(DBI)
library(odbc)
library(dplyr)
library(dbplyr)

dbFileName <- "1BOW.00.101 Prattsville (10-2014).accdb"

connect_to_access_dbi <- function(db_file_path)  {
  require(DBI)
  # make sure that the file exists before attempting to connect
  if (!file.exists(db_file_path)) {
    stop("DB file does not exist at ", db_file_path)
  }
  # Assemble connection strings
  dbq_string <- paste0("DBQ=", db_file_path)
  driver_string <-
    "driver={SQL Server};"
  dsn_string <- "dsn = {Microsoft Access Driver (*.mdb, *.accdb)};"
  db_connect_string <- paste0(driver_string, dsn_string, dbq_string)
  myconn <- dbConnect(odbc::odbc(),
                      .connection_string = db_connect_string)
  return(myconn)
}

test <- connect_to_access_dbi(dbFileName)


