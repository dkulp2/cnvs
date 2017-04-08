sql2R <- function(data.type) {
  # return an R type given the SQL type
  if (grepl('^text|^varchar|^char', data.type, ignore.case = TRUE)) {
    'character'
  } else if (grepl('^int', data.type, ignore.case = TRUE)) {
    'integer'
  } else if (grepl('^double|^float', data.type, ignore.case = TRUE)) {
    'numeric'
  } else if (grepl('^bool', data.type, ignore.case = TRUE)) {
    'logical'
  } else {
    warning("Unknown data type '",data.type,"'")
    'character'
  }
}

#' read a delimited file with metadata into a tibble
#' 
#' This function reads a tab-delimited file that contains
#' special metadata in the first line and returns a tibble.
#' @param filename the name of a file to read or an open connection. Can be a compressed file ending in .gz.
#' @export
#' @examples 
#' my.tbl <- readTbl("file.txt")
readTbl <- function(filename) {
  require(tibble)
  
  # file.conn <- gzfile(filename,"r")  # will read uncompressed, too. But I'm consistently getting "seek on a gzfile connection returned an internal error"
  if (inherits(filename,"connection")) {
    file.conn <- filename
    skip <- 1  # bug? if connection passed in, then read.table seems to reset. otherwise, doesn't.
  } else {
    if (grepl('.gz$', filename)) {
      # warning messages are raised regarding seek problems, but still seems to work
      file.conn <- gzfile(filename,"r")
    } else {
      file.conn <- file(filename, "r")
    }
    skip <- 0
  }
  
  header <- readLines(file.conn, 1)   # read the first line of file and parse the metadata

  header <- sub('^#\\s*','', header, perl=TRUE)
  parts <- unlist(strsplit(header, "\\s*:\\s*", perl=TRUE))
  
  tblName <- parts[1]
  columns <- parts[2]
  
  # extract any key data to store later as an attribute
  keys.pos <- regexpr('\\s*;\\s*(.*)', columns, perl=TRUE)
  keys <- substr(columns, attr(keys.pos, "capture.start"), nchar(columns))
  
  # ignore key data and split by column
  columns <- sub('\\s*;.*','', columns, perl=TRUE)
  columns <- unlist(strsplit(columns, "\\s*,\\s*", perl=TRUE))
  
  # iterate over each column, extracting the data type, if present
  colDefs <-
    lapply(columns, function(col) {
      matches <- regexpr("(\\S+)\\s*\\{\\s*(.+)\\}", col, perl=TRUE)
      match.pos <- attr(matches, "capture.start")
      match.len <- attr(matches, "capture.length")
      
      colName <- substr(col, match.pos[1], match.pos[1] + match.len[1] - 1)
      
      if (nchar(colName)==0) {
        colName <- col; colClass <- 'character'
      } else {
        colClass <- sql2R(substr(col, match.pos[2], match.pos[2] + match.len[2] - 1))
      }
      return(c(colName, colClass))
    })

  t <- as_tibble(read.delim(file.conn, quote="", col.names=sapply(colDefs, function(x) x[1]), 
                            colClasses=sapply(colDefs, function(x) x[2]), fill=FALSE, header=FALSE, skip=skip))
  close(file.conn)
  
  if (nchar(keys) > 0) {
    attr(t, 'keys') <- keys
  }
  attr(t, 'tablename') <- tblName
  
  t
}


#' read a delimited file with metadata and create an object in the current (or specified) environment
#' 
#' This function reads a tab-delimited file that contains
#' special metadata in the first line. It creates a tibble with the name of the
#' table as named in the metadata. There is return value.
#' @param filename the name of a file to read or an open connection. Can be a compressed file ending in .gz.
#' @param environ the environment to load the object into. Default is the current environment.
#' @param append add data to an existing tibble, if already present. Default is FALSE.
#' @export
#' @examples 
#' loadTbl("file.txt")
#' loadTbl("file.txt", append=TRUE)
loadTbl <- function(filename, environ=parent.frame(), append=FALSE) {
  tmp.tbl <- readTbl(filename)
  tbl.name <- attr(tmp.tbl,'tablename')
  if (append && exists(tbl.name, envir=environ)) {
    assign(tbl.name, rbind(get(tbl.name, envir = environ),tmp.tbl), envir=environ)
  } else {
    assign(attr(tmp.tbl, 'tablename'), tmp.tbl, envir=environ)
  }
}  
