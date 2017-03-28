# Routines for creating a "flat" set of environment variables that can be quickly read and set within R
# without calling the shell. This could be done without R, but this code shows how to read shell configuration
# files and interpret them directly in R.

# read export name=value lines from fname
# evaluating them with the shell so that other symbols and shell commands are properly evaluated.
do.conf <- function(fname) {
  conf <- grep('^export',readLines(fname), ignore.case=TRUE, value=TRUE)
  eq <- regexpr('=',conf)
  sapply(1:length(conf), function(i) {
    name <- substr(conf[i], 8, eq[i]-1)
    value.str <- substr(conf[i], eq[i]+1, nchar(conf[i]))
    value <- system(sprintf("%s -lc 'echo %s'", Sys.getenv("SHELL"), value.str), intern=T)
    names(value) <- name
    do.call(Sys.setenv, as.list(value))
    return(value)
  })
}

# write the values to fname. These are all expanded values that can be read and set without using the shell
write.conf <- function(fname, vals) {
  lines <- sprintf("export %s=%s", names(vals), vals)
  writeLines(lines, fname)
}

# just read and set, don't evaluate.
read.conf <- function(fname) {
  conf <- grep('^export',readLines(fname), ignore.case=TRUE, value=TRUE)
  eq <- regexpr('=',conf)
  env.vars <-
    sapply(1:length(conf), function(i) {
      name <- substr(conf[i], 8, eq[i]-1)
      value <- substr(conf[i], eq[i]+1, nchar(conf[i]))
      names(value) <- name
      return(value)
  })
  do.call(Sys.setenv, as.list(env.vars))
  
  cat(sprintf("%s=%s\n", names(env.vars), Sys.getenv(names(env.vars))))  
}


# returns TRUE if the named environment variable has a true-ish value
check.conf <- function(name) {
  grepl('^TRUE\\b|^T\\b|^1\\b|^YES\\b|^Y\\b', Sys.getenv(name), ignore.case = TRUE, perl=TRUE)
}
