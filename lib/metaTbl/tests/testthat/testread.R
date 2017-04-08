# metaTbl tests. This is just a placeholder for possibly a suite of tests.
# test with devtools::test() when working dir in library.

library(metaTbl)
context("Reading")

tbl.fn <- 'tbl.txt'

test_that("readTbl basic facts", {
  loadTbl(tbl.fn)
  expect_equal(nrow(x),2)
  expect_equal(length(x),2)
  expect_equal(x$foo.baz,c('AB','CD'))
  expect_equal(x$bar,c(334,445))
})

test_that("readTbl and loadTbl have the same results", {
  z <- readTbl(tbl.fn)
  loadTbl(tbl.fn)  # loads into x

  expect_equal(z, x)
})

test_that("loadTbl appends data correctly", {
  z <- x <- readTbl(tbl.fn)
  loadTbl(tbl.fn, append=TRUE)
  
  expect_equal(nrow(z)*2, nrow(x))
})