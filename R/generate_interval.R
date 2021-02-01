#############################################

#' Title generate_interval_integer
#' function to generate all intervals by random combination of two number from the given two number
#' the intervals must has order, the start number must smaller than the end
#' generate all non-duplicated intervals combination
#'
#' @param i the start number, must be integer, such as 1
#' @param j the end number, must be integer, such as 10
#'
#' @return interval
#' @export generate_interval_integer
#'
#' @examples
generate_interval_integer <- function(i, j){
  start_list = ""  #generate an empty list to store the output from for loop
  end_list = "" #generate an empty list to store the output from for loop
  for(i in i:j-1){
    start_i = rep(i,j-i)
    start_list <- c(start_list, start_i)
    end_i = seq(i+1, j, 1)
    end_list <- c(end_list, end_i)
  }
  interval <- unique(data.table(start_list, end_list)) #convert two vectors into data table
  interval <- interval[-c(0:j+1),] #remove the intervals not whthin the given number, it caused by the for loop
  return(interval)
}

#############################################
#' Title generate_interval_numeric
#' this function is used to generate interval with decimal
#' function to generate all intervals by random combination of two number from the given two number
#' the intervals must has order, the start number must smaller than the end
#' generate all non-duplicated intervals combination
#' @param i the start value, it must be A Decimal number, such as 0.1
#' @param j the end value, it must be A Decimal number, such sa 0.2
#'
#' @return interval
#' @export generate_interval_decimal
#'
#' @examples
generate_interval_decimal <- function(i, j){
  start_list = ""
  end_list = ""
  distance = j-i + 0.1
  i = i * 10
  j = j * 10
  for(i in i:j-1){
    start_i = rep(i,j-i)
    start_list <- c(start_list, start_i)
    end_i = seq(i+1, j, 1)
    end_list <- c(end_list, end_i)
  }
  interval <- unique(data.table(start_list, end_list))
  interval <- interval[-c(0:(distance*10+1)),] #here is the distace between two value time 10, plus 1 empty
  interval$start_list <-  as.numeric(interval$start_list)/10
  interval$end_list <- as.numeric(interval$end_list)/10
  return(interval)
}

