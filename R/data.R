#' Senate Roll Call Voting Data
#'
#'
#' This dataset contains the voting records from 100 senators to 675 roll calls
#' in years 2003 and 2004.
#'
#'
#' @format A list with 4 objects:
#' \describe{
#'   \item{data}{voting matrix. 'Yea', 'Nay' and 'Not Voting', treated as 1,0 and NA, respectively.}
#'   \item{senate}{names of senators}
#'   \item{party}{parties senators belong to}
#'   \item{state}{states senators represent}
#' }
#' @source \url{https://legacy.voteview.com/dwnl.htm.}
"vote108"
