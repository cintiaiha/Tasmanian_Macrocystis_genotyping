#' Computes similarity coefficient for two binary variables
#'
#' @param i,j Boolean (logical) or numeric vector with 0s and 1s. i and j must be of the same length
#' @param method character string indicating which coefficient to be computed. Default to "phi". Can be abbreviated.
#'
#' @Details
#'
#' Given 2 binary (presence/TRUE/positive/1 or absence/FALSE/negative/0) variables i and j,
#' with a 2x2 contingency table or confusion matrix:
#'
#'          j                     predicted
#'        1   0            a        P   N
#'       -------           c      ---------
#'   1  | a | b |    OR    t   P | TP | FN |
#' i     -------           u      ---------
#'   0  | c | d |          a   N | FP | TN |
#'       -------           l      ---------
#'
#' Where:
#' a = number of cases on which i and j are 1 (positive matches / true positive)
#' b = number of cases where i is 1 and j is 0 (mismatches - i occurrence / false negative)
#' c = number of cases where i is 0 and j is 1 (mismatches - j occurrence / false positive)
#' d = number of cases where both i and j are 0 (negative matches / true negative)
#' a+b+c+d = n, the number of variables
#'
#' Compreensive lists of similarity metrics:
#' https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0247751
#' http://www.iiisci.org/journal/PDV/sci/pdfs/GS315JG.pdf


similarity <- function(i, j, method = c("phi", "ari", "baroni-urbani", "consonni-todeschini",
                                        "dennis", "dispersion", "hamann", "maxwell-pilliner",
                                        "michael", "peirce1", "peirce2", "peirce3", "scott",
                                        "tarantula", "yuleQ", "yuleW") ) {

  if (anyNA(i)) warning("i contains NA")
  if (anyNA(j)) warning("j contains NA")
  if (length(i) != length(j)) stop("i and j must have the same length")
  if (is.logical(i)) i <- as.integer(i)
  if (is.logical(j)) i <- as.integer(j)
  if (!all(i %in% c(0:1, NA)) | !all(j %in% c(0:1, NA))) stop("i and j must be logical or 0/1")
  if (!is.numeric(method)) method <- match.arg(method)

  tbl <- table(i,j)
    a <- tbl[2,2]
    b <- tbl[2,1]
    c <- tbl[1,2]
    d <- tbl[1,1]
    n <- sum(tbl) # a+b+c+d

  # ARI
  if (method == "ari") return( (n*(a+d)-((a+b)*(a+c)+(c+d)*(b+d)))/(n^2-((a+b)*(a+c)+(c+d)*(b+d))) )

  # Baroni-Urbani and Buser II
  if (method == "baroni-urbani-buser2") return( (sqrt(a*d)+a-b-c)/(sqrt(a*d)+a+b+c) )

  # Consonni and Todeschini
  if (method == "consonni-todeschini") return( (log(1+a*d)-log(1+b*c))/log(1+n^2/4) )

  # Dennis
  if (method == "dennis") return( (a*d-b*c)/sqrt(n*(a+b)*(a+c)) )

  # Dispersion
  if (method == "dispersion") return( (a*d-b*c)/n^2 )

  # Hamann
  if (method == "hamann") return( (a+d-b-c)/n )

  # Maxwell and Pilliner
  if (method == "maxwell-pilliner") return( 2*(a*d-b*c)/((a+b)*(c+d)+(a+c)*(b+d)) )

  # Michael
  if (method == "michael") return( 4*(a*d-b*c)/((a+d)^2+(b+c)^2) )

  # Peirce I
  if (method == "peirce1") return( (a*d-b*c)/((a+b)*(c+d)) )

  # Peirce II
  if (method == "peirce2") return( (a*d-b*c)/((a+c)*(b+d)) )

  # Peirce III (0:1)
  if (method == "peirce3") return( (a*b+b*c)/(a*b+2*b*c+c*d) )

  # Phi / Matthews Correlation Coefficient (MCC)
  if (method == "phi" | method == "mcc") return( (a*d-b*c)/sqrt((a+b)*(a+c)*(b+d)*(c+d)) )

  # Scott
  if (method == "scott") return( (4*a*d-(b+c)^2)/((2*a+b+c)*(2*d+b+c)) )

  # Tarantula (0:1)
  if (method == "tarantula") return( (a*(c+d))/(c*(a+b)) )

  # Yule’s Q
  if (method == "yuleQ") return( (a*d-b*c)/(a*d+b*c) )

  # Yule’s W
  if (method == "yuleW") return( (sqrt(a*d)-sqrt(b*c))/(sqrt(a*d)+sqrt(b*c)) )
}
