

#' plot posterior bars and compare to correct or incorrect
#' 
#' Makes a nice plot. 
#' @param SS vector of sibsizes
#' @param post  vector of posteriors associated with each sibship (between 0 and 1)
#' @param false_positive A logical vector indicating whether the sibship included
#' any non-sibs
#' @param lscale how big should the lengths be on the y-scale. i.e. one sibling = lscale
#' @param correct  logical vector saying whether each inferred sibship is correct or not
#' @param add TRUE means add to existing plot
#' @export
sib_result_plot <- function(SS, post, false_positive, lscale = .01, correct, add = FALSE) {
  if(length(SS) != length(post)) stop("length mismath between SS and post")
  L <- length(SS)
  
  # tops and bottoms of our segments
  bots <- post - lscale * SS / 2
  tops <- post + lscale * SS / 2
  
  # set up plot area
  if(add == FALSE) plot(c(1,L), c(0,1.2), type="n", 
                        ylab = "Estimated Posterior Probability of Inferred Sibling Group", 
                        xlab = "Index of Inferred Sibling Group (Ordered by Posterior, Size)")
  
  # plot the segments 
  segments(x0 = 1:L, y0 = bots, x1 = 1:L, y1 = tops, col=c("gray", "blue", "red")[(!correct) + 1 + false_positive])
  lines(1:L, post)
  
  # put the legend in the lower right
  legend("bottomleft", 
         legend = c("Correctly-inferred sibling group", "Incomplete sibling group", "Grouping of non-siblings"), 
         col = c("gray", "blue", "red"), lty = "solid"
         )
}

