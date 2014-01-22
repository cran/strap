\name{geoscalePhylo}
\alias{geoscalePhylo}
\title{Adds the geological timescale to an phylogeny plot.}
\description{Plots a time-scaled phylogeny and includes the geological time-scale with a number of possible scales}
\usage{geoscalePhylo(tree,ages, cex.age=0.3, cex.ts=0.3, cex.tip=0.3, width=1,
  ranges=FALSE, scale=boxes, units=c("Age", "Epoch", "Period"), ts.col=TRUE,
  boxes="Age", nscale, x.lim, vers="ICS2013",...)}
\arguments{
  \item{tree}{A tree as a phylo object.}
  \item{ages}{The second dataset, containing the upper and lower ages of taxa.}
  \item{cex.age}{Size of the text on the scale bar.}
  \item{cex.ts}{Size of the text on the geological time scale.}
  \item{cex.tip}{Size of the tip labels on the phylogeny}
  \item{width}{Width of the edges of the phylogeny.}
  \item{ranges}{Option to add geological ranges to the phylogeny, requires an ages file with first (FAD) and last (LAD)  appearence datums.}
  \item{scale}{The resolution of the tick marks at the base of the timescale, the default is the same as boxes.}
  \item{units}{The different levels to be used in the timescale, options include: "Eon", "Era", "Period", "Epoch", "Age" and "Other".}
  \item{ts.col}{Option for using standard ICS colours on the time scale.}
  \item{boxes}{Option for including grey boxes at a certain temporal resolution, options are the same as for units.}
  \item{nscale}{If the option "Other" is selected in source then the new scale can be added here, requires "Start", "End", "Midpoint" and "Name" fields, see data(biozones) as an example.}
  \item{x.lim}{A two item statement for the x limit, in millions of years, of the plot i.e. (0,65), if only one value is used it will be used as the upper limit.}
  \item{vers}{The version of the time scale to use. Options include: "ICS2013", "ICS2012", "ICS2010", "ICS2009" or "ICS2008".}
  \item{...}{All other arguments passed to plot.phylo}
  
}

\author{Mark A. Bell \email{mark.bell521@gmail.com}}

\examples{

# Example lungfish data
data(Dipnoi)

tree_l <- DatePhylo(Dipnoi$tree, Dipnoi$ages, method="equal", rlen=1)
   		
geoscalePhylo(tree_l, Dipnoi$ages, ranges=TRUE, boxes="Age", cex.tip=0.4)

# Example trilobite data
data(Asaphidae)

tree_a <- DatePhylo(Asaphidae$trees[[1]], Asaphidae$ages, method="equal", rlen=1)

geoscalePhylo(ladderize(tree_a), Asaphidae$ages, ranges=TRUE, boxes="Age", 
	x.lim=c(504,440), cex.tip=0.5, cex.ts=0.5) 
}

\keyword{time-scale, geological, plot, phylogeny, tree}
