\name{geoscalePhylo}
\alias{geoscalePhylo}
\title{Adds the geological timescale to an phylogeny plot.}
\description{Plots a time-scaled phylogeny and includes the geological time-scale with a number of possible scales}
\usage{geoscalePhylo(tree, ages, upwards=FALSE, units=c("Period", "Epoch", "Age"),
  boxes="Age",tick.scale=boxes, user.scale, cex.age=0.3, cex.ts=0.3, cex.tip=0.3, 
  width=1, ts.col=TRUE, ranges=FALSE, vers="ICS2013", x.lim, quat.rm=FALSE,...)}
\arguments{
  \item{tree}{A tree as a phylo object.}
  \item{ages}{A dataset containing the first and last appearence datums,"FAD" and "LAD" respectively, of all taxa in the phylogeny. See the object $ages in data(Dipnoi) for an example.}
  \item{upwards}{If set to TRUE the tree will be plotted vertically rather than horizontal, default set to FALSE}
  \item{units}{The temporal unit(s) to be included in the timescale, options include: "Eon", "Era", "Period", "Epoch", "Age" and "User". The option "User" is used when a user-defined timescale is included. This requires an object to be assigned to user.scale.}
  \item{boxes}{Option for including grey boxes at a certain temporal resolution, options are the same as for units.}
  \item{tick.scale}{The resolution of the tick marks at the base of the timescale, the default is the same as units. The resolution of the scale can also be chosen by specifiying a value or removed entirely by using "no".} 
  \item{user.scale}{The data object to be used when including a user-defined time-scale, requires the option "User" to be included in units. See data(biozones) as an example of the required data format.}
  \item{cex.age}{Size of the text on the scale bar.}
  \item{cex.ts}{Size of the text on the geological time scale.}
  \item{cex.tip}{Size of the tip labels on the phylogeny}
  \item{width}{Width of the edges of the phylogeny.}
  \item{ts.col}{Option for using standard ICS colours on the time scale.}
  \item{ranges}{Option to add geological ranges of taxa to the phylogeny, requires an ages file with first (FAD) and last (LAD) appearence datums.}
  \item{vers}{The version of the geological time scale to use. Options include: "ICS2013", "ICS2012", "ICS2010", "ICS2009" or "ICS2008".}
  \item{x.lim}{A two item statement for the geological range, in millions of years, of the plot i.e. (0,65). If only one value is used it will be used as the upper limit.}
  \item{quat.rm}{Option to remove the names from Quaternary time bins, useful when plotting clades with long durations that range through to the recent.}
  \item{...}{All other arguments passed to plot.phylo} 
}

\author{Mark A. Bell \email{mark.bell521@gmail.com}}

\examples{

### Example lungfish data
data(Dipnoi)

tree_l <- DatePhylo(Dipnoi$tree, Dipnoi$ages, method="equal", rlen=1)
   		
geoscalePhylo(tree_l, Dipnoi$ages, ranges=TRUE, boxes="Age", cex.tip=0.4)

# Including all temporal units into the stratigraphic column

geoscalePhylo(tree_l, Dipnoi$ages, units = c("Eon","Era","Period","Epoch","Age"),
  ranges=TRUE, boxes="Age", cex.tip=0.4)

### Example trilobite data
data(Asaphidae)

tree_a <- DatePhylo(Asaphidae$trees[[1]], Asaphidae$ages, method="equal", rlen=1)

geoscalePhylo(ladderize(tree_a), Asaphidae$ages, ranges=TRUE, boxes="Age", 
	x.lim=c(504,435), cex.tip=0.5, cex.ts=0.5,vers="ICS2009") 

# Plotting the tree vertically

geoscalePhylo(ladderize(tree_a), Asaphidae$ages, ranges=TRUE, boxes="Age", 
  x.lim=c(504,435), cex.tip=0.5, cex.ts=0.5,upwards=TRUE,vers="ICS2009") 


### Including a user-defined time-scale
data(UKzones)
data(Asaphidae)

tree_a <- DatePhylo(Asaphidae$trees[[1]], Asaphidae$ages, method="equal", rlen=1)

geoscalePhylo(ladderize(tree_a), Asaphidae$ages, units = c("Eon","Era","Period",
  "Epoch","Age","User"),ranges=TRUE, boxes="Age", cex.tip=0.4,user.scale=UKzones,
  vers="ICS2009",cex.ts=0.5,x.lim=c(440),upwards=TRUE)

 
}

\keyword{time-scale, geological, plot, phylogeny, tree}
