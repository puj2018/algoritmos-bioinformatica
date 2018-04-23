#!/usr/bin/Rscript

#--------------------------------------------------
# A Hidden Markov Model of DNA sequence evolution 
#--------------------------------------------------
# The transition matrix and emission matrix for a HMM
nucleotides   <- c("A", "C", "G", "T")   # Define the alphabet of nucleotides
states        <- c("H", "L")             # Define the names of the states Hight (CG-rich) and Low (AT-rich)

#--------------------------------------------------------------
# Transitions
#--------------------------------------------------------------
GCHight       <- c(0.9, 0.1)             # Set the probabilities of switching states, where the previous state was "AT-rich"
GCLow         <- c(0.5, 0.5)             # Set the probabilities of switching states, where the previous state was "GC-rich"
transitions <- matrix(c(GCHight, GCLow), 2, 2, byrow = TRUE) # Create a 2 x 2 matrix
rownames(transitions) <- states
colnames(transitions) <- states

cat (">>> The transition matriz for the High GC and Low GC:\n")
print (transitions)

#--------------------------------------------------------------
# Emissions
#--------------------------------------------------------------
highEmissions    <- c(0.1, 0.4, 0.4, 0.1) # Set the values of the probabilities, for the GC-rich state
lowEmissions    <- c(0.4, 0.1, 0.1, 0.4) # Set the values of the probabilities, for the AT-rich state
emissions <- matrix(c(highEmissions, lowEmissions), 2, 4, byrow = T) # Create a 2 x 4 matrix
rownames(emissions) <- states
colnames(emissions) <- nucleotides
 
cat ("\n>>> The HMM emission matrix:\n")
print (emissions)

#-----------------------------------------------------------------------
# This fill the viterbi graph using a matrix vkMatrix 
# This function is called by the viterbi function (below)
#-----------------------------------------------------------------------
fillViterbiGraph <- function(sequence, transitions, emissions) {
	# Find out how many states are in the HMM
	numstates       <- dim(transitions)[1]
	numObservations <- length (sequence)
	# Make a matrix with as many rows as positions in the sequence, and as many
	# columns as states in the HMM
	vkMatrix <- matrix(NA, nrow=numstates, ncol=numObservations )
	# Set the values in the first column of matrix vkMatrix (representing the first position of the sequence) to 0
	vkMatrix [,1] <- 0
	# Set the value in the first row of matrix vkMatrix, first column to 1
	vkMatrix [1,1] <- 1
	# Fill in the matrix vkMatrix:

	for (i in 2:numObservations) { # For each position in the DNA sequence: 
		for (l in 1:numstates) {  # For each of the states of in the HMM:
		   # Find the probabilility, if we are in state l, of choosing the nucleotide at position in the sequence
		   emission_li <- emissions[l,sequence[i]]

		   vkMatrix[l,i] <-  emission_li * max(vkMatrix[,(i-1)] * transitions[,l])
		}
	}
	return(vkMatrix)
}

	
#-----------------------------------------------------------------------
# This carries out the Viterbi algorithm by calling the fillViterbiGraph function
#-----------------------------------------------------------------------
viterbi <- function(sequence, transitions, emissions) {
     theMostProbPath = character()
     # Get the names of the states in the HMM:
     states <- rownames(emissions)
	 numObservations <- length (sequence)

     # Make the Viterbi matrix vkMatrix:
     vkMatrix <- fillViterbiGraph(sequence, transitions, emissions)

     # Go through each of the rows of the matrix vkMatrix (where each row represents
     # a position in the DNA sequence), and find out which column has the
     # maximum value for that row (where each column represents one state of
     # the HMM):
     mostprobablestatepath <- apply(vkMatrix, 2, function(x) which.max(x))

     # Print out the most probable state path:
     prevnucleotide <- sequence[1]
     prevmostprobablestate <- mostprobablestatepath[1]
     prevmostprobablestatename <- states[prevmostprobablestate]
	 theMostProbPath = c(theMostProbPath, prevmostprobablestatename)
     startpos <- 1
     for (i in 2:numObservations) {
        nucleotide <- sequence[i]
        mostprobablestate <- mostprobablestatepath[i]
        mostprobablestatename <- states[mostprobablestate]
        if (mostprobablestatename != prevmostprobablestatename) {
           print(paste("Positions",startpos,"-",(i-1), "Most probable state = ", prevmostprobablestatename))
           startpos <- i
        }
        prevnucleotide <- nucleotide
        prevmostprobablestatename <- mostprobablestatename
		theMostProbPath = c(theMostProbPath, mostprobablestatename)
     }
     print(paste("Positions",startpos,"-",i, "Most probable state = ", prevmostprobablestatename))
	 return (theMostProbPath)
}


#--------------------------------------------------------------
# Main Test
#--------------------------------------------------------------
stringSequence = "CCCGGGGTTTCCC"
sequence     = strsplit (stringSequence, split="")[[1]]

# Call to Viterbi to compute the most probable path
theMostProbPath = viterbi (sequence, transitions, emissions)

cat ("\n>>> The most probable path is:\n")
cat (sequence)
cat ("\n")
cat (theMostProbPath)
cat ("\n")
