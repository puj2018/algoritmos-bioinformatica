#!/usr/bin/Rscript

# Code that shows how to create a HMM and define a function to
# generate a sequence of observationa and states based on the
# model parameters (transitions, emissions, and initials)

nucleotides <- c ("A","C","G","T")
states = c("H","L")
initials = c(0.5, 0.5)
names (initials) = states

# Create Transitions Matrix
GCHight = c(0.9,0.1)
GCLow = c(0.5, 0.5)
transitions = matrix (c (GCHight, GCLow), 2,2, byrow =TRUE)
rownames (transitions) = states
colnames (transitions) = states

# Create Emissions Matrix
highEmissions = c (0.1, 0.4, 0.4, 0.1)
lowEmissions = c (0.4,0.1,0.1,0.4)
emissions = matrix (c (highEmissions, lowEmissions), 2, 4, byrow=T)
rownames (emissions) = states
colnames (emissions) = nucleotides

# Test sequence 
stringSequence = "CGCAAAGTTCG"
sequence = strsplit (stringSequence,split="")[[1]]
nRows = length (states)   # Number of rows
nCols = length (sequence) # Number of columns

vkMatrix  = matrix (NA, nRows, nCols) = pathMatrix
rownames (vkMatrix) = states
colnames (vkMatrix) = sequence

#----------------------------------------------------------------------
# Recursive  Viterbi Function for calculating the most probable path
# ending at "l" at position "i"
#----------------------------------------------------------------------
Vk <- function (l, i) {
	cat (">>> ",l,i,">>>\n")
	if (i==1) {      #V_0(0) = 1
		value = emissions [l,i]*initials [l]
		vkMatrix [l,i] <<- value
		return (value)
	} else {          #V_l(i) = E_l(i)*max (V_k(i-1)*a_li)
		emission_li = emissions [l, sequence [i]]

		# Get the Vk max calling the Vk(i-1)
		vkMax =0
		maxK = 0
		n = length (states)
		for (k in states) {
			vkTemp = Vk (k, i-1)*transitions[k,l]
			cat ("vtTemp:",vkTemp,">\n")
			cat  ("\n>",k,l,"<\n")
			cat ("vtMax:",vkMax,">\n")
			if (vkTemp > vkMax) {
				vkMax = vkTemp
				vkMatrix [l,i] <<- round (vkMax,7)
			}
		}
		return (emission_li*vkMax)
	}
}
viterbiRecursive <- function (sequence) {
	n = length (sequence)
	cat ("n:")
	print (n)
	print (sequence)
	for (s in states) {
		valVk = Vk (s, n, sequence)
		cat ("\n")
		print ("Recursive: ")
		cat ("State ", s, valVk,"\n")
		vkMatrix [s, n] <<- round (valVk,7)
	}

	mostProbPathValues = rep (0, n)
	mostProbPathStates = rep (0, n)
	for (i in n:1) {
		mostProbPathValues [i] <- max (round (vkMatrix [,i],6))
		st = which.max (vkMatrix [,i])
		cat ("pm:",format (vkMatrix),">\n")
		print (vkMatrix)
		cat ("st:",st,">\n")
		cat ("i:",i,">\n")
		cat ("mostP:",mostProbPathStates,">\n")
		mostProbPathStates [i] <- states[st]
	}
	return (list (mostProbPathValues, mostProbPathStates))
}


cat ("\nTransitions:\n")
print (transitions)
cat ("\nEmissions:\n")
print (emissions)

#N = length (sequence)
#probVkH = Vk ("H", 2)

#cat ("\nPath Matrix\n")
#print (vkMatrix)
#cat ("\n Prob of ending in High", probVkH)
#cat ("\n", sequence)
#cat ("\n")
#print (vkMatrix)

#probVkL = Vk ("L", 2)
#cat ("\n", probVkL)


probPathStates = viterbiRecursive (sequence)
cat ("\nMaxVk:\n")
print (vkMatrix)
cat ("\n")
print (sequence)
cat ("\n")
print (probPathStates[[1]])
cat ("\n")
print (probPathStates [[2]])
