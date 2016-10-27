#!/usr/bin/env python

'''
Andras Palfi
10/23/2016
Program 3 Rosalind
A solution to the Rosalind problem: 
Inferring Protein from Spectrum
'''
 
import numpy as np

ProteinMassInput=[]
#import the Rosalind data file
with open("/Users/Andras/Documents/UConn/2016-Fall/CSE-3800/rosalind_spec.txt") as input:
	pml = input.readlines()


#Clean up the data, convert to floats,  and append the values to the mass list that 
# will be used for matching:
for lines in pml:
	strip = lines.strip()
	floatConv = float(strip)
	ProteinMassInput.append(floatConv)


#Build the list of Monoisotopic masses provided by Rosalind
ProteinMassData=np.array([('A',71.03711),('C',103.00919), ('D',115.02694) , 
	('E',129.04259), ('F',147.06841), ('G',57.02146), ('H',137.05891) , 
	('I',113.08406), ('K',128.09496), ('L',113.08406) , ('M',131.04049) , 
	('N',114.04293), ('P',97.05276), ('Q',128.05858), ('R',156.10111) , 
	('S',87.03203), ('T',101.04768), ('V',99.06841), ('W',186.07931) , 
	('Y',163.06333)])

#Get the length of the data, useful for while loop later
lengthPMD = len(ProteinMassData)

#Here I take the ProteinMassData array and simply strip the characters.
#This leaves just the weight values which makes it a little bit easier 
#later on
weights=[]
i=0
while(i <= lengthPMD-1):
	weightPOS = ProteinMassData[i][1]
	weights.append(float(weightPOS))
	i= i+1



lengthPMI = len(ProteinMassInput)

j=0
#Initilaze a final string that is empty at first
FINAL_STRING=''

#The main while loop where all the magic happens
while (j <= lengthPMI -2 ):
	#We extract the first protein in position j
	PROTEIN_1 = ProteinMassInput[j]
	#We extract the second protein in position j+1 
	PROTEIN_2 = ProteinMassInput[j+1]
	#The prefix weight is the difference between the 2nd protein and the 1st one
	prefix = PROTEIN_2 - PROTEIN_1
	
	#Create an empty list to caluclate the BEST fit for each one
	compare_weights=[]
	'''
	We do this by using a for loop that goes through each protein weight and subtracting
	it from the prefix we found in the last step. All of these values are saved to the 
	compare_weight list. After this list is compiled, we find the closest match which by definitiom
	is the smallest value in this array. 

	Once we find this "smallest value" we know that its index in the weights list 
	corresponds to the same index in the ProteinMassData array. Once we have that, it is 
	trivial to find the correct amino acid residue character and append it to FINAL_STRING.

	'''
	for eachWeight in weights:
		diff = abs(prefix - eachWeight)
		compare_weights.append(diff)
	closest_match = min(compare_weights)
	closest_match_index = compare_weights.index(closest_match)
	protein_char = str(ProteinMassData[closest_match_index][0])
	FINAL_STRING = FINAL_STRING + protein_char
	j = j+1
print(FINAL_STRING)
