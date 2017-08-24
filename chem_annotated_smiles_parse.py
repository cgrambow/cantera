# -*- coding: utf-8 -*-


import numpy as np
import re



'''

S.T. Stober 2017

##### IMPORTANT ######

If in the chem_annotated.inp there are "non-SMILES" entries, e.g.:

SPECIES
    dodecane            ! dodecane
    Ar                  ! Ar
    He                  ! He
    Ne                  ! Ne
    N2                  ! N2
    HO2(1)              ! hydroperoxyradical(1)
    HO(2)               ! hydroxylradical(2)
    oxygen(3)           ! oxygen(3)
    C16H32(4)           ! reactant_a(4)
    C16H31(14)          ! CCCCCCCCCC[C]1CCCCC1(14)
    
you should first fill in any SMILES strings so that the dictionary
is built correctly.  E.g.:

SPECIES
    dodecane            ! CCCCCCCCCCCC
    Ar                  ! Ar
    He                  ! He
    Ne                  ! Ne
    N2                  ! N#N
    HO2(1)              ! [O]O(1)
    HO(2)               ! [OH](2)
    oxygen(3)           ! O=O(3)
    C16H32(4)           ! CCCCCCCCCCC1CCCCC1(4)
    C16H31(14)          ! CCCCCCCCCC[C]1CCCCC1(14)


To open the npy file:
dict_loaded = np.load('chem_annotated.inp.dict.npy').item()




###################### How it works

1) Read lines until you get to "SPECIES"
2) Parse, a line that looks like this:

S(194)              ! CC(CCCCCCCCC1CCCCC1)OO(194)

3) Remove trailing "(194)"
4) Save the output as a dictionary that converts the species name to SMILES

'''



# Find the line number with the tag "SPECIES"
f = open('chem_annotated.inp','r')

# create the dictionary
species_d = {'species_name':'SMILES_string'}

raw = f.readline()
print 'Line read: ' + raw

# keep reading lines until we find SPECIES
while raw[-1]== '\n':
    if raw == 'SPECIES\n' :
        print 'Found species'

        ### SPECIES found... here goes the parser ###

        # read the next line     
        raw = f.readline()
        print 'Found species: ' + raw
        # if it is not "END", parse it        
        while raw != 'END\n' :
            # we found a valid line ########### main parser
            print 'Parsing line: ' + raw

            # split into a list of two strings according to '!'
            out1 = re.split('!', raw)
            
            # strip all the spaces from the species key
            out_species_key = out1[0].strip()
            print 'Species key parsed is: ' + out_species_key            

            # strip all the spaces from the SMILES
            out1_strip = out1[1].strip()
            
            # match any number of numbers within ()
            li = re.split('\([0-9]*\)', out1_strip)
            out_smiles_value = li[0]
            print 'Species smiles value is: ' + out_smiles_value
            
            # add it to the dictionary
            species_d[out_species_key] =  out_smiles_value
            
            print 'Parsed... reading new line'
            raw = f.readline()
            print 'Read line: ' + raw
        

        print 'End of species reached... saving dictionary'
        np.save('chem_annotated.inp.dict',species_d)
        break

    else : 
        raw = f.readline()
        print 'Line read: ' + raw































