# -*- coding: utf-8 -*-
"""
MD Analysis Functions Compatible with NAMD MD Engine Output

2015 Spencer T. Stober

Version 0.1.2 ---> when this changes, change also in def initialize_data_output_file():

Version 0.1.2: changed output CSV format, "|" delimeter, better org. for post
processing

"""

import numpy as np
#import string
#import scipy.integrate
#from scipy.optimize import curve_fit
#import matplotlib.pyplot as plt
#from matplotlib.ticker import AutoMinorLocator
#from matplotlib.ticker import MultipleLocator
import csv
import datetime
import time
import os
from rdkit.Chem import AllChem as Chem
from rdkit.Chem import Descriptors as ChemDesc

def initialize_data_output_file(path):
    
    # open/append the data file and write results
    
    print('STS/MD-Processor>  Initializing Data File')    
    
    st = datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d %H:%M:%S')    
        
    txt1 = ('Ouptu> Version 0.1.2, 2015-09-15')
    txt2 = ('STS/MD-Processor> Timestamp ' + st)
    txt3 = ('STS/MD-Processor> Running in directory: ' + os.getcwd())    
    txt4 = ('UNITS | Vis(Pa) | T(K) | P(Pa) | Rho(kg/m**3) | fit_time(fs)')
    txt5 = ('UNITS | c_p(J/mol*K) | beta_t(1/Pa) | alpha_p(1/K) | P(Pa) | T(K) | Vmolar(m**3/mol) | density(kg/m**3)')
       
    print(txt1)
    print(txt2)
    print(txt3)
    print(txt4)
    print(txt5)

    f = open(path, 'a')    
    out_writer = csv.writer(f, delimiter='|',quoting=csv.QUOTE_MINIMAL)
    out_writer.writerow([txt1])
    out_writer.writerow([txt2])    
    out_writer.writerow([txt3])
    out_writer.writerow([txt4])
    out_writer.writerow([txt5])
    out_writer.writerow('')
    f.close()

    return
    
    
def output_writer(value, path):

    strngs = []
    
   
    for i in range(len(value)):
        strngs.append( str(value[i]) )
    
    # commented out as hack to make parsing easier with excel
    #strngs.append(title)

    # write the file
    f = open(path, 'a')    
    out_writer = csv.writer(f, delimiter='|',quoting=csv.QUOTE_NONE)
    out_writer.writerow(strngs)
    f.close()  

    return




def label_to_smile(out_lbl):
    # take a list of labels and use the dictionary to parse them into SMILES
    # strings    
    
    dict_loaded = np.load('chem_annotated.inp.dict.npy').item()
    
    #allocate
    out_smiles = list(np.zeros(len(out_lbl)))
    
    for i in range(len(out_lbl)):
        out_smiles[i] = dict_loaded[out_lbl[i]]
        
    return out_smiles


def rd_kit_exact_mass(smiles_list):
    # use RDKit to convert a list of SMILES strings into exact mass for
    # comparison to mass spec. analysis
    
    out_data = list(np.zeros(len(smiles_list)))
    
    for i in range(len(out_data)):
    
        # create a mol object from the smiles string
        m = Chem.MolFromSmiles(smiles_list[i])     
        # use the ExactMolWt descriptor and save to list        
        out_data[i] =   str(ChemDesc.ExactMolWt(m))
        
    return out_data
    
def rd_kit_num_rad_electrons(smiles_list):
    # use RDKit to convert a list of SMILES strings into exact mass for
    # comparison to mass spec. analysis
    
    out_data = list(np.zeros(len(smiles_list)))
    
    for i in range(len(out_data)):
    
        # create a mol object from the smiles string
        m = Chem.MolFromSmiles(smiles_list[i])     
        # use the ExactMolWt descriptor and save to list        
        out_data[i] =   str(ChemDesc.NumRadicalElectrons(m))
        
    return out_data    
    






