from mp_api import MPRester
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import math
from math import gcd
from functools import reduce
import chemparse as cp


def download_oxides_from_matproj(MAPI_KEY, elements_to_avoid, file_name):
    mpr = MPRester(MAPI_KEY)

    data = mpr.summary.search(exclude_elements = elements_to_avoid, 
                              elements = ['O'], 
                              nelements = [3],
                              crystal_system = ['Cubic'],
                              fields =['material_id',
                                       'formula_pretty',
                                       'theoretical',
                                       'nsites',
                                       'decomposes_to',
                                       'is_stable',
                                       'n',
                                       'band_gap',
                                       'density',
                                       'k_vrh',
                                       'homogeneous_poisson',
                                       'formula_anonymous'])

    data2 = mpr.summary.search(exclude_elements = elements_to_avoid, 
                              elements = ['O'], 
                              nelements = [2],
                              crystal_system = ['Cubic'],
                              fields =['material_id',
                                       'formula_pretty',
                                       'theoretical',
                                       'nsites',
                                       'decomposes_to',
                                       'is_stable',
                                       'n',
                                       'band_gap',
                                       'density',
                                       'k_vrh',
                                       'homogeneous_poisson',
                                       'formula_anonymous'])

    myfile = 'misc_data.txt'
    pdf = pd.read_csv(myfile,sep='\t')
    k_B = 1.38064852*1e-23

    d = {'film':[],'film ID':[],'ICSD?':[],'band gap (eV)':[],'n':[], 'stability':[],'B':[],'mu':[], 'nsites':[],'elements':[],'stoich':[], 'MW':[], 'density':[],'CTE':[], 'x-tal':[]}
    for ii in range(0,len(data)):
      # organize films names 
      d['film'].append(data[ii].formula_pretty)
      d['film ID'].append(data[ii].material_id[:])

      # Is the material in the ICSD database?
      if data[ii].theoretical:
        d['ICSD?'].append('false')
      if not data[ii].theoretical:
        d['ICSD?'].append('true')

      # bandgap
      d['band gap (eV)'].append(data[ii].band_gap)

      # refractive index
      if data[ii].n:
        d['n'].append(round(data[ii].n,3))
      if not data[ii].n:
        d['n'].append(0)

      # bulk modulus
      if data[ii].k_vrh:
        d['B'].append(data[ii].k_vrh)
      if not data[ii].k_vrh:
        d['B'].append(0)

      # Poisson ratio
      if data[ii].homogeneous_poisson:
        d['mu'].append(data[ii].homogeneous_poisson)
      if not data[ii].homogeneous_poisson:
        d['mu'].append(0)

      # stoichiometry, molecular weight, denisty and CTE
      element_dict = cp.parse_formula(data[ii].formula_pretty)
      d['elements'].append(list(element_dict.keys()))
      stoich = []
      weight = []
      for key in element_dict.keys():
        stoich.append(element_dict[key])
        elmt_idx = pdf.loc[pdf['Element'] == key].index[0]
        weight.append(float(pdf['Weight'][elmt_idx]))
      d['stoich'].append(stoich)
      d['nsites'].append(sum(stoich))
      d['MW'].append(np.sum(np.array(weight)*np.array(stoich)))
      d['density'].append(data[ii].density)
      N = d['nsites'][ii]*6.02214076*1e23
      if data[ii].k_vrh:
        d['CTE'].append(((3*N*k_B*(1+d['mu'][ii])/(2 - 3*d['mu'][ii]))/(2*(d['MW'][ii]/d['density'][ii])*d['B'][ii])/1000)*1e6)
      if not data[ii].k_vrh:
        d['CTE'].append(0)
    
      

      # crystal system
      if np.array_equal(np.sort(stoich),[1, 1, 3]) is True:
        d['x-tal'].append('perovskite')
      elif np.array_equal(np.sort(stoich),[1, 3, 1]) is True:
        d['x-tal'].append('spinel')
      elif np.array_equal(np.sort(stoich),[1, 2, 4]) is True:
        d['x-tal'].append('spinel')
      elif np.array_equal(np.sort(stoich),[2, 2, 7]) is True:
        d['x-tal'].append('pyrochlore')
      else:
        d['x-tal'].append('unknown')

      # stability
      if data[ii].decomposes_to is None:
        d['stability'].append('stable')
      else:
        d['stability'].append('metastable')


    for ii in range(0,len(data2)):
      # organize films names 
      d['film'].append(data2[ii].formula_pretty)
      d['film ID'].append(data2[ii].material_id[:])

      # Is the material in the ICSD database?
      if data2[ii].theoretical:
        d['ICSD?'].append('false')
      if not data2[ii].theoretical:
        d['ICSD?'].append('true')

      # bandgap
      d['band gap (eV)'].append(data2[ii].band_gap)

      # refractive index
      if data2[ii].n:
        d['n'].append(round(data2[ii].n,3))
      if not data2[ii].n:
        d['n'].append(0)

      # bulk modulus
      if data2[ii].k_vrh:
        d['B'].append(data2[ii].k_vrh)
      if not data2[ii].k_vrh:
        d['B'].append(0)

      # Poisson ratio
      if data2[ii].homogeneous_poisson:
        d['mu'].append(data2[ii].homogeneous_poisson)
      if not data2[ii].homogeneous_poisson:
        d['mu'].append(0)

      # stoichiometry, molecular weight, denisty and CTE
      element_dict = cp.parse_formula(data2[ii].formula_pretty)
      d['elements'].append(list(element_dict.keys()))
      stoich = []
      weight = []
      for key in element_dict.keys():
        stoich.append(element_dict[key])
        elmt_idx = pdf.loc[pdf['Element'] == key].index[0]
        weight.append(float(pdf['Weight'][elmt_idx]))
      d['stoich'].append(stoich)
      d['nsites'].append(sum(stoich))
      d['MW'].append(np.sum(np.array(weight)*np.array(stoich)))
      d['density'].append(data2[ii].density)
      N = d['nsites'][ii]*6.02214076*1e23 
      if data2[ii].k_vrh:
        d['CTE'].append(((3*N*k_B*(1+d['mu'][ii])/(2 - 3*d['mu'][ii]))/(2*(d['MW'][ii]/d['density'][ii])*d['B'][ii])/1000)*1e6)
      if not data2[ii].k_vrh:
        d['CTE'].append(0)
        
      # crystal system
      if np.array_equal(np.sort(stoich),[1, 1]) is True:
        d['x-tal'].append('rocksalt')
      elif np.array_equal(np.sort(stoich),[1, 2]) is True:
        d['x-tal'].append('fluorite')
      else:
        d['x-tal'].append('unknown')

      # stability
      if data2[ii].decomposes_to is None:
        d['stability'].append('stable')
      else:
        d['stability'].append('metastable')

    df = pd.DataFrame(d)

    d = {'a':[], 'b':[],'c':[]}
    for ii in range(0, len(df)):
      struct = mpr.get_structure_by_material_id(df['film ID'][ii], final=True, conventional_unit_cell=True)
      d['a'].append(round(struct.lattice.abc[0],3))
      d['b'].append(round(struct.lattice.abc[1],3))
      d['c'].append(round(struct.lattice.abc[2],3))

    df['a'] = d['a']
    df['b'] = d['b']
    df['c'] = d['c']

    df.to_pickle(file_name)
    return

def find_gcd(list):
    x = reduce(gcd, list)
    return x

def find_films(file_name, substrate, sub_xtal_type,lattice_parameter, orientation, strain, CTE):
    df = pd.read_pickle(file_name)
    df.drop('b', inplace=True, axis=1)
    df.drop('c', inplace=True, axis=1)
    df = df[~(df['x-tal']==sub_xtal_type)]
    df=df.reset_index(drop=True)

    # 5 --- calculate all possible bonding types and minimum strain
    d = {'film orientation':[], 'film lattice parameter':[],'% strain':[]}
    for ii in range(0, len(df)):
      lattice_matching = find_xtal_match(sub_xtal_type, lattice_parameter,orientation, df['x-tal'][ii], df['a'][ii])
      d['film orientation'].append(lattice_matching[2])
      d['film lattice parameter'].append(lattice_matching[0])
      d['% strain'].append(lattice_matching[1])
    df['film orientation'] = d['film orientation']
    df['film lattice parameter'] = d['film lattice parameter']
    df['% strain'] = d['% strain']
    df = df.loc[(df['% strain'] < strain)]
    df = df.reset_index(drop=True)
    df = df.sort_values(by=['ICSD?','% strain', 'stability','n','band gap (eV)'] , ascending=[False,True,False, False, False])
    df = df.reset_index(drop=True)

    df.drop('film ID', inplace=True, axis=1)
    df.drop('stoich', inplace=True, axis=1)
    df = df[~(df['x-tal']== 'unknown')]
    df = df.reset_index(drop=True)
    #df.to_csv(r'sub_'+substrate+'_'+orientation+'.csv', index = False, header=True)
    pd.options.display.max_columns = 260
    pd.options.display.max_rows = 260
    return df
    
def find_xtal_match(sub_type, sub_param, sub_orient, film_type, film_param):
    if sub_orient == '(100)':
        
         if sub_type == 'rocksalt':
            if film_type == 'fluorite':
                 orient = ['(110)']
            elif film_type == 'perovskite':
                orient = ['(100)', '(110)', '(111)']
            elif film_type == 'spinel':
                orient = ['(100)', '(110)']
            elif film_type == 'pyrochlore':
                orient = ['(110)', '(111)']
            else: 
                orient = ['(100)', '(110)', '(111)']

         elif sub_type == 'fluorite':
            if film_type == 'rocksalt':
                 orient = ['(111)']
            elif film_type == 'perovskite':
                orient = ['(110)', '(111)']
            elif film_type == 'spinel':
                orient = ['(100)', '(111)']
            elif film_type == 'pyrochlore':
                orient = ['(100)', '(111)']
            else: 
                orient = ['(100)', '(110)', '(111)']

         elif sub_type == 'perovskite':
            if film_type == 'rocksalt':
                 orient = ['(100)', '(110)']
            elif film_type == 'fluorite':
                orient = ['(110)']
            elif film_type == 'spinel':
                orient = ['(100)', '(110)']
            elif film_type == 'pyrochlore':
                orient = ['(110)', '(111)']
            else: 
                orient = ['(100)', '(110)', '(111)']

         elif sub_type == 'spinel':   
            if film_type == 'rocksalt':
                 orient = ['(100)', '(110)', '(111)']
            elif film_type == 'fluorite':
                orient = ['(100)', '(110)', '(111)']
            elif film_type == 'perovskite':
                orient = ['(100)', '(110)', '(111)']
            elif film_type == 'pyrochlore':
                orient = ['(100)', '(110)', '(111)']
            else: 
                orient = ['(100)', '(110)', '(111)']

         elif sub_type == 'pyrochlore':
            if film_type == 'rocksalt':
                 orient = ['(111)']
            elif film_type == 'fluorite':
                orient = ['(100)', '(111)']
            elif film_type == 'perovskite':
                orient = ['(110)', '(111)']
            elif film_type == 'spinel':
                orient = ['(100)', '(111)']
            else: 
                orient = ['(100)', '(110)', '(111)']
                
    elif sub_orient == '(110)':
        
         if sub_type == 'rocksalt':
            if film_type == 'fluorite':
                 orient = ['(110)']
            elif film_type == 'perovskite':
                orient = ['(100)', '(110)', '(111)']
            elif film_type == 'spinel':
                orient = ['(100)', '(110)']
            elif film_type == 'pyrochlore':
                orient = ['(110)', '(111)']
            else: 
                orient = ['(100)', '(110)', '(111)']
                
         elif sub_type == 'fluorite':
            if film_type == 'rocksalt':
                 orient = ['(100)','(110)']
            elif film_type == 'perovskite':
                orient = ['(100)','(110)', '(111)']
            elif film_type == 'spinel':
                orient = ['(100)', '(110)']
            elif film_type == 'pyrochlore':
                orient = ['(110)', '(111)']
            else: 
                orient = ['(100)', '(110)', '(111)']
                
         elif sub_type == 'perovskite':
            if film_type == 'rocksalt':
                 orient = ['(100)', '(110)','(111)']
            elif film_type == 'fluorite':
                orient = ['(100)', '(110)','(111)']
            elif film_type == 'spinel':
                orient = ['(100)', '(110)','(111)']
            elif film_type == 'pyrochlore':
                orient = ['(100)', '(110)','(111)']
            else: 
                orient = ['(100)', '(110)', '(111)']

         elif sub_type == 'spinel':   
            if film_type == 'rocksalt':
                 orient = ['(100)', '(110)']
            elif film_type == 'fluorite':
                orient = ['(110)']
            elif film_type == 'perovskite':
                orient = ['(100)', '(110)', '(111)']
            elif film_type == 'pyrochlore':
                orient = ['(110)', '(111)']
            else: 
                orient = ['(100)', '(110)', '(111)']

         elif sub_type == 'pyrochlore':
            if film_type == 'rocksalt':
                 orient = ['(100)','(110)']
            elif film_type == 'fluorite':
                orient = ['(110)']
            elif film_type == 'perovskite':
                orient = ['(100)','(110)', '(111)']
            elif film_type == 'spinel':
                orient = ['(100)','(110)', '(111)']
            else: 
                orient = ['(100)', '(110)', '(111)']   

    elif sub_orient == '(111)':
         if sub_type == 'rocksalt':
            if film_type == 'fluorite':
                 orient = ['(100)', '(111)']
            elif film_type == 'perovskite':
                orient = ['(110)', '(111)']
            elif film_type == 'spinel':
                orient = ['(100)', '(111)']
            elif film_type == 'pyrochlore':
                orient = ['(100)', '(111)']
            else: 
                orient = ['(100)', '(110)', '(111)']

         elif sub_type == 'fluorite':
            if film_type == 'rocksalt':
                 orient = ['(111)']
            elif film_type == 'perovskite':
                orient = ['(110)', '(111)']
            elif film_type == 'spinel':
                orient = ['(100)', '(111)']
            elif film_type == 'pyrochlore':
                orient = ['(100)', '(111)']
            else: 
                orient = ['(100)', '(110)', '(111)']

         elif sub_type == 'perovskite':
            if film_type == 'rocksalt':
                 orient = ['(100)', '(110)','(111)']
            elif film_type == 'fluorite':
                orient = ['(100)', '(110)','(111)']
            elif film_type == 'spinel':
                orient = ['(100)', '(110)','(111)']
            elif film_type == 'pyrochlore':
                orient = ['(100)', '(110)','(111)']
            else: 
                orient = ['(100)', '(110)', '(111)']

         elif sub_type == 'spinel':   
            if film_type == 'rocksalt':
                 orient = ['(111)']
            elif film_type == 'fluorite':
                orient = ['(100)', '(111)']
            elif film_type == 'perovskite':
                orient = ['(110)', '(111)']
            elif film_type == 'pyrochlore':
                orient = ['(100)', '(110)', '(111)']
            else: 
                orient = ['(100)', '(110)', '(111)']

         elif sub_type == 'pyrochlore':
            if film_type == 'rocksalt':
                 orient = ['(100)', '(110)', '(111)']
            elif film_type == 'fluorite':
                orient = ['(100)', '(110)', '(111)']
            elif film_type == 'perovskite':
                orient = ['(100)', '(110)', '(111)']
            elif film_type == 'spinel':
                orient = ['(100)', '(110)', '(111)']
            else: 
                orient = ['(100)', '(110)', '(111)'] 
    
    lp = np.array([])
    final_orient = []
    if '(100)' in orient:
        lp = np.append(lp,[0.5, 1, 2])
        final_orient.extend(['(100)/2','(100)','(100)*2'])
    if '(110)' in orient:
        lp = np.append(lp,[0.5*math.sqrt(2), 1*math.sqrt(2), 2*math.sqrt(2)])
        final_orient.extend(['(110)/2','(110)','(110)*2'])
    if '(111)' in orient:
        lp = np.append(lp,[0.5*math.sqrt(3), 1*math.sqrt(3), 2*math.sqrt(3)])
        final_orient.extend(['(111)/2','(111)','(111)*2'])
    lp = film_param*lp


    lattice_strain = (abs(sub_param - lp)/sub_param)*100
    lattice_strain2 = ((sub_param - lp)/sub_param)*100
    idx = np.where(lattice_strain == np.amin(lattice_strain))
    min_strain = np.amin(lattice_strain)
    orientation = final_orient[idx[0][0]]
    return [lp[idx[0][0]],min_strain, orientation, lattice_strain2[idx[0][0]]]

