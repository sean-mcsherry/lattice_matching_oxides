import sys
from pymatgen.ext.matproj import MPRester
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import math
import os
from math import gcd
from functools import reduce

def download_oxides_from_matproj(MAPI_KEY, elements_to_avoid, bandgap_min, file_name):

    mpr = MPRester(MAPI_KEY)
    data = mpr.query(criteria={"nelements": 3, "elements": {"$nin": elements_to_avoid,"$all": ["O"]}}, 
                     properties=['material_id','pretty_formula','spacegroup','theoretical','band_gap','elasticity','diel','unit_cell_formula', 'structure','decomposes_to'])
    data2 = mpr.query(criteria={"nelements": 2, "elements": {"$nin": elements_to_avoid,"$all": ["O"]}}, 
                     properties=['material_id','pretty_formula','spacegroup','theoretical','band_gap','elasticity','diel','unit_cell_formula','structure','decomposes_to'])


    # 1 --- all binary and ternary oxides (1)
    # 2 --- exclude impractical elements (3)

    d = {'film':[],'film ID':[], 'x-tal':[], 'ICSD?':[],'band gap (eV)':[],'n':[],'stoich':[],'x-tal type':[], 'stability':[]}
    for ii in range(0,len(data)):
      d['film'].append(data[ii]['pretty_formula'])
      d['x-tal'].append(data[ii]['spacegroup']['crystal_system'])
      d['film ID'].append(data[ii]['material_id'])
      if data[ii]['theoretical'] == True:
        d['ICSD?'].append('false')
      elif data[ii]['theoretical'] == False:
        d['ICSD?'].append('true')
      d['band gap (eV)'].append(data[ii]['band_gap'])
      if data[ii]['diel'] == None:
        d['n'].append(0)
      else:
        d['n'].append(round(data[ii]['diel']['n'],3))
      new_list = []
      lattice_stoich = [] 
      for jj in range(0, len(list(data[ii]['unit_cell_formula']))):
        elmt = list(data[ii]['unit_cell_formula'])[jj]
        new_list.append(elmt)
        lattice_stoich.append(int(data[ii]['unit_cell_formula'][elmt]))
      GCD = find_gcd(lattice_stoich)
      lattice_stoich_2 = np.array(lattice_stoich)/GCD
      d['stoich'].append(np.sort(lattice_stoich_2))
      if np.array_equal(np.sort(lattice_stoich_2),[1, 1, 3]) is True:
        d['x-tal type'].append('perovskite')
      elif np.array_equal(np.sort(lattice_stoich_2),[1, 2, 4]) is True:
        d['x-tal type'].append('spinel')
      elif np.array_equal(np.sort(lattice_stoich_2),[2, 2, 7]) is True:
        d['x-tal type'].append('pyrochlore')
      else:
        d['x-tal type'].append('unknown')
      if data[ii]['decomposes_to'] is None:
        d['stability'].append('stable')
      else:
        d['stability'].append('metastable')


    for ii in range(0,len(data2)):
      d['film'].append(data2[ii]['pretty_formula'])
      d['x-tal'].append(data2[ii]['spacegroup']['crystal_system'])
      d['film ID'].append(data2[ii]['material_id'])
      if data2[ii]['theoretical'] == True:
        d['ICSD?'].append('false')
      elif data2[ii]['theoretical'] == False:
        d['ICSD?'].append('true')
      d['band gap (eV)'].append(data2[ii]['band_gap'])
      if data2[ii]['diel'] == None:
        d['n'].append(0)
      else:
        d['n'].append(round(data2[ii]['diel']['n'],3))
      new_list = []
      lattice_stoich = [] 
      for jj in range(0, len(list(data2[ii]['unit_cell_formula']))):
        elmt = list(data2[ii]['unit_cell_formula'])[jj]
        new_list.append(elmt)
        lattice_stoich.append(int(data2[ii]['unit_cell_formula'][elmt]))
      GCD = find_gcd(lattice_stoich)
      lattice_stoich_2 = np.array(lattice_stoich)/GCD
      d['stoich'].append(np.sort(lattice_stoich_2))
      if np.array_equal(np.sort(lattice_stoich_2),[1, 1]) is True:
        d['x-tal type'].append('rocksalt')
      elif np.array_equal(np.sort(lattice_stoich_2),[1, 2]) is True:
        d['x-tal type'].append('fluorite')
      else:
        d['x-tal type'].append('unknown')
      if data2[ii]['decomposes_to'] is None:
        d['stability'].append('stable')
      else:
        d['stability'].append('metastable')

    df = pd.DataFrame(d)

    # 3 --- cubic crystal system
    df = df[(df['x-tal']=='cubic')]
    df=df.reset_index(drop=True)


    # 4 --- limit by bandgap
    df = df[(df['band gap (eV)'] == 0) | (df['band gap (eV)'] >= bandgap_min)]
    df=df.reset_index(drop=True)

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


def find_films(file_name, substrate, sub_xtal_type,lattice_parameter, orientation, strain):

    df = pd.read_pickle(file_name)
    df.drop('b', inplace=True, axis=1)
    df.drop('c', inplace=True, axis=1)
    df = df[~(df['x-tal type']==sub_xtal_type)]
    df=df.reset_index(drop=True)

    # 5 --- calculate all possible bonding types and minimum strain
    d = {'film orientation':[], 'film lattice parameter':[],'% strain':[]}
    for ii in range(0, len(df)):
      lattice_matching = find_xtal_match(sub_xtal_type, lattice_parameter,orientation, df['x-tal type'][ii], df['a'][ii])
      d['film orientation'].append(lattice_matching[2])
      d['film lattice parameter'].append(lattice_matching[0])
      d['% strain'].append(lattice_matching[1])
    df['film orientation'] = d['film orientation']
    df['film lattice parameter'] = d['film lattice parameter']
    df['% strain'] = d['% strain']
    df = df.loc[(df['% strain'] < strain)]
    df = df.reset_index(drop=True)
    df = df.sort_values(by=['stability','% strain','n','band gap (eV)'] , ascending=[False, True, False, False])
    df = df.reset_index(drop=True)

    df.drop('film ID', inplace=True, axis=1)
    df.drop('x-tal', inplace=True, axis=1)
    df.drop('stoich', inplace=True, axis=1)
    df = df[~(df['x-tal type']== 'unknown')]
    df = df.reset_index(drop=True)
    df.to_csv(r'sub_'+substrate+'_'+orientation+'.csv', index = False, header=True)
    pd.options.display.max_columns = 150
    pd.options.display.max_rows = 150
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
    idx = np.where(lattice_strain == np.amin(lattice_strain))
    min_strain = np.amin(lattice_strain)
    orientation = final_orient[idx[0][0]]
    return [lp[idx[0][0]],min_strain, orientation]

