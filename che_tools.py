
import re
import requests
import string
from collections import deque
import numpy as np
from itertools import combinations

extract_single_props = {'Molecular Weight' : 'Mw',
                        'Critical Temperature' : 'Tc',
                        'Critical Pressure' : 'Pc',
                        'Critical Volume' : 'Vc',
                        'Acentric factor' : 'w',
                        'Normal boiling point' : 'Tb',
                        'Heat of vaporization' : 'Hvap'}

extract_coeff_props={'Vapor Pressure' : 'Pvap',
                     'Ideal Gas Heat Capacity' : 'CpIG',
                     'Liquid Heat Capacity' : 'CpL',
                     'Solid Heat Capacity' : 'CpS',
                     'Heat of Vaporization' : 'Hvap'}
base_url = 'https://raw.githubusercontent.com/profteachkids/CHE2064/master/data/'

def import_ChemCad(comps, base_url=base_url, extract_single_props=extract_single_props,
                   extract_coeff_props=extract_coeff_props, suffix='Props.txt'):

    N_comps = len(comps)

    id_pat = re.compile(r'ID\s+(\d+)')
    formula_pat = re.compile(r'Formula:\s+([A-Z1-9]+)')
    single_props_pat = re.compile('^\s+([\w\s]+?)\s+:\s+([-.0-9e+]+)\s+[\w\s/]*$', re.MULTILINE)
    coeffs_name_pat = re.compile("([\w ]+)\s[^\n]*?Equation.*?Coeffs:([- e\d.+]+)+?", re.DOTALL)
    coeffs_pat = re.compile('([-\de.+]+)')

    props_deque=deque()
    for comp in comps:
        text = requests.get(base_url+comp + suffix).text
        props={'Name': comp}
        props['ID']=id_pat.search(text).groups(1)[0]
        props['Formula']=formula_pat.search(text).groups(1)[0]
        single_props = dict(single_props_pat.findall(text))
        for k,v in extract_single_props.items():
            props[v]=float(single_props.pop(k))

        coeffs_name_strings = dict(coeffs_name_pat.findall(text))
        for k,v in extract_coeff_props.items():
            coeffs = coeffs_pat.findall(coeffs_name_strings[k])
            for letter, value in zip(string.ascii_uppercase,coeffs):
                props[v+letter]=float(value)
        props_deque.append(props)
    props={}
    for prop in props_deque[0].keys():
        if N_comps>1:
            values = np.array([comp[prop] for comp in props_deque])
        else:
            values = props_deque[0][prop]
        props[prop]=values
    return props

data_file = 'https://raw.githubusercontent.com/profteachkids/CHE2064/master/data/BinaryNRTL.txt'
def import_NRTL(ids, data_file=data_file):

    N_comps = len(ids)
    text = requests.get(data_file).text

    comps_string = '|'.join(ids)
    id_name_pat = re.compile(r'^\s+(\d+)[ ]+(' + comps_string +')[ ]+[A-Za-z]',re.MULTILINE)
    id_str = id_name_pat.findall(text)
    #maintain order of components
    id_dict = {v:k for k,v in id_str}
    id_str = [id_dict[id] for id in ids]
    comb_strs = combinations(id_str,2)
    comb_indices = combinations(range(N_comps),2)
    NRTL_A, NRTL_B, NRTL_C, NRTL_D, NRTL_alpha = np.zeros((5, N_comps,N_comps))
    start=re.search(r'Dij\s+Dji',text).span()[0]

    for comb_str, comb_index in zip(comb_strs, comb_indices):
        comb_str = '|'.join(comb_str)
        comb_values_pat = re.compile(r'^[ ]+(' + comb_str +
                                     r')[ ]+(?:' + comb_str + r')(.*)$', re.MULTILINE)

        first_id, values = comb_values_pat.search(text[start:]).groups(1)
        #if matched order is flipped, also flip indices
        if first_id != comb_index[0]:
            comb_index = (comb_index[1],comb_index[0])
        bij, bji, alpha, aij, aji, cij, cji, dij, dji  = [float(val) for val in values.split()]
        np.add.at(NRTL_B, comb_index, bij)
        np.add.at(NRTL_B, (comb_index[1],comb_index[0]), bji)
        np.add.at(NRTL_A, comb_index, aij)
        np.add.at(NRTL_A, (comb_index[1],comb_index[0]), aji)
        np.add.at(NRTL_C, comb_index, cij)
        np.add.at(NRTL_C, (comb_index[1],comb_index[0]), cji)
        np.add.at(NRTL_D, comb_index, dij)
        np.add.at(NRTL_D, (comb_index[1],comb_index[0]), dji)
        np.add.at(NRTL_alpha, comb_index, alpha)
        np.add.at(NRTL_alpha, (comb_index[1],comb_index[0]), alpha)
    return dict(A=NRTL_A, B=NRTL_B, C=NRTL_C, D=NRTL_D, alpha=NRTL_alpha)