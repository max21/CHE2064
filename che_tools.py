import re
import requests
import string
from collections import deque
import numpy as np

def import_ChemCad(comps, base_url=None, extract_single_props=None,
                   extract_coeff_props=None, suffix='Props.txt'):

  N_comps = len(comps)
  if base_url is None:
    base_url = 'https://raw.githubusercontent.com/profteachkids/CHE2064/master/data/'

  if extract_single_props is None:
    extract_single_props = {'Molecular Weight' : 'Mw',
                            'Critical Temperature' : 'Tc',
                            'Critical Pressure' : 'Pc',
                            'Critical Volume' : 'Vc',
                            'Acentric factor' : 'w',
                            'Normal boiling point' : 'Tb',
                            'Heat of vaporization' : 'Hvap'}
  if extract_coeff_props is None:
    extract_coeff_props={'Vapor Pressure' : 'Pvap'}

  single_props_pat = re.compile('^\s+([\w\s]+?)\s+:\s+([-.0-9e+]+)\s+[\w\s/]*$', re.MULTILINE)
  coeffs_name_pat = re.compile("([\w ]+)\s[^\n]*?Equation.*?Coeffs:([- e\d.+]+)+?", re.DOTALL)
  coeffs_pat = re.compile('([-\de.+]+)')

  props_deque=deque()
  for comp in comps:
    text = requests.get(base_url+comp + suffix).text
    single_props = dict(single_props_pat.findall(text))
    props={'Name': comp}
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
    values = np.array([comp[prop] for comp in props_deque])
    props[prop]=values
  return props