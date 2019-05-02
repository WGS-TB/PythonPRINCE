from scipy import stats
import numpy as np
import csv

# -*- coding: utf-8 -*-
# Trying to fix testing bug ^ https://stackoverflow.com/questions/21639275/python-syntaxerror-non-ascii-character-xe2-in-file

def get_X_and_Y(data,template):
    X,Y = [],[]
    for genome,templateName,tnum,cn,score in data:
        if score == 0: continue
        if template == int(tnum):
            X.append(score)
            Y.append(cn)
    return(X,Y)

def get_equations(data, number_of_equations):
    equations = []
    for t in range(number_of_equations):
        X,Y = get_X_and_Y(data,t)
        X = np.array(X, dtype=np.float64)
        Y = np.array(Y, dtype=np.float64)
        slope, intercept, r_value, p_value, std_err = stats.linregress(X,Y)
        equations.append((slope,intercept))
    return(equations)

def get_copy_number(x,slope,intercept):
    return(round(slope*x+intercept,2))

def get_data(file):
    with open(file, 'r') as f:
        reader = csv.reader(f)
        data = list(reader)
        return(data)
