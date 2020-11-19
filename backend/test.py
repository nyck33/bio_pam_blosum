#!/usr/bin/env python3
import numpy as np
import pandas as pd
 

def main():
    #matrix= np.genfromtxt('./sub_matrix/PAM30', dtype=None, skip_header=True, encoding=None, usecols=np.arange(1,25))
    df = pd.read_csv('./sub_matrix/PAM30', comment="#", header=1, index_col=0)
    print(df['A'].loc['U'])


if __name__ == "__main__":
    main()