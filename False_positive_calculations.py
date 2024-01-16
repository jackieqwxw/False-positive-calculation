#!/usr/bin/env python
# coding: utf-8

# # Using all targets of decoys to calculate the false positive

# In[1]:


import pandas as pd

csv_decoy_GB = 'decoys-MMGBSA.csv'
csv_active_GB = 'active-MMGBSA.csv'
df1 = pd.read_csv(csv_decoy_GB)
df2 = pd.read_csv(csv_active_GB)
#print(df1)
#print(df2)

csv_decoy_PB = 'decoys-MMPBSA.csv'
csv_active_PB = 'active-MMPBSA.csv'
df3 = pd.read_csv(csv_decoy_PB)
df4 = pd.read_csv(csv_active_PB)
#print(df3)
#print(df4)

csv_decoy_FEP = 'decoys-FEP.csv'
csv_active_FEP = 'active-FEP.csv'
df5 = pd.read_csv(csv_decoy_FEP)
df6 = pd.read_csv(csv_active_FEP)
#print(df5)
#print(df6)


# In[2]:


threshold = [i for i in range(-56,-10)]


# In[3]:


def false_positive_cal(df1,df2,method):

    decoy_decoy = []
    decoy_active = []
    active_decoy = []
    active_active = [] 

    for x in threshold:
        decoy0 = []
        active0 = []
        for i in df1:
            if i > x:
                decoy0.append(i)
            elif i <= x:
                active0.append(i)
        decoy_decoy.append(decoy0)
        decoy_active.append(active0)
    
        decoy1 = []
        active1 = []
        for i in df2:
            if i <= x:
                active1.append(i)
            elif i > x:
                decoy1.append(i)
        active_active.append(active1)
        active_decoy.append(decoy1)

    ad = [len(i) for i in active_decoy]
    aa = [len(i) for i in active_active]
    dd = [len(i) for i in decoy_decoy]
    da = [len(i) for i in decoy_active]

    #print('Confusion Table list')
    #print(f'//threshold {threshold} kcal/mol//')
    #print(f'Active_Active: {aa}\nActive_Decoy: {ad}\nDecoy_Decoy: {dd}\nDecoy_Active: {da}\n')

    #print(f"false positive for {method}:")
    for i in range(len(threshold)):
        false_positive = da[i] / (da[i] + aa[i])
        values = f'{false_positive*100:.2f}% for the cutoff {threshold[i]} kcal/mol'
        print(values)
    return 


# In[4]:


def false_positive_cal(df1,df2,method):

    decoy_decoy = []
    decoy_active = []
    active_decoy = []
    active_active = [] 

    for x in threshold:
        decoy0 = []
        active0 = []
        for i in df1:
            if i > x:
                decoy0.append(i)
            elif i <= x:
                active0.append(i)
        decoy_decoy.append(decoy0)
        decoy_active.append(active0)
    
        decoy1 = []
        active1 = []
        for i in df2:
            if i <= x:
                active1.append(i)
            elif i > x:
                decoy1.append(i)
        active_active.append(active1)
        active_decoy.append(decoy1)

    ad = [len(i) for i in active_decoy]
    aa = [len(i) for i in active_active]
    dd = [len(i) for i in decoy_decoy]
    da = [len(i) for i in decoy_active]

    #print('Confusion Table list')
    #print(f'//threshold {threshold} kcal/mol//')
    #print(f'Active_Active: {aa}\nActive_Decoy: {ad}\nDecoy_Decoy: {dd}\nDecoy_Active: {da}\n')

    #print(f"false positive for {method}:")
    false_positive = []
    cutoff = []
    for i in range(len(threshold)):
        false_positive.append(da[i] / (da[i] + aa[i])) 
        cutoff.append(threshold[i])
    values = list(zip(cutoff,false_positive))
    return values


# In[5]:


false_positive_cal(df1['Energy'],df2['Energy'],'MM/GBSA')


# In[6]:


false_positive_cal(df3['Energy'],df4['Energy'],'MM/PBSA')


# In[7]:


false_positive_cal(df5['Energy'],df6['Energy'],'FEP')


# # Using 7 out of total 9 targets of decoys to calculate the false positive

# In[8]:


from itertools import combinations
import numpy as np
import os

target_titles = ['ACE', 'ADRB1', 'FAK1','GRIK1','HMDH','MCR','PGH2','PRGR','TRYB1']

num_groups_to_extract = 7
combinations_of_groups = list(combinations(target_titles, num_groups_to_extract))


# In[9]:


def sort_target(df1,df2):
    a = []
    for i in df1:
        a.append(i)
    b = []
    for i in df2:
        b.append(i)
    data = list(zip(a,b))
    data.sort()
    return data


# In[10]:


def finding_minimum_cutoff(data):
    min_tuple = min(data, key=lambda x: x[1])
    corresponding_var = min_tuple[0]
    #minum = f"The minimum cutoff is: {corresponding_var}"   
    return corresponding_var


# In[11]:


def decoy_combination_std_dev(combinations_of_groups,data_method,df2,method):
    all_data = []
    for combination in combinations_of_groups:
        selected_groups = [x for x in data_method if x[0] in combination]
        all_data.append(selected_groups)

    cutoff_min = []
    fp_value = []
    for index,i in enumerate(all_data): # all_data has 36 groups data; each group has 137 or 140 data
        each_group = []
        for j in i:
            each_group.append(j[1])
        #print('\n')
        value = false_positive_cal(each_group,df2,method)
        fp_value.append(value)
        minum = finding_minimum_cutoff(value)
        cutoff_min.append(minum)

    second_elements_list = [[item[1] for item in sublist] for sublist in fp_value]
    transposed_data = np.array(second_elements_list).T

    np.savetxt(f'false_positive_{method}.txt', transposed_data, fmt='%f', delimiter='\t')
    
    average_value = np.mean(cutoff_min)
    std_dev = np.std(cutoff_min)
    output = f'The threshold values on false positive with standard deviation by {method}: {average_value:.2f} +/ {std_dev:.2f} kcal/mol'
    return output


# In[12]:


data_GBSA = sort_target(df1['Name'],df1['Energy'])
data_PBSA = sort_target(df3['Name'],df3['Energy'])
data_FEP = sort_target(df5['Name'],df5['Energy'])


# In[13]:


decoy_combination_std_dev(combinations_of_groups,data_GBSA,df2['Energy'],'MM-GBSA')


# In[14]:


decoy_combination_std_dev(combinations_of_groups,data_PBSA,df4['Energy'],'MM-PBSA')


# In[15]:


decoy_combination_std_dev(combinations_of_groups,data_FEP,df6['Energy'],'FEP')

