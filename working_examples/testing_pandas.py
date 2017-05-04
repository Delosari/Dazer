import pandas as pd

table_address   = 'example_table.txt'

#The two dimensional table is called dataframe in pandas. In contrast each column and row would be a series
#The parameters are:
#table_address (the location of your file)
#header, the number of the column you want to be the table column names (must be unique names)
#index_col, the number of the column you want the name of the rows (must be unique) (if you do not give one it will be 0,1,2,3,4,5...)
#delim_whitespace, true if the columns are separated by whitespace
df = pd.read_csv(table_address, header=0, delim_whitespace = True)

#If you print df you will get something like this:
print df

'''
   indexR  Col1  Col2  Col3 Col4
0       1     4    10   100  AAA
1       1     5    20    50  BBB
2       2     6    30   -30  AAA
3       3     7    40   -50  CCC
'''

#If you want a column you can just use the column name (the .values gives you the output in numpy array):
print 'One column', df['Col1'].values
print 'One column', df.Col1.values

#if you want a row you use .loc with the row name (the .values gives you the output in numpy array):
print 'One row', df.loc[0].values

#If you want a given cell you have several approaches
print 'One cell', df['Col1'][0]
print 'One cell', df.loc[0, 'Col1']
print 'One cell', df.loc[0].Col1

#In order to select rows which meet several Conditions
idx = (df.indexR == 1)
print  df.loc[idx].values

idx = (df.indexR == 1) & (df['Col2'] > 10)
print df.loc[idx].values

idx = (df.indexR == 1) | (df.indexR == 2)
print df.loc[idx].values

indexR_we_want = [1,2]
idx = (df.indexR.isin([1,2]))
print  df.loc[idx].values

#If we just one one column which meets the conditions condition
indexR_we_want = [1,2]
idx = (df.indexR.isin([1,2]))
print  df.loc[idx, 'Col1'].values
print  df.loc[idx].Col1.values







