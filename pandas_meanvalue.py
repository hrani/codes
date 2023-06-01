#From the file get max value for the groupby Name and divide each value for the specific name's value with max number
import pandas as pd
import io
##csv file

# Time,Value,Name,Type
# 0.0,6.952e-05,BDNF,Stimulus
# 3990.0,6.952e-05,BDNF,Stimulus
# 4020.0,0.0072999999999999,BDNF,Stimulus
# 4980.0,0.0072999999999999,BDNF,Stimulus
# 5010.0,0.0072999999999999,BDNF,Stimulus
# 5040.0,0.0072999999999999,BDNF,Stimulus
# 5070.0,0.0072999999999999,BDNF,Stimulus
# 7590.0,0.0072999999999999,BDNF,Stimulus
# 0.0,0.01,glu,Stimulus
# 3990.0,0.01,glu,Stimulus
# 4020.0,0.01,glu,Stimulus
# 4980.0,0.01,glu,Stimulus
# 5010.0,0.01,glu,Stimulus
# 5040.0,0.01,glu,Stimulus
# 5070.0,0.01,glu,Stimulus
# 7590.0,0.01,glu,Stimulus
#df = pd.read_csv("new.csv")
data = '''Time Value Name Type
0.0 6.952e-05 BDNF Stimulus
3990.0 6.952e-05 BDNF Stimulus
4020.0 0.0072999999999999 BDNF Stimulus
4980.0 0.0072999999999999 BDNF Stimulus
5010.0 0.0072999999999999 BDNF Stimulus
5040.0 0.0072999999999999 BDNF Stimulus
5070.0 0.0072999999999999 BDNF Stimulus
7590.0 0.0072999999999999 BDNF Stimulus
0.0 0.01 glu Stimulus
3990.0 0.01 glu Stimulus
4020.0 0.01 glu Stimulus
4980.0 0.01 glu Stimulus
5010.0 0.01 glu Stimulus
5040.0 0.01 glu Stimulus
5070.0 0.01 glu Stimulus
7590.0 0.01 glu Stimulus
'''
df = pd.read_csv(io.StringIO(data), delim_whitespace=True)
newvalue = pd.DataFrame(df.groupby('Name')['Value'].apply(lambda x: (x/ x.max())))
result = pd.merge(df["Time"], newvalue, left_index=True, right_index=True)
result1 = pd.merge(result,df["Name"], left_index=True, right_index=True)
result2 = pd.merge(result1, df["Type"], left_index=True, right_index=True)
resultF = pd.DataFrame(result2.iloc[0: len(result2.index)])
frame = pd.DataFrame(resultF, columns=["Time", "Value", "Name", "Type"])
frame.to_csv("/tmp/out.csv",index=False)
