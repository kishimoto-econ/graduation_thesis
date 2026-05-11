import pandas as pd
import numpy as np
from statsmodels.tsa.stattools import adfuller 

data = pd.read_csv(
    r"C:\Users\Kohsu\OneDrive\Desktop\graduation_thesis\data\data.csv",
    header=[0]
)
# print(data)

# log-transform
df = pd.DataFrame({
    "b": np.log(data["b"]).astype(float),
    "k": np.log(data["k"]).astype(float),
    "q": np.log(data["q"]).astype(float),
    "Y": np.log(data["Y"]).astype(float)
})
df = df.dropna()

dates = pd.period_range(
    start="1960Q1",
    periods=len(df),
    freq="Q"
)

df.index = dates
df = df["1960Q1":"2004Q1"]

# ADF test
variables_type = ['b', 'k', 'q', 'Y']
reg_types = ['n', 'c', 'ct']

for varia in variables_type:
    print(f"\n===== {varia} =====")
    for reg in reg_types:
        result = adfuller(df[varia], regression=reg)
        print(f"regression = '{reg}'")
        print(f"ADF statistic : {result[0]:.4f}")
        print(f"p-value       : {result[1]:.4f}")
        print("-" * 30)

# get the difference
df_diff = pd.DataFrame({
    "b_diff": df["b"].diff(),
    "k_diff": df["k"].diff(),
    "q_diff": df["q"].diff(),
    "Y_diff": df["Y"].diff()
})
df_diff = df_diff.dropna()
# print(df_diff)

variables_type_diff = ['b_diff', 'k_diff', 'q_diff', 'Y_diff']
for varia in variables_type_diff:
    print(f"\n===== {varia} =====")
    for reg in reg_types:
        result = adfuller(df_diff[varia], regression=reg)
        print(f"regression = '{reg}'")
        print(f"ADF statistic : {result[0]:.4f}")
        print(f"p-value       : {result[1]:.4f}")
        print("-" * 30)

