import pandas as pd
import numpy as np
from statsmodels.tsa.stattools import adfuller
from statsmodels.tsa.vector_ar.var_model import VAR
import matplotlib.pyplot as plt

data = pd.read_csv(r"C:\Users\Kohsu\OneDrive\Desktop\graduation_thesis\data\data.csv", header=[0])
data = data.rename(columns = {'Unnamed: 0': 'date'})
print(data)

# log-transform
df = pd.DataFrame({
    "b": np.log(data["b"].astype(float)),
    "k": np.log(data["k"].astype(float)),
    "q": np.log(data["q"].astype(float)),
    "Y": np.log(data["Y"].astype(float)),
    "i": np.log(data["i"].astype(float))
})
df = df.dropna()
dates = pd.period_range(
    start="1960Q1",
    periods=len(df),
    freq="2Q"
)
df.index = dates

# get the difference
df_diff = pd.DataFrame({
    "K": df["k"].diff(),
    "b": df["b"].diff(),
    "Y": df["Y"].diff(),
    "q": df["q"].diff(),
    "i": df["i"].diff()
})
df_diff = df_diff.dropna()
df = df["1980Q1":"2007Q2"]

# VAR analysis
model = VAR(df_diff)
maxlags = 5
lag = model.select_order(maxlags).selected_orders
print("optimal lag：",lag['aic'],'\n')

results = model.fit(lag['aic'])
period = 12
irf = results.irf(period)
irf.plot(orth=True)
plt.title("Impulse Response")
plt.show()


