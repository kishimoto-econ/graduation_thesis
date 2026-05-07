import pandas as pd
import numpy as np
from statsmodels.tsa.stattools import adfuller
from statsmodels.tsa.vector_ar.var_model import VAR
import matplotlib.pyplot as plt

data = pd.read_csv(r"add_path", header=[0])
data = data.rename(columns = {'Unnamed: 0': 'date'})
print(data)

# log-transform
df = pd.DataFrame({
    "year": data["year"],
    "b": np.log(data["b"].astype(float)),
    "K": np.log(data["k"].astype(float)),
    "q": np.log(data["q"].astype(float)),
    "Y": np.log(data["Y"].astype(float)),
    "i": np.log(data["i"].astype(float))
})
df = df.set_index("year")


# get the difference
df_diff = pd.DataFrame({
    "K": data["k"].diff(),
    "b": data["b"].diff(),
    "Y": data["Y"].diff(),
    "q": data["q"].diff(),
    "i": data["i"].diff()
})
df_diff = df_diff.dropna()

# VAR analysis
model = VAR(df_diff)
maxlags = 10
lag = model.select_order(maxlags).selected_orders
print("optimal lag：",lag['aic'],'\n')

results = model.fit(lag['aic'])
period = 12
irf = results.irf(period)
irf.plot(orth=True)
plt.title("Impulse Response")
plt.show()


