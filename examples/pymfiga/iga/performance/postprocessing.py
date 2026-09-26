from matplotlib import pyplot as plt
from matplotlib import ticker as ticker
import seaborn as sns
import pandas as pd
import os

FILEPATH = os.path.dirname(os.path.realpath(__file__))
FOLDER2RESU = f"{FILEPATH}/results"
FOLDER2DATA = f"{FILEPATH}/data"
if not os.path.isdir(FOLDER2RESU):
    os.mkdir(FOLDER2RESU)
if not os.path.isdir(FOLDER2DATA):
    os.mkdir(FOLDER2DATA)
sns.set_style("whitegrid")


def set_axis(ax):
    ax.set_xlabel("Number of elements")
    ax.set_ylabel("CPU time (s)")
    ax.set_axis_on()
    ax.set(xscale="log")
    ax.set(yscale="log")
    ax.set_xlim((8, 64))
    ax.set_ylim(bottom=1e-3, top=1e1)
    ax.xaxis.set_minor_locator(ticker.NullLocator())
    ax.set_xticks([8, 16, 32, 64], ["8", "16", "32", "64"])


filename = f"{FOLDER2DATA}/cpu_time_elliptic_fastdiag"
df = pd.read_csv(f"{filename}.csv")

fig, ax = plt.subplots(figsize=(4, 4))
sns.lineplot(data=df, x="nbel", y="time", color="tab:blue", errorbar="sd", ax=ax)

set_axis(ax)
fig.tight_layout()
fig.savefig(f"{filename}.pdf")
plt.close()

################################################################################

filename = f"{FOLDER2DATA}/cpu_time_elliptic_matvec"
df = pd.read_csv(f"{filename}.csv")

fig, ax = plt.subplots(figsize=(4, 4))
sns.lineplot(
    data=df,
    x="nbel",
    y="time",
    hue="degree",
    ax=ax,
)

set_axis(ax)
fig.tight_layout()
fig.savefig(f"{filename}.pdf")
plt.close()

################################################################################

fig, ax = plt.subplots(figsize=(4, 4))
filename = f"{FOLDER2DATA}/cpu_time_sptm_fastdiag"

df = pd.read_csv(f"{filename}_legacy.csv")
sns.lineplot(
    data=df,
    x="nbel",
    y="time",
    color="tab:orange",
    errorbar="sd",
    ax=ax,
    label="Legacy",
)

df = pd.read_csv(f"{filename}_arrowhead.csv")
sns.lineplot(
    data=df,
    x="nbel",
    y="time",
    color="tab:blue",
    errorbar="sd",
    ax=ax,
    label="Arrowhead",
)

set_axis(ax)
fig.tight_layout()
fig.savefig(f"{filename}.pdf")
plt.close()

################################################################################

filename = f"{FOLDER2DATA}/cpu_time_sptm_matvec"
df = pd.read_csv(f"{filename}.csv")

fig, ax = plt.subplots(figsize=(4, 4))
sns.lineplot(
    data=df,
    x="nbel",
    y="time",
    hue="degree",
    ax=ax,
)

set_axis(ax)
fig.tight_layout()
fig.savefig(f"{filename}.pdf")
plt.close()


################################################################################

fig, ax = plt.subplots(figsize=(4, 4))
filename = f"{FOLDER2DATA}/cpu_time_sptm_fastdiag_3D"

df = pd.read_csv(f"{filename}.csv")
sns.lineplot(
    data=df,
    x="nbel",
    y="time",
    color="tab:blue",
    errorbar="sd",
    ax=ax,
)

set_axis(ax)
fig.tight_layout()
fig.savefig(f"{filename}.pdf")
plt.close()

################################################################################

filename = f"{FOLDER2DATA}/cpu_time_sptm_matvec_3D"
df = pd.read_csv(f"{filename}.csv")

fig, ax = plt.subplots(figsize=(4, 4))
sns.lineplot(
    data=df,
    x="nbel",
    y="time",
    hue="degree",
    ax=ax,
)

set_axis(ax)
fig.tight_layout()
fig.savefig(f"{filename}.pdf")
plt.close()
