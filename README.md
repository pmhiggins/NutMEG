# NutMEG

![Project Logo](https://nutmeg-astrobiology.readthedocs.io/en/nm_v2/_images/NutMEG_logo_2026.png)

This is the beta release of NutMEG v2 (v2.0.0b1), which is more flexible for biokinetics than v1 but may still show some unexpected bugs.

## Nutrients, Maintenance, Energy &amp; Growth.  A flexible model for predicting habitability and productivity.

NutMEG is designed for problems in biogeochemistry and astrobiology where we know very little about the biology present---or if it even is present. If you want to get a taste of what NutMEG can do, check out the [documentation](https://nutmeg-astrobiology.readthedocs.io)!

## Installation
The NutMEG v2 beta can be installed via `pip`. However, this should be performed in a `conda` (or mamba) environment due to NutMEG's dependence on [`reaktoro`](https://reaktoro.org/). The minimal setup can be achieved in a terminal. Below we create an environment called `nm`:

```
conda config --add channels conda-forge
conda create --name nm reaktoro
conda activate nm
pip install NutMEG-astrobiology
```

From there, you can `import NutMEG` or its submodules as needed from within the
`nm` virtual environment. A complete setup guide is available in NutMEG's [documentation](https://nutmeg-astrobiology.readthedocs.io/en/nm_v2/guides/installation.html).

Tested and working on Python 3.12.x, with no specific hardware requirements.
