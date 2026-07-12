# NutMEG

![Project Screenshot](docs/source/_static/NutMEG_logo_dark.png)

### Nutrients, Maintenance, Energy &amp; Growth.  A flexible model for predicting habitability and productivity.

NutMEG is designed for problems in biogeochemistry and astrobiology, where we know very little about the biology present---or if it even is present. If you want to get a taste of what NutMEG can do, check out the [documentation](https://nutmeg-astrobiology.readthedocs.io)!

This branch is an early build of NutMEG v2, which is more flexible for biokinetics than v1 but may still show some unexpected bugs.

While v2 is in development, after cloning this repository, use the NutMEG folder as if a local python package e.g. `import NutMEG`. A complete setup guide is available in NutMEG's [documentation](https://nutmeg-astrobiology.readthedocs.io/en/nm_v2/guides/installation.html).

Tested and working on Python 3.7.x, with no specific hardware requirements. Dependencies: numpy, pandas, reaktoro, uncertainties. Due to the reaktoro dependence, we strongly advise to use NutMEG from within a conda envrionment to enhance ease of use. Advice on getting this working is [here](https://nutmeg-astrobiology.readthedocs.io/en/latest/guides/installation.html).


More examples for using NutMEG can be found in the [NutMEG-Implementations respository](http://github.com/pmhiggins/NutMEG-Implementations), but please note these implementations used NutMEG v1.
