
# shim module for handling ufloats in math throughout NutMEG.
try:
    from uncertainties import umath as math
except ImportError: # fallback to regular math if uncertainties is not installed.
    import math
