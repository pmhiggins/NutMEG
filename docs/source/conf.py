import os
import sys
sys.path.insert(0, os.path.abspath(".."))


# -- Project information -----------------------------------------------------

project = 'NutMEG'
copyright = '2026, P. M. Higgins'
author = 'P. M. Higgins'
# html_logo = 'NutMEG_logo_2026.png'

# The full version, including alpha/beta/rc tags
release = '2.0'

master_doc = 'index'



extensions = [
    "sphinx.ext.napoleon",
    "sphinx.ext.intersphinx",
    "sphinx.ext.viewcode",
    'sphinx_rtd_theme',
    "autoapi.extension",
    "sphinx.ext.mathjax",
    # "sphinx_copybutton",
    # "myst_parser",
    # "sphinx_nested_apidoc",
]

napoleon_use_ivar = True

myst_enable_extensions = [
    "colon_fence",
    "deflist",
]

autoapi_type = "python"
autoapi_dirs = ["../../NutMEG"]  # adjust if src/ layout
# autoapi_dirs = ["../../NutMEG/core", "../../NutMEG/models"]  # adjust if src/ layout
autoapi_root = "api"
autoapi_keep_files = True
autoapi_generate_api_docs = True
autoapi_python_use_implicit_namespaces = True
autoapi_python_class_content = 'class'

autoapi_options = [
    # "members",
    "undoc-members",
    "show-inheritance",
    "show-module-summary",
]

# add to this files or directories to be ignored in the API
autoapi_ignore = [
    "*/_v1_to_incorp/*",
    "*/util/*",
]

html_theme = 'furo'#'sphinx_rtd_theme'
# html_theme_options = {
#     "navigation_depth": 6,
#     "collapse_navigation": False,
#     "titles_only": True, # supress sphinx's urge to add 'submodules' subheadings in the TOC
#     "logo_only":True,
# }
html_static_path = ["_static"]
html_theme_options = {
    'sidebar_hide_name':True,
    "light_logo": "NutMEG_logo_2026.png",
    "dark_logo": "NutMEG_logo_dark.png",
}



# functions below stop autoAPI from adding the attribute
# definitions in code to the documentation (ie, creating duplicates).
def skip_util_members(app, what, name, obj, skip, options):
    # Napoleon handles the 'Attributes' section in the docstring.
    # AutoAPI finds the actual assignments in the code as 'attribute' or 'data'.
    if what in ("attribute", "data"):
        return True  # Skip the code-discovered version
    return skip

def setup(app):
    app.connect("autoapi-skip-member", skip_util_members)
