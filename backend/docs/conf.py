import os
import sys
sys.path.insert(0, os.path.abspath('..'))
project = 'MyFastAPIDocs'
author = 'Author'
extensions = ['sphinx.ext.autodoc', 'sphinx.ext.viewcode']
html_theme = 'alabaster'