from distutils.core import setup

from pathlib import Path
this_directory = Path(__file__).parent
long_description = (this_directory / "docs/source/usage.rst").read_text()

setup(
  name = 'AlgDiff',
  packages = ['AlgDiff'],
  version = '2.4',
  license='bsd-3-clause',
  description = 'AlgDiff is a Python class implementing all necessary tools for the design, analysis, and discretization of algebraic differentiators. An interface to Matlab is also provided.',
  
  long_description = long_description,
  long_description_content_type='text/x-rst',
  
  author = 'Amine Othmane',
  author_email = 'amine.othmane@uni-saarland.de',
  url = 'https://github.com/aothmane-control/Algebraic-differentiators',
  download_url = 'https://github.com/aothmane-control/Algebraic-differentiators/releases/tag/v2.4',
  keywords = ['numerical-differentiation ', 'fir-filters', 'orthogonal-polynomials', 'numerical-methods '],
  install_requires=[
          'scipy',
          'mpmath',
          'numpy',
          'sphinx_rtd_theme'
      ],
  classifiers=[
    'License :: OSI Approved :: BSD License',
    'Programming Language :: Python :: 3',
  ],
)
