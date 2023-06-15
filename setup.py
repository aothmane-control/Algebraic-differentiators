from distutils.core import setup
setup(
  name = 'AlgDiff',
  packages = ['AlgDiff'],
  version = '2.0',
  license='bsd-3-clause',
  description = 'AlgDiff is a Python class implementing all necessary tools for the design, analysis, and discretization of algebraic differentiators. An interface to Matlab is also provided.',
  author = 'Amine Othmane',
  author_email = 'amine.othmane@uni-saarland.de',
  url = 'https://github.com/aothmane-control/Algebraic-differentiators',
  download_url = 'https://github.com/aothmane-control/Algebraic-differentiators/releases/tag/v1.1.3',
  keywords = ['numerical-differentiation ', 'fir-filters', 'orthogonal-polynomials', 'numerical-methods '],
  install_requires=[            # I get to this in a second
          'scipy',
          'mpmath',
          'numpy',
      ],
  classifiers=[
    'License :: OSI Approved :: BSD License',
    'Programming Language :: Python :: 3',
  ],
)
