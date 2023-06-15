Usage
=====


.. _installation:

Installation
------------



To use AlgDiff using pip run the following in the command line: 

.. code-block:: console

   $ pip install AlgDiff


.. _examples:

Some examples
----------------

To initialize an algebraic differentiator with parameter :math:`T=0.1`, :math:`\alpha=\beta=4`, :math:`N=0` and a sampling period of the measured signal equal to :math:`ts=0.01` do the following:

>>> from AlgDiff import *
>>> g = AlgebraicDifferentiator(N=0,alpha=4.,beta=4,T=0.1, ts=0.01)

The cutoff frequency in rad/s of this differentiator can be computed by

>>> wc = g.get_cutoffFreq()


The first order derivative of a signal :math:`y` given as a numpy array and sampled by :math:`ts=0.01` can be achieved by

>>> dydt = g.estimateDer(1,y) # y is a numpy array

Further detailed examples are provided in the `jupyter notebooks <https://github.com/aothmane-control/Algebraic-differentiators/tree/master/examples>`_ given the repository.

