Welcome to the documentation of AlgDiff
===================================

**AlgDiff** is a Python class that provides all necessary tools for the
design, analysis, and discretization of algebraic differentiators. An
interface to Matlab is also provided. This implementation was released
as part of the survey `[1] <#1>`__.

The toolbox is licensed under the BSD-3-Clause License, which is
suitable for both academic and industrial/commercial purposes.

This code has been created for research purposes at the `Chair of
Systems Theory and Control
Engineering <https://www.uni-saarland.de/en/chair/rudolph.html>`__ of
Saarland University, Germany. We apply algebraic differentiators to
solve different problems related to control theory and signal
processing: Parameter estimation, feedback control, fault detection and
fault tolerant control, model-free control … 


Motivation
==========

Estimating the derivatives of noisy signals is of paramount importance in many fields of engineering and applied mathematics. It is, however, a longstanding ill-posed and challenging problem, in the sense that a small error in measurement data can induce a significant error in the estimated derivatives.

Figure 1 shows the results of the numerical estimation of the first time derivative of a noisy signal based on an algebraic differentiator on the one hand and the simple difference quotient rule on the other. This simulation shows the excellent performance of this numerical differentiation approach.

.. figure:: motivationAlgDiff.png
   :scale: 50 %
   :alt: Numerical differentiation of a noisy signal

   Figure 1. Numerical differentiation of a noisy signal using a simple difference quotient on the one hand and an algebraic differentiator on the other.


Check out the :doc:`usage` section for further information, including
how to :ref:`installation` the project.


Contents
--------

.. toctree::
   
   onAlgebraicDiff
   gui
   usage
   api
   references
