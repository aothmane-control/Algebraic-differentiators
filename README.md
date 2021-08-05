# Algebraic differentiators
A Python class implementing all necessary tools for the design, analysis, and discretization of algebraic differentiators. An interface to Matlab is also provided.

The toolbox is licensed unter the BSD-3-Clause License, which is suitabe for both academic and industrial/commercial purposes.
# Motivation 
Estimating the derivatives of noisy signals is of paramount importance in many
fields of engineering and applied mathematics. It is, however, a longstanding ill-posed
and challenging problem, in the sense that a small error in measurement data can
induce a significant error in the estimated derivatives.

Algebraic differentiators have been derived and discussed in the systems and control theory community. The initial works based on differential-algebraic methods have been developped by M. Mboup, C. Join, and M. Fliess. These numerical non-asymptotic approximation approaches
for higher-order derivatives of noisy signals are suited for real-time embedded systems. 

The following figure shows the results of the numerical estimation of the first time derivative of a noisy signal with an algebraic differentiator and the simple quotient rule. This simulation shows the good estimation results using these approaches.
![Motivation example](https://github.com/aothmane-control/Algebraic-differentiators-mirror/blob/master/data/motivationAlgDiff.png)

This code has been created for the research at the [Chair of Systems Theory and Control Engineering](https://www.uni-saarland.de/lehrstuhl/rudolph.html) of Saarland University, Germany.
 We apply algebraic differentiators to solve different problems related to control theory an signal processing: Parameter estimation, feedback control, fault detection and fault tolerant control, modell-free control ...

# Prerequisites
The code runs in Python 3. To use all functionalities the following packages are required: scipy, numpy, mpmath, and math. The provided examples are written in jupyter notebooks which require the packages jupyter_latex_envs for the generation of usefull documentations and matplotlib for the creation of plots. The functions can also be used in Matlab. Different examples in Matlab are provided. Check the Matlab [documentation](https://de.mathworks.com/help/matlab/matlab_external/install-supported-python-implementation.html) for more details on the compatibility of your Matlab version with python.

# Structure
The file algebraicDifferentiator.py implements the class AlgebraicDifferentiator. This class contins all necessary functions for the design, analysis, and discretization of the differentiators.

# How to use
Algebraic differentiators are linear time-invariant filters with finite-duration impulse response. These filters can be approximated as lowpass filters with a computable cutoff frequency and stopband slope. The following figure presents the amplitude and phase spectra of two example filters. The lowpass approximation is also shown. 
![filter_characteristics](https://github.com/aothmane-control/Algebraic-differentiators-mirror/blob/master/data/filterSpectrum.png)

The contribution of this implementation is an easy to use framework for the design and discretization of algebraic differentiators to achieve desired filter characteristics, i.e., specify the cutoff frequency and the stopband slope. Different examples are provided as jupyter notebooks and Matlab code in the folder examples. In the folder documentation a sphinxs Makefile for the automatic generation of a documentation is given. The pre-compiled html documentation can be found in documentation/_build/html/index.html

# Questions  & Contact
Feel free to contact [Amine](https://www.uni-saarland.de/lehrstuhl/rudolph/personen/aothmane.html) in case of suggestions or questions.

# License
BSD 3-Clause "New" or "Revised" License, see [License-file](https://github.com/aothmane-control/Algebraic-differentiators/blob/master/LICENSE).
