[![PyPI version](https://badge.fury.io/py/AlgDiff.svg)](https://badge.fury.io/py/AlgDiff)
# AlgDiff
AlgDiff: A Python class that provides all necessary tools for the design, analysis, and discretization of algebraic differentiators. An interface to Matlab is also provided.
This implementation was released as part of the survey [[1]](#1).  

The toolbox is licensed under the BSD-3-Clause License, which is suitable for both academic and industrial/commercial purposes.

This code has been created for research purposes at the [Chair of Systems Theory and Control Engineering](https://www.uni-saarland.de/en/chair/rudolph.html) of Saarland University, Germany.
 We apply algebraic differentiators to solve different problems related to control theory and signal processing: Parameter estimation, feedback control, fault detection and fault tolerant control, model-free control ...

Table of Contents
=================
* [Motivation](#motivation)
* [On algebraic differentiators](#on-algebraic-differentiators)
* [GUI](#gui)
* [Prerequisites for the implementation](#prerequisites-for-the-implementation)
* [Installation](#installation)
* [How to use  the implementation](#how-to-use--the-implementation)
* [Troubleshooting](#troubleshooting)
* [Documentation](#documentation)
* [Questions  &amp; Contact](#questions---contact)
* [License](#license)
* [References](#references)


# Motivation 
Estimating the derivatives of noisy signals is of paramount importance in many
fields of engineering and applied mathematics. It is, however, a longstanding ill-posed
and challenging problem, in the sense that a small error in measurement data can
induce a significant error in the estimated derivatives.



Figure 1 shows the results of the numerical estimation of the first time derivative of a noisy signal based on an algebraic differentiator on the one hand and the simple difference quotient rule on the other. This simulation shows the excellent performance of this numerical differentiation approach. 
| ![Motivation example](https://github.com/aothmane-control/Algebraic-differentiators/blob/master/data/motivationAlgDiff.png) |
|:--:| 
| Figure 1. Numerical differentiation of a noisy signal using a simple difference quotient on the one hand and an algebraic differentiator on the other |




# On algebraic differentiators
Algebraic differentiators have been derived and discussed in the systems and control theory community. The initial works based on differential-algebraic methods have been developed by Mboup,  Join, and Fliess in [[2]](#2). These numerical, non-asymptotic approximation approaches
for higher-order derivatives of noisy signals are well suited for real-time embedded systems. A historical overview and a detailed discussion of these differentiators and their time-domain and frequency-domain properties are given in the survey [[1]](#1).  

The approximation-theoretic derivation recalled in the survey [[1]](#1) permits the interpretation of the estimation process by the following three steps illustrated in the figure below stemming from [[1]](#1):

1. Projection: At time <img src="https://render.githubusercontent.com/render/math?math=t">, the sough <img src="https://render.githubusercontent.com/render/math?math=n">-th order time derivative <img src="https://render.githubusercontent.com/render/math?math=y^{(n)}"> over the interval <img src="https://render.githubusercontent.com/render/math?math=I_{T}(t)"> is projected onto the space of polynomials of degree <img src="https://render.githubusercontent.com/render/math?math=\mathrm{N}">. This yields the polynomial <img src="https://render.githubusercontent.com/render/math?math=p_\mathrm{N}"> depicted in the left and middle part of Figure 2.
2. Evaluation: The polynomial <img src="https://render.githubusercontent.com/render/math?math=p_\mathrm{N}"> is evaluated at <img src="https://render.githubusercontent.com/render/math?math=t-\delta_t">, which gives an estimate <img src="https://render.githubusercontent.com/render/math?math={\hat{y}^{(n)}(t)=p_{\N}(t-\delta_t)}"> for the derivative <img src="https://render.githubusercontent.com/render/math?math=y^{(n)}(t)"> as depicted in the central part of Figure 2. Choosing the delay to be the largest root of a special Jacobi polynomial increases the approximation order by 1 with a minimal delay. Alternatively, a delay-free estimation or even a prediction of the future derivative might be
    selected, at the cost of a reduced accuracy.
3. Repetition: The first two steps are repeated at each discrete time instant <img src="https://render.githubusercontent.com/render/math?math=t_i"> while keeping the parameters of the differentiator constant. This yields  the estimate <img src="https://render.githubusercontent.com/render/math?math=\hat{y}^{(n)}"> depicted in the right part of the Figure 2.

| ![filter_characteristics](https://github.com/aothmane-control/Algebraic-differentiators/blob/master/data/interpretationDifferentiators.png) |
|:--:| 
| Figure 2. Three-step process of the estimation of the  <img src="https://render.githubusercontent.com/render/math?math=n">-th order derivative <img src="https://render.githubusercontent.com/render/math?math={y^{(n)}:t\mapsto y^{(n)}(t)}"> of a signal <img src="https://render.githubusercontent.com/render/math?math=y:t\mapsto y(t)"> using  algebraic differentiators (figure from [[1]](#1)) |


Algebraic differentiators can be interpreted as linear time-invariant filters with a finite-duration impulse response. Figure 3 visualizes the online estimation process of the first derivative of a noisy signal. The filter window, the buffered signal, and the filter kernel can be clearly seen. 

|  <img src="https://github.com/aothmane-control/Algebraic-differentiators/blob/master/data/animationEstimation.gif" height="500">|
|:--:| 
| Figure 3. Visualization of the online estimation of the first derivative a noisy signal using an algebraic differentiator.|

These filters can be approximated as lowpass filters with a known cutoff frequency and a stopband slope. Figure 4 presents the amplitude and phase spectra of two exemplary filters. The lowpass approximation is also shown. 
| ![filter_characteristics](https://github.com/aothmane-control/Algebraic-differentiators/blob/master/data/filterSpectrum.png) |
|:--:| 
| Figure 4. Amplitude and phase spectra of two different filters and the corresponding lowpass approximation of the amplitude spectrum |

See [[1]](#1), [[3]](#3), [[4]](#4), and [[5]](#5) for more details on the parametrization of these differentiators.

# GUI
Since Version 1.1 a GUI is provided. Executable files for Linux and Windows operating systems are provided and do not require the installation of additional software. The binary files of the GUIs for different operating systems can be downloaded from the latest release [page](https://github.com/aothmane-control/Algebraic-differentiators/releases).
 Neither Python not Matlab have to be installed to start designing algebraic differentiators, get discrete filter coefficients, and estimate derivatives. This GUI can be used to plot relevant data (impulse and step responses, amplitude and phase spectra, estimated derivatives, ...), display relevant properties of the differentiators (estimation delay, cutoff frequency, window length, discretization effects, ...), and load measured signals for the estimation of their derivatives without a single line of code. Relevant properties, signals, spectra, and discrete filter coefficients can be exported for further processing. For testing the import of measurement data, a [file](https://github.com/aothmane-control/Algebraic-differentiators/blob/master/examples/QuickStart.ipynb) has been provided in the folder DataForGUI.

| ![GUI](https://github.com/aothmane-control/Algebraic-differentiators/blob/master/data/figureGUI.png) |
|:--:| 
| Figure 5. GUI for the interactive design, analysis, and use of algebraic differentiators |

# Prerequisites for the implementation
The code is implemented in Python 3. To use all functionalities the required packages are given in the requirements.txt file. The examples implemented in Python are written in [jupyter notebooks](https://jupyter.org/) and require the packages [jupyter_latex_envs](https://github.com/jfbercher/jupyter_latex_envs) for the generation of useful documentations and [matplotlib](https://matplotlib.org/) for the creation of plots. The functions in the toolbox can also be used in Matlab for which different examples are also included. Check the Matlab [documentation](https://de.mathworks.com/help/matlab/matlab_external/install-supported-python-implementation.html) for more details on the compatibility of your Matlab version with Python.

# Installation
To use AlgDiff using pip run the following in the command line: 

```
   $ pip install AlgDiff
```

# How to use  the implementation
The contribution of this implementation is an easy to use framework for the design and discretization of algebraic differentiators to achieve desired filter characteristics, i.e., to specify the cutoff frequency and the stopband slope. The file [algebraicDifferentiator.py](https://github.com/aothmane-control/Algebraic-differentiators/blob/master/algebraicDifferentiator.py) implements the class AlgebraicDifferentiator. This class contains all necessary functions for the design, analysis, and discretization of the differentiators.

Different examples are provided as jupyter notebooks and Matlab code in the following:
* A quick start in a jupyter [notebook](https://github.com/aothmane-control/Algebraic-differentiators/blob/master/examples/QuickStart.ipynb) available also as an [HTML file](https://htmlpreview.github.io/?https://github.com/aothmane-control/Algebraic-differentiators/blob/master/examples/QuickStart.html)
* A detailed jupyter [notebook](https://github.com/aothmane-control/Algebraic-differentiators/blob/master/examples/DetailedExamples.ipynb) available also as an [HTML file](https://htmlpreview.github.io/?https://github.com/aothmane-control/Algebraic-differentiators/blob/master/examples/DetailedExamples.html)
* The simultaneous elimination of a harmonic disturbance and approximation of derivatives is demonstrated in the jupyter [notebook](https://github.com/aothmane-control/Algebraic-differentiators/blob/master/examples/EliminationDisturbancesExample.ipynb) available also as an [HTML file](https://htmlpreview.github.io/?https://github.com/aothmane-control/Algebraic-differentiators/blob/master/examples/EliminationDisturbancesExample.html)
* A quick start in [Matlab](https://github.com/aothmane-control/Algebraic-differentiators/blob/master/examples/QuickStart.mlx)
* A [Matlab](https://github.com/aothmane-control/Algebraic-differentiators/blob/master/examples/DetailedExamples.mlx) code with several examples

# Troubleshooting
A list of known issues and fixes is also [provided](https://algebraic-differentiators.readthedocs.io/en/latest/troubleshooting.html).

# Documentation
 A detailed [documentation](https://algebraic-differentiators.readthedocs.io/en/latest/documentation.html) for all the functions is also provided in the webpage of the [project](https://algebraic-differentiators.readthedocs.io/en/latest/index.html).

# Questions  & Contact
Feel free to contact [Amine](https://www.uni-saarland.de/en/chair/rudolph/staff/aothmane.html) in case of suggestions or questions.

# License
BSD 3-Clause "New" or "Revised" License, see [License-file](https://github.com/aothmane-control/Algebraic-differentiators/blob/master/LICENSE).

# References
<a id="5">[1]</a> A. Othmane, L. Kiltz, and J. Rudolph, "Survey on algebraic numerical differentiation: historical developments, parametrization, examples, and applications", Int. J. Syst. Sci. https://doi.org/10.1080/00207721.2022.2025948

<a id="1">[2]</a> M. Mboup,  C. Join, and M. Fliess, "Numerical differentiation with annihilators in noisy environment", Numerical Algorithms, 50 (4), 439–467, 2009, https://doi.org/10.1007/s11075-008-9236-1


<a id="2">[3]</a> L. Kiltz and J. Rudolph, “Parametrization of algebraic numerical
differentiators to achieve desired filter characteristics,” in Proc. 52nd
IEEE Conf. on Decision and Control, Firenze, Italy, 2013, pp. 7010–
7015, https://doi.org/10.1109/CDC.2013.6761000

<a id="3">[4]</a> M. Mboup and S. Riachy, "Frequency-domain analysis and tuning of the algebraic differentiators," Int. J. Control , 91 (9), 2073–2081, 2018, https://doi.org/10.1080/00207179.2017.1421776 

<a id="4">[5]</a> A. Othmane, J. Rudolph, and H. Mounier, "Systematic comparison of numerical differentiators and an application to model-free control", Eur. J. Control. https://doi.org/10.1016/j.ejcon.2021.06.020

