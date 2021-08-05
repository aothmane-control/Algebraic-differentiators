# Algebraic differentiators
A Python class implementing all necessary tools for the design, analysis, and discretization of algebraic differentiators. An interface to Matlab is also provided.

The toolbox is licensed unter the BSD-3-Clause License, which is suitabe for both academic and industrial/commercial purposes.
# Motivation 
Estimating the derivatives of noisy signals is of paramount importance in many
fields of engineering and applied mathematics. It is, however, a longstanding ill-posed
and challenging problem, in the sense that a small error in measurement data can
induce a significant error in the estimated derivatives.

Algebraic differentiators have been derived and discussed in the systems and control theory community. The initial works based on differential-algebraic methods have been developped by M. Fliess, M. Mboup, C. Join, and H. Sira-Ramirez.  

The following figure shows the results of the numerical estimation of the first time derivative of a noisy signal with an algebraic differentiator and the simple quotient rule.
![Motivation example](https://github.com/aothmane-control/Algebraic-differentiators-mirror/blob/master/data/motivationAlgDiff.png)
