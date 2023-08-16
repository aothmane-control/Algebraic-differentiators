"""
algebraicDifferentiator.py
=========================================
The implementation of the class AlgebraicDifferentiator
"""

from scipy import special
from mpmath import *
import numpy as np
import warnings
import math


class AlgebraicDifferentiator(object):
    """
    This class implements an algebraic differentiator and all the \
    methods necessary for its use, analysis, tuning, and discretization. The 
    differentiators which are LTI filters can be parametrized to achieve desired filter
    characteristics: Cutoff frequency and stopband slope. A cutoff frequency
    :math:`\omega_{\mathrm{c}}` or a filter window length :math:`T` can be
    specified. If the filter window length is specified the cutoff frequency
    is computed. For frequencies lower than the cutoff frequency the amplitude 
    of the Fourier transform is :math:`0` dB. For frequencies higher the
    amplitude falls with :math:`\\mu=\min\{\\alpha,\\beta\}+1` dB per decade.

    :param ts: Sampling period.
    :type ts: float
    :param alpha: Parameter :math:`\\alpha` of the weight function of the
        Jacobi polynomials. This parameter should satisfy :math:`\\alpha>n-1`,
        with :math:`n` the highest derivative to be estimated.
    :type alpha: float
    :param beta: Parameter :math:`\\beta` of the weight function of the
        Jacobi polynomials. This parameter should satisfy :math:`\\beta>n-1`,
        with :math:`n` the highest derivative to be estimated. The stopband
        slope is given by :math:`\\mu=\min\{\\alpha,\\beta\}+1`, i.e., frequencies 
        higher than the cutoff frequency are attenuated by :math:`20\\mu` dB
        per decade. Frequencies lower than the cutoff frequency are attenuated
        by :math:`0\,` dB.
    :type beta: float
    :param N: Truncation order of the generalized Fourier series. A delay-free
        derivative approximation is only possible for :math:`N\geq 1`. The 
        differentiator is parametrized by default such that the non-zero delay is
        minimized for the choice :math:`N\geq 1`. See the
        method set_theta for more details.
    :type N: int
    :param T: Filter window length. Takes a positive values if the length has
        to be specified. The cutoff frequency is then computed automatically.
        It should take the value None if the cutoff frequency is
        specified.
    :type T: float
    :param wc: Cutoff frequency of the differentiator. Takes a positive values
        if the cutoff frequency has to be specified. The filter window length
        is then computed automatically. It should take the value None 
        if the filter window length is specified.
    :type wc: float
    :type T: float
    :param display: Boolean variable indicating if all the information
        characterizing the differentiator should be printed.
    :type display: bool
    :param corr: Boolean variable that indicates if errors in the DC component of the 
        approximated signal stemming from the discretization should be corrected.
    :type corr: bool

    """

    ## The constructor    
    def __init__(self, ts=0.01, alpha=1., beta=1., N=0, T=1., wc =
                 None,display=True,corr=True):
        # Sampling period
        self.__ts = ts
        # Parameter alpha of weight function
        self.__alpha = alpha
        # Parameter beta of weight function
        self.__beta = beta
        # Order of polynomial approximating the derivative
        self.__N = N
        # Correction of discretization error
        self.correction = corr
        # Compute delay parameter
        if self.__N==0:
            self.__theta = 0
            self.__thetaBool = False
        else:
            r,_ = special.roots_jacobi(self.__N+1,self.__alpha,self.__beta)
            self.__theta = np.max(r)
            self.__thetaBool = True

        # Continuous window length
        if T is not None:
            self.__T = T
            self.checkParameters()
            self.__wc = self.get_cutoffFreq()
            self.__setUp = 'T'
        elif wc is not None:
            self.__wc = wc
            self.computeTfromWc(self.__wc)
            self.__setUp = 'wc'
        else:
            raise Exception('Cutoff frequency wx or window length T have to be given.')

        # Create dict to save the discretized filter values for each der.
        self.__w = {}
        # Create dict for the estimation delay when filter is discretized
        self.__delayDisc = {}
        # Print all properties of differentiator
        if display:
            self.printParam()

        self.warningRelativeAttenuation = -10

    def checkParameters(self, der=None):
        """
            Check if window length is too small. Function raises an error if T<ts and gives a warning
            if aliasing effects are to be expected.
        """
        if self.__T < self.__ts:
            e = "Window length T of the filter is smaller than the given sampling period ts" + \
                    " Choose a larger window, i.e., increase the parameter T."
            raise Exception(e)
        if der is not None:
            # Check if parameters alpha and beta satisfy conditions
            if not (np.min((self.__alpha, self.__beta)) > der-1):
                print('OK')
                txt = "The parameters \u03B1 have  \u03B2 have to satisfy min(\u03B1,\u03B2)>n-1," \
                      " with n the order of the sought " \
                      "derivative. Increase min(\u03B1,\u03B2)."
                raise Exception(txt)
            # Check for aliasing effects
            wc = self.__wc
            wN = np.pi/self.__ts
            mu = np.min((self.__alpha, self.__beta))+1
            attNum = (wc/wN)**(mu-der)
            if attNum > self.warningRelativeAttenuation:
                txt = "Attention: \n\tAliasing effects are to be expected. Compare for example the amplitude" \
                      " spectra of the continuous-time and the discrete-time filters. You can use the implemented" \
                      " functions for these computations." \
                      "\nPossible solutions:\
                            \n\tincrease min(\u03B1,\u03B2)\
                            \n\tdecrease the cutoff frequency (if in cutoff frequency mode)\
                            \n\tincrease window length (if in window length mode)"
                warnings.warn(txt)

    def discretize(self,der,method="mid-point",redFilLength=False,\
                   redTol=0.01,discreteSpectrum=False):
        """
        This function performs the discretization of the filter for the 
        estimation of the derivative of order :math:`der`.
        Three discretization schemes are implemented: midpoint, trapezoidal
        and using the analytical integration. The mid-point rule uses one
        filter coefficient less than the trapezoidal rule. It also reduces
        the estimation delay by half a sampling period. The analytical
        integration rule is recommended for small filter window lengths,
        i.e., in general less than 20 filter coefficients. This discretization
        method can only be used for the first order or any higher derivative.
        The error stemming from
        the discretization is corrected using a correction factor, if the 
        differentiator has been initialized to do so.
        When the parameters :math:`\\alpha` and :math:`\\beta` are large,
        the filter kernel has very low values at the beginning and the at
        the end of the filter window. The length of the window can be
        reduced to save computation time and memory. This can be performed
        by enabling the redFilterLength parameter and giving a tolerance
        redTolerance. This tolerance automatically computes the size of
        the new window such that the truncation is performed symmetrically.
        If this is enabled, the functions also returns the times at which
        the filter was truncated.

        :param der: Order of the derivative to be estimated.
        :type der: int
        :param method: Discretization scheme: "mid-point", "trapezoidal",\
            "analytic", "simpson rule".
        :type method: string
        :param reduceFilLength: Reduce or not the filter window length.
        :type reduceFilLength: bool
        :param redTol: Tolerance to be used when the filter length is reduced.
        :type redTol: float.
        :param discreteSpectrum: If it is set, then the parameter \
            :math:`\\theta` used in the discretization is returned. See survey paper
            for more details.
        :type discreteSpectrum: bool
        
        :returns:
        	- coeff (dictionary) - Discretized filter in a dict where the keys are the derivative \
            		for which a filter has been discretized. Each element is a dict. with \
            		keys the used discretization methods. If the correction of the DC \
            		component has been enabled with the parameter corr of the class \
            		initialization, this output contains the corrected filter coefficients.
        	- tau_1 (:py:class:`float`) - If redFilLength is set, tau_1 is the time where the filter \
        		window is reduced \
            		before at the left side of the interval. The estimation delay is reduced by tau_1.
        	- tau_2 (:py:class:`float`) - If redFilLength is set, tau_2 is the time where the filter \
        	 	window is reduced\
            		before at the right side of the interval. This value does not affect the delay.
        	- theta (:py:class:`float`) - If discreteSpectrum is set then the parameter\
        		 :math:`\\theta` is also returned.
        """
        self.checkParameters(der)
        theta0 = 0
        theta = theta0
        L0 = int(self.__T/self.__ts)
        red = ""
        tau1 = 0
        
        def newton_cotes_rules(p,order,L):
            order -= 1
            out = np.zeros((L,))
            out[0] = p[0]
            for i in range(1,L-1):
                if i%order==0:
                    out[i] = p[0]+p[-1]
                else:
                    out[i] = p[i%order]
            out[-1] = p[-1]
            return out

                
        if redFilLength:
            tau1,tau2,_,_,_ = self.reduceFilterLength(der,tol=redTol)
            theta0 = tau1/self.__ts
            L0 = int((tau2-tau1)/self.__ts)+1
            red = "-red"
        if method=="mid-point":
            # Discretized window length
            self.__L = L0
            theta = 0.5+theta0
            self.__delayDisc[method+red] =\
            self.get_delay()-self.__ts/2-tau1
            k = (np.array(range(self.__L))+theta)*self.__ts
            if der==0:
                if der in self.__w.keys():
                    self.__w[der][method+red] = self.__ts*self.evalKernel(k)
                else:
                    self.__w[der]= {method+red: self.__ts*self.evalKernel(k)}
            else:
                if der in self.__w.keys():
                    self.__w[der][method+red] = self.__ts*self.evalKernelDer(k,der)
                else:
                    self.__w[der] =\
                    {method+red:self.__ts*self.evalKernelDer(k,der)}
        if method=="euler":
            # Discretized window length
            self.__L = L0+1
            theta = theta0
            self.__delayDisc[method+red] =\
            self.get_delay()-self.__ts/2-tau1
            k = (np.array(range(self.__L))+theta)*self.__ts
            if der==0:
                if der in self.__w.keys():
                    self.__w[der][method+red] = self.__ts*self.evalKernel(k)
                else:
                    self.__w[der]= {method+red: self.__ts*self.evalKernel(k)}
            else:
                if der in self.__w.keys():
                    self.__w[der][method+red] = self.__ts*self.evalKernelDer(k,der)
                else:
                    self.__w[der] =\
                    {method+red:self.__ts*self.evalKernelDer(k,der)}
        elif method=="trapezoidal":
            self.__L = 1+L0
            theta = theta0
            self.__delayDisc[method+red] = self.get_delay()-tau1*self.__ts
            k = (np.array(range(self.__L))+theta)*self.__ts
            if der==0:
                if der in self.__w.keys():
                    self.__w[der][method+red] = self.__ts*self.evalKernel(k)
                else:
                    self.__w[der]= {method+red: self.__ts*self.evalKernel(k)}
            else:
                if der in self.__w.keys():
                    self.__w[der][method+red] = self.__ts*self.evalKernelDer(k,der)
                else:
                    self.__w[der] =\
                    {method+red:self.__ts*self.evalKernelDer(k,der)}
            self.__w[der][method+red][0] = self.__w[der][method+red][0]/2
            self.__w[der][method+red][-1] = self.__w[der][method+red][-1]/2
        elif method=="simpson rule":
            order = 3
            self.__L = L0
            if self.__L%(order)!=0:
                self.__L += order-self.__L%order
            theta = theta0
            self.__delayDisc[method+red] = self.get_delay()-tau1*self.__ts
            k = np.linspace(0,self.__T,self.__L)
            p = [1/3,4/3,1/3]
            weight = newton_cotes_rules(p,order,self.__L)
            if der==0:
                if der in self.__w.keys():
                    self.__w[der][method+red] = weight*self.__ts*self.evalKernel(k)
                else:
                    self.__w[der]= {method+red: weight*self.__ts*self.evalKernel(k)}
            else:
                if der in self.__w.keys():
                    self.__w[der][method+red] = weight*self.__ts*self.evalKernelDer(k,der)
                else:
                    self.__w[der] =\
                    {method+red:weight*self.__ts*self.evalKernelDer(k,der)}
        elif method=="simpson 3/8 rule":
            order = 4
            self.__L = L0
            if self.__L%(order)!=0:
                self.__L += order-self.__L%order
            theta = theta0
            self.__delayDisc[method+red] = self.get_delay()-tau1*self.__ts
            k = np.linspace(0,self.__T,self.__L)
            p = np.array([1/8,3/8,3/8,1/8])*3
            weight = newton_cotes_rules(p,order,self.__L)
            if der==0:
                if der in self.__w.keys():
                    self.__w[der][method+red] = weight*self.__ts*self.evalKernel(k)
                else:
                    self.__w[der]= {method+red: weight*self.__ts*self.evalKernel(k)}
            else:
                if der in self.__w.keys():
                    self.__w[der][method+red] = weight*self.__ts*self.evalKernelDer(k,der)
                else:
                    self.__w[der] =\
                    {method+red:weight*self.__ts*self.evalKernelDer(k,der)}
        elif method=="boole rule":
            order = 5
            self.__L = L0
            if self.__L%(order)!=0:
                self.__L += order-self.__L%order
            theta = theta0
            self.__delayDisc[method+red] = self.get_delay()-tau1*self.__ts
            k = np.linspace(0,self.__T,self.__L)
            p = 1/90*np.array([7,32,12,32,7])*(order-1)
            weight = newton_cotes_rules(p,order,self.__L)
            if der==0:
                if der in self.__w.keys():
                    self.__w[der][method+red] = weight*self.__ts*self.evalKernel(k)
                else:
                    self.__w[der]= {method+red: weight*self.__ts*self.evalKernel(k)}
            else:
                if der in self.__w.keys():
                    self.__w[der][method+red] = weight*self.__ts*self.evalKernelDer(k,der)
                else:
                    self.__w[der] =\
                    {method+red:weight*self.__ts*self.evalKernelDer(k,der)}
        elif method=="newton-cotes order 6":
            order = 6
            self.__L = L0
            if self.__L%(order)!=0:
                self.__L += order-self.__L%order
            theta = theta0
            self.__delayDisc[method+red] = self.get_delay()-tau1*self.__ts
            k = np.linspace(0,self.__T,self.__L)
            p = np.array([41/840, 9/35, 9/280, 34/105, 9/280, 9/35, 41/840])*(order-1)
            weight = newton_cotes_rules(p,order,self.__L)
            if der==0:
                if der in self.__w.keys():
                    self.__w[der][method+red] = weight*self.__ts*self.evalKernel(k)
                else:
                    self.__w[der]= {method+red: weight*self.__ts*self.evalKernel(k)}
            else:
                if der in self.__w.keys():
                    self.__w[der][method+red] = weight*self.__ts*self.evalKernelDer(k,der)
                else:
                    self.__w[der] =\
                    {method+red:weight*self.__ts*self.evalKernelDer(k,der)}
        elif method=="newton-cotes order 7":
            order = 7
            self.__L = L0
            if self.__L%(order)!=0:
                self.__L += order-self.__L%order
            theta = theta0
            self.__delayDisc[method+red] = self.get_delay()-tau1*self.__ts
            k = np.linspace(0,self.__T,self.__L)
            p = np.array([751/17280, 3577/17280, 49/640, 2989/17280,\
                          2989/17280, 49/640, 3577/17280, 751/17280])*(order-1)
            weight = newton_cotes_rules(p,order,self.__L)
            if der==0:
                if der in self.__w.keys():
                    self.__w[der][method+red] = weight*self.__ts*self.evalKernel(k)
                else:
                    self.__w[der]= {method+red: weight*self.__ts*self.evalKernel(k)}
            else:
                if der in self.__w.keys():
                    self.__w[der][method+red] = weight*self.__ts*self.evalKernelDer(k,der)
                else:
                    self.__w[der] =\
                    {method+red:weight*self.__ts*self.evalKernelDer(k,der)}
        elif method=="analytic":
            self.__L = L0
            self.__delayDisc[method+red] = self.get_delay()-self.__ts/2\
                                            -tau1*self.__ts
            if der==0:
                print("Error: Method not available for this derivative order")
            elif der>1:
                theta = 0.5+theta0
                k = (np.array(range(self.__L))+theta)*self.__ts
                p1 = self.evalKernelDer(k,der-1)
                p2 = self.evalKernelDer(k+self.__ts,der-1)
                if der in self.__w.keys():
                    self.__w[der][method+red] = p2-p1
                else:
                    self.__w[der] = {method+red:p2-p1}
            elif der==1:
                theta = 0.5+theta0
                k = (np.array(range(self.__L))+theta)*self.__ts
                p1 = self.evalKernel(k)
                p2 = self.evalKernel(k+self.__ts)
                if der in self.__w.keys():
                    self.__w[der][method+red] = p2-p1
                else:
                    self.__w[der] = {method+red:p2-p1}

        # Correction discretization error
        if self.correction:
            g = 0
            for ik in range(len(self.__w[der][method+red])):
                g += self.__w[der][method+red][ik]*(-ik)**(der)
            g = g*self.__ts**der/math.factorial(der)
            self.__w[der][method+red] = self.__w[der][method+red]/g
        if redFilLength:
            if discreteSpectrum:
                return self.__w, tau1, tau2, theta
            else:
                return self.__w, tau1, tau2
        else:
            if discreteSpectrum:
                return self.__w,theta
            else:
                return self.__w

    def estimateDer(self,k,x,method="mid-point",conv='same',\
                    redFilLength=False,redTol=0.01):
        """
        This function estimates the derivative of order :math:`k\geq0`\
        using the discretization method given. If the filter was not\
        discretized before the function performs the discretization.\
        See the method discretize for more details on the discretization.

        :param k: Order of the derivative to be estimated.
        :type k: int
        :param x: Signal whose derivative has to be estimated.
        :type x: numpy array
        :param conv: Parameter to specify the type of convolution
                the numpy conv function does. See the numpy documentation
                for more details.
        :param method: Discretization method: "mid-point",
                "trapezoidal", "analytic"

        :type method: string
        :type conv: string
        :param reduceFilterLength: Specify whether or not the filter window
                length should be reduced or not. See the discretize method for more
                details.
        :type reduceFilterLength: bool
        :param redTol: Tolerance to be used when the filter length is reduced
        :type redTol: float.
        :return: Estimated derivative in a numpy array with the same dimensions
            as the variable x
        """
        red = ""
        if redFilLength:
            red = "-red"

        # If filter was not discretized yet, discretize it
        if k not in self.__w:
            self.discretize(k,method,redFilLength=redFilLength,redTol=redTol)
        elif method+red not in self.__w[k]:
            self.discretize(k,method,redFilLength=redFilLength,redTol=redTol)

        # Estimate derivative of x
        
        N = int(len(self.__w[k][method+red]))-1
        if conv == 'valid':
            return np.hstack((np.zeros((N,)),\
                                       np.convolve(x,self.__w[k][method+red],'valid')))
        else:
            N = int(len(self.__w[k][method+red])/2)
            dx2 = np.convolve(x,self.__w[k][method+red],'same')
            return np.concatenate((np.zeros((N,)),dx2[:-N]))


    def get_alpha(self):
        """
        This function returns the parameter :math:`\\alpha`.

        :return: :math:`\\alpha` as a float.
        """

        return self.__alpha

    def get_beta(self):
        """
        This function returns the parameter :math:`\\beta`.
        
        :return: :math:`\\beta` as a float.
        """

        return self.__beta
    
    def get_N(self):
        """
        This function returns the parameter :math:`N`.

        :return: :math:`N` as a float.
        """

        return self.__N
    
    def get_theta(self):
        """
        This function returns the parameter :math:`\\vartheta`.

        :return: :math:`\\vartheta` as a float.
        """

        return self.__theta
    
    def get_ts(self):
        """
        This function returns the sampling period.

        :return: Sampling period as a float.
        """

        return self.__ts

    def get_T(self):
        """
        This function returns the filter window length of the algebraic differentiator.

        :return: Filter window length as a float.
        """

        return self.__T

    def get_degreeExactness(self,n):
        """
        This function returns the degree of exactness :math:`\gamma`\
                when the :math:`n`-th derivative is approximated.

        :return: Degree of exactness as an integer.
        """

        if (self.__N==0) or (self.__thetaBool):
            gamma = n+self.__N+1
        else:
            gamma = n+self.__N

        return gamma

    def get_cutoffFreq(self):
        """
        This function returns the cutoff frequency :math:`\omega_c`\
        of the algebraic differentiator.

        :return: Cutoff frequency as a float.
        """

        a,b,c = self.get_asymptotesAmpFilter(np.array([1.0]))

        return self.__wc

    def timeShift(self,t):
        """
        This function computes the time transformation :math:`\\theta(t)=1-2/T`, with\
        :math:`T` the window length of the filter. 

        :param t: Time instants where the transformation should be evaluated.
        :type t: numpy array
        :return: Evaluated transformation with the same type as the variable t.
        """

        return 1-2/self.__T*t
       
    def weightFcn(self,a,b,t):
        """
        This function evaluates the weight function of the algebraic differentiator at\
                the time instants :math:`t`. 

        :param a: Parameter :math:`\\alpha>-1` of the weight function
        :type a: float
        :param b: Parameter :math:`\\beta>-1` of the weight function
        :type b: float
        :param t: Time instants where weight function should be evaluated.
        :type t: numpy array
        :return: Evaluated weight function in a numpy array.
        """
        if a<0 or b<0:
            t[t>=1] = 1
            t[t<=-1] = -1
            w = np.multiply(np.power(1.-t,a),np.power(t+1., b))
            w[t>=1] = 0
            w[t<=-1] = 0
        else:
            t[t>1] = 1.1
            t[t<-1] = -1.1
            w = np.multiply(np.power(1.-t,a),np.power(t+1., b))
            w[t>1] = 0
            w[t<-1] = 0

        return w

    def get_stepResponse(self,t):
        """
        This function returns the step response of the differentiator \
        evaluated at the time instants in t.

        :param t: Time instants where the step response should be evaluated.
        :type t: numpy array
        :return: Evaluated step response in a numpy array.
        """
        response = np.zeros(t.shape)
        a = self.__alpha + 1
        b = self.__beta + 1
        for it in range(len(t)):
            if t[it]<=0:
                response[it] = 0
            elif t[it]>self.__T:
                response[it] = 1
            else:
                response[it] = special.betainc(a,b,t[it]/self.__T)
                for i in range(1,self.__N+1):
                    h = math.pow(2.,a+b+1)*special.gamma(i+a+1)*special.gamma(i+b+1)\
                    /(math.factorial(i)*(2*i+a+b+1)*special.gamma(i+a+b+1))
                    w = self.weightFcn(a,b,self.timeShift(np.array([t[it]])))
                    P = special.eval_jacobi(i-1,a,b,self.timeShift(np.array([t[it]])))
                    response[it]+=w/2*special.eval_jacobi(i,a,b,self.__theta)/i/h*P
        return response

    def evalKernel(self,t):
        """
        This function evaluates the kernel of the algebraic differentiator at\
        times t. It corresponds to the impulse response of the filter. 

        :param t: Time instants where the kernel should be evaluated.
        :type t: numpy array
        :return: Evaluated kernel in a numpy array with same dimensions as t.
        """

        a = self.__alpha
        b = self.__beta
        g = 0
        # Evaluate the filter kernel
        for i in range(self.__N+1):
            h = math.pow(2.,a+b+1)*special.gamma(i+a+1)*special.gamma(i+b+1)\
                    /(math.factorial(i)*(2*i+a+b+1)*special.gamma(i+a+b+1))
            w = self.weightFcn(a,b,self.timeShift(t))
            P = special.eval_jacobi(i,a,b,self.timeShift(t))
            g +=1/h*special.eval_jacobi(i,a,b,self.__theta)*np.multiply(w,P)

        return 2./self.__T*g
    
    def evalKernelDer(self,t,k):
        """
        This function evaluates the derivative of order :math:`k` of the\
        algebraic differentiator at times t. If :math:`k=0`\
        the kernel itself is evaluated.

        :param t: Time instants where the derivative should be evaluated.
        :type t: numpy array
        :param k: Order of the derivative to be evaluated.
        :type k: int
        :return: Evaluated derivative in a numpy array with same dimensions as
            t.
        """

        a = self.__alpha
        b = self.__beta
        dg = 0
        
        # Evaluate the filter kernel
        w = self.weightFcn(a-k,b-k,self.timeShift(t))
        for i in range(self.__N+1):
            h = pow(2,a+b+1)*special.gamma(i+a+1)*special.gamma(i+b+1)\
                    /(math.factorial(i)*(2*i+a+b+1)*special.gamma(i+a+b+1))
            P = special.eval_jacobi(i+k,a-k,b-k,self.timeShift(t))
            dg +=1/h*special.eval_jacobi(i,a,b,self.__theta)*P*special.gamma(i+1+k)\
                    /special.gamma(i+1)

        return 2**(2*k+1)/(self.__T**(k+1))*dg*w

    def get_ampAndPhaseFilter(self,omega):
        """
        This function evaluates the amplitude and the phase of the \
        Fourier transform of the algebraic differentiator\
        at the frequencies omega.

        :param omega: Frequencies where the amplitude and phase should be \
                        evaluated
        :type omega: numpy array
        
        :returns:
        	- amp (:py:class:`numpy array`) - Amplitude of the Fourier transform of the filter.
        	- phase (:py:class:`numpy array`) - Phase of the Fourier transform of the filter.
        """

        mp.dps = 150; mp.pretty = True
        a = self.__alpha
        b = self.__beta
        theta = self.__theta
        F = 0
        M = np.frompyfunc(hyp1f1, 3, 1)
        absVec = np.frompyfunc(fabs,1,1)
        argVec = np.frompyfunc(arg,1,1)
        for i in range(self.__N+1):
            bi = (a+b+2*i+1)*special.eval_jacobi(i,a,b,theta)/(a+b+i+1)
            tmp = 0
            for k in range(i+1):
                cik = (-1)**(i-k)*special.binom(i, k)
                tmp += cik*M(b+k+1,a+b+i+2,1j*omega*self.__T)
            F += bi*tmp
        F = F*np.exp(-1j*omega*self.__T)
        amp = absVec(F)
        phase = argVec(F)
        ampF = np.array(amp.tolist(),dtype=np.float64)
        phaseF = np.unwrap(np.array(phase.tolist(),dtype=np.float64))

        return ampF, phaseF

    def get_ampSpectrumDiscreteFilter(self,omega,n,method='mid-point'):
        """
        This function computes the amplitude spectrum of the discretized filter
        for the estimation of the :math:`n` -th derivative at the frequencies
        :math:`\\omega`. The filter is discretized using the given method. See
        the discretize method for more details.

        :param omega: Frequencies where the amplitude should be evaluated.
        :type omega: numpy array
    
        :param n: Derivative to be estimated.
        :type n: integer

        :param method: Method used for the discretization: "mid-point",
            "trapezoidal", "analytic".
        :type method: string 
        
        :returns:
        	- amplitude (:py:class:`numpy array`) - Amplitude of the Fourier transform of the discretized filter evaluated at omega. 
        	- phase (:py:class:`numpy array`) - Phase of the Fourier transform of the discretized filter evaluated at omega.
        """
        w,theta = self.discretize(n,method,discreteSpectrum=True)
        k = np.arange(0,len(w[n][method]))
        Gdis =  np.array([np.exp(-1j*omegai*(k+theta)*self.__ts).dot(w[n][method]) for omegai in omega])
        return np.abs(Gdis), np.unwrap(np.angle(Gdis))

    def get_asymptotesAmpFilter(self,omega):
        """
        This function evaluates the upper and lower bounds of the \
        amplitude of the Fourier transform of the algebraic differentiator\
        at the frequencies omega. Note that these bounds hold only for\
        frequencies greater than the cutoff frequency.

        The function also return the approximation as a low pass filter.

        :param omega: Frequencies where the bounds should be evaluated
        :type omega: numpy array

        :returns:
	        - ub - (:py:class:`numpy array`) - Upper bound for the amplitude spectrum evaluated at omega..
        	- lb - (:py:class:`numpy array`) -  Lower bound for the amplitude spectrum evaluated at omega..
        	- mb - (:py:class:`numpy array`) - Lowpass approximation of the filter evaluated at omega..
        """

        kappa = np.abs(self.__beta-self.__alpha)
        mu = 1+np.minimum(self.__beta,self.__alpha)
        if self.__beta>=self.__alpha:
            sigma = 1
        else:
            sigma = -1
        r = 0
        s = 0
        for i in range(self.__N+1):
            ci = (2*mu+kappa+2*i-1)*special.gamma(2*mu+kappa+i-1)
            P = special.eval_jacobi(i,mu-1,mu+kappa-1,sigma*self.__theta)
            r += ci/special.gamma(mu+kappa+i)*P
            s +=(-1)**i*ci/special.gamma(mu+i)*P
        Omega = omega*self.__T
        T1 = np.abs(r)/np.abs(np.power(Omega,mu))
        T2 = np.abs(s)/np.abs(np.power(Omega,mu+kappa))

        # Compute upper bound
        uB = T1+T2

        # Compute lower bound
        lB = np.abs(T1-T2)

        # Compute asymptote
        if kappa==0:
            q = special.gamma(mu)*np.maximum(r,s)
        else:
            q = special.gamma(mu+kappa)*np.abs(r)

        # Compute cutoff frequency
        self.__wc = 1./self.__T*np.power(q/special.gamma(mu+kappa),1/mu)
        mB = q/special.gamma(mu+kappa)/np.power(np.abs(Omega),mu)
        #mB = np.power(self.__wc/omega,mu)
        mB[omega<self.__wc] = 1

        return uB, lB, mB

    def computeTfromWc(self,wc):
        """
        This function computes the filter window length for a desired
        cutoff frequency.

        :param wc: Cutoff frequency.
        :type wc: float

        :returns:
        	- T (:py:class:`float`) -  filter window length
        """

        kappa = np.abs(self.__beta-self.__alpha)
        mu = 1+np.minimum(self.__beta,self.__alpha)
        if self.__beta>=self.__alpha:
            sigma = 1
        else:
            sigma = -1
        r = 0
        s = 0
        for i in range(self.__N+1):
            ci = (2*mu+kappa+2*i-1)*special.gamma(2*mu+kappa+i-1)
            P = special.eval_jacobi(i,mu-1,mu+kappa-1,sigma*self.__theta)
            r += ci/special.gamma(mu+kappa+i)*P
            s +=(-1)**i*ci/special.gamma(mu+i)*P

        # Compute asymptote
        if kappa==0:
            q = special.gamma(mu)*np.maximum(np.abs(r),np.abs(s))
        else:
            q = special.gamma(mu+kappa)*np.abs(r)

        # Compute filter winfow length and convert to a muliple of ts
        self.__T = int(1./wc*np.power(q/special.gamma(mu+kappa),1/mu)\
                       /self.__ts)*self.__ts
        self.checkParameters()


    def get_delay(self):
        """
        This function returns the estimation delay :math:`\delta_t` of the\
        the algebraic differentiator. 

        :return:  :math:`\delta_t` as a float 
        """
        if self.__N==0:
            delay = (self.__alpha+1)/(self.__alpha+self.__beta+2)*self.__T
        else:
            delay = (1-self.__theta)/2*self.__T
        return delay

    def get_delayDiscrete(self):
        """
        This function returns the estimation delay :math:`\delta_{t,\mathrm{d}}` of the\
                discrete algebraic differentiator. The delay can vary depending on the\
                used discretization method.

        :return:  Dictionary with keys the used discretization schemes.\
                for each key a value of the delay is provided.
        """
        return self.__delayDisc

    def set_T(self,T):
        """
        This function sets the parameter :math:`T` of the window length of\
        the algebraic differentiator to the value T. 

        :param T: Filter window length.
        :type T: float
        
        """
        self.__T = T
        self.checkParameters()
        self.__L = np.round(self.__T/self.__ts)
        for der in self.__w.keys():
            for m in self.__w[der].keys():
                tmp = self.discretize(der,m)
    
    def set_samplingPeriod(self,ts):
        """
        This function sets the sampling period :math:`t_s` of the discrete
        measurements. 

        :param ts: Sampling period.
        :type ts: float
        
        """

        self.__ts = ts

    def set_alpha(self,alpha):

        """
        This function sets the parameter :math:`\\alpha` of the weight\
        function of the algebraic\
        differentiator to the value alpha. The value must satisfy\
        :math:`\\alpha>n-1`, with :math:`n` the highest derivative to be estimated.

        :param alpha: Parameter of the weight function.
        :type alpha: float
        
        """

        self.__alpha = alpha
        for der in self.__w.keys():
            for m in self.__w[der].keys():
                tmp = self.discretize(der,m)


    def set_beta(self,beta):
        """
        This function sets the parameter :math:`\\beta` of the weight\
        function of the algebraic\
        differentiator to the value beta. The value must satisfy\
        :math:`\\beta>n-1`, with :math:`n` the derivative to be estimated.\ 
        If a filter was discretized earlier, the functions computes the\
        discretized filter with this new parameter.


        :param beta: Parameter of the weight function.
        :type beta: float
        """

        self.__beta = beta
        for der in self.__w.keys():
            for m in self.__w[der].keys():
                tmp = self.discretize(der,m)


    def set_N(self,N): 
        """
        This function sets the parameter :math:`N` of the algebraic\
        differentiator to the value N. If a filter was discretized\
        earlier, the functions computes the discretized filter with this new\
        parameter.


        :param N: Order of the approximating polynomial.
        :type N: int
        """

        self.__N = N
        for der in self.__w.keys():
            for m in self.__w[der].keys():
                tmp = self.discretize(der,m)

    
    def set_theta(self,theta,rootJacobiPol):
        """
        This function sets the parameter :math:`\\vartheta` of the algebraic\
        differentiator to the value theta. The differentiator is parametrized 
        by default such that :math:`\\vartheta` is the largest zero of the
        Jacobi polynomial of degree :math:`N+1`. Thus, the delay is 
        minimized. For a delay-free approximation 
        the choice :math:`\\vartheta=1` is required. If a filter was discretized\
        earlier, the functions computes the discretized filter with this new\
        parameter.


        :param theta: Parameter :math:`\\vartheta` of the algebraic
            differentiator.
        :type theta: float
        :param rootJacobiPol: True if the value theta is equal to a zero of\
        the jacobi polynomial of order :math:`N+1`. For the delay-free
            approximation use False, since :math:`1` is not a zero.
        :type rootJacobiPol: bool
        """

        self.__thetaBool = rootJacobiPol
        self.__theta = theta
        for der in self.__w.keys():
            for m in self.__w[der].keys():
                tmp = self.discretize(der,m)
    
        # Recompute window length from cut off frequency
        if self.__setUp=='wc':
            self.computeTfromWc(self.__wc)
        elif self.__setUp=='T':
            self.__wc = self.get_cutoffFreq()


    def printParam(self):
        """This function prints all the parameters of the algebraic
        differentiator.
        """

        print("The differentiator has the parameters:")
        print("Alpha: %2.6f"%(self.__alpha))
        print("Beta: %2.6f"%(self.__beta))
        print("Window length in s: %2.6f"%(self.__T))
        print("Sampling period in s: %2.6f"%(self.__ts))
        print("Polynomial degree: %d"%(self.__N))
        print("Estimation delay in s: %2.6f"%(self.get_delay()))
        print("Cutoff Frequency in rad/s: %3.6f"%(self.get_cutoffFreq()))
        print("Cutoff Frequency in Hz: %3.6f"%(self.get_cutoffFreq()/2/np.pi))
        print("Discrete window length: %d"%(int(self.__T/self.__ts)))

    def get_ratioNyquistCutoff(self,k):
        """This function computes the ratio
        :math:`k_N=20\\log_{10}\\left(\\frac{\omega_N^k\\big|\mathcal{G}(\omega_N)\\big|}{\omega_c^k\Big|\mathcal{G}(\omega_c)\Big|}\\right)`, with :math:`\omega_N=\pi/t_s`, the Nyquist frequency and :math:`\omega_c` the                cutoff frequency.

        :param k: Orders of derivatives to be estimated.
        :type k: list of natural numbers
        :return:  :math:`k_N` in a numpy array with same dimensions as the
            variable k 
        """

        wc = self.get_cutoffFreq()
        wN = np.pi/self.__ts
        G, tmp = self.get_ampAndPhaseFilter(np.array([wc,wN]))

        return 20*np.log10(np.power(wN,k)*G[1]/np.power(wc,k)/G[0])
    
    def get_discretizationError(self,k,Omega,n=1000,method='mid-point'):
        """This function computes the error
        :math:`\\mathcal{J}=\\frac{\\int_{0}^{\\Omega}\\left|{\\mathcal{F}\\left\\{g^{(n)}-\\hat{g} ^{(n)}\\right\\}(\\omega)}\\right|^2\\mathrm{d}\\omega}
        {\\int_{0}^{\\Omega}\\left|{\\mathcal{F}\\left\\{g^{(n)}\\right\\}(\\omega)}\\right|^2\\mathrm{d}\\omega}`. 
        Therein :math:`g` and :math:`\\hat{g}` are the continuous and discrete time differentiators, respectively. 
        A reasonable choice for :math:`\\Omega_N` is the Nyquist frequency :math:`\\pi/t_s`. 
        The integral is approximated using the trapezoidal rule.

        :param k: Order of derivative to be estimated.
        :type k: integer
        :param Omega: Upper bound of integration interval.
        :type Omega: float
        :param n: The interval :math:`[0,\\Omega]` is divided into n parts. 
        :type n: integer
        :return:  cost function :math:`\\mathcal{J}`
        """

        omega = np.linspace(0,Omega,n)
        amp,phase = self.get_ampAndPhaseFilter(omega)
        F = amp*np.exp(1j*phase)*(1j*omega)**k
        ampDisc,phaseDisc = self.get_ampSpectrumDiscreteFilter(omega,k,\
                                        method)
        Fd = ampDisc*np.exp(1j*phaseDisc)
        e = np.trapz(np.abs(Fd-F)**2,x=omega)/np.trapz(np.abs(F)**2,x=omega)
        return e

    def reduceFilterLength(self,der,tol=0.01):
        """ This function returns a distribution function that can
        be used to reduce the filter window length. The length is reduced such
        that the truncation is of the same amount on both sides of the filter
        window.

        :param der: The sought derivative.
        :type der: int
        
        :param tol: The tolerance that for reducing the filter length
        :type tol: float
        
        :returns:
        	- tau1 (:py:class:`float`) - New starting point of the window. The value 0, i.e., the old
            		starting point, is taken as the reference.
        	- tau2 (:py:class:`float`) - New end point of the window. The value 0, i.e., the old starting
           	 	point, is taken as the reference.
        	- d (:py:class:`numpy array`) - Distribution function used for the truncation.
        	- t (:py:class:`numpy array`) - Times instants where distribution function has been evaluated.
        	- n (:py:class:`int`) - Factor relating the sampling rate of the measurements and the
			    sampling rate used for the evaluation of the distribution
			    function.
        """
        # Compute the first moment of the n-th derivative of the kernel
        n = 10
        t = np.arange(0,self.__T,self.__ts/n)
        g = self.evalKernelDer(t,der)
        m0 = np.trapz(np.abs(g),x=t)

        # Compute the distribution function
        disFcn =\
        np.array([np.trapz(np.abs(self.evalKernelDer(t[:i],der)),x=t[:i])\
                  for i in range(1,len(t)+1)])/m0

        # Find values where to reduce window length
        tau1 = np.argmax(disFcn>=tol)
        tau2 = np.argmax(disFcn>=1-disFcn[tau1])
        return t[tau1], t[tau2], disFcn, t, n
