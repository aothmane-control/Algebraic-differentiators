classdef AlgDiff <  matlab.mixin.CustomDisplay
    %ALGDIFF Wrapper for py.AlgDiff.AlgebraicDifferentiator.
    %
    %   This class implements an algebraic differentiator and all the
    %   methods necessary for its use, analysis, tuning, and discretization.
    %   The differentiators which are LTI filters can be parametrized to
    %   achieve desired filter characteristics: Cutoff frequency and
    %   stopband slope.
    %
    %   For more information and installation instructions, see <a href=
    %   "matlab:web('https://github.com/aothmane-control/Algebraic-differentiators')">the source repository</a>.
    %
    %   See also py.AlgDiff.algebraicDifferentiator

    properties (Constant, GetAccess=protected)
        % Error Identifier
        ErrID = "AlgDiff"
    end

    properties(SetAccess=immutable, GetAccess=public)
        % Instance of py.AlgDiff.AlgebraicDifferentiator
        inst
    end
    
    %% Construction
    methods
        function obj = AlgDiff(ts, alpha, beta, N, opts)
            %ALGDIFF Constructs an algebraic differentiator.
            %   Either a cutoff frequency, omega_c ('CutoffFrequency'), or a
            %   filter window length, T ('FilterWindowLength'), can be
            %   specified. If the filter window length is specified the
            %   cutoff frequency is computed. For frequencies lower than
            %   the cutoff frequency the amplitude of the Fourier transform
            %   is 0 dB. For frequencies higher than the cutoff frequency
            %   the amplitude falls with min(alpha, beta)+1 dB per decade.
            %
            %   Arguments:
            %     - ts: Sampling period.
            %     - alpha: Parameter alpha of the weight function of the
            %         Jacobi polynomials. Has to satisfy alpha > n - 1, 
            %         where n is highest derivative to be estimated.
            %     - beta: Parameter beta of the weight function of the
            %         Jacobi polynomials. Has to satisfy beta > n - 1, 
            %         where n: highest derivative to be estimated. The
            %         stopband slope is given by mu=min(alpha, beta)+1, i.e.
            %         frequencies higher than the cutoff frequency are
            %         attenuated by 20*mu dB per decade. Frequencies lower
            %         than the cutoff frequency are attenuated by 0 dB.
            %     - N: Truncation order (int) of the generalized Fourier 
            %         series. A delay-free derivative approximation is only
            %         possible for N >= 1. The differentiator is
            %         parametrized by default such that the non-zero delay
            %         is minimized for the choice N >= 1. See the method
            %         set_theta for more details.
            %
            %   Key-value arguments:
            %     - FilterWindowLength: Filter window length (float). Takes
            %         a positive value if the length has to be specified.
            %         The cutoff frequency is then computed automatically.
            %         It should take the value [] if the cutoff frequency
            %         is specified.
            %     - CutoffFrequency: Cutoff frequency w_c of the 
            %         differentiator in rad (float). Takes a positive value 
            %         if the cutoff frequency has to be specified. The 
            %         filter window length is then computed automatically. 
            %         It should take the value [] if the filter window 
            %         length is specified.
            %     - Correction: Boolean variable that indicates if errors
            %         in the DC component of the approximated signal
            %         stemming from the discretization should be corrected.
            %
            %   See also py.AlgDiff.algebraicDifferentiator.AlgebraicDifferentiator

            arguments
                ts (1,1) {mustBeReal, mustBePositive} = 0.01
                alpha (1,1) {mustBeReal} = 1
                beta (1,1) {mustBeReal} = 1
                N (1,1) {mustBeInteger, mustBeNonnegative} = 0

                opts.Correction (1,1) matlab.lang.OnOffSwitchState = true
                opts.FilterWindowLength {mustBeScalarOrEmpty} = []
                opts.CutoffFrequency {mustBeScalarOrEmpty} = []

            end

            if isempty(opts.FilterWindowLength) ...
                && isempty(opts.CutoffFrequency)

                % Default to FilterWindowLength=1 if neither is provided
                opts.FilterWindowLength = 1;

            end

            % Either FilterWindowLength or CutoffFrequency has to be
            % specified.
            if ~isempty(opts.FilterWindowLength)
                if ~isempty(opts.CutoffFrequency)
                    error(AlgDiff.ErrID + ":algdiff:tooManyArguments", ...
                        "Either FilterWindowLength or CutoffFrequency " ...
                        + "can be specified, but not both");
                end

                if ~isreal(opts.FilterWindowLength) ...
                        || opts.FilterWindowLength <= 0
                    error(AlgDiff.ErrID + ":algdiff:FilterWindowLength", ...
                        "FilterWindowLength has to be positive and real");
                end

                opts.CutoffFrequency = string(missing);

            else
                if isempty(opts.CutoffFrequency)
                    error(AlgDiff.ErrID + ":algdiff:tooFewArguments", ...
                        "Either FilterWindowLength or CutoffFrequency " ...
                        + "has to be specified");
                end

                if ~isreal(opts.CutoffFrequency) ...
                        || opts.CutoffFrequency <= 0
                    error(AlgDiff.ErrID + ":algdiff:CutoffFrequency", ...
                        "CutoffFrequency has to be positive and real");
                end

                opts.FilterWindowLength = string(missing);

            end

            % Instantiate python class
            try
                obj.inst = py.AlgDiff.AlgebraicDifferentiator(...
                    pyargs( ...
                        "ts", ts, ...
                        "alpha", alpha, ...
                        "beta", beta, ...
                        "N", int32(N), ...
                        "T", opts.FilterWindowLength, ...
                        "wc", opts.CutoffFrequency, ...
                        "corr", logical(opts.Correction), ...
                        "display", false ...
                    ));
            catch ME
                if strcmp(ME.identifier, "MATLAB:undefinedVarOrClass")
                    % Missing library
                    missingLibException = MException( ...
                        AlgDiff.ErrID + ":algdiff:MissingPythonLibrary", ...
                        "Could not load the AlgDiff library from python. " ...
                        + "See https://github.com/aothmane-control/Algebraic-differentiators" ...
                        + " for installation instructions");

                    throw(missingLibException.addCause(ME));
                end

                rethrow(ME);
            end

        end
    end

    %% Estimation
    methods
        function der = estimateDer(obj, k, x, opts)
            %ESTIMATEDER Estimates the derivative of order 'k' >= 0 using
            %the given discretization method. If the filter was not yet
            %discretized prior this call, this functions performs the
            %discretization.
            %   Arguments:
            %     - k: Order of the derivative to be estimated.
            %     - x: Signal (vector) whose derivative has to be
            %       estimated.
            %
            %   Key-value arguments:
            %     - Conv: Parameter specifying the type of convolution the
            %       numpy conv function performs. See the numpy docs for
            %       more details.
            %
            %   Returns:
            %     - der: Estimated derivative as a vector of same length as
            %       the signal x.
            %
            %   See also AlgDiff.discretize

            arguments
                obj (1,1) AlgDiff
                k (1,1) {mustBeInteger,mustBeNonnegative}
                x {mustBeVector}

                opts.Method {mustBeMember(opts.Method, ["mid-point",...
                    "trapezoidal", "analytic", "simpson rule"])} = "mid-point";
                opts.Conv {mustBeTextScalar} = "same";
                opts.ReduceFilterLength (1,1) matlab.lang.OnOffSwitchState ...
                    = "off";
                opts.RedTol (1,1) {mustBeNonnegative} = 0.01;
            end
            
            % Always returns a row vector
            der = obj.inst.estimateDer(int32(k), x, ...
                pyargs("method", opts.Method, ...
                       "conv", opts.Conv, ...
                       "redFilLength", logical(opts.ReduceFilterLength), ...
                       "redTol", opts.RedTol));

            der = double(der);
        end

        function [coeffs, tau1, tau2, theta] = discretize(obj, der, opts)
            %DISCRETIZE Performs the discretization of the filter for the
            %estimation of the derivative of order 'der'. 
            %   The following discretization schemes are implemented: 
            %   midpoint, trapezoidal and using analytical integration. 
            %   They may be set using the key-value parameter 'Method'. The 
            %   midpoint rule uses one filter coefficient less than the 
            %   trapezoidal rule. It also reduces the estimation delay by 
            %   half a sampling period. The analytical integration rule is 
            %   recommended for small filter window lengths, i.e. in 
            %   general less than 20 filter coefficients. This 
            %   discretization method can only be used for the first order 
            %   or any higher derivative. The error stemming from the 
            %   discretization is corrected using a correction factor, if 
            %   the differentiator has been initialized to do so.
            %
            %   When the parameters 'alpha' and 'beta' are large, the
            %   filter kernel has very low values at the beginning and at
            %   the end of the filter window. The length of the window can
            %   be reduced to save computation time and memory. This can be
            %   performed by enabling the redFilterLength parameter and
            %   giving a tolerance by enabling the key-value parameter
            %   'ReduceFilterLength' and giving a tolerance by setting
            %   'RedTol'. This tolerance automatically computes the size of
            %   the new window such that the truncation is performed
            %   symmetrically. If this is enabled, the functions also
            %   returns the times at which the filter was truncated.
            %
            %   Arguments:
            %     - der: Order of the derivative to be estimated (int >= 0)
            %
            %   Key-value arguments:
            %     - Method: Discretization scheme. One of "mid-point",
            %       "trapezoidal", "analytic" or "simpson rule". Defaults
            %       to "mid-point".
            %     - ReduceFilterLength: true/false or "on"/"off". Wether or
            %       not to reduce the filter window length. Defaults to
            %       false.
            %     - RedTol: Tolerance to be used if ReduceFilterLength is
            %       set to true. Defaults to 0.01.
            %     - DiscreteSpectrum: true/false or "on"/"off". If set the
            %       parameter 'theta' used in the discretization is
            %       returned. See the survey paper for more details.
            %
            %   Returns:
            %     - coeffs: Discretized filter for the given derivative and
            %       discretization method. If the correction of the DC
            %       component has been enabled with the parameter corr of
            %       the class construction, the output are the corrected
            %       filter coefficients.
            %     - tau1: If ReduceFilterLength is set contains the time
            %       where the filter window is reduced before at the left
            %       side of the interval. The estimation delay is reduced
            %       to tau1.
            %     - tau2: If ReduceFilterLength is set contains the time
            %       where the filter window is reduced before at the right
            %       side of the interval. This value does not affect the
            %       delay.
            %     - theta: Returned if DiscreteSpectrum is set.

            arguments
                obj (1,1) AlgDiff
                der (1,1) {mustBeInteger,mustBeNonnegative}

                opts.Method {mustBeMember(opts.Method, ["mid-point", ...
                    "trapezoidal", "simpson rule",  "analytic"])} ...
                    = "mid-point";
                opts.ReduceFilterLength (1,1) matlab.lang.OnOffSwitchState ...
                    = "off";
                opts.RedTol (1,1) {mustBeNonnegative} = 0.01;
                opts.DiscreteSpectrum (1,1) matlab.lang.OnOffSwitchState ...
                    = "off";
            end

            % Tuple of varying length
            data = obj.inst.discretize(int32(der), ...
                pyargs("method", opts.Method, ...
                       "redFilLength", logical(opts.ReduceFilterLength), ...
                       "redTol", opts.RedTol, ...
                       "discreteSpectrum", logical(opts.DiscreteSpectrum)));
            
            % Only return sought-after coefficients
            % Check if a tuple is returned or just the dict
            if opts.ReduceFilterLength || opts.DiscreteSpectrum
                coeffs = data{1};
            else
                coeffs = data;
            end
            
            % Get coefficients
            if opts.ReduceFilterLength
                coeffs = double(coeffs{der}{opts.Method + "-red"});
            else
                coeffs = double(coeffs{der}{opts.Method});
            end

            if opts.ReduceFilterLength
                tau1 = double(data{2});
                tau2 = double(data{3});

            else
                tau1 = [];
                tau2 = [];

            end

            if opts.DiscreteSpectrum
                theta = double(data{end});

            else
                theta = [];

            end

        end
    end

    %% Special functions
    methods
        function [amp, phase] = get_ampAndPhaseFilter(obj, omega)
            %GET_AMPANDPHASEFILTER Evaluates the amplitude and the phase of
            %the Fourier transform of the algebraic differentiator at the
            %frequencies 'omega'.
            %   Arguments:
            %     - omega: Vector of frequencies where the amplitude and
            %       the phase should be evaluated.
            %
            %   Returns:
            %     - amp: Amplitude of the Fourier transform of the filter
            %       as a vector of same length as 'omega'.
            %     - phase: Phase of the Fourier transform of the filter
            %       as a vector of same length as 'omega'.

            arguments
                obj (1,1) AlgDiff
                omega {mustBeVector,mustBeReal}
            end

            % Returns tuple
            data = obj.inst.get_ampAndPhaseFilter(py.numpy.array(omega));
            
            amp = double(data{1});
            phase = double(data{2});

        end

        function [amp, phase] = get_ampSpectrumDiscreteFilter(obj, omega, n, opts)
            arguments
                obj (1,1) AlgDiff
                omega {mustBeVector,mustBeReal}
                n (1,1) {mustBeInteger,mustBeNonnegative}

                opts.Method {mustBeMember(opts.Method, ["mid-point",...
                    "trapezoidal", "analytic"])} = "mid-point";
            end

            % Returns tuple
            data = obj.inst.get_ampSpectrumDiscreteFilter(py.numpy.array(omega), ...
                int32(n), pyargs("method", opts.Method));
            
            amp = double(data{1});
            phase = double(data{2});

        end

        function response = get_stepResponse(obj, t)
            %GET_STEPRESPONSE Returns the step response of the
            %differentiator evaluated at the time instants in 't'.
            %   Arguments:
            %     - t: Time instants (as a vector) where the step response
            %       should be evaluated at.
            %
            %   Returns:
            %     - response: The response as a vector of same length as
            %       't'.

            arguments
                obj (1,1) AlgDiff
                t {mustBeVector,mustBeReal}
            end

            response = double(obj.inst.get_stepResponse(py.numpy.array(t)));
        end

        function kern = evalKernel(obj, t)
            %EVALKERNEL Evaluates the kernel of the algebraic
            %differentiator at time instants 't'. It corresponds to the
            %impulse response of the filter.
            %   Arguments:
            %     - t: Time instants (as a vector) where the kernel should
            %       be evaluated at.
            %
            %   Returns:
            %     - kern: Evaluated kernel as a vector of same length as
            %       't'.

            arguments
                obj (1,1) AlgDiff
                t {mustBeVector,mustBeReal}
            end

            kern = double(obj.inst.evalKernel(py.numpy.array(t)));
        end

        function kern = evalKernelDer(obj, t, k)
            %EVALKERNELDER Evaluates the derivative of order 'k' of the
            %algebraic differentiator at time instants 't'. If k is zero,
            %the kernel itself is evaluated.
            %   Arguments:
            %     - t: Time instants (as a vector) where the kernel should
            %       be evaluated at.
            %     - k: Order of the derivative to be evaluated (int).
            %
            %   Returns:
            %     - kern: Evaluated derivative as a vector of same length 
            %       as 't'.

            arguments
                obj (1,1) AlgDiff
                t {mustBeVector,mustBeReal}
                k (1,1) {mustBeInteger,mustBeNonnegative}
            end

            kern = double(obj.inst.evalKernelDer(py.numpy.array(t), ...
                int32(k)));
        end

        function [uB, lB, mB] = get_asymptotesAmpFilter(obj, omega)
            arguments
                obj (1,1) AlgDiff
                omega {mustBeVector,mustBeReal}
            end

            % Returns tuple
            data = obj.inst.get_asymptotesAmpFilter(py.numpy.array(omega));
            
            uB = double(data{1});
            lB = double(data{2});
            mB = double(data{3});

        end
    
        function t = timeShift(obj, t)
            %TIMESHIFT Computes the time transformation theta(t)=1-2/T,
            %with the window length T of the filter.
            %   Arguments:
            %     - t: Time instants (as a vector) where the transformation
            %       should be evaluated at.
            %
            %   Returns:
            %     - t: Evaluated transformation as a vector of same length
            %       as given t.

            arguments
                obj (1,1) AlgDiff
                t {mustBeVector,mustBeReal}
            end

            t = double(obj.inst.timeShift(py.numpy.array(t)));
        end

        function w = weightFcn(obj, a, b, t)
            %WEIGHTFCN Evaluates the weight function of the algebraic
            %differentiator at the time instants 't'.
            %   Arguments:
            %     - a: Parameter alpha (> -1) of the weight function.
            %     - b: Parameter beta (> -1) of the weight function.
            %     - t: Time instants (as a vector) where the weight
            %       function should be evaluated at.
            %
            %   Returns:
            %     - w: Evaluated weight function as a vector of same length
            %       as t.
            
            arguments
                obj (1,1) AlgDiff
                a (1,1) {mustBeGreaterThan(a, -1)}
                b (1,1) {mustBeGreaterThan(b, -1)}
                t {mustBeVector,mustBeReal}
            end

            w = double(obj.inst.weightFcn(a, b, py.numpy.array(t)));
        end

    end

    %% Getter
    methods
        function alpha = get_alpha(obj)
            %GET_ALPHA Returns the paramater 'alpha'.

            alpha = double(obj.inst.get_alpha());
        end

        function beta = get_beta(obj)
            %GET_BETA Returns the parameter 'beta'.

            beta = double(obj.inst.get_beta());
        end

        function N = get_N(obj)
            %GET_N Returns the parameter 'N'.

            N = double(obj.inst.get_N());
        end

        function theta = get_theta(obj)
            %GET_THETA Returns the parameter 'theta'.

            theta = double(obj.inst.get_theta());
        end

        function ts = get_ts(obj)
            %GET_TS Returns the sampling period 'ts'.

            ts = double(obj.inst.get_ts());
        end

        function T = get_T(obj)
            %GET_T Returns the filter window length of the algebraic
            %differentiator.

            T = double(obj.inst.get_T());
        end

        function wc = get_cutoffFreq(obj)
            %GET_CUTOFFFREQ Returns the cutoff frequency, omega_c, of the
            %algebraic differentiator.

            wc = double(obj.inst.get_cutoffFreq());
        end

        function delay = get_delay(obj)
            %GET_DELAY Returns the estimation delay, delta_t, of the
            %algebraic differentiator.

            delay = double(obj.inst.get_delay());
        end

        function gamma = get_degreeExactness(obj, n)
            %GET_DEGREEEXACTNESS Returns the degree of exactness, gamma,
            %when the 'n'-th derivative is approximated.
            %   Arguments:
            %     -n: Order of derivative (>= 0).

            arguments
                obj (1,1) AlgDiff
                n (1,1) {mustBeInteger,mustBeNonnegative}
            end

            gamma = double(obj.inst.get_degreeExactness(int32(n)));
        end

        function k_N = get_ratioNyquistCutoff(obj, k)
            %GET_RATIONYQUISTCUTOFF Computes the ratio
            %k_N=20\\log_{10}\\left(\\frac{\omega_N^k\\big|\mathcal{G}(\omega_N)\\big|}{\omega_c^k\Big|\mathcal{G}(\omega_c)\Big|}\\right)
            %with the Nyquist frequency omega_N=\pi/t_s, and the cutoff
            %frequency \omega_c.
            %   Arguments:
            %     - k: Vector of orders of derivatives to be estimated.
            %
            %   Returns:
            %     - k_N: Ratios as a vector of same length as 'k'.

            arguments
                obj (1,1) AlgDiff
                k {mustBeVector,mustBeInteger,mustBeNonnegative}
            end

            k_N = double(obj.inst.get_ratioNyquistCutoff(py.numpy.array(k)));
        end

        function J = get_discretizationError(obj, k, omega, opts)
            %GET_DISCRETIZATIONERROR Computes the error
            %\\mathcal{J}=\\frac{\\int_{0}^{\\Omega}\\left|{\\mathcal{F}\\left\\{g^{(n)}-\\hat{g} ^{(n)}\\right\\}(\\omega)}\\right|^2\\mathrm{d}\\omega}{\\int_{0}^{\\Omega}\\left|{\\mathcal{F}\\left\\{g^{(n)}\\right\\}(\\omega)}\\right|^2\\mathrm{d}\\omega}
            %Therein g and \\hat{g} are the continous and discrete time
            %differentiators, respectively.
            %   A reasonable choice for \\Omega_N is the Nyquist frequency
            %   \\pi/t_s.
            %   The integral is approximated using the trapezoidal rule.
            %
            %   Arguments:
            %     - k: Order of derivative to be estimated.
            %     - omega: Upper bound of integration interval.
            %
            %   Key-value arguments:
            %     - NumberOfParts: The interval [0, 'omega'] is divided
            %       into 'n' parts. Defaults to 1000.
            %
            %   Returns:
            %     - J: Evaluated cost function \\mathcal{J}.
            %
            %   See also AlgDiff.discretize

            arguments
                obj (1,1) AlgDiff
                k (1,1) {mustBeInteger,mustBeNonnegative}
                omega (1,1) {mustBeFloat,mustBeNonnegative}

                opts.NumerOfParts (1,1) {mustBeInteger,mustBePositive} = 1000; 
                opts.Method {mustBeTextScalar,mustBeMember(opts.Method, ["mid-point",...
                    "trapezoidal", "analytic", "simpson rule"])} = "mid-point";
            end

            J = double(obj.inst.get_discretizationError( ...
                int32(k), omega, int32(opts.NumerOfParts), opts.Method ...
                ));
        end

        function delta_t_d = get_delayDiscrete(obj, method, opts)
            %GET_DELAYDISCRETE Returns the estimation delay delta_td of the
            %discrete algebraic differentiator. The delay can vary
            %depending on the chosen discretization method.
            %   
            %   See also AlgDiff.discretize

            arguments
                obj (1,1) AlgDiff
                method {mustBeTextScalar,mustBeMember(method, ["mid-point",...
                    "trapezoidal", "analytic", "simpson rule"])}
                
                opts.ReduceFilterLength (1,1) matlab.lang.OnOffSwitchState ...
                    = "off";
            end

            % Dictionary key
            key = method;
            if opts.ReduceFilterLength
                key = key + "-red";
            end

            % Access dict
            data = obj.get_property_delayDisc();
            value = data.get(key);

            if isnumeric(value)
                delta_t_d = value;

            else
                % Not in dictionary
                error(AlgDiff.ErrID + ":algdiff:MissingDiscretization", ...
                      "Method '%s' with ReduceFilterLength='%s'" ...
                      + " was not yet discretized", ...
                      method, opts.ReduceFilterLength);

            end
            
        end
    end

    %% Setter
    methods
        function set_theta(obj, theta, rootJacobiPol)
            %SET_THETA Sets the parameter theta of the algebraic
            %differentiator to the value 'theta'.
            %   The differentiator is parametrized by default such that
            %   theta is the largest zero of the Jacobi polynomial of
            %   degree N+1. Thus, the delay is minimized.
            %   For a delay-free approximation the choice theta=1 is
            %   required.
            %   If a filter was discretized earlier, this function also
            %   computes the discretized filter with the new parameter.
            %
            %   Arguments:
            %     - theta: Parameter theta of the algebraic differentiator.
            %     - rootJacobiPol: True if the value theta is equal to a
            %       zero of the Jacobi polynomial of order N+1. For a
            %       delay-free approximation use False, as theta=1 is not a
            %       zero.

            arguments
                obj (1,1) AlgDiff
                theta (1,1) {mustBeReal,mustBeFinite,mustBeFloat}
                rootJacobiPol (1,1) matlab.lang.OnOffSwitchState
            end

            obj.inst.set_theta(theta, logical(rootJacobiPol));

        end

        function set_T(obj, T)
            %SET_T Sets the window length T of the algebraic differentiator
            %to the value 'T'.
            %   Arguments:
            %     - T: Filter window length

            arguments
                obj (1,1) AlgDiff
                T {mustBePositive,mustBeReal}
            end

            obj.inst.set_T(double(T));
        end

        function set_samplingPeriod(obj, ts)
            %SET_SAMPLINGPERIOD Sets the sampling period ts of the discrete 
            %measurements to the value 'ts'.
            %   Arguments:
            %     - ts: Sampling period.

            arguments
                obj (1,1) AlgDiff
                ts {mustBePositive,mustBeReal}
            end

            obj.inst.set_samplingPeriod(double(ts));
        end

        function set_alpha(obj, alpha)
            %SET_ALPHA Sets the parameter alpha of the weight function of
            %the algebraic differentiator to the value 'alpha'.
            %   The value must satisfy alpha > n-1, where n it the highest
            %   derivative to be estimated

            arguments
                obj (1,1) AlgDiff
                alpha {mustBeReal}
            end

            obj.inst.set_alpha(double(alpha));
        end

        function set_beta(obj, beta)
            %SET_BETA Sets the parameter beta of the weight function of
            %the algebraic differentiator to the value 'beta'.
            %   The value must satisfy beta > n-1, where n it the highest
            %   derivative to be estimated.
            %   If a filter was discretized earlier, the function also
            %   computes the discretized filter with this new parameter.

            arguments
                obj (1,1) AlgDiff
                beta {mustBeReal}
            end

            obj.inst.set_beta(double(beta));
        end

        function set_N(obj, N)
            %SET_N Sets the parameter N of the weight function of
            %the algebraic differentiator to the value 'N'.
            %   If a filter was discretized earlier, the function also
            %   computes the discretized filter with this new parameter.

            arguments
                obj (1,1) AlgDiff
                N {mustBeInteger,mustBeNonnegative}
            end

            obj.inst.set_N(int32(N));
        end
    end

    %% Custom methods
    methods      
        function coeff = get_filter_coefficients(obj, der, opts)
            %GET_FILTER_COEFFICIENTS Returns the discrete filter
            %coefficients for the 'der'-th. derivative.
            %   A call with the similar parameters to discretize had to be
            %   done before.
            %
            %   See also AlgDiff.discretize

            arguments
                obj (1,1) AlgDiff
                der (1,1) {mustBeInteger,mustBeNonnegative}

                opts.Method {mustBeMember(opts.Method, ["mid-point", ...
                    "trapezoidal", "simpson rule",  "analytic"])} ...
                    = "mid-point";
                opts.ReduceFilterLength (1,1) matlab.lang.OnOffSwitchState ...
                    = "off";
            end

            % Query dictionary
            w = obj.get_property_w();
            
            % First key
            data = w.get(int32(der));
            if isa(data, "py.NoneType")
                % Not in dictionary
                error(AlgDiff.ErrID + ":algdiff:MissingDiscretization", ...
                      "%d-th derivative was not yet discretized", ...
                      der);

            end

            % Second key
            key = opts.Method;
            if opts.ReduceFilterLength
                key = key + "-red";
            end

            data = data.get(key);
            if isa(data, "py.NoneType")
                % Not in dictionary
                error(AlgDiff.ErrID + ":algdiff:MissingDiscretization", ...
                      "Method '%s' with ReduceFilterLength='%s'" ...
                      + " was not yet discretized for %d-th derivative", ...
                      opts.Method, opts.ReduceFilterLength, der);

            end

            coeff = double(data);
        end

    end

    %% Display (i.e. disp(obj))
    methods (Access = protected)

        function s = getFooter(obj)
            %GETFOOTER Prints properties of the implemented differentiator
            %into a string

            s = sprintf(['\tParameters of the differentiator:\n', ...
                            '\t\tAlpha: %.6f\n', ...
                            '\t\tBeta: %.6f\n', ...
                            '\t\tWindow length in s: %.6f\n', ...
                            '\t\tSampling period in s: %.6f\n', ...
                            '\t\tPolynomial degree: %d\n', ...
                            '\t\tEstimation delay in s: %.6f\n', ...
                            '\t\tCutoff Frequency in rad/s: %.6f\n', ...
                            '\t\tCutoff Frequency in Hz: %.6f\n', ...
                            '\t\tDiscrete window length: %d\n'], ...
                        obj.get_alpha(), ...
                        obj.get_beta(), ...
                        obj.get_T(), ...
                        obj.get_ts(), ...
                        obj.get_N(), ...
                        obj.get_delay(), ...
                        obj.get_cutoffFreq(), ...
                        obj.get_cutoffFreq() / (2*pi), ...
                        int32(obj.get_T() / obj.get_ts()));
        end

    end
    
    %% Property access helpers
    methods (Access = protected)

        function w = get_property_w(obj)
            %GET_PROPERTY_W Returns the property '__w' as a py.dict

            w = py.getattr(obj.inst, "_AlgebraicDifferentiator__w");
        end

        function delayDisc = get_property_delayDisc(obj)
            %GET_PROPERTY_DELAYDISC Returns the property '__delayDisc' as a
            %py.dict
            
            delayDisc = py.getattr(obj.inst, "_AlgebraicDifferentiator__delayDisc");
        end
    end
end

