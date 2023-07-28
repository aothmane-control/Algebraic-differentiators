classdef AlgDiff <  matlab.mixin.CustomDisplay
    %ALGDIFF Wrapper for py.AlgDiff.AlgebraicDifferentiator.
    %   For more information and installation instructions, see <a href=
    %   "matlab:web('https://github.com/aothmane-control/Algebraic-differentiators')">
    %   the source repository</a>.

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
            %ALGDIFF Instantiates an object of the python class
            %AlgebraicDifferentiator
            %   For more information see AlgebraicDifferentiator.__init__
            %
            %   Arguments:
            %       - ts: sampling period in seconds
            %       - alpha: parameter alpha of Jacobi polynomials. Has to
            %       satisfy alpha > n - 1, where n: highest derivative to
            %       be estimated
            %       - beta: parameter beta of Jacobi polynomials. Has to
            %       satisfy beta > n - 1, where n: highest derivative to be
            %       estimated
            %       - N: truncation order of generalized Fourier series
            %
            %   Key-value arguments:
            %       - Correction: See corr@AlgebraicDifferentiator.__init__ 
            %       - FilterWindowLength: See T@AlgebraicDifferentiator.__init__ 
            %       - CutoffFrequency: See wc@AlgebraicDifferentiator.__init__ 
            %
            arguments
                ts (1,1) {mustBeReal, mustBePositive} = 0.01
                alpha (1,1) {mustBeReal} = 1
                beta (1,1) {mustBeReal} = 1
                N (1,1) {mustBeInteger, mustBeNonnegative} = 0

                opts.Correction (1,1) matlab.lang.OnOffSwitchState = true
                opts.FilterWindowLength {mustBeScalarOrEmpty} = 1
                opts.CutoffFrequency {mustBeScalarOrEmpty} = []

            end

            % Either FilterWindowLength or CutoffFrequency has to be
            % specified
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
                if contains(ME.message, "Unable to resolve the name")
                    fprintf("Did not find AlgDiff library on python's path.\n"...
                          + "See: https://github.com/aothmane-control/Algebraic-differentiators"...
                          + " for installation instructions.");
                end

                rethrow(ME);
            end

        end
    end

    %% Estimation
    methods
        function der = estimateDer(obj, k, x, opts)
            %ESTIMATEDER    Estimates the derivative of order <k> of the
            %signal <x>.
            %
            %
            arguments
                obj (1,1) AlgDiff
                k (1,1) {mustBeInteger,mustBeNonnegative}
                x {mustBeVector}

                opts.Method {mustBeMember(opts.Method, ["mid-point",...
                    "trapezoidal", "analytic"])} = "mid-point";
                opts.Conv {mustBeMember(opts.Conv, ["same", "valid"])} ...
                    = "same"; % TODO: check numpy
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
            arguments
                obj (1,1) AlgDiff
                t {mustBeVector}
            end

            response = double(obj.inst.get_stepResponse(py.numpy.array(t)));
        end

        function kern = evalKernel(obj, t)
            arguments
                obj (1,1) AlgDiff
                t {mustBeVector}
            end

            kern = double(obj.inst.evalKernel(py.numpy.array(t)));
        end

        function kern = evalKernelDer(obj, t, k)
            arguments
                obj (1,1) AlgDiff
                t {mustBeVector}
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
            t = double(obj.inst.timeShift(py.numpy.array(t)));
        end
    end

    %% Getter
    methods
        function alpha = get_alpha(obj)
            %GET_ALPHA Returns the paramater <alpha> of the differentiator.
            %
            alpha = double(obj.inst.get_alpha());
        end

        function beta = get_beta(obj)
            %GET_BETA Returns the paramater <beta> of the differentiator.
            %
            beta = double(obj.inst.get_beta());
        end

        function N = get_N(obj)
            %GET_N Returns the paramater <N> of the differentiator.
            %
            N = double(obj.inst.get_N());
        end

        function theta = get_theta(obj)
            %GET_THETA Returns the paramater <theta> of the differentiator.
            %
            theta = double(obj.inst.get_theta());
        end

        function ts = get_ts(obj)
            %GET_TS Returns the sampling period.
            %
            ts = double(obj.inst.get_ts());
        end

        function T = get_T(obj)
            %GET_T Returns the filter window length of the algebraic
            %differentiator
            %
            T = double(obj.inst.get_T());
        end

        function wc = get_cutoffFreq(obj)
            %GET_CUTOFFFREQ Returns the cutoff frequency, omega_c, of the
            %algebraic differentiator
            %
            wc = double(obj.inst.get_cutoffFreq());
        end

        function delay = get_delay(obj)
            %GET_DELAY Returns the estimation delay, delta_t, of the
            %algebraic differentiator
            %
            delay = double(obj.inst.get_delay());
        end

        function gamma = get_degreeExactness(obj, n)
            %GET_DEGREEEXACTNESS Returns the degree of exactness, gamma,
            %when the <n>-th derivative is approximated
            %
            arguments
                obj (1,1) AlgDiff
                n (1,1) {mustBeInteger,mustBeNonnegative}
            end

            gamma = double(obj.inst.get_degreeExactness(int32(n)));
        end

        function k_N = get_ratioNyquistCutoff(obj, k)
            arguments
                obj (1,1) AlgDiff
                k {mustBeVector,mustBeInteger,mustBeNonnegative}
            end

            k_N = double(obj.inst.get_ratioNyquistCutoff(py.numpy.array(k)));
        end

        function J = get_discretizationError(obj, k, omega, opts)
            arguments
                obj (1,1) AlgDiff
                k (1,1) {mustBeInteger,mustBeNonnegative}
                omega (1,1) {mustBeFloat,mustBeNonnegative}

                opts.NumerOfParts (1,1) {mustBeInteger,mustBePositive} = 1000; 
                opts.Method {mustBeTextScalar,mustBeMember(opts.Method, ["mid-point",...
                    "trapezoidal", "analytic"])} = "mid-point";
            end

            J = double(obj.inst.get_discretizationError( ...
                int32(k), omega, int32(opts.NumerOfParts), opts.Method ...
                ));
        end

        function delta_t_d = get_delayDiscrete(obj, method, opts)
            arguments
                obj (1,1) AlgDiff
                method {mustBeTextScalar,mustBeMember(method, ["mid-point",...
                    "trapezoidal", "analytic"])}
                
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

        function coeff = get_filter_coefficients(obj, der, opts)
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

    %% Setter
    methods
        function set_theta(obj, theta, rootJacobiPol)
            arguments
                obj (1,1) AlgDiff
                theta (1,1) {mustBeReal,mustBeFinite,mustBeFloat}
                rootJacobiPol (1,1) matlab.lang.OnOffSwitchState
            end

            obj.inst.set_theta(theta, logical(rootJacobiPol));

        end

        function set_T(obj, T)
            %SET_T Sets the window length of the algebraic differentiator
            %to the value <T>
            %
            arguments
                obj (1,1) AlgDiff
                T {mustBePositive,mustBeReal}
            end

            obj.inst.set_T(double(T));
        end

        function set_samplingPeriod(obj, ts)
            %SET_SAMPLINGPERIOD Sets the sampling period of the discrete 
            %algebraic differentiator to the value <ts>
            %
            arguments
                obj (1,1) AlgDiff
                ts {mustBePositive,mustBeReal}
            end

            obj.inst.set_samplingPeriod(double(ts));
        end

        function set_alpha(obj, alpha)
            %SET_ALPHA Sets the parameter alpha of the weight function of
            %the algebraic differentiator to the value <alpha>.
            %   The value must satisfy alpha > n-1, where n it the highest
            %   derivative to be estimated
            %
            arguments
                obj (1,1) AlgDiff
                alpha {mustBeReal}
            end

            obj.inst.set_alpha(double(alpha));
        end

        function set_beta(obj, beta)
            %SET_BETA Sets the parameter beta of the weight function of
            %the algebraic differentiator to the value <beta>.
            %   The value must satisfy beta > n-1, where n it the highest
            %   derivative to be estimated.
            %   If a filter was discretized earlier, the function also
            %   computes the discretized filter with this new parameter.
            %
            arguments
                obj (1,1) AlgDiff
                beta {mustBeReal}
            end

            obj.inst.set_beta(double(beta));
        end

        function set_N(obj, N)
            %SET_N Sets the parameter N of the weight function of
            %the algebraic differentiator to the value <N>.
            %   If a filter was discretized earlier, the function also
            %   computes the discretized filter with this new parameter.
            %
            arguments
                obj (1,1) AlgDiff
                N {mustBeInteger,mustBeNonnegative}
            end

            obj.inst.set_N(int32(N));
        end
    end

    %% Display (i.e. disp(obj))
    methods (Access = protected)
        function s = getFooter(obj)
            %GETFOOTER
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
            w = py.getattr(obj.inst, "_AlgebraicDifferentiator__w");
        end

        function delayDisc = get_property_delayDisc(obj)
            delayDisc = py.getattr(obj.inst, "_AlgebraicDifferentiator__delayDisc");
        end
    end
end

