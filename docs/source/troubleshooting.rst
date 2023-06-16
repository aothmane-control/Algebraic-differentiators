Troubleshooting
===============

A list of known issues and fixes is presented here.

- **Python Error: ValueError: illegal value in argument 11 of internal sbevd**


	**Error message**

		When calling the algebraic differentiators in Matlab the error

		>>> Python Error: ValueError: illegal value in argument 11 of internal sbevd

		can appread. In the Terminal where matlab has been called the following error message can be found: 

		>>> Intel MKL ERROR: Parameter 11 was incorrect on entry to DSBEVD.

	**Fix**

		This error is because the Python code is incompatible with MATLAB's MKL.

		This can be avoided by calling

		>>> py.sys.setdlopenflags(int32(bitor(int64(py.os.RTLD_NOW), int64(py.os.RTLD_DEEPBIND))));

		or 	

		>>> py.sys.setdlopenflags(int32(bitor(int64(py.os.RTLD_LAZY), int64(py.os.RTLD_DEEPBIND))));

		directly before any calls to py and after calls to pyenv. See this `support answer <https://de.mathworks.com/matlabcentral/answers/358233-matlab-python-interface-broken?s_tid=email_ans_new_ans_ans_h>`_. The simples fix is to add the lines to the startup script of matlab.
		This fixs has first been published on the `archlinux wiki <https://wiki.archlinux.org/title/MATLAB>`_ in section 4.25.