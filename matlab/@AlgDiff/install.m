function success = install()
%INSTALL Tries to install AlgDiff python library if not already installed.
%   Detailed explanation goes here

% Try to simply instantiate python class
if try_instantiate()
    success = true;
end

% Check for python 
info = pyenv;
if isempty(info.Version)
    %% Python is not installed
    % Query supported python versions
    data = webread('https://www.mathworks.com/support/requirements/python-compatibility.html');
    [startIndex, endIndex] = regexp(data, '<table.*?>.*?</table>');
    tbl = data(startIndex:endIndex);
    
    % Find supported versions (3.xy)
    v = version('-release');
    [startIndex, endIndex] = regexp(tbl, ['<tr.*?>(?!</tr>).*?R' v '.*?</tr>']);
    matches = regexp(tbl(startIndex:endIndex), '3\.\d{1,2}', 'match');
    py_versions = containers.Map('KeyType', 'char', 'ValueType', 'double');
    for idx=1:length(matches)
        if py_versions.isKey(matches{idx})
            py_versions(matches{idx}) = py_versions(matches{idx}) + 1;
        else
            py_versions(matches{idx}) = 1;
        end
    end

    

end

end

function versions = get_supported_python_versions()
    % Query supported python versions
    data = webread('https://www.mathworks.com/support/requirements/python-compatibility.html');
    [startIndex, endIndex] = regexp(data, '<table.*?>.*?</table>');
    tbl = data(startIndex:endIndex);
    
    % Find supported versions (3.xy)
    v = version('-release');
    [startIndex, endIndex] = regexp(tbl, ['<tr.*?>(?!</tr>).*?R' v '.*?</tr>']);
    matches = regexp(tbl(startIndex:endIndex), '3\.\d{1,2}', 'match');
    py_versions = containers.Map('KeyType', 'char', 'ValueType', 'double');
    for idx=1:length(matches)
        if py_versions.isKey(matches{idx})
            py_versions(matches{idx}) = py_versions(matches{idx}) + 1;
        else
            py_versions(matches{idx}) = 1;
        end
    end

    % None found (?)
    if py_versions.Count == 0
        versions = [];
    end

    % Sort by the following rules
    %   1. versions where count == 3    (then by decreasing version)
    %   2. versions where count == 2    (then by decreasing version)
    %   3. Decreasing version
    function y = sort_by_version(x)
        
    end

end

function success = try_instantiate()
    % Instantiate python class
    try
        py.AlgDiff.AlgebraicDifferentiator();
        success = true;

    catch ME
        success = false;
    end
end