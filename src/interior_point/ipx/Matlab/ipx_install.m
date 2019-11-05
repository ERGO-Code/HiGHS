%IPX_INSTALL compiles BASICLU and IPX for MATLAB.
%
% The current working directory must be ipx/Matlab. If the BASICLU root
% directory is located side by side with the IPX root directory, then the
% function can be called without argument. Otherwise the path to the BASICLU
% root directory must be given as argument.
%
% Example:
%   ipx_install();  % assumes that ../../basiclu is the BASICLU root dir
%   ipx_install('/path/to/basiclu');
function ipx_install(basicluroot)
    if nargin == 0
        basicluroot = fullfile('..', '..', 'basiclu');
    end
    if exist(basicluroot, 'dir') ~= 7 || ...
       exist(fullfile(basicluroot, 'include', 'basiclu.h'), 'file') ~= 2
        msg = sprintf(['Directory %s does not exist or is not the BASICLU ' ...
                       'root directory.'], basicluroot);
        error(msg)
    end

    % Compile each BASICLU source file into its object file and store them in
    % build/basiclu.
    fprintf('Compiling BASICLU for Matlab ...\n');
    basiclu_include_flag = ['-I' fullfile(basicluroot, 'include')];
    basiclu_source_files = fullfile(basicluroot, 'src', '*.c');
    basiclu_build_dir = fullfile('build', 'basiclu');
    mex('-largeArrayDims', '-c', '-outdir', basiclu_build_dir, ...
        basiclu_include_flag,  basiclu_source_files);

    basiclu_obj_files = dir(fullfile(basiclu_build_dir, '*.o'));
    basiclu_obj_files = { basiclu_obj_files(:).name };
    basiclu_obj_files = fullfile('build', 'basiclu', basiclu_obj_files);

    % Compile each IPX source file into its object file and store them in
    % build/ipx.
    fprintf('Compiling IPX for Matlab ...\n');
    ipx_include_flag = ['-I' fullfile('..', 'include')];
    ipx_source_files = fullfile('..', 'src', '*.cc');
    ipx_build_dir = fullfile('build', 'ipx');
    mex('-largeArrayDims', '-c', '-outdir', ipx_build_dir, ...
        ipx_include_flag, basiclu_include_flag, ipx_source_files);

    ipx_obj_files = dir(fullfile(ipx_build_dir, '*.o'));
    ipx_obj_files = { ipx_obj_files(:).name };
    ipx_obj_files = fullfile('build', 'ipx', ipx_obj_files);

    % Compile each C source file from src/ into a Matlab executable and link
    % with IPX and BASICLU.
    fprintf('Compiling MEX interface ...\n');
    mex_source_files = dir(fullfile('src', '*.c'));
    for k = 1:length(mex_source_files)
        mex_source_file = fullfile('src', mex_source_files(k).name);
        mex('-largeArrayDims', ipx_include_flag, ...
            mex_source_file, ipx_obj_files{:}, basiclu_obj_files{:});
    end
end
