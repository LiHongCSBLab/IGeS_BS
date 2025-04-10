function p = mortarpath
% MORTARPATH Get location of Mortar library

if exist('MORTARPATH','var')
	% disp('using existing MORTARPATH');
    p = MORTARPATH;
else
	% disp('NOT using existing MORTARPATH');
    p = strrep(which(mfilename), sprintf('%slib%s%s.m',filesep, filesep, mfilename), '');    
end
