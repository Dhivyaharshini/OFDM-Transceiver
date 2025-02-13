function indices = helperLTEIntrlvrIndices(blkLen)
%HELPERLTEINTRLVRINDICES Generate turbo code internal interleaver indices. 
%
%   INDICES = helperLTEIntrlvrIndices(BLKLEN) returns the turbo code
%   internal interleaver indices for a specified block length, BLKLEN, as
%   per 3GPP TS 36.212 v9.0.0, Section 5.1.3.2.3.
%
%   The format of the returned INDICES column vector is such that it can be
%   directly used as the parameter value in the general block interleaver.
%
%   See also getf1f2.

%   Reference:
%   3GPP TS 36.212 v9.0.0, "3rd Generation partnershiop project;
%   Technical specification group radio access network; Evolved Universal
%   Terrestrial Radio Acess (E-UTRA); Multiplexing and channel coding
%   (release 9)", 2009-12.

%   Copyright 2011-2020 The MathWorks, Inc.

%#codegen

[f1, f2] = getf1f2(blkLen);
Idx      = (0:blkLen-1).';
indices  =  mod(f1*Idx + f2*Idx.^2, blkLen) + 1;

% [EOF]
