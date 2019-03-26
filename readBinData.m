% READBINDATA - read individual channel RF binary data from antares
%
% X = READBINDATA(FILENAME) reads the binary file FILENAME and outputs
% the image data into X.  The image data in X is sorted first by sample,
% then by element number, and lastly by vector.  
%
% [X,NUMVEC,NUMELEM,NUMSAMP] = READBINDATA(FILENAME) will output the
% number of vectors, elements, and samples.
%
% Ex:
%   [x,numVectors,numElements,numSamples] = readBinData('imageData.bin');
%
% by Jeremy Dahl 11/06/02

function [x,numVectors,numElements,numSamples] = readBinData(filename)

% Read header data from binary file.
fid = fopen(filename);
header = fread(fid,7,'uint32');
fclose(fid);

numVectors = header(1);
numElements = header(2);
numSamples = header(3);

% Get data from binary file and sort it into a 3-D matrix.
fid = fopen(filename);
x = fread(fid,inf,'uint16');
fclose(fid);
x = reshape(x(15:end),[numSamples numElements numVectors]);

x = x-512;

return
