function [ MatOutput1 ] = thresholdstatistics( MatInput1, Column1, Threshold1 )
%THRESHOLDSTATISTICS takes a Matrix, column, and threshold, and produces an
%output matrix with that threshold applied
%   takes the TrackSummary cell array produced from 
%   trackstatistics and returns a cell array less the members below the
%   threshold (use negative column for less than and positive column for
%   greater than)

if Column1>0
    goodIndx1 = find(MatInput1(:,abs(Column1))>Threshold1);
else
    goodIndx1 = find(MatInput1(:,abs(Column1))<Threshold1);
end

MatOutput1 = MatInput1(goodIndx1,:);