%   Extended phase graph for CPMG sequences
%   
%   Usage: S = epgMEX(T1,T2,esp,alpha)
%   Author: RM Lebel
%   Date: 06/2011
%   
%   Input:
%   T1: T1 relaxation time (s)
%   T2: T2 relaxation time (s)
%   esp: echo spacing (s)
%   alpha: array of flip angles during the echo train (rad)
%       This must be of size [ETL x nPoints] where ETL is the echo train
%       length and nPoints is the number of points across the slice
%
%   Output:
%   S: Relative echo amplitude. Size: [ETL x nPoints]
