clear all

% positions: [w_n,l_n,w_c,l_c,h,a_y]
% For unit test to pass the structure must be
% w_n = [1,7,13,19,25]
% l_n = [2,8,14,20,26]
% w_c = [3,9,15,21,27]
% l_c = [4,10,16,22,28]
% h = [5,11,17,23,29]
% a_y = [6,12,18,24,30]
%
% function passes unit test

n = 5;
testgeo = 1:1:30;
unpacked = unpackgeometry(testgeo,n);