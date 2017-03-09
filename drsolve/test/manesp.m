function [man, esp] = manesp(x)
%MANESP	get mantissa and exponent of a number.
%
%   Called by test functions, to format LaTeX tables.
%   Syntax: [man,eps] = manesp(x).

%   Antonio Arico' & Giuseppe Rodriguez, University of Cagliari, Italy
%   Email: {arico,rodriguez}@unica.it
%
%   Last revised Feb 26, 2010

esp = floor(log10(x));
man = x / 10.^esp;

