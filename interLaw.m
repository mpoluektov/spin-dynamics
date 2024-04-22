function [ J ] = interLaw( x, decCoef )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

J = sin( pi*x - pi/2 ).*exp( -decCoef*(x-1) ).*( x.^(-3) );

end

