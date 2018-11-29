function [output1,output2] = eg_func(input1,input2)
% This is a test funciton for miranda
%>>Input>>
%input 1 - A number between 1 and 10
%input 2 - A number between 50 and 100
%<<Ouput<<
%output1 - the product of all the numbers up to input 1
%output2 - the sum of all the numbers up to input 2

vec1 = 1;
for ii = 1:input1
    vec1 = vec1*ii;
end
output1 = vec1;

vec2 = [1:input2];
output2 = sum(vec2);