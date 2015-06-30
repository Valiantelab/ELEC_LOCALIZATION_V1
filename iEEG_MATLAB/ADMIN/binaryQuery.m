function resp=binary_query()
%function resp=binary_query()
%
% Asks user to "Answer y or n:" and continues to ask until a 'y' or 'n'
% is input.  Function is NOT case sensitive
%
% D. Groppe

resp=[];
while ~strcmpi(resp,'y') && ~strcmpi(resp,'n')
    resp=input('Answer y or n: ','s');
end