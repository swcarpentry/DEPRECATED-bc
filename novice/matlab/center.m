
function out = center(data, desired)
    %   Center data around a desired value.
    %
    %       center(DATA, DESIRED) 
    %
    %   Returns a new array containing the values in
    %   DATA centered around the value.

    out = (data  - mean(data)) + desired;
