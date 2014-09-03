function inPlaceExample()

% Create an array occupying about 2% of available memory.
m = memory;
mSize = round(sqrt(0.02*m.MemAvailableAllArrays/8));
A = rand(mSize);

disp('Memory available after creating large array (GB):')
m = memory;
disp(m.MemAvailableAllArrays/1073741824)

A = iNotInPlace(A);
disp('Memory available after calling the not in-place function (GB):')
m = memory;
disp(m.MemAvailableAllArrays/1073741824)

A = iInPlace(A); %#ok<NASGU>
disp('Memory available after calling the in-place function (GB):')
m = memory;
disp(m.MemAvailableAllArrays/1073741824)

end

function A = iInPlace(A)

A = 2*A + 1;
disp('Memory available inside the in-place function (GB):')
m = memory;
disp(m.MaxPossibleArrayBytes/1073741824)


end

function B = iNotInPlace(A)

B = 2*A + 1;
disp('Memory available inside the not in-place function (GB):')
m = memory;
disp(m.MaxPossibleArrayBytes/1073741824)

end

