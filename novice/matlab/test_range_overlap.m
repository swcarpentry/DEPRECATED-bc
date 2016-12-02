%script test_range_overlap.m


% assert(range_overlap([0.0, 1.0], [5.0, 6.0]) == NaN);
% assert(range_overlap([0.0, 1.0], [1.0, 2.0]) == NaN);
assert(range_overlap([0, 1.0]) == [0, 1.0]);
assert(range_overlap([2.0, 3.0], [2.0, 4.0]) == [2.0, 3.0]);
assert(range_overlap([0.0, 1.0], [0.0, 2.0], [-1.0, 1.0]) == [0.0, 1.0]);

