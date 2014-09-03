%% Application code for GETWORDS.
addpath('Literature')
[W, T] = getwords('sherlock_holmes.txt');
[W2, T2] = getwords('don_quixote.txt');
rmpath('Literature')
