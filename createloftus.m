clear

ncond=5; % number of conditions
ncondobs=50000; % number of observations in each condition (n)
nobs=ncond*ncondobs; % total number of observations (N)


%%%%%%%%
% Create variables of n=50000 with the same mean and variance as the Loftus
% example in Abdi et al. 2010.

Contact=normrnd(30,10.7909427,ncondobs);
Hit=normrnd(35,9.297550454,ncondobs);
Bump=normrnd(38,11.05541597,ncondobs);
Collide=normrnd(41,6.446359869,ncondobs);
Smash=normrnd(46,5.773502692,ncondobs);

% Now concatenate the data
R=[Contact; Hit; Bump; Collide; Smash];

% create a vector with the condition numbers (1=Contact; 2=Hit; 3=Bump;
% 4=Collide; 5=Smash
r=ones(ncondobs,1);
Group=[];
for a=1:ncond;
    Group=[Group; r*a];
end
clear r

% Create a vector of condition names
n=0;
n=n+1; nomvar{n}='Contact';
n=n+1; nomvar{n}='Hit';
n=n+1; nomvar{n}='Bump';
n=n+1; nomvar{n}='Collide';
n=n+1; nomvar{n}='Smash';

% Create vector of condition names for each observation
m=0;
for i= 1:n
     for j=1:ncondobs     
          m=m+1;
          nomobs{m}=nomvar{i};
     end
end



% Now create a variable for reaction time
RT=normrnd(300,140,nobs);

% Now create a variable for age (need some variability in the data.
Age=nan(nobs,1);
Age(1:50000,1)=normrnd(42.25,16,ncondobs);
Age(50001:100000,1)=normrnd(42,16,ncondobs);
Age(100001:150000,1)=normrnd(41.75,17,ncondobs);
Age(150001:200000,1)=normrnd(41.5,15,ncondobs);
Age(200001:250000,1)=normrnd(41.25,18,ncondobs);
Age=round(Age);

% Now create a variable for IQ
IQ=nan(nobs,1);
IQ(1:50000,1)=normrnd(100.5,16,ncondobs);
IQ(50001:100000,1)=normrnd(100,15,ncondobs);
IQ(100001:150000,1)=normrnd(99.8,14,ncondobs);
IQ(150001:200000,1)=normrnd(99,16,ncondobs);
IQ(200001:250000,1)=normrnd(100.5,14,ncondobs);
IQ=round(IQ);

% Create a vector of participant numbers
Px=[1:nobs]';



% Concatenate the data so far

 toto =    num2cell(Px);
 toto1 = nomobs';
 toto2 = num2cell([Group, Age, IQ, R, RT]);
 
 R2=[toto toto1 toto2];
 clear toto toto1 toto2

k=0;
k=k+1; colnames{k}='Participant';
k=k+1; colnames{k}='Condition Name';
k=k+1; colnames{k}='Condition Number';
k=k+1; colnames{k}='Age';
k=k+1; colnames{k}='IQ';
k=k+1; colnames{k}='Estimated Speed (mph)';
k=k+1; colnames{k}='Reaction Time (ms)';

R2=[colnames;R2];

% write line-by-line to csv file
fid = fopen('loftus.csv','wt');
for i=1:size(R2,1)
    fprintf(fid, '%s,%d,%d\n', R2{i,:});
end
fclose(fid);

 




