clear

ncond=5;
ncondobs=50000;
nobs=ncond*ncondobs;

% Create variables of n=50000 with the same mean and variance as the Loftus
% example in Abdi et al. 2010.

Contact=normrnd(30,10.7909427,50000);

Hit=normrnd(35,9.297550454,50000);

Bump=normrnd(38,11.05541597,50000);

Collide=normrnd(41,6.446359869,50000);

Smash=normrnd(46,5.773502692,50000);

% Now concatenate the data
R=[Contact; Hit; Bump; Collide; Smash];

% create a vector with the condition numbers (1=Contact; 2=Hit; 3=Bump;
% 4=Collide; 5=Smash
r=ones(50000,1);
Group=[];
for a=1:5;
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
     for j=1:50000     
          m=m+1;
          nomobs{m}=nomvar{i};
     end
end



% Now create a variable for reaction time
RT=normrnd(300,140,250000);

% Now create a variable for age (need some variability in the data.
Age=nan(250000,1);
Age(1:50000,1)=normrnd(42.25,16,50000);
Age(50001:100000,1)=normrnd(42,16,50000);
Age(100001:150000,1)=normrnd(41.75,17,50000);
Age(150001:200000,1)=normrnd(41.5,15,50000);
Age(200001:250000,1)=normrnd(41.25,18,50000);
Age=round(Age);

% Now create a variable for IQ
IQ=nan(250000,1);
IQ(1:50000,1)=normrnd(100.5,16,50000);
IQ(50001:100000,1)=normrnd(100,15,50000);
IQ(100001:150000,1)=normrnd(99.8,14,50000);
IQ(150001:200000,1)=normrnd(99,16,50000);
IQ(200001:250000,1)=normrnd(100.5,14,50000);
IQ=round(IQ);

% Create a vector of participant numbers
Px=[1:250000]';



% Concatenate the data so far
R2=[Px, Group, Age, IQ, R, RT];
dlmwrite('loftus.csv',R2)



