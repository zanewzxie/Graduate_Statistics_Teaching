


% ============== 1. Data tranform.
%genderate a skew population, such as the gamma distribution
randomsample =gamrnd(3,4,10000,1);
figure(1);subplot(1,3,1);
hist(randomsample,30); %plot it to check the distribution
axis([0,40,0,1500])
title('original skewed random sample');

transform1 = log(randomsample);
figure(1);subplot(1,3,2);
hist(transform1,30); %plot it to check the distribution
axis([0,5,0,1500])
title('log transformed sample');

transform2 = sqrt(randomsample);
figure(1);subplot(1,3,3);
hist(transform2,30); %plot it to check the distribution
axis([1,8,0,1500])
title('Squared root transformed sample');


%% ============== 2. Simple resampling example for correlation

%Data
gpa = [3.4,4.0,3.1,2.9,4.0,3.3,3.5,2.8,2.7,3.5,3.8,3.9,3.7,2.9];
gre = [300,320,310,290,300,320,325,310,280,330,310,320,330,280];

%scatter plot
figure(2);clf;subplot(1,2,1);
plot(gpa,gre,'ko','MarkerFaceColor','b');
xlabel('GPA');
ylabel('GER');

%do the resampling for one time.
%get the index for the pairs
nReps = 10000;
Npairs=length(gpa);
PairIndex = 1:Npairs;
resampledIndex = datasample(PairIndex,Npairs);% randomly sample the dat with replacement.
tmp = corrcoef(gpa(resampledIndex),gre(resampledIndex));
r = tmp(1,2);

%repeat the process for 10000 times
for i = 1:nReps
resampledIndex = datasample(PairIndex,Npairs);% randomly sample the dat with replacement.
tmp = corrcoef(gpa(resampledIndex),gre(resampledIndex));
r(i) = tmp(1,2);
end
    
% calculate 95% CI
CIrange = 95;  %alpha =.05 (two-tailed)
sortedr = sort(r);
CI(1) = sortedr(single(nReps*(1-CIrange/100)/2));
CI(2) = sortedr(nReps-single(nReps*(1-CIrange/100)/2));

% See the r distribution 
figure(2);subplot(1,2,2);hist(r,100);
xlabel('boostrapped correlation');
hold on;plot(mean(r)*[1,1],ylim,'g-','LineWidth',2.5);
hold on;plot(CI(1)*[1,1],ylim,'r-','LineWidth',2.5);
hold on;plot(CI(2)*[1,1],ylim,'r-','LineWidth',2.5);
title(sprintf('r = %5.2f,  %2.1f%% CI: [%5.2f, %5.2f]',mean(r),CIrange,CI(1),CI(2)));

