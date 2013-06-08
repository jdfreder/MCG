experiment = importdata('radioactivedecay.dat');
t = experiment.data(:,1);
N = experiment.data(:,2);
figure(2)
plot(t,N,'.b');
hold on;
%%
%function to fit
fun=@(a,b,t) a*exp(-b*t);
%Find a starting point for the parameters a, b, and c.
guess = fun(12,.35,t); %fun(a,b,t)
plot(t,guess,'r-')
%fit the data
fittedmodel=fit(t',N',fun,'StartPoint',[12 .35])
%plot the result
plot(fittedmodel,'r-');
legend('fit')