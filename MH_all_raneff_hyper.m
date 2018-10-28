clc
clear all

%Trials for sample specific R2 with hyper parameter u and tau^2, R2i ~
%N(u,tau^2)

%different measurements number among samples
dirvar='xxx';

fillist=dir(strcat(dirvar,'*.txt'));
dattab=[];
%A20: 7 samples with same time point 96 hours
samp_num=length(fillist);
measinsamp=zeros(samp_num,1);
for i=1:samp_num
    tmptab=readtable(strcat(dirvar,fillist(i).name),'Delimiter','\t','ReadVariableNames',false);
    measinsamp(i)=length(table2array(tmptab(1:end,1)));
    dattab=vertcat(dattab,tmptab);
end

max_iter=22000;

%circuit model shared parameters: R1,R3,C1,C2,alpha1,alpha2
%population shared parameters: sigma^2,u,tau^2
%individual parameter: R2
r1_alliter=zeros(max_iter,1);
r1_alliter(1,1)=300;

r2_alliter=zeros(max_iter,samp_num);
r2_alliter(1,1)=2e+4;
r2_alliter(1,2)=4e+4;
r2_alliter(1,3)=2e+4;
r2_alliter(1,4)=3e+4;

r3_alliter=zeros(max_iter,1);
r3_alliter(1,1)= 1e+6;

c1_alliter=zeros(max_iter,1);
c1_alliter(1,1)=1e-5;  

c2_alliter=zeros(max_iter,1);
c2_alliter(1,1)=0.0005;

alpha1_alliter=zeros(max_iter,1);
alpha1_alliter(1,1)=0.8;       

alpha2_alliter=zeros(max_iter,1);
alpha2_alliter(1,1)=0.4;   

%conjugate prior for sigmasqr, u, tau
mu_alliter=zeros(max_iter,1);
mu_alliter(1,1)=3e+4;

tausqr_alliter=zeros(max_iter,1);
tausqr_alliter(1,1)=1e4;

sigmasqr_alliter=zeros(max_iter,1);
sigmasqr_alliter(1,1)=1;

prior_mean=0;
prior_sigmasqr=1e20;

step_sigma1=1;
%step_sigma2(1:7)=[100 10 10 100 100 100 100];
step_sigma2=100;
step_sigma3=100;
step_sigma4=1e-6;
step_sigma5=1e-6;%1e-5

prior_a=0.01;    %non-informative IG prior param for sigma^2 
prior_b=0.01;

prior_a1=0.01;  %non-informative IG prior param for tau^2
prior_b1=0.01;

%prior_uR=3e+4;   
%prior_sigsqrR=1e+4;
prior_uR=0;
prior_sigsqrR=1e+20;

theta_alliter=horzcat(r1_alliter,r2_alliter,r3_alliter,c1_alliter,c2_alliter,...
    alpha1_alliter,alpha2_alliter,mu_alliter,tausqr_alliter,sigmasqr_alliter);

startind=6;
%In each iteration, Gibbs sampling for sigmasqr,mu,tausqr and MH for others
for i=2:max_iter
    thetavec=theta_alliter(i-1,1:end);
    
    %draw individual R2
    for j=1:samp_num
        thetavectmp=thetavec([1 j+1 startind:end]);
        r2_cur=normrnd(thetavectmp(end-2),thetavectmp(end-1));%MH_sampling_impl_rand(thetavectmp,dattab(1+sum(measinsamp(1:j-1)):measinsamp(j)+sum(measinsamp(1:j-1)),1:end),2,step_sigma2,[thetavec(end-2),thetavec(end-1)],0,samp_num,startind,measinsamp);
        thetavec(j+1)=r2_cur;
    end
    %draw shared R3
    r3_cur=MH_sampling_impl_rand(thetavec,dattab,startind,step_sigma3,[prior_mean,prior_sigmasqr],1,samp_num,startind,measinsamp);
    thetavec(startind)=r3_cur;
    %draw shared C1
    c1_cur=MH_sampling_impl_rand(thetavec,dattab,startind+1,step_sigma4,[prior_mean,prior_sigmasqr],1,samp_num,startind,measinsamp);
    thetavec(startind+1)=c1_cur;
    %draw shared C2
    c2_cur=MH_sampling_impl_rand(thetavec,dattab,startind+2,step_sigma5,[prior_mean,prior_sigmasqr],1,samp_num,startind,measinsamp);
    thetavec(startind+2)=c2_cur;
    %draw shared R1
    r1_cur=MH_sampling_impl_rand(thetavec,dattab,1,step_sigma1,[prior_mean,prior_sigmasqr],1,samp_num,startind,measinsamp);
    thetavec(1)=r1_cur;    
    %draw shared alpha1
    alpha1_cur=MH_sampling_impl_rand(thetavec,dattab,startind+3,-1,[0;1],1,samp_num,startind,measinsamp);
    thetavec(startind+3)=alpha1_cur;
    %draw shared alpha2
    alpha2_cur=MH_sampling_impl_rand(thetavec,dattab,startind+4,-1,[0;1],1,samp_num,startind,measinsamp);
    thetavec(startind+4)=alpha2_cur;

    %calculate new rss and draw shared sigma square
    diffvec=[];
    for j=1:samp_num
        thetatmp=thetavec([1 j+1 startind:end]);
        H2=Calc_transfer_func(thetatmp);
        resp2=freqresp(H2,table2array(dattab(1+sum(measinsamp(1:j-1)):measinsamp(j)+sum(measinsamp(1:j-1)),3))'*2*pi);
        resp_new_re=squeeze(real(resp2)');
        resp_new_img=squeeze(abs(imag(resp2)'));
        resp_new_mean=[resp_new_re;resp_new_img];
        tmpdiff=[table2array(dattab(1+sum(measinsamp(1:j-1)):measinsamp(j)+sum(measinsamp(1:j-1)),1));table2array(dattab(1+sum(measinsamp(1:j-1)):measinsamp(j)+sum(measinsamp(1:j-1)),2))]-resp_new_mean;
        diffvec=vertcat(diffvec,tmpdiff);
    end
    
    %draw sample which has conjugate prior    
    %draw sample u
    post_uR=(sum(thetavec(2:samp_num+1))*prior_sigsqrR+prior_uR*thetavec(end-1))/(samp_num*prior_sigsqrR+thetavec(end-1));
    post_sigsqrR=thetavec(end-1)*prior_sigsqrR/(samp_num*prior_sigsqrR+thetavec(end-1));
    pd1=makedist('Normal','mu',post_uR,'sigma',sqrt(post_sigsqrR));
    distt=truncate(pd1,1e4,Inf);
    mu_cur=random(distt,1,1);
    thetavec(end-2)=mu_cur;
    
    %draw sample tau square
    post_a1=prior_a1+samp_num/2;
    post_b1=prior_b1+1/2*sum((thetavec(2:samp_num+1)-thetavec(end-2)).^2);
    tausqr_cur=1/gamrnd(post_a1,1/post_b1);
    thetavec(end-1)=tausqr_cur;
    
    %draw sample sigma square
    post_a=prior_a+sum(measinsamp(1:end))/2;
    post_b=prior_b+1/2*dot(diffvec,diffvec);
    sigmasqr_cur=1/gamrnd(post_a,1/post_b);
    thetavec(end)=sigmasqr_cur;
    
    theta_alliter(i,1:end)=thetavec;
end

save('A20-K7-60-all_22000it_hyper.mat');
save('A20-K7-60-all_22000it_hyper.txt','-ascii');
save('A20-K7-60-all_theta_22000it_hyper.txt','theta_alliter','-ascii');
xlswrite('A20-K7-60-all_hyper.xls',theta_alliter);

%derive posterior mode for sample specific parameters and shared parameters
thetatmp=theta_alliter(end,1:end);
for i=1:length(thetatmp);
    [tmph tmpx]=hist(theta_alliter(2000:end,i));
    [maxfreq maxval]=max(tmph);
    thetatmp(i)=tmpx(maxval);
end

%fitting of 4 samples with individual Rcoating
resps1=[];
respother=[];
titarr=['Nyquist plot of Sample 1';'Nyquist plot of Sample 2';'Nyquist plot of Sample 3';'Nyquist plot of Sample 4'];
for j=1:samp_num
wfreq=table2array(dattab(1+sum(measinsamp(1:j-1)):measinsamp(j)+sum(measinsamp(1:j-1)),3))'*2*pi;
Gtmp=Calc_transfer_func(thetatmp([1 j+1 6:end]));
resptmp=freqresp(Gtmp,wfreq);
realtmp=squeeze(real(resptmp)');
imagtmp=squeeze(abs(imag(resptmp)'));
if j==1
    resps1=horzcat(resps1,realtmp,imagtmp);
else
   respother=horzcat(respother,realtmp,imagtmp);
end
end

save('A20_22000it_hyper_resp_S1.txt','resps1','-ascii');
save('A20_22000it_hyper_resp_other.txt','respother','-ascii');

%%%%%%save template
%save('step100.mat')
%save('step100.txt','-ascii')
%save('step100_theta.txt','theta_alliter','-ascii')
%xlswrite('step100_theta.xls',theta_alliter)