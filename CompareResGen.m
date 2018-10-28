clc
clear all

%3 components model sampling results
load('A20-K7-60-all_22000it_3comp.mat');
theta_alliter_3comp=theta_alliter;      %converge at beginning
%5 components model sampling results
load('A20-K8-36-all_12000it.mat');
theta_alliter_5comp=theta_alliter;      %converge at 2000 iteration

%check bode plot/nyquist plot point, check visualized fitting
thetatmp1=theta_alliter_3comp(end,1:end);
for i=1:length(thetatmp1);
    [tmph tmpx]=hist(theta_alliter_3comp(500:end,i));
    [maxfreq maxval]=max(tmph);
    thetatmp1(i)=tmpx(maxval);
end
thetatmp2=theta_alliter_5comp(end,1:end);
for i=1:length(thetatmp2)
    [tmph tmpx]=hist(theta_alliter_5comp(2000:end,i));
    [maxfreq maxval]=max(tmph);
    thetatmp2(i)=tmpx(maxval);
end

%nyquist plot for 3 components and 5 model estimation
wfreq=table2array(dattab(1:36,3))'*2*pi;
Gtmp1=Cal_transfer_func_3comp(thetatmp1);
% figure(1)
% nyquist(Gtmp1,wfreq);
% hold on
% plot(table2array(dattab(1:end,1)),table2array(dattab(1:end,2)),'s');
 Gtmp2=Calc_transfer_func(thetatmp2);
% figure(2)
% nyquist(Gtmp2,wfreq);
% hold on
% plot(table2array(dattab(1:end,1)),table2array(dattab(1:end,2)),'s');

%nyquist plot together: 3 components, 5 components, real data
resp1=freqresp(Gtmp1,wfreq);
resp1_re=squeeze(real(resp1)');
resp1_img=squeeze(abs(imag(resp1)'));
resp2=freqresp(Gtmp2,wfreq);
resp2_re=squeeze(real(resp2)');
resp2_img=squeeze(abs(imag(resp2)'));

figure(3)
realplot=plot(table2array(dattab(1:end,1)),table2array(dattab(1:end,2)),'r-s');
hold on
plot1=plot(resp1_re,resp1_img,'b-o');       %3 components
hold on
plot2=plot(resp2_re,resp2_img,'k-*','LineWidth',2);   %5 components

%nyquist plot of initial guess 
Htmp=Calc_transfer_func(thetatmp111);
wfreqtmp=table2array(dattab(1:36,3))'*2*pi;
resptmp=freqresp(Htmp,wfreqtmp);
resptmpre1=squeeze(real(resptmp)');
resptmpimg1=squeeze(abs(imag(resptmp)'));
refplot=plot(resptmpre1,resptmpimg1,'b--o');
resptmp_mean=[resptmpre1;resptmpimg1];
resptmp_diff=resptmp_mean-[table2array(dattab(1:36,1));table2array(dattab(1:36,2))];
rss=dot(resptmp_diff,resptmp_diff);

thetatmplower2=theta_alliter_5comp(end,1:end);
thetatmpupper2=theta_alliter_5comp(end,1:end);
for i=1:length(thetatmplower2)
    [tmpval tmpind]=min(theta_alliter_5comp(2000:end,i));
    thetatmplower2(i)=tmpval;
    [tmpval tmpind]=max(theta_alliter_5comp(2000:end,i));
    thetatmpupper2(i)=tmpval;    
end
%generate all parameter combinations
%tmpset1={[thetatmplower1(1) thetatmpupper1(1)],[thetatmplower1(2) thetatmpupper1(2)],...
%    [thetatmplower1(3) thetatmpupper1(3)],[thetatmplower1(4) thetatmpupper1(4)]};
%[tmp1 tmp2 tmp3 tmp4]=ndgrid(tmpset1{:});
%combparam1=[tmp1(:) tmp2(:) tmp3(:) tmp4(:)];
tmpset2={[thetatmplower2(1) thetatmpupper2(1)],[thetatmplower2(2) thetatmpupper2(2)],...
    [thetatmplower2(3) thetatmpupper2(3)],[thetatmplower2(4) thetatmpupper2(4)],[thetatmplower2(5) thetatmpupper2(5)],...
    [thetatmplower2(6) thetatmpupper2(6)],[thetatmplower2(7) thetatmpupper2(7)]};
[tmp1 tmp2 tmp3 tmp4 tmp5 tmp6 tmp7]=ndgrid(tmpset2{:});
combparam2=[tmp1(:) tmp2(:) tmp3(:) tmp4(:) tmp5(:) tmp6(:) tmp7(:)];

candlen=length(combparam2(:,1));
respre5comp=zeros(candlen,36);
respimg5comp=zeros(candlen,36);
for i=1:candlen
    G5tmp=Calc_transfer_func(combparam2(i,:));
    resp5comp=freqresp(G5tmp,wfreq);
    respre5comp(i,:)=squeeze(real(resp5comp));
    respimg5comp(i,:)=squeeze(abs(imag(resp5comp)));
end
final_lower_re2=zeros(36,1);
final_lower_img2=zeros(36,1);
final_upper_re2=zeros(36,1);
final_upper_img2=zeros(36,1);
for i=1:25
    if i==1
        respref_img=respimg5comp(:,i);
        respref_re=respre5comp(:,i);
    else
        tmpidx=find(respre5comp(:,i)>=final_lower_re2(i-1));
        respref_re=respre5comp(tmpidx,i);
        respref_img=respimg5comp(tmpidx,i);
    end
    
    [tmpval tmpind]=max(respref_re);
    if length(tmpind) > 1
        final_lower_re2(i)=tmpval(1);
        tmpv=respref_img(tmpind);
        [tmpval2 tmpind2]=min(tmpv);
        final_lower_img2(i)=tmpval2;
    else
        final_lower_img2(i)=respref_img(tmpind);
        final_lower_re2(i)=respref_re(tmpind);
    end
    
    if i==1
        respref_img=respimg5comp(:,i);
        respref_re=respre5comp(:,i);
    else
        tmpidx=find(respre5comp(:,i)>=final_upper_re2(i-1));
        respref_re=respre5comp(tmpidx,i);
        respref_img=respimg5comp(tmpidx,i);
    end
    [tmpval tmpind]=min(respref_re);
    if length(tmpind) > 1
        final_upper_re2(i)=tmpval(1);
        tmpv=respref_img(tmpind);
        [tmpval2 tmpind2]=max(tmpv);
        final_upper_img2(i)=tmpval2;
    else
        final_upper_img2(i)=respref_img(tmpind);
        final_upper_re2(i)=respref_re(tmpind);
    end
end

for i=26:36
    if i==1
        respref_img=respimg5comp(:,i);
        respref_re=respre5comp(:,i);
    else
        tmpidx=find(respre5comp(:,i)>=final_lower_re2(i-1));
        respref_re=respre5comp(tmpidx,i);
        respref_img=respimg5comp(tmpidx,i);
    end
    
    [tmpval tmpind]=min(respref_img);
    if length(tmpind) > 1
        final_lower_img2(i)=tmpval(1);
        tmpv=respref_re(tmpind);
        [tmpval2 tmpind2]=max(tmpv);
        final_lower_re2(i)=tmpval2;
    else
        final_lower_img2(i)=respref_img(tmpind);
        final_lower_re2(i)=respref_re(tmpind);
    end
    
    if i==1
        respref_img=respimg5comp(:,i);
        respref_re=respre5comp(:,i);
    else
        tmpidx=find(respre5comp(:,i)>=final_upper_re2(i-1));
        respref_re=respre5comp(tmpidx,i);
        respref_img=respimg5comp(tmpidx,i);
    end
    [tmpval tmpind]=max(respref_img);
    if length(tmpind) > 1
        final_upper_img2(i)=tmpval(1);
        tmpv=respref_re(tmpind);
        [tmpval2 tmpind2]=min(tmpv);
        final_upper_re2(i)=tmpval2;
    else
        final_upper_img2(i)=respref_img(tmpind);
        final_upper_re2(i)=respref_re(tmpind);
    end
end

%check plot
%plot(table2array(dattab(1:end,1)),table2array(dattab(1:end,2)),'r-s');
%hold on
lbplot=plot(final_lower_re2,final_lower_img2,'k--');
hold on
ubplot=plot(final_upper_re2,final_upper_img2,'k-.');


legend([realplot plot1 plot2 refplot lbplot ubplot],'Real Response','Simplified Randles Cell Model','Failed Coating Equivalent Circuit Model','Gamry Software','Lower Bound','Upper Bound');

load('A20-K8-36-all_22000it_const.mat');
%check nyquist plot
thetatmp=theta_alliter(end,1:end);
for i=1:length(thetatmp)
    [hval xval]=hist(theta_alliter(5000:end,i));
    [tmpval tmpind]=max(hval);
    thetatmp(i)=xval(tmpind);
end
Htmp=Calc_transfer_func(thetatmp);
resptmp=freqresp(Htmp,table2array(dattab(1:end,3))'*2*pi);
resp_const_re=squeeze(real(resptmp)');
resp_const_img=squeeze(abs(imag(resptmp)'));
plot(table2array(dattab(1:end,1)),table2array(dattab(1:end,2)),'r-o');
hold on
plot(resp_const_re,resp_const_img,'b-s');

%save for R plot
compares=horzcat(table2array(dattab(1:end,1:2)),resp1_re,resp1_img,resp2_re,resp2_img,resptmpre1,...
    resptmpimg1,final_lower_re2,final_lower_img2,final_upper_re2,final_upper_img2,resp_const_re,resp_const_img);

save('comparison_plot.txt','compares','-ascii');