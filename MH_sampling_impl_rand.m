function param_newiter=MH_sampling_impl_rand(theta_vec,dat,upd_ind,stepsiz,prior_param,shared_ind,sampnum,sh_startind,measnum)
old_diff=[];
if shared_ind==0
    %individual specific parameter sampling: RSS calculate on individual data
    Hold=Calc_transfer_func(theta_vec);
    resp1=freqresp(Hold,table2array(dat(1:end,3))'*2*pi);
    old_re=squeeze(real(resp1)');
    old_img=squeeze(abs(imag(resp1)'));
    old_mean=[old_re;old_img];
    old_diff=[table2array(dat(1:end,1));table2array(dat(1:end,2))]-old_mean;
else
    for j=1:sampnum
        thetatmp1=theta_vec([1 j+1 sh_startind:end]);
        Holdtmp=Calc_transfer_func(thetatmp1);
        resp1tmp=freqresp(Holdtmp,table2array(dat(1+sum(measnum(1:j-1)):measnum(j)+sum(measnum(1:j-1)),3))'*2*pi);
        tmpold_re=squeeze(real(resp1tmp)');
        tmpold_img=squeeze(abs(imag(resp1tmp)'));
        tmpold_mean=[tmpold_re;tmpold_img];
        tmpold_diff=[table2array(dat(1+sum(measnum(1:j-1)):measnum(j)+sum(measnum(1:j-1)),1));table2array(dat(1+sum(measnum(1:j-1)):measnum(j)+sum(measnum(1:j-1)),2))]-tmpold_mean;
        old_diff=vertcat(old_diff,tmpold_diff);
    end
end
old_rss=dot(old_diff,old_diff);
oldval=theta_vec(upd_ind);

tmpval=oldval;
if (stepsiz ~=-1)    %draw sample from normal
    rand_dist=normrnd(0,stepsiz);
    tmpval=max(0,oldval+rand_dist);
    if tmpval == 0
       tmpval=oldval; 
    end
else          %draw sample from uniform
    rand_dist=(rand(1,1) - 0.5*ones(1,1))*0.1;
    tmpval=min(1,max(0.01,oldval+rand_dist));
end

if stepsiz ~=0
    %calculate r to determine the update
    theta_tmp=theta_vec;
    theta_tmp(upd_ind)=tmpval;
    new_diff=[];
    if shared_ind==0
        %individual parameter new RSS calculation
        Hnew=Calc_transfer_func(theta_tmp);
        resp2=freqresp(Hnew,table2array(dat(1:end,3))'*2*pi);
        new_re=squeeze(real(resp2)');
        new_img=squeeze(abs(imag(resp2)'));
        new_mean=[new_re;new_img];
        new_diff=[table2array(dat(1:end,1));table2array(dat(1:end,2))]-new_mean;
    else
        %shared parameter new RSS calculation
        for j=1:sampnum
            thetatmp2=theta_tmp([1 j+1 sh_startind:end]);
            Hnewtmp=Calc_transfer_func(thetatmp2);
            resp2tmp=freqresp(Hnewtmp,table2array(dat(1+sum(measnum(1:j-1)):measnum(j)+sum(measnum(1:j-1)),3))'*2*pi);
            tmpnew_re=squeeze(real(resp2tmp)');
            tmpnew_img=squeeze(abs(imag(resp2tmp)'));
            tmpnew_mean=[tmpnew_re;tmpnew_img];
            tmpnew_diff=[table2array(dat(1+sum(measnum(1:j-1)):measnum(j)+sum(measnum(1:j-1)),1));table2array(dat(1+sum(measnum(1:j-1)):measnum(j)+sum(measnum(1:j-1)),2))]-tmpnew_mean;
            new_diff=vertcat(new_diff,tmpnew_diff);
        end
    end
    new_rss=dot(new_diff,new_diff);

    if (stepsiz ~=-1)  
        logrval=old_rss/2/theta_vec(end)+log(normpdf(oldval,prior_param(1),prior_param(2)))...
        -new_rss/2/theta_vec(end)-log(normpdf(tmpval,prior_param(1),prior_param(2)));
    else
        logrval=old_rss/2/theta_vec(end)+log(unifpdf(oldval,prior_param(1),prior_param(2)))...
        -new_rss/2/theta_vec(end)-log(unifpdf(tmpval,prior_param(1),prior_param(2)));
    end

    critrval=log(rand);
    if critrval < logrval
        param_newiter=tmpval;
    else
        param_newiter=oldval;
    end
else
    param_newiter=oldval;
end

end


