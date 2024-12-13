function [Sr,sys,flag]=RatAppLocal(S,sys)
% [Sr,sys,flag]=RatAppLocal(S,sys)
% Computes a rational model of the data present in S by using the bootstrap
% structure sys - the optimisation is done ine l2 terms directly on the data 

DVC=GetDefaultValuesAndConstants(S);
sign_in_zero_trans_12=DVC.CDAFS.sign_in_zero_trans_12;
adjust=DVC.RA.shift_adjustment & IsFilter(S);
adjust_post=DVC.RA.shift_adjustment_post & IsFilter(S);
stop_tresh=DVC.RA.func_stop_crit_optim;
stop_tresh_post=DVC.RA.func_stop_crit_post_optim;
solver_flag=DVC.RA.solver_flag;
n=length(sys.a);
flag='true';

thetas=L2C(S.freq);
try
    switch(solver_flag)
        case 'MRARL2'
            tic;
            [p,q,beta]=Run_MRarl2_l2(thetas,S.value,sys,n,adjust_post,stop_tresh_post);
            fprintf('Time to compute post-optim  in rat. app. %f\n',toc);
            if (adjust_post)
                [S,sys,p]=ApplyAdditionalBeta(S,beta,sign_in_zero_trans_12,sys,p);
            end
        otherwise
            tic;
            [p,q,beta,sys]=Run_Rarl2_l2(thetas,S.value,sys,n,adjust_post,stop_tresh_post);
            fprintf('Time to compute post-optim  in rat. app. %f\n',toc);
            if (adjust_post)
                [S,sys,p]=ApplyAdditionalBeta(S,beta,sign_in_zero_trans_12,sys,p);
            end
    end
catch err
    err
    flag='false';
    warning('Local Optimsation failed');
    Sr=S;
end
% Passing to plane
[p,q]=MoebSubsInTrans(p,q);
Sr=S;
Sr.p=p;
Sr.q=q;
Sr=AdjustValueAtInfinityFromRattApp(Sr);