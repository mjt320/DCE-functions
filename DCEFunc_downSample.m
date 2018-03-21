function [t_ds,C_ds]=DCEFunc_downSample(t,C,t_res_ds,mode)
% down-sample the curve to a lower temporal resolution; used mainly for simulations
% t_ds: down-sampled time points
% C_ds: down-sampled curve values
% t: original time points
% t_res_ds: new temporal resolution
% mode: method for sampling: 'nearest point' or 'mean'

t_ds=((t_res_ds/2):t_res_ds:max(t)).';
N_ds=size(t_ds,1);

C_ds=nan(N_ds,1);

switch mode
    case 'nearestPoint' %down sample by finding the nearest high-res time-point for each time
        for n=1:N_ds
            [temp,idx]=min(abs(t-t_ds(n)));
            C_ds(n)=C(idx);
        end
    case 'mean' %down sample by taking average of high-res time points within range
        for n=1:N_ds
            C_ds(n)=mean(C((t>(t_ds(n)-0.5*t_res_ds)) & (t<(t_ds(n)+0.5*t_res_ds))));
        end
end

end