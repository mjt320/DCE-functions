function [t_ds,C_ds]=DCEFunc_downSample(t,C,t_res_ds,mode)
% down-sample a curve to a lower temporal resolution; used mainly for DCE simulations
% OUTPUT:
% t_ds: down-sampled time points
% C_ds: down-sampled curve values
% INPUT:
% t: original time points
% t_res_ds: new temporal resolution required
% mode: method for downsampling: base new curve value on the 'nearest
% point' or 'mean' value

t_ds=((t_res_ds/2):t_res_ds:max(t)).'; %new time points for d/s data
N_ds=size(t_ds,1); %number of d/s data points

C_ds=nan(N_ds,1); %curve values at new time points

switch mode
    case 'nearestPoint' %down sample by finding the nearest high-res time point for each down sampled time point
        for n=1:N_ds
            [temp,idx]=min(abs(t-t_ds(n)));
            C_ds(n)=C(idx);
        end
    case 'mean' %down sample by taking average of high-res time points within range
        for n=1:N_ds
            C_ds(n)=mean(C((t>(t_ds(n)-0.5*t_res_ds)) & (t<=(t_ds(n)+0.5*t_res_ds))));
        end
end

end