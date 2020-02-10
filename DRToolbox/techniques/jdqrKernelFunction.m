function y=jdqrKernelFunction(v,flag)


%make defaults dynamic
global dataX verbose widthparam
if verbose; fprintf('In jdqrKernelFunction\n'); end
if nargin < 2
y=kernel_function_fast(v,dataX,true,'gauss',widthparam,3,'Normal');
else
    if strcmp(flag,'dimension')
        y=size(dataX,2);
        if verbose; fprintf('Size was requested\n'); end
    end
end


end 