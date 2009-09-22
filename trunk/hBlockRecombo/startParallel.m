function [] = startParallel(n)
    isOpen = matlabpool('size') > 0;
    
    if ~isOpen
        matlabpool(n);
    end
end