function [] = startParallel(n)
    try
        isOpen = matlabpool('size') > 0;
        if ~isOpen
            matlabpool(n);
        end
    catch exception
        disp 'can not run in parallele'
        disp(exception.message);
    end
end