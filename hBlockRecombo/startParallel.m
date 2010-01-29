function [] = startParallel(n)
    try
        isOpen = matlabpool('size') > 0;
        if ~isOpen
            if nargin == 1
                matlabpool(n);
            else
                matlabpool();
            end
        end
    catch exception
        disp 'can not run in parallele'
        disp(exception.message);
    end
end