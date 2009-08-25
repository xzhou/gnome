function [] = dataAnalysis()
    cd '/home/xzhou/research_linux/gnome/workspace/data/realdata'

    [double single] = readSNPplotter('caseonly.freq.txt');
    Sr = double.r;
    Sp = double.p;
    [double single] = readSNPplotter('controlonly.freq.txt');
    Rr = double.r;
    Rp = single.p;

    save;
    load;

    sampleSign = sign(Sr);
    refSign = sign(Rr);

    signDiff = sampleSign .* refSign;

    nrow = length(Sr);
    x = true(nrow);
    index = triu(x);
    mixR = Rr;
    mixR(index) = Sr(index);

    plotLdWithSign(mixR, signDiff, 0);
end