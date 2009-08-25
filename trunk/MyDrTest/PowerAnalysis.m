function [StatS, StatR, StatT, Truth] = PowerAnalysis(rev, replaceSign)
	loadFromFile = 0;
	if loadFromFile == 0
		precision = 4;
		Trials = 100;
		estimateStd = 0;
		Trialsonce = 1000;
		stdFromTest = 1; % esitimate variance from the Test values, with some power loss.
		muFromTest = 0; % estimate mu from Test values
		if estimateStd == 0
			Trials = Trialsonce;
		end
		
		% for quick debug
		mini = 0;
		nSeq = 100;
		nSnps = 30;

		if nargin == 0
			rev = 1;
		end
		%% read the haplotype sequences
		cd '/home/xzhou/research_linux/gnome/workspace/data/dist100x77'
		if rev
			fastafile = 'SIM_100x77_smp.fasta';
			fastaRef = 'SIM_100x77_ctl.fasta';
			fastafile_testgroup = 'SIM_100x77_ctl.fasta';
			groupname = 'sample+case100Sign0.1';
			outgroupname = 'samegroup100';
		else
			fastafile = '80SNP_CEU_sim_4000seq_control1000.fasta';
			fastaRef = '80SNP_CEU_sim_4000seq_sample1000.fasta';
			fastafile_testgroup = '80SNP_CEU_sim_4000seq_sample1000.fasta';
			groupname = '';
			outgroupname = '';
		end
		
		seqS = fastaread(fastafile);
		nSample = length(seqS);
		Len = length(seqS(1).Sequence);
		
		seqR = fastaread(fastaRef);
		seqT = fastaread(fastafile_testgroup); % test group
		seqT = seqT(1:nSample);

		%seqT = [seqOther; seqT];  % test group
		nT = length(seqT);

		int4S = zeros(nSample,Len);
		int4R = int4S;
		for i = 1:length(seqS)
			int4S(i,:)=nt2int(seqS(i).Sequence) - 1;
			int4R(i,:)=nt2int(seqR(i).Sequence) - 1;
		end

		int4T = zeros(nT,Len); % test group
		for i = 1:length(seqT)
			int4T(i,:)=nt2int(seqT(i).Sequence) - 1;
		end

		%% convert SNP nucleotide to 0/1
		% define allele as allele 1 for each of the SNP location, thus that SNP
		% with non 0/1 label can be converted to 0/1 label
		% int4* is the 4 integer labeling of SNPs ,and int2* is the 2 integer
		% labeling of SNPs
		allele1 = int4T(end,:);
		int2S = (int4S == repmat(allele1,nSample,1)) + 0; 
		int2R = (int4R == repmat(allele1,nSample,1)) + 0;
		int2T = (int4T == repmat(allele1,nT,1)) + 0;
		% samesign = sum(sign(sum(int2S,1)-nS/2).* sign(sum(int2R,1)-nR/2))
		% Len

		%% compute the r-values of the sample group, and reference group
		all_r_S = corrcoef(int2S);
		all_r_S(isnan(all_r_S))= 0; % deal with invariant sites

		%TODO comment this line
		%alterSR = replaceSignByT(all_r_S);
		%all_r_S = alterSR;

		all_r_R = corrcoef(int2R);
		all_r_R(isnan(all_r_R))= 0; % deal with invariant sites

		%% compute the Tr test statistic values
		%% simulation using reference group and obtain the estimated variance
		%% using z-test to optain the p-values for each sample/reference individual
		StatS.Tr = zeros(nSample,1);
		StatS.std = StatS.Tr;
		StatS.p = StatS.Tr;
		StatS.signRate = StatS.Tr;
		
		StatR.Tr = zeros(nSample,1);
		StatR.std = StatR.Tr;
		StatR.p = StatR.Tr;
		StatR.signRate = StatR.Tr;
		
		StatT.Tr = zeros(nT,1);
		StatT.std = StatT.Tr;
		StatT.p = StatT.Tr;
		StatT.signRate = StatT.Tr;

		iMCmodel_R = iMCmodelBuild(seqR, 0); % use refernence to build a model for simulation
		nS = nSample;
		nR = nSample;
		Truth = zeros(nS + nR + nT,1);
		std = getstdTr([], iMCmodel_R, allele1, nS, nR, Trialsonce, precision);
		save();
	else
		load();
	end
		TvsPower = zeros(0, 6);
		k = 0;
		step = 0.01;
		for T = 0:step:1
			k = k + 1;
			%calculate the r value for sample and reference
			for i = 1:nSample
				[StatS.Tr(i), StatS.signRate(i)] = getTr(int2S(i,:), all_r_S, all_r_R, T, step);
				[StatR.Tr(i), StatR.signRate(i)] = getTr(int2R(i,:), all_r_S, all_r_R, T, step);
				
				if estimateStd ==1 
					StatS.std(i) = getstdTr(int4S(i,:), iMCmodel_R, allele1, nS, nR, Trials, precision);
					StatR.std(i) = getstdTr(int4R(i,:), iMCmodel_R, allele1, nS, nR, Trials, precision);
				end
				%???
				sim = sum(repmat(int4S(i,:),nS,1)==int4S,2);
				Truthmax(i) = 1;
				Truthmean(i) = mean(sim)/Len;
				sim = sum(repmat(int4R(i,:),nS,1)==int4S,2);
				Truthmax(i + nS) = max(sim)/Len;
				Truthmean(i + nS) = mean(sim)/Len;
			end

			for i = 1:nT
				[StatT.Tr(i), StatT.signRate(i)] = getTr(int2T(i,:), all_r_S, all_r_R, T, step);
				if estimateStd ==1 
					StatT.std(i) = getstdTr(int4T(i,:), iMCmodel_R, allele1, nS, nR, Trials, precision);
				end
				sim = sum(repmat(int4T(i,:),nS,1)==int4S,2);
				Truthmax(i + nS + nR) = max(sim)/Len;
				Truthmean(i + nS + nR) = mean(sim)/Len;
			end

			StatS.Tr = StatS.Tr/sqrt(Len*(Len-1)/2);
			StatR.Tr = StatR.Tr/sqrt(Len*(Len-1)/2);
			StatT.Tr = StatT.Tr/sqrt(Len*(Len-1)/2);

			TvsPower(k, 1) = T;
			%test the average Tr difference, larger means strong power
			%(depending on the variance)
			TvsPower(k, 2) = mean(StatS.Tr) - mean(StatR.Tr);
			TvsPower(k, 3) = mean(StatS.Tr);
			TvsPower(k, 4) = mean(StatR.Tr);
			TvsPower(k, 5) = mean(StatR.signRate);
			TvsPower(k, 6) = mean(StatS.signRate);
		end

    save();
    figure;
    hold on;
    plot(TvsPower(:,1), TvsPower(:,2), 'marker', '*');
	plot(TvsPower(:,1), TvsPower(:,3), 'marker', '+');
	plot(TvsPower(:,1), TvsPower(:,4), 'marker', 'o');
	legend({['Tr diff'] ['Sample.Tr'] ['Ref.Tr']});
    xlabel('r threshold');
    ylabel('T_r mean diff');
    saveas(gcf, 'TvsPower.pdf');
	
	figure;
	hold on;
	
	plot(TvsPower(1:end-2,1), TvsPower(1:end-2,5), 'marker', '*');
	plot(TvsPower(1:end-2,1), TvsPower(1:end-2,6), 'marker', '+');
	legend('Ref', 'Sample');
	xlabel('r threshold');
	ylabel('signMatchRate');
	saveas(gcf, 'TvsSignRate.pdf');
    %end execution here
    %return;
    
    if stdFromTest ==1
        std = sqrt(var(StatT.Tr));
    end
    mu = 0;
    if muFromTest == 1
        mu = mean(StatT.Tr);
    end
    if estimateStd == 0 
        StatS.std = StatS.std+std;
        StatR.std = StatR.std+std;
        StatT.std = StatT.std+std;
    end
    StatS.p = normcdf(-StatS.Tr,0,StatS.std);
    StatR.p = normcdf(-StatR.Tr,0,StatR.std);
    StatT.p = normcdf(-StatT.Tr,0,StatT.std);
    % plot the p-values

    indxSim95 = find(Truthmax>0.95 & Truthmax<1); 
    indxSim95 = indxSim95(find(indxSim95>nS));
    indxSim = find(Truthmax==1); 
    indxSim = indxSim(find(indxSim>nS));
    Trall = [StatS.Tr; StatR.Tr; StatT.Tr];
    Stdall = [StatS.std; StatR.std; StatT.std];
    pall = [StatS.p; StatR.p; StatT.p];
    indxS = 1:nS;
    indxR = nS+1:nS+nR;
    indxOther = nS+nR+1:nSample;
    %indxT = nSample+1:nS+nR+nT;
	indxT = nS+nR+nSample+1:nS+nR+nSample+nSample;

    %index adjustment
    indxR = setdiff(indxR,union(indxSim95, indxSim));
    indxOther = setdiff(indxOther,union(indxSim95, indxSim));
    indxT = setdiff(indxT,union(indxSim95, indxSim));
    %cd 'E:\1.2.2.Research Proteomics\GWAS information security\manuscripts\fig 2';
    figure;
    hold on;
    % title(['Tr values; ' 'precision: ' num2str(precision) '; Trials: ' num2str(Trials)]);
    nspace = 8;
    indxSim95old = indxSim95;
    I1 = find(indxSim95>nSample);
    indxSim95(I1) = indxSim95(I1)+nspace;
    I1 = find(indxSim95>nS+nR);
    indxSim95(I1) = indxSim95(I1)+nspace;
    I1 = find(indxSim95>nS);
    indxSim95(I1) = indxSim95(I1)+nspace;

    indxSimold = indxSim;
    I1 = find(indxSim>nSample);
    indxSim(I1) = indxSim(I1)+nspace;
    I1 = find(indxSim>nS+nR);
    indxSim(I1) = indxSim(I1)+nspace;
    I1 = find(indxSim>nS);
    indxSim(I1) = indxSim(I1)+nspace;

    plot(indxS, Trall(indxS), '.r');
    plot(nspace+indxR, Trall(indxR), '.g');
    %plot(nspace*2+indxOther,  Trall(indxOther), '.b');
    plot(nspace*3+indxT,  Trall(indxT), '.k');
    plot(indxSim95, Trall(indxSim95old), 'xr');
    plot(indxSim, Trall(indxSimold), '*r');
    a = axis;
    %plot([nS+0.5+nspace nS+0.5+nspace], [a(3) a(4)], '--b', [nS+nR+0.5+2*nspace nS+nR+0.5+2*nspace], [a(3) a(4)], '--b',[nSample+0.5+3*nspace nSample+0.5+3*nspace], [a(3) a(4)], '--b');
    y0 = std*norminv(0.95)+mu;
    plot([a(1) a(2)], [y0 y0], '--k')
    xlabel('indx of individuals');
    ylabel('T_r values');
    legend({['Case'] ['Ref'] ['Test1'] ['Test2'] 'S95' 'S100'});
    batch = floor(rand*10000);
    saveas(gcf, ['Tr_' groupname '_vs_' outgroupname '_' num2str(batch) '.fig']);
    saveas(gcf, ['Tr_' groupname '_vs_' outgroupname '_' num2str(batch) '.pdf']);

    figure;
    hold off;
    semilogy(Truthmax(1:nS), StatS.p, '.r');
    hold on;
    semilogy(Truthmax(nS+1:nS+nR), StatR.p, '.g');
    semilogy(Truthmax(nS+nR+1:nSample), StatT.p(1:nSample-(nS+nR)), '.b');
    semilogy(Truthmax(nSample+1:nS+nR+nT), StatT.p(nSample-(nS+nR)+1:nT), '.k');
    xlabel('indx of individuals');
    ylabel('p value');
    xlabel('max similarity to the sample individuals');
    % title(['max similarity vs. p-value values; ' 'precision: ' num2str(precision) '; Trials: ' num2str(Trials)]);
    legend({['Case'] ['Reference'] ['Test'] ['Test outgroup']  'Simmax>0.95'});
    a = axis;
    a(2) = 1.05;
    axis(a);
    saveas(gcf, ['simmax_ ' groupname '_vs_' outgroupname '_' num2str(batch) '.fig']);
    saveas(gcf, ['simmax_ ' groupname '_vs_' outgroupname '_' num2str(batch) '.png']);

    nIDS_e2 = sum(StatS.p<0.01)
    nIDR_e2 = sum(StatR.p<0.01)
    nIDT_e2 =sum(StatT.p<0.01)

    nIDS_e10 = sum(StatS.p<10e-10)
    nIDR_e10 = sum(StatR.p<10e-10)

    nIDT_e10 = sum(StatT.p<10e-10)
end

function [std_Dr, mean_Dr] = getstdTr(int4, iMCmodel, refhaplotype, nS, nR, Trials, precision)
    Len = length(int4);
    if Len ~= 0
        allele1 = int4;
    else
        allele1 = iMCgenerate(iMCmodel,1);
    end
    Len = length(allele1);
    % Sample one individual form the model
    A = (allele1 == refhaplotype) + 0;
    A2 = (2*A'-1)*(2*A-1);

    D_r_M0 = zeros(Trials,1);

    for i = 1: Trials
        sample_R = iMCgenerate(iMCmodel,nR);
        sample_S = iMCgenerate(iMCmodel,nS);
        sample_R = (sample_R == repmat(refhaplotype,nR,1)) + 0;
        sample_S = (sample_S == repmat(refhaplotype,nS,1)) + 0;

        %sample_S(end,:)=A;
        all_r_R = corrcoef(sample_R);
        all_r_R(isnan(all_r_R))= 0; % deal with invariant sites
        all_r_S = corrcoef(sample_S);
        all_r_S(isnan(all_r_S))= 0; % deal with invariant sites

        % Test statistic for r
        D_r_S(i) = sum(sum((all_r_S - all_r_R).* A2))/2;
    end
    D_r_S = D_r_S/sqrt(Len*(Len-1)/2);
    % figure
    % hist(D_r_S);
    mean_Dr = mean(D_r_S);
    std_Dr = sqrt(var(D_r_S));
    return
end

% T is the threshold
function [Tr signRate] = getTr(Y, r_S, r_R, T, step)
    A2 = (2*Y'-1)*(2*Y-1);
    if nargin == 3
        T = 1;
	end
	Total = (r_S - r_R).*A2;
	%filter the value
	Total(abs(r_S) > T) = 0;
	Tr = sum(sum(Total))/2;
	
	index = (abs(r_S) < T)|(abs(r_S) >= (T+step));
	%check the sign
	rSign = sign(r_S - r_R);
	signMatch = rSign .* A2;
	signMatch(index) = 0;
	
	matched = signMatch;
	matched(matched ~= 1) = 0;
	signRate = sum(sum(matched))/sum(sum(abs(signMatch)));
end