function [double single]= readSNPplotter(filename)
	% read pairwise frequency information and compute the 
	% i	j	pA	nameA	pA'	namea	cAMax	cAMin	cANA	pB	nameB	pB'	nameb	cBMax	cBMin	cBNA	pAB	r^2
	% 1	2	0.5901408	2	0.4098592	1	1257	873	NA	0.7744361	2	0.2255639	1	1648	480	2	0.3661234	0.1955739

	if nargin == 0
		filename = '80SNP_CEU_sim_4000seq_control1000.freq.txt'
	end
	
	headerpattern = '%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s';
	datapattern = '%f32 %f32%f32%s%f32%s%f32%f32%s%f32%s%f32%s%f32%f32%s%f32%f32%f%d%d%d%d%d%d%d%d%d';

	fid = fopen(filename, 'r');
	x = textscan(fid,headerpattern,1);
	All = textscan(fid,datapattern);
	% [I J pA0 nameA0 pA1 nameA1 cAMax cAMin cANA pB0 nameB0 pB1 nameB1 cBMax cBMin cBNA pAB r2] =
	%  1 2 3    4      5   6      7     8     9    10  11    12   13     14   15   16   17   18  

	%i j pA nameA pA' namea cAMax cAMin cANA pB nameB pB' nameb cBMax cBMin cBNA pAB r r^2 n11 n12 n13 n21 n22 n23 n31 n32 n33
	%1 2 3   4     5  6     7     8     9    10  11   12  13    14    15    16   17  18 19 20  21  22  23  24  25  26  27  28 

	Len = max(All{2});
	r = zeros(Len,Len);
	p1 = zeros(Len,1);
	p11 = zeros(Len,Len);


	p1 = [All{3}(1);All{10}(1:Len-1)];


	indx = All{1}+ Len*(All{2}-1);

	%read p11
	p11(indx) = All{17};

	%
	indx = All{2}+ Len*(All{1}-1);


	p11(indx) = All{17};


	p11((1:Len)*(Len+1)-Len) = p1;


	r = (p11 - p1*p1')./sqrt((p1 * p1') .* ((1-p1) * (1-p1')));

	single.name = {All{4}{1} All{11}{1:Len-1}; All{6}{1} All{13}{1:Len-1}};

	single.p = p1';

	double.p = p11;

	double.r = r;

