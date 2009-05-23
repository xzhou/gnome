



public class Test {
	class Item { //
		int iAC;
		int iAD;
		int iBC;
		int iBD;
	}
	
	//calculate rsquare
	private static double calculateRSQAURE(int iAC, int iAD, int iBC, int iBD) {
		if (iAC < 0 || iAD < 0 || iBC < 0 || iBD < 0) {
			System.err.println("error calculateRSQUARE("+iAC+","+iAD+","+iBC+","+iBD+")");
			System.exit(0);
		}
		
		long total = iAC + iAD + iBC + iBD;
		double x11,x12,x21,x22,p1,q1,p2,q2,rsquare;
		x11 = iAC * 1.0 / total;
		x12 = iAD * 1.0 / total;
		x21 = iBC * 1.0 / total;
		x22 = iBD * 1.0 / total;
		
		p1 = x11 + x12;
		q1 = x11 + x21;
		p2 = x21 + x22;
		q2 = x12 + x22;
		
		rsquare = (x11-p1*q1)*(x11-p1*q1)/(p1*q1*p2*q2);
		
		return rsquare;
	}
	
	//alternative of the method calculateRSQAURE
	private static double calculateRSQAURE_a(int iAC, int iAD, int iBC, int iBD) {
		if (iAC < 0 || iAD < 0 || iBC < 0 || iBD < 0) {
			System.err.println("error calculateRSQUARE_a("+iAC+","+iAD+","+iBC+","+iBD+")");
			System.exit(0);
		}
		
		int total = iAC + iAD + iBC + iBD;
		long x11,x12,x21,x22,p1,q1,p2,q2;
		x11 = iAC;
		x12 = iAD;
		x21 = iBC;
		x22 = iBD;
		
		p1 = x11 + x12;
		q1 = x11 + x21;
		p2 = x21 + x22;
		q2 = x12 + x22;
		
		//System.out.println(""+Math.pow(x11*total-p1*q1,2));
		//System.out.println(""+(p1*q1*p2*q2));
		double rsquare = Math.pow(x11*total-p1*q1,2) / (p1*q1*p2*q2);
		
		return rsquare;
	}
	
	private static double calculateDPrime(int iAC, int iAD, int iBC, int iBD) {
		if (iAC < 0 || iAD < 0 || iBC < 0 || iBD < 0) {
			System.err.println("error calculateDPrime("+iAC+","+iAD+","+iBC+","+iBD+")");
			System.exit(0);
		}
		
		long total = iAC + iAD + iBC + iBD;
		long x11,x12,x21,x22,p1,q1,p2,q2;
		x11 = iAC;
		x12 = iAD;
		x21 = iBC;
		x22 = iBD;
		
		p1 = x11 + x12;
		q1 = x11 + x21;
		p2 = x21 + x22;
		q2 = x12 + x22;
		
		double D = (x11*total-p1*q1)*1.0;
		double Dmax = 0;
		if (D >= 0) 
			Dmax = Math.min(p1*q2, p2*q1);
		else
			Dmax = Math.max(-p1*q1, -p2*q2);
		
		//System.out.println("D="+D+",Dmax="+Dmax);
		return D/Dmax;
	}
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		
		
		int total = 100;
		int iAC = 0, iAD = 0, iBC = 0, iBD = 0;
		
//		iAC = 23; iAD = 336; iBC = 233; iBD = 80;
//		iAC = 642; iAD = 358; iBC = 0; iBD = 0;
//		//int a[] = {23,336,233,80};
//		int a[] = {5,294,544,157};
//		int count = 0;
//		for (int i=0;i<4;i++) {
//			for (int j=0;j<4;j++) {
//				if (i == j) continue;
//				for (int k=0;k<4;k++) {
//					if (k == i || k == j) continue;
//					
//					for (int l=0;l<4;l++) {
//						if (l == k || l == i || l == j) continue;
//						
//						//System.out.println(""+a[i]+" "+a[j]+" "+a[k]+" "+a[l]+" "+calculateRSQAURE(a[i],a[j],a[k],a[l]));
//						System.out.println(""+(++count)+":"+i+" "+j+" "+k+" "+l+" "+calculateRSQAURE_a(a[i],a[j],a[k],a[l]));
//						//System.exit(0);
//					}
//				}
//			}
//		}
//		
//		System.exit(0);
		
		double x11,x12,x21,x22,p1,q1,p2,q2,rsquare,chisquare,dprime;
		int k = 0;
		for (iAC=0;iAC<=total;iAC++) {
			for (iAD=0;iAD<=total;iAD++) {
				for (iBC=0;iBC<=total;iBC++) {
					iBD = total - iAC - iAD - iBC;
					if (iBD < 0) continue;
					
					rsquare = calculateRSQAURE(iAC,iAD,iBC,iBD);
					chisquare = calculateChiSquare(iAC,iAD,iBC,iBD);
					dprime = calculateDPrime(iAC,iAD,iBC,iBD);

//					System.out.println(Integer.toString(iAC) + " "
//							+ Integer.toString(iAD) + " "
//							+ Integer.toString(iBC) + " "
//							+ Integer.toString(iBD) + " "
//							+ rsquare + " "
//							+ chisquare + " "
//							+ dprime + " ");
//							//+ new PrintfFormat("%10.6f").sprintf(rsquare));
					
					
					if (chisquare >= 11.4 && chisquare < 12.8) {
						if (iAC+iAD != 50) continue;
						k++;
						System.out.println(Integer.toString(iAC) + " "
							+ Integer.toString(iAD) + " "
							+ Integer.toString(iBC) + " "
							+ Integer.toString(iBD) + " "
							+ rsquare + " "
							+ chisquare + " "
							+ dprime + " ");
							//+ new PrintfFormat("%10.6f").sprintf(rsquare));
					}
				}
			}
		}
		System.out.println("k="+k);
	}

	private static double calculateChiSquare(int iAC, int iAD, int iBC, int iBD) {
		if (iAC < 0 || iAD < 0 || iBC < 0 || iBD < 0) {
			System.err.println("error calculateRSQUARE("+iAC+","+iAD+","+iBC+","+iBD+")");
			System.exit(0);
		}
		
		long total = iAC + iAD + iBC + iBD;
		double x11,x12,x21,x22,p1,q1,p2,q2,chisquare;
		x11 = iAC * 1.0;
		x12 = iAD * 1.0;
		x21 = iBC * 1.0;
		x22 = iBD * 1.0;
		
		p1 = x11 + x12;
		q1 = x11 + x21;
		p2 = x21 + x22;
		q2 = x12 + x22;
		
		chisquare = Math.pow((x11*x22-x12*x21), 2)*(p1+p2)/(p1*p2*q1*q2);
		
		return chisquare;
	}
}
