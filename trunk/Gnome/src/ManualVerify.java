import java.io.BufferedReader;
import java.io.DataInputStream;
import java.io.FileInputStream;
import java.io.InputStreamReader;
import java.util.Vector;

import lPSolver.LPSolver;



/**
 * 
 */

/**
 * @author xzhou
 *
 */

public class ManualVerify {

	Allele[] firstLine = null;	//standard line
	Allele[] secondLine = null; //alternative line

	public SNP[] readFastaFile(String fileName, int cutSNP, int cutRec)
	{
		Vector<Allele[]> alleles = new Vector<Allele[]>();
		int numberOfRecord = cutRec;	//back up ot cutRec
		try {
			// Open the file that is the first
			// command line parameter
			FileInputStream fstream = new FileInputStream(fileName);
			// Get the object of DataInputStream
			DataInputStream in = new DataInputStream(fstream);
			BufferedReader br = new BufferedReader(new InputStreamReader(in));
			String strLine;
			// Read File Line By Line
			Allele[] firstLine_new = null;
			firstLine = null;
			secondLine = null;
			while ((strLine = br.readLine()) != null) {
				// Print the content on the console
				strLine = strLine.trim();
				if (!strLine.isEmpty() && strLine.charAt(0) == '>') { // snp
					//now strLine is the neucleotide sequence
					strLine = br.readLine();
					if (cutSNP > 0 && cutSNP < strLine.length())
						strLine = strLine.substring(0, cutSNP);
					
					//convert the string to alleles
					Allele[] temp = GWASAttack.getAlleleFromStr(strLine);
					Allele[] rewrite;
					
					//first line is the 		
					if (firstLine == null) {
						firstLine = temp;
						secondLine = new Allele[temp.length];
						for (int m = 0; m < temp.length; m++)
							secondLine[m] = firstLine[m];		//copy first line

						firstLine_new = new Allele[firstLine.length];
						for (int i = 0; i < firstLine.length; i++)
							firstLine_new[i] = Allele.A;	
						rewrite = firstLine_new;	//the first line is the standard
					} else {
						// rewrite the string
						rewrite = new Allele[temp.length];
						for (int i = 0; i < temp.length; i++) {
							if (temp[i] == firstLine[i]) {
								rewrite[i] = Allele.A;
							} else {
								secondLine[i] = temp[i];	//this is a different type of allele, we name it as T
								rewrite[i] = Allele.T;
							}
						}
					}
					alleles.add(rewrite);
					// System.out.println(strLine);

					cutRec--;
					if (cutRec == 0)
						break;
				}
			}
			// Close the input stream
			in.close();
		} catch (Exception e) {// Catch exception if any
			System.err.println("Error: " + e.getMessage());
		}
		
		
		Vector<Allele[]> all = new Vector<Allele[]>();
		for (int i = 0; i < alleles.get(0).length; i++) {
			Allele[] a = new Allele[alleles.size()];
			for (int j = 0; j < alleles.size(); j++) {
				a[j] = alleles.get(j)[i];
			}
			all.add(a);
		}
		
		SNP[] all_snps = new SNP[all.size()]; // 
		for (int i = 0; i < all.size(); i++) {
			all_snps[i] = new SNP("snp" + i, i, all.get(i), numberOfRecord);
		}
		
		computeRSquares(all_snps);
		return all_snps;
	}
	
	private void computeRSquares(SNP[] all_snps) {
		GWASAttack.computeRSquares(all_snps);
	}

	public void testFalseNegative()
	{
		SNP[] snps = readFastaFile("./data/casectrl_77SNP.fasta", 77, 100);
		
		
		long totalPair = 0;
		long numberOfFalse = 0;
		//verify the with some correct data;
		for(int index1 = 0; index1 < snps.length - 1; index1 ++)
			for(int index2 = 0; index2 < snps.length - 1; index2 ++)
				for(int index3 = 0; index3 < snps.length -1; index3 ++)
					for(int index4 = 0; index4 < snps.length - 1; index4 ++)
					{
						
						if(index1 == index2 || index1 == index3 || index1 == index4 || index2 == index3 || index2 == index4 || index3 == index4)
							continue;
						RSquare rs12 = UtilityFunctions.getRSquare(snps, index1, index2);
						RSquare rs13 = UtilityFunctions.getRSquare(snps, index1, index3);
						RSquare rs14 = UtilityFunctions.getRSquare(snps, index1, index4);
						RSquare rs23 = UtilityFunctions.getRSquare(snps, index2, index3);
						RSquare rs24 = UtilityFunctions.getRSquare(snps, index2, index4);
						RSquare rs34 = UtilityFunctions.getRSquare(snps, index3, index4);
						
						Double p1A = snps[index1].getpA();
						Double p2A = snps[index2].getpA();
						Double p3A = snps[index3].getpA();
						Double p4A = snps[index4].getpA();
						
						Double p12 = rs12.getpAB();
						Double p13 = rs13.getpAB();
						Double p14 = rs14.getpAB();
						Double p23 = rs23.getpAB();
						Double p24 = rs24.getpAB();
						Double p34 = rs34.getpAB();
						//p12 = 1.0;
						boolean isValid = SignRecover.isConsistent4(false, p1A, p2A, p3A, p4A, p12, p13, p14, p23, p24, p34);
						totalPair ++;
						if(!isValid)
							numberOfFalse ++;
						
						System.out.println("(" + index1 + "," + index2 + "," + index3 + "," + index4 + ")" + isValid);
						
					}
		
		double fp = numberOfFalse*1.0/totalPair;
		System.err.println(fp);
	}
}


