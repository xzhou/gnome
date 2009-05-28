import java.io.FileOutputStream;
import java.io.IOException;
import java.io.ObjectOutputStream;
import java.util.Vector;

import matcher.Matcher;

import com.mathworks.toolbox.javabuilder.MWArray;
import com.mathworks.toolbox.javabuilder.MWClassID;
import com.mathworks.toolbox.javabuilder.MWComplexity;
import com.mathworks.toolbox.javabuilder.MWNumericArray;

public class UtilityFunctions {
	
	static Matcher	bnb;
	
	public UtilityFunctions(){
		try{
			bnb = new Matcher();
		}
		catch(Exception e)
		{
			bnb = null;
			System.out.println("Exception: " + e.toString());
			System.exit(0);
		}
	}
	
	public static RSquare getRSquare(SNP[] snps, int i, int j) {
		return snps[Math.max(i, j)].rsquares.get(Math.min(i, j));
	}
	
	public static boolean isMatchByMatlab(boolean print, Double p1, Double p2,
			Double p3, Double p23, Double p12, Double p13) {
		if (print)
			System.out.println("" + p1 + "," + p2 + "," + p3 + "," + p23 + ","
					+ p12 + "," + p13);
		MWNumericArray B = null; /* Stores input values a */
		Object[] result = null; /* Stores the result */
		boolean ret = false;
		try {
			
			if(bnb == null)
			{
				return false;
			}
			
			int[] dims = { 1, 8 };
			B = MWNumericArray.newInstance(dims, MWClassID.DOUBLE,
					MWComplexity.REAL);
			int[] index = { 1, 1 };
			double devi = -0.01;
			// double devi = 0;
			index[0] = 1;
			index[1] = 1;
			B.set(index, devi); // devi
			double pr = 0;
			if (print)
				pr = 1;
			index[0] = 1;
			index[1] = 2;
			B.set(index, pr); // print
			index[0] = 1;
			index[1] = 3;
			B.set(index, p1);
			index[0] = 1;
			index[1] = 4;
			B.set(index, p2);
			index[0] = 1;
			index[1] = 5;
			B.set(index, p3);
			index[0] = 1;
			index[1] = 6;
			B.set(index, p23);
			index[0] = 1;
			index[1] = 7;
			B.set(index, p12);
			index[0] = 1;
			index[1] = 8;
			B.set(index, p13);

			
			// System.out.println(B);
			result = bnb.isMatch(1, B);
			// System.out.println("Result:");
			MWNumericArray re = (MWNumericArray) result[0];
			if (re.getDouble(1) != 0)
				ret = true;
			else
				ret = false;
		} catch (Exception e) {
			System.out.println("Exception: " + e.toString());
			System.exit(0);
		}

		finally {
			/* Free native resources */
			MWArray.disposeArray(B);
			MWArray.disposeArray(result);
		}
		return ret;
	}
	
	
	public static double computeConfidenceValue(SNP[] snps, RSquare rs,
			boolean pos)
	{
		int index1 = rs.snp1.getID();
		int index2 = rs.snp2.getID();
		double p23 = rs.v2;
		if (pos)
			p23 = rs.v1;
		double ret = 0;
		int total = 0;
		for (int i = 0; i < snps.length; i++) {
			if (i == index1)
				continue;
			if (i == index2)
				continue;

			Vector<Double> pAB1 = new Vector<Double>();
			Vector<Double> pAB2 = new Vector<Double>();
			RSquare rs1 = UtilityFunctions.getRSquare(snps, i, index1);
			RSquare rs2 = UtilityFunctions.getRSquare(snps, i, index2);

			// if (!rs1.isSignRecovered()) continue;
			// if (!rs2.isSignRecovered()) continue;

			if (rs1.isSignRecovered()) {
				pAB1.add(rs1.pAB);
			} else {
				pAB1.add(rs1.v1);
				pAB1.add(rs1.v2);
			}

			if (rs2.isSignRecovered()) {
				pAB2.add(rs2.pAB);
			} else {
				pAB2.add(rs2.v1);
				pAB2.add(rs2.v2);
			}

			int count = 0, count1 = 0;
			for (int i1 = 0; i1 < pAB1.size(); i1++) {
				for (int i2 = 0; i2 < pAB2.size(); i2++) {
					double p12 = pAB1.get(i1), p13 = pAB2.get(i2);

					boolean angel = false;
					
					angel = isMatchByMatlab(false, snps[i].getpA(),
							snps[index1].getpA(), snps[index2].getpA(), p23,
							p12, p13);
					if (angel) {
						count1++;
					}
					count++;
				}
			}
			ret += count1 * 1.0 / count;
			total++;
		}
		System.out.println("pair (" + index1 + "," + index2 + ") pos=" + pos
				+ ",ConfidenceInterval=" + (ret / total) + ",r=" + rs.r
				+ ",total=" + total);
		return ret / total;
	}
	
	public static void saveSNPsToFile(String fileName, SNP[] snps)
		throws IOException {
			// Use a FileOutputStream to send data to a file
			// called myobject.data.
			FileOutputStream f_out = new FileOutputStream(fileName);
			
			// Use an ObjectOutputStream to send object data to the
			// FileOutputStream for writing to disk.
			ObjectOutputStream obj_out = new ObjectOutputStream(f_out);
			
			// Pass our object to the ObjectOutputStream's
			// writeObject() method to cause it to be written out
			// to disk.
			obj_out.writeObject(snps);
		}

	public static int signRecoverRate(SNP[] snps) {
		int count = 0, err = 0;
		for (int index1 = 0; index1 < snps.length - 1; index1++) {
			for (int index2 = index1 + 1; index2 < snps.length; index2++) {
				RSquare rs = UtilityFunctions.getRSquare(snps, index1, index2);

				if (rs.isSignRecovered()) {
					count++;
					if ((rs.r > 0 && rs.pAB == rs.v2)
							|| (rs.r < 0 && rs.pAB == rs.v1)) {
						err++;
					}
				}
			}
		}
		int total = snps.length * (snps.length - 1) / 2;
		System.out.println("" + count + " out of " + total + "("
				+ (count * 1.0 / total) + ") are recovered,err=" + err);
		return count;
	}

}
