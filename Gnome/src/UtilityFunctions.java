import java.io.FileOutputStream;
import java.io.IOException;
import java.io.ObjectOutputStream;
import java.util.Vector;

public class UtilityFunctions {
	public static RSquare getRSquare(SNP[] snps, int i, int j) {
		return snps[Math.max(i, j)].rsquares.get(Math.min(i, j));
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
					
//					angel = isMatchByMatlab(false, snps[i].getpA(),
//							snps[index1].getpA(), snps[index2].getpA(), p23,
//							p12, p13);
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
