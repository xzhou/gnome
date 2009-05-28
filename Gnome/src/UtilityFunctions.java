import java.io.FileOutputStream;
import java.io.IOException;
import java.io.ObjectOutputStream;
import java.util.Vector;

public class UtilityFunctions {
	public static RSquare getRSquare(SNP[] snps, int i, int j) {
		return snps[Math.max(i, j)].rsquares.get(Math.min(i, j));
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
