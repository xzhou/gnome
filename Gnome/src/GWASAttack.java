import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.DataInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.PrintStream;
import java.util.Collection;
import java.util.Enumeration;
import java.util.Random;
import java.util.StringTokenizer;
import java.util.Vector;

import matcher.Matcher;

//import org.gnu.glpk.GlpkSolver;

//import m_j.m_jclass;


import com.mathworks.toolbox.javabuilder.MWArray;
import com.mathworks.toolbox.javabuilder.MWClassID;
import com.mathworks.toolbox.javabuilder.MWComplexity;
import com.mathworks.toolbox.javabuilder.MWNumericArray;

import jp.ac.kobe_u.cs.cream.DefaultSolver;
import jp.ac.kobe_u.cs.cream.IBBSearch;
import jp.ac.kobe_u.cs.cream.IntVariable;
import jp.ac.kobe_u.cs.cream.Network;
import jp.ac.kobe_u.cs.cream.SASearch;
import jp.ac.kobe_u.cs.cream.Solution;
import jp.ac.kobe_u.cs.cream.Solver;
import jp.ac.kobe_u.cs.cream.TabooSearch;

public class GWASAttack {
	private Network net;
	private Solver solver;
	private Solution solution;
	private Vector<Record> records = new Vector<Record>();
	public static int RSQUARE_PRECISION = 1000;
	public static int PVALUE_PRECISION = 10000;
	public static int DEVIATION = 4;
	public static int CASE;
	public static Integer BLOCKS;
	public static boolean allGoodSNP = false;

	private int[] counts = null;
	private static char sign = '*'; // '*' means the sign is unknown; '-' means
									// negative;
	// '+' means positive; '0' means zero
	// 'u' means unknown
	private static Matcher bnb = null;

	private static int iSolution = 0;
	private static IntVariable mm = null;
	private static double MAXDEVIATION = 1e-4;
	private static int SOLUTIONS = -1;
	private static int CUTOFF;
	private static double gap = 0.001; // gap for r2
	private static Allele[] firstLine = null;
	private static Allele[] secondLine = null;

	public GWASAttack() {
		net = new Network();
		solver = new DefaultSolver(net);
		// solver = new SASearch(net);
	}

	// set the constraints for RSQUARE
	public void setConstraintForRSquare(RSquare rsquare) {
		if (rsquare.rsquare < 0) { // rsquare is not presented
			rsquare.constraint = null;
			return;
		}

		IntVariable vAC = new IntVariable(net);
		IntVariable vAD = new IntVariable(net);
		IntVariable vBC = new IntVariable(net);
		IntVariable vBD = new IntVariable(net);

		vAC.ge(0);
		vAD.ge(0);
		vBC.ge(0);
		vBD.ge(0); // all frequencies should be larger than or equal to 0

		IntVariable p1 = vAC.add(vAD);
		IntVariable q1 = vAC.add(vBC);
		IntVariable p2 = vBC.add(vBD);
		IntVariable q2 = vAD.add(vBD);

		int total1 = rsquare.snp1.getCasePopulation();
		int total2 = rsquare.snp2.getCasePopulation();
		IntVariable vTotal = vAC.add(vAD).add(vBC).add(vBD);
		if (total1 < total2) {
			vTotal.ge(total1 - PVALUE_PRECISION);
			vTotal.le(total2 + PVALUE_PRECISION);
		} else {
			vTotal.ge(total2 - PVALUE_PRECISION);
			vTotal.le(total1 + PVALUE_PRECISION);
		}

		if (allGoodSNP) {
			rsquare.snp1.setBadSNP(false);
			rsquare.snp2.setBadSNP(false);
		}

		if (rsquare.snp1.isBadSNP() || rsquare.snp2.isBadSNP()) {
			if (!rsquare.snp1.isInvolved()) {
				if (rsquare.snp1.isBadSNP()) {
					Vector<Entry> entries = rsquare.snp1.entries.entries;
					IntVariable zeroA = null, zeroT = null, zeroM = null;
					for (int i = 0; i < entries.size(); i++) {
						Entry entry = entries.get(i);
						if (zeroA == null)
							zeroA = p1.subtract((int) entry.getCaseA());
						else
							zeroA = zeroA.multiply(p1.subtract((int) entry
									.getCaseA()));
						if (zeroT == null)
							zeroT = p2.subtract((int) entry.getCaseT());
						else
							zeroT = zeroT.multiply(p2.subtract((int) entry
									.getCaseT()));
						if (zeroM == null)
							zeroM = p1.multiply(p2)
									.subtract(
											(int) (entry.getCaseA() * entry
													.getCaseT()));
						else
							zeroM = zeroM.multiply(p1.multiply(p2)
									.subtract(
											(int) (entry.getCaseA() * entry
													.getCaseT())));
					}
					zeroA.equals(0);
					zeroT.equals(0);
					zeroM.equals(0);
				} else {
					p1.equals(rsquare.snp1.getCaseAFreq());
					p2.equals(rsquare.snp1.getCaseTFreq());
				}
			}

			if (!rsquare.snp2.isInvolved()) {
				if (rsquare.snp2.isBadSNP()) {
					Vector<Entry> entries = rsquare.snp2.entries.entries;
					IntVariable zeroA = null, zeroT = null, zeroM = null;
					for (int i = 0; i < entries.size(); i++) {
						Entry entry = entries.get(i);
						if (zeroA == null)
							zeroA = q1.subtract((int) entry.getCaseA());
						else
							zeroA = zeroA.multiply(q1.subtract((int) entry
									.getCaseA()));
						if (zeroT == null)
							zeroT = q2.subtract((int) entry.getCaseT());
						else
							zeroT = zeroT.multiply(q2.subtract((int) entry
									.getCaseT()));
						if (zeroM == null)
							zeroM = q1.multiply(q2)
									.subtract(
											(int) (entry.getCaseA() * entry
													.getCaseT()));
						else
							zeroM = zeroM.multiply(q1.multiply(q2)
									.subtract(
											(int) (entry.getCaseA() * entry
													.getCaseT())));
					}
					zeroA.equals(0);
					zeroT.equals(0);
					zeroM.equals(0);
				} else {
					q1.equals(rsquare.snp2.getCaseAFreq());
					q2.equals(rsquare.snp2.getCaseTFreq());
				}
			}
		} else {
			int ip1 = rsquare.snp1.getCaseAFreq();
			int ip2 = rsquare.snp1.getCaseTFreq();
			int iq1 = rsquare.snp2.getCaseAFreq();
			int iq2 = rsquare.snp2.getCaseTFreq();

			double minRSQUARE = (int) (rsquare.rsquare * RSQUARE_PRECISION)
					* 1.0 / RSQUARE_PRECISION;
			double maxRSQUARE = ((int) (rsquare.rsquare * RSQUARE_PRECISION) + 1)
					* 1.0 / RSQUARE_PRECISION;
			int min = (int) Math.sqrt(minRSQUARE * ip1 * iq1 * ip2 * iq2);
			int max = (int) Math.sqrt(maxRSQUARE * ip1 * iq1 * ip2 * iq2) + 1;

			System.out.println("min=" + min + ",max=" + max + ",minR="
					+ minRSQUARE + ",maxR=" + maxRSQUARE);

			p1.equals(ip1);
			p2.equals(ip2);
			q1.equals(iq1);
			q2.equals(iq2);

			IntVariable absAC = vAC.multiply(vTotal).subtract(ip1 * iq1).abs();
			IntVariable absAD = vAD.multiply(vTotal).subtract(ip1 * iq2).abs();
			IntVariable absBC = vBC.multiply(vTotal).subtract(ip2 * iq1).abs();
			IntVariable absBD = vBD.multiply(vTotal).subtract(ip2 * iq2).abs();
			absAC.ge(min);
			absAC.le(max);
			absAD.ge(min);
			absAD.le(max);
			absBC.ge(min);
			absBC.le(max);
			absBD.ge(min);
			absBD.le(max);
		}

		rsquare.constraint = new ConstraintForRSquare(vAC, vAD, vBC, vBD);
	}

	// set the constraints for PVALUE
	public void setConstraintForPvalue(RSquare pvalue) {
		if (pvalue.rsquare < 0) { // pvalue is not presented
			pvalue.constraint = null;
			return;
		}

		int total = pvalue.getTotal();

		IntVariable vAC = new IntVariable(net);
		IntVariable vAD = new IntVariable(net);
		IntVariable vBC = new IntVariable(net);
		IntVariable vBD = new IntVariable(net);

		vAC.ge(0);
		vAD.ge(0);
		vBC.ge(0);
		vBD.ge(0); // all frequencies should be larger than or equal to 0

		vAC.add(vAD).add(vBC).add(vBD).equals(total);

		IntVariable p1 = vAC.add(vAD);
		IntVariable q1 = vAC.add(vBC);
		IntVariable p2 = vBC.add(vBD);
		IntVariable q2 = vAD.add(vBD);

		int ip1 = pvalue.iAC + pvalue.iAD;
		int ip2 = pvalue.iBC + pvalue.iBD;
		System.out.println("" + ip1 + " " + ip2);
		p1.equals(ip1);
		p2.equals(ip2);

		// System.out.println(""+Math.pow(x11*total-p1*q1,2));
		// System.out.println(""+(p1*q1*p2*q2));
		int iRSQUARE = (int) Math.round(pvalue.rsquare * PVALUE_PRECISION);
		IntVariable div = p1.multiply(q1).multiply(p2).multiply(q2);
		IntVariable left = div.multiply(iRSQUARE);
		IntVariable right = vAC.multiply(total).subtract(p1.multiply(q1));
		right = right.multiply(right).multiply(PVALUE_PRECISION);

		div.notEquals(0);
		IntVariable e = left.subtract(right);
		e.ge(-DEVIATION);
		e.le(DEVIATION);

		pvalue.constraint = new ConstraintForRSquare(vAC, vAD, vBC, vBD);
	}

	public void solve() {
		solution = solver.findFirst();
	}

	public int getIntValue(IntVariable x) {
		return solution.getIntValue(x);
	}

	/**
	 * @param args
	 * @throws Exception
	 */
	public static void main(String[] args) throws Exception {
		
		if (args.length < 1 || args.length > 15) {
			System.out.println("usage:java -jar test.jar filename [cutSNP] [cutRec] [allGoodSNP]"
							+ " [rsquare_precision] [pvalue_precision] [case] "
							+ "[sigSNPIndex1] [sigSNPIndex2]");
			System.out.println("'filename': the full path of the input fasta file");
			System.out.println("'cutSNP': the limit of number of SNPs");
			System.out.println("'cutRec': the limit of number of individual records");
			System.out.println("'allGoodSNP': if true, single SNP frequencies will not be computed");
			System.out.println("'rsquare_precision': the precision of rsquares");
			System.out.println("'pvalue_precision': the precision of p-values");
			System.out.println("'case': # of individuals in case group");
			System.out.println("'sigSNPIndex1': the index of sig SNP 1"); // sig SNP is the one whose frequency is known
			System.out.println("'sigSNPIndex2': the index of sig SNP 2");
			System.out.println("'approachChoice': approach choice"); // there are two available choices: 1-first approach;2-second approach 
			System.exit(0);
		}
		int cutSNP = 0, cutRec = 0;
		if (args.length >= 2) { // cutSNP is provided
			cutSNP = Integer.valueOf(args[1]);
		}
		if (args.length >= 3) { // cutRec is provided
			cutRec = Integer.valueOf(args[2]);
		}
		if (args.length >= 4) { // deviation is provided
			GWASAttack.allGoodSNP = Boolean.valueOf(args[3]);
		}
		if (args.length >= 5) { // rsquare_precision is provided
			GWASAttack.RSQUARE_PRECISION = Integer.valueOf(args[4]);
		}
		if (args.length >= 6) { // pvalue_precision is provided
			GWASAttack.PVALUE_PRECISION = Integer.valueOf(args[5]);
		}
		if (args.length >= 7) { // case is provided
			GWASAttack.CASE = Integer.valueOf(args[6]);
		} else {
			GWASAttack.CASE = cutRec;
		}
		int sigSNPIndex1 = 0, sigSNPIndex2 = 15;
		if (args.length >= 8) { // sigSNPIndex1 is provided
			sigSNPIndex1 = Integer.valueOf(args[7]);
		}
		if (args.length >= 9) { // sigSNPIndex2 is provided
			sigSNPIndex2 = Integer.valueOf(args[8]);
		}
		int approachChoice = 1;
		if (args.length >= 10) { // approachChoice is provided
			approachChoice = Integer.valueOf(args[9]);
		}

		SNP[] all_snps = readSNPsFromFastaFile(args[0], cutSNP, cutRec); // read SNPs from fasta file

		updateRSquaresFromDefaultFile(all_snps);

		all_snps[sigSNPIndex1].setSigSNP(true);
		all_snps[sigSNPIndex2].setSigSNP(true);
		if (!allGoodSNP) computeSingleSNPFrequencies(all_snps, sigSNPIndex1); //compute single SNP frequencies

		if (approachChoice == 1) {   //triangle approach
			computeSignOfRSquare(all_snps);
			propagateSignOfRSquare(all_snps);
		} else if (approachChoice == 2) { //Hapmap approach
			//computeSignOfRSquare(all_snps);
			SNP[] hapmap_snps = readSNPsFromFastaFile("./data/hapmap_chr7_80SNP_CEU_haplotype.fasta", cutSNP,cutRec);
			computeSignOfRSquareFromHapmap(all_snps, hapmap_snps);
			//removeHapmapNoise(all_snps);
		} else {
			throw new Exception("unsupported approach choice "+approachChoice);
		}
	}

	private static void removeHapmapNoise(SNP[] snps) throws Exception {
		double devi = 0.001;
		boolean propagating = false;
		int undecidable = 0;
		while (true) {
			for (int index1 = 0; index1 < snps.length - 1; index1++) {
				undecidable = 0;
				for (int index2 = index1 + 1; index2 < snps.length; index2++) {
					RSquare rs = getRSquare(snps, index1, index2);

					if (rs.rsquare < gap)
						continue;

					if (rs.isSignRecovered()) {
						if (rs.getRecoveredSign() * rs.r < 0)
							System.out.println("noise");
						else {
							// System.out.println("hapmap");
							continue;
						}
					} else {
						continue;
					}

					double i1 = computeConfidenceInterval(snps, rs, true);
					double i2 = computeConfidenceInterval(snps, rs, false);
					double d = i1 - i2;

					if (d >= devi) {
						rs.setSignRecovered(true, true);
						propagating = true;
						if (rs.r > 0)
							System.out.println("right\n");
						else
							System.out.println("wrong\n");
					} else if (d <= -devi) {
						rs.setSignRecovered(false, true);
						propagating = true;
						if (rs.r < 0)
							System.out.println("right\n");
						else
							System.out.println("wrong\n");
					} else {
						System.out.println("undecidable");
						undecidable++;
					}
				}
				signRecoverRate(snps);
				System.out.println("undecidable=" + undecidable);
			}

			if (!propagating)
				break;
			propagating = false;
		}
		System.exit(0);
	}

	private static void propagateSignOfRSquare(SNP[] snps) throws Exception {
		boolean propagating = false;
		double c_i = 0.5;
		while (true) {
			for (int index1 = 0; index1 < snps.length - 1; index1++) {
				for (int index2 = index1 + 1; index2 < snps.length; index2++) {
					RSquare rs = getRSquare(snps, index1, index2);

					if (rs.isSignRecovered()) {
						for (int i = 0; i < snps.length; i++) {
							if (i == index1 || i == index2)
								continue;

							RSquare rs1 = getRSquare(snps, index1, i);
							RSquare rs2 = getRSquare(snps, index2, i);
							if (rs1.isSignRecovered() && rs2.isSignRecovered())
								continue;

							Vector<Double> pAB1 = new Vector<Double>();
							Vector<Double> pAB2 = new Vector<Double>();

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

							double devi = 1e-2;
							for (int i1 = 0; i1 < pAB1.size(); i1++) {
								for (int i2 = 0; i2 < pAB2.size(); i2++) {
									boolean god = false;
									double p23 = rs.pAB, p12 = pAB1.get(i1), p13 = pAB2
											.get(i2);
									if (p12 > 1)
										p12 -= 1;
									if (p13 > 1)
										p13 -= 1;

									if (Math.abs(p12 - rs1.pAB) < devi
											&& Math.abs(p13 - rs2.pAB) < devi) {
										god = true;
									}

									boolean angel = isMatchByMatlab(false,
											snps[i].getpA(), snps[index1]
													.getpA(), snps[index2]
													.getpA(), p23, p12, p13);
									if (angel) {
										if (pAB1.get(i1) < 1)
											pAB1.set(i1, p12 + 1);
										if (pAB2.get(i2) < 1)
											pAB2.set(i2, p13 + 1);
									} else if (god) {
										isMatchByMatlab(true, snps[i].getpA(),
												snps[index1].getpA(),
												snps[index2].getpA(), p23, p12,
												p13);
										// throw new
										// Exception("angel should be true because god is");
										System.err
												.println("angel should be true because god is");
									}
									// System.out.println("angel="+angel);
								}
							}

							for (int i1 = 0; i1 < pAB1.size(); i1++) {
								double p12 = pAB1.get(i1);
								if (p12 < 1) {
									if (p12 == rs1.v1) {
										double interval = computeConfidenceInterval(
												snps, rs1, false);
										if (interval >= c_i) {
											rs1.setSignRecovered(false, true);
											propagating = true;
										}
									} else if (p12 == rs1.v2) {
										double interval = computeConfidenceInterval(
												snps, rs1, true);
										if (interval >= c_i) {
											rs1.setSignRecovered(true, true);
											propagating = true;
										}
									} else {
										throw new Exception(
												"error detected: p12=" + p12
														+ ",v1=" + rs1.v1
														+ ",v2=" + rs1.v2);
									}
								}
							}
							for (int i2 = 0; i2 < pAB2.size(); i2++) {
								double p13 = pAB2.get(i2);
								if (p13 < 1) {
									if (p13 == rs2.v1) {
										double interval = computeConfidenceInterval(
												snps, rs2, false);
										if (interval >= c_i) {
											rs2.setSignRecovered(false, true);
											propagating = true;
										}
									} else if (p13 == rs2.v2) {
										double interval = computeConfidenceInterval(
												snps, rs2, true);
										if (interval >= c_i) {
											rs2.setSignRecovered(true, true);
											propagating = true;
										}
									} else {
										throw new Exception(
												"error detected: p13=" + p13
														+ ",v1=" + rs2.v1
														+ ",v2=" + rs2.v2);
									}
								}
							}
						}
						signRecoverRate(snps);
						saveSNPsToFile("snps.dat", snps);
					}
				}
			}
			if (!propagating)
				break;
			saveSNPsToFile("snps.dat", snps);
			System.out.println("propagating");
			propagating = false;
		}
	}

	//
	private static double computeConfidenceInterval(SNP[] snps, RSquare rs,
			boolean pos) {
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
			RSquare rs1 = getRSquare(snps, i, index1);
			RSquare rs2 = getRSquare(snps, i, index2);

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

					boolean angel = isMatchByMatlab(false, snps[i].getpA(),
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

	private static int signRecoverRate(SNP[] snps) {
		int count = 0, err = 0;
		for (int index1 = 0; index1 < snps.length - 1; index1++) {
			for (int index2 = index1 + 1; index2 < snps.length; index2++) {
				RSquare rs = getRSquare(snps, index1, index2);

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

	private static boolean isMatchByMatlab(boolean print, Double p1, Double p2,
			Double p3, Double p23, Double p12, Double p13) {
		if (print)
			System.out.println("" + p1 + "," + p2 + "," + p3 + "," + p23 + ","
					+ p12 + "," + p13);
		MWNumericArray B = null; /* Stores input values a */
		Object[] result = null; /* Stores the result */
		boolean ret = false;
		try {
			if (bnb == null) {
				bnb = new Matcher();
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

	//Simply use the signs of Hapmap as the signs for the real data
	private static void computeSignOfRSquareFromHapmap(SNP[] snps,
			SNP[] hapmap_snps) {
		int rec = 0, rev_rec = 0, true_rec = 0, true_rev_rec = 0;
		int hap_signs = 0, zeros = 0;
		double meanR = 0, meanNR = 0;
		System.out.println("SNPi SNPj r(recovered) r(real) r2");
		for (int index1 = 0; index1 < snps.length - 1; index1++) {
			for (int index2 = index1 + 1; index2 < snps.length; index2++) {
				RSquare rs = getRSquare(snps, index1, index2);
				RSquare hap_rs = getRSquare(hapmap_snps, index1, index2);

				double hap_sign = hap_rs.getSign();
				if (hap_sign > 0) {
					System.out.println("" + rs.snp1.getID() + " "
							+ rs.snp2.getID() + " " + Math.abs(rs.r) + " "
							+ rs.r + " " + rs.rsquare);
				} else if (hap_sign < 0) {
					System.out.println("" + rs.snp1.getID() + " "
							+ rs.snp2.getID() + " " + (-Math.abs(rs.r)) + " "
							+ rs.r + " " + rs.rsquare);
				}

				if (rs.rsquare == 0) {
					zeros++;
					continue;
				}
				
				double value = rs.rsquare;
				if (value >= gap) {
					rec++;
					if (rs.r * hap_rs.getSign() > 0)
						true_rec++;
				} else {
					rev_rec++;
					if (rs.r * hap_rs.getSign() > 0) {
						true_rev_rec++;
					}
				}
			}
		}

		int total = snps.length * (snps.length - 1) / 2;
		System.out.println("rec=" + true_rec + "/" + rec + "("
				+ (true_rec * 1.0 / rec) + "),rev_rec=" + true_rev_rec + "/"
				+ rev_rec + "(" + (true_rev_rec * 1.0 / rev_rec) + "),total="
				+ total + ",zeros=" + zeros);
		System.exit(0);
	}

	private static void computeSignOfRSquare(SNP[] snps) throws Exception {
		int i = 0, count1 = 0, count2 = 0;
		for (int index1 = 0; index1 < snps.length - 1; index1++) {
			for (int index2 = index1 + 1; index2 < snps.length; index2++) {
				RSquare rs = getRSquare(snps, index1, index2);
				double p1 = rs.getpA();
				double q1 = rs.getpB();
				double p2 = 1 - p1;
				double q2 = 1 - q1;
				double r = Math.sqrt(rs.rsquare * p1 * p2 * q1 * q2);
				double v1 = r + p1 * q1;
				double v2 = -r + p1 * q1;
				rs.v1 = v1;
				rs.v2 = v2;

				double v = -1;
				double devi = 0.01;
				if (((v1 - p1) > devi || (v1 - q1) > devi)) {
					i++;
					v = v2;
				} else if ((v2 - p1) > devi || (v2 - q1) > devi) {
					i++;
					v = v1;
				} else {
					if (q2 - p1 + v1 < -devi) {
						i++;
						v = v2;
					} else if (q2 - p1 + v2 < -devi) {
						i++;
						v = v1;
					}
				}

				if (v < 0)	continue;
				
				if ((v == v1 && rs.r < 0) || (v == v2 && rs.r > 0)) {
					System.err.println("incorect sign");
				}

				rs.pAB = v;
				rs.setSignRecovered(true);
			}
		}
		signRecoverRate(snps);
		System.out.println("i=" + i + ",count1=" + count1 + ",count2=" + count2);
	}

	private static void updateRSquaresFromDefaultFile(SNP[] snps)
			throws Exception {
		
		//TODO change filename
		String fileName = "./data/us.freq.txt";
		
		// Open the file
		FileInputStream fstream = new FileInputStream(fileName);
		// Get the object of DataInputStream
		DataInputStream in = new DataInputStream(fstream);
		BufferedReader br = new BufferedReader(new InputStreamReader(in));
		String strLine;
		// Read File Line By Line
		int lines = 0;
		int cor = 0, corSign = 0;
		while ((strLine = br.readLine()) != null) {
			// Print the content on the console
			strLine = strLine.trim();
			if (strLine.charAt(0) == '#') // ignore comments
				continue;
			StringTokenizer st = new StringTokenizer(strLine);
			int k = 0;
			int i = 0, j = 0;
			double r2 = 0, r = 0, pAB = 0, pA = 0, pB = 0;
			int countA1 = 0, countT1 = 0, countA2 = 0, countT2 = 0;
			while (st.hasMoreTokens()) {
				String token = st.nextToken();
				if (k == 0)
					i = Integer.valueOf(token) - 1; // index1
				else if (k == 1)
					j = Integer.valueOf(token) - 1; // index2
				else if (k == 17)
					r = Double.valueOf(token); // r
				else if (k == 18)
					r2 = Double.valueOf(token); // r2
				else if (k == 6)
					countA1 = Integer.valueOf(token);
				else if (k == 7)
					countT1 = Integer.valueOf(token);
				else if (k == 13)
					countA2 = Integer.valueOf(token);
				else if (k == 14)
					countT2 = Integer.valueOf(token);
				else if (k == 16)
					pAB = Double.valueOf(token);
				else if (k == 2) {
					pA = Double.valueOf(token);
				} else if (k == 9) {
					pB = Double.valueOf(token);
				}

				k++;
			}

			if (i == 73 || i == 49 || i == 64)
				continue;
			if (j == 73 || j == 49 || j == 64)
				continue;
			if (i > 49 && i < 64)
				i--;
			else if (i > 64 && i < 73)
				i = i - 2;
			else if (i > 73)
				i = i - 3;
			if (j > 49 && j < 64)
				j--;
			else if (j > 64 && j < 73)
				j = j - 2;
			else if (j > 73)
				j = j - 3;

			RSquare rsquare = getRSquare(snps, i, j);

			rsquare.realr2 = rsquare.rsquare;
			rsquare.rsquare = getrSQUAREByPrecision(r2);
			int s1 = 1, s2 = 1;
			if (rsquare.snp1.getCaseAFreq() > rsquare.snp1.getCaseTFreq())
				s1 = -1;
			if (rsquare.snp2.getCaseAFreq() > rsquare.snp2.getCaseTFreq())
				s2 = -1;
			if (r * rsquare.r * s1 * s2 > 0)
				corSign++;
			rsquare.r = r;
			rsquare.sign = getSignFromDouble(r);

			rsquare.pAB = pAB;
			rsquare.pA = pA;
			rsquare.pB = pB;

			lines++;
		}
		// Close the input stream
		in.close();
		System.out.println("corSign=" + corSign);
	}

	private static double getrSQUAREByPrecision(double v) {
		return ((int) (v * RSQUARE_PRECISION)) * 1.0 / RSQUARE_PRECISION;
	}

	private static RSquare getRSquare(SNP[] snps, int i, int j) {
		return snps[Math.max(i, j)].rsquares.get(Math.min(i, j));
	}

	private static void printRecords(Vector<Record> records) {
		int k = 0;
		for (Enumeration<Record> rds = records.elements(); rds
				.hasMoreElements();) {
			Record rd = rds.nextElement();
			System.out.println("" + k + ":" + rd.toString());
			k++;
		}
	}

	private static void computeRangesForRecords(Vector<Double[]> M,
			Vector<Double[]> N, Vector<Record> records, SNP[] snps)
			throws Exception {
		// N = null; //test
		int count = 0;
		boolean[] freeVars = null;
		boolean removed = true;
		while (removed) {
			int lead = 0;
			freeVars = new boolean[records.size()];
			count = 0;
			Vector<Integer> zeroRows = new Vector<Integer>();
			for (int i = 0; i < freeVars.length; i++)
				freeVars[i] = false;
			int row = 0;
			for (Enumeration<Double[]> m = M.elements(); m.hasMoreElements();) {
				Double[] mm = m.nextElement();
				for (int i = lead; i < mm.length; i++) {
					if (mm[i] == 1) {
						lead = i;
						break;
					}
				}
				if (mm[lead] != 1) {
					zeroRows.add(row);
					row++;
					continue;
				}
				int freeCount = 0;
				int lastFree = -1;
				boolean gt = false;
				for (int k = lead + 1; k < mm.length - 1; k++) {
					if (mm[k] != 0.0) {
						if (!freeVars[k]) {
							freeVars[k] = true;
							count++;
							records.get(k).setFree(true);
						}

						freeCount++;
						lastFree = k;
						gt = (mm[k] > 0);
					}
				}
				if (freeCount == 1 && gt && mm[mm.length - 1] == 0) {
					records.get(lastFree).setRemovable(true);
					records.get(lead).setRemovable(true);
					System.out.println("remove lead " + lead);
					System.out.println("remove free variable " + lastFree);
					printArray(mm);
				}
				if (freeCount == 0 && mm[mm.length - 1] == 0) {
					records.get(lead).setRemovable(true);
					System.out.println("remove lead " + lead);
					printArray(mm);
				}

				row++;
			}

			removed = false;
			for (int i = records.size() - 1; i >= 0; i--) {
				Record record = records.get(i);
				if (record.isRemovable()) {
					int ccc = calculateRecordCount(snps, record);
					if (ccc > 0) {
						System.out.println("error detected: record " + i + " "
								+ record.toString()
								+ " should not be deleted, count=" + ccc);
						// throw new
						// Exception("error detected: record "+i+" "+record.toString()+" should not be deleted, count="+ccc);
					}
					records.remove(i);
					removeColFromMatrix(M, i);
					if (N != null)
						removeColFromMatrix(N, i);
					System.out.println("record " + i + " removed");
					removed = true;
				}
			}
			for (int i = zeroRows.size() - 1; i >= 0; i--) {
				int kk = zeroRows.get(i);
				M.remove(kk);
				removed = true;
				System.out.println("Row " + kk + " of M is removed");
			}
		}
		// System.exit(0);
		System.out.println("after removing columns");
		if (!isRealSolutionValidBeforeReduce(records, M, snps)) {
			throw new Exception(
					"error detected: real solution is inconsistent with M");
		}
		printMatrix(M);

		if (N == null)
			N = new Vector<Double[]>();

		Double[][] A = new Double[M.size() + N.size()][count];
		Double[] b = new Double[M.size() + N.size()];
		Double[] c = new Double[count];
		Integer[] lb = new Integer[count];
		Integer[] ub = new Integer[count];
		for (int i = 0; i < c.length; i++)
			c[i] = 0.0;
		if (c.length >= 1)
			c[0] = 1.0;
		int n = 0;
		for (Enumeration<Double[]> m = M.elements(); m.hasMoreElements();) {
			Double[] mm = m.nextElement();

			count = 0;
			for (int k = 0; k < mm.length - 1; k++) {
				if (freeVars[k]) {
					A[n][count] = mm[k];

					count++;
				}
			}
			b[n] = mm[mm.length - 1];

			n++;
		}

		// System.out.println("N=");
		// printMatrix(N);
		for (Enumeration<Double[]> u = N.elements(); u.hasMoreElements();) {
			Double[] uu = u.nextElement();
			Double[] cc = new Double[uu.length - 1];
			for (int k = 0; k < cc.length; k++)
				cc[k] = 0.0;
			b[n] = uu[uu.length - 1];
			for (int k = 0; k < cc.length; k++) {
				if (uu[k] != 0) {
					int index = getRowByLeadIndex(M, k);
					if (index >= 0) {
						Double[] tmp = M.get(index);
						for (int j = k + 1; j < cc.length; j++) {
							cc[j] -= uu[k] * tmp[j];
						}
						b[n] -= uu[k] * tmp[tmp.length - 1];
					} else {
						cc[k] += uu[k];
					}
				}
			}
			count = 0;
			for (int k = 0; k < cc.length; k++) {
				if (freeVars[k]) {
					A[n][count] = cc[k];

					count++;
				}
			}

			n++;
		}

		count = 0;
		for (int i = 0; i < records.size(); i++) {
			if (freeVars[i]) {
				lb[count] = records.get(i).getLB();
				ub[count] = records.get(i).getUB();
				count++;
			}
		}

		// find Xs of 0
		for (int i = 0; i < A.length; i++) {
			if (A[i].length > 0) {
				boolean gt = true;
				for (int k = 0; k < A[i].length; k++) {
					if (A[i][k] < 0) {
						gt = false;
						break;
					}
				}
				if (gt && (b[i] == 0)) {
					for (int k = 0; k < A[i].length; k++) {
						if (A[i][k] > 0) {
							ub[k] = 0;
						}
					}
				}
			}
		}
		// remove entries of 0s in A
		// if (A[0].length > 0) {
		// int cc = 0;
		// for (int i=0;i<A.length;i++) {
		// if (A[i].length > 0) {
		// boolean z = true;
		// for (int k=0;k<A[i].length;k++) {
		// if (A[i][k] != 0) {
		// z = false;
		// break;
		// }
		// }
		// if (z) {
		// b[i] = Double.MIN_VALUE;
		// cc++;
		// }
		// }
		// }
		// if (cc != 0) {
		// cc = b.length - cc;
		// System.out.println("cc="+cc);
		// Double[][] AA = new Double[cc][A[0].length];
		// Double[] bb = new Double[cc];
		// cc = 0;
		// for (int kk=0;kk<b.length;kk++) {
		// if (b[kk] != Double.MIN_VALUE) {
		// bb[cc] = b[kk].intValue();
		// for (int jj=0;jj<A[kk].length;jj++) {
		// AA[cc][jj] = A[cc][jj].intValue();
		// }
		// cc++;
		// }
		// }
		// A = AA;
		// b = bb;
		// }
		// }

		System.out.println("A=");
		for (int i = 0; i < A.length; i++) {
			for (int k = 0; k < A[i].length; k++) {
				System.out.print("" + A[i][k] + " ");
			}
			System.out.println(";");
		}
		System.out.println("b=");
		for (int i = 0; i < b.length; i++) {
			System.out.print("" + b[i] + ";");
		}
		System.out.println();

		System.out.println("c=");
		for (int i = 0; i < c.length; i++) {
			System.out.print("" + c[i] + ",");
		}
		System.out.println();

		System.out.println("lb=");
		for (int i = 0; i < lb.length; i++) {
			System.out.print("" + lb[i] + ",");
		}
		System.out.println();
		System.out.println("ub=");
		for (int i = 0; i < ub.length; i++) {
			System.out.print("" + ub[i] + ",");
		}
		System.out.println();

		// if (true) System.exit(0);
		Integer[] min = new Integer[c.length];
		Integer[] max = new Integer[c.length];
		if (c.length > 0) {
			Integer[] result = cpaMatlab(c, A, b, lb, ub, records);
			for (int k = 0; k < c.length; k++) {
				min[k] = result[k];
				max[k] = result[k + c.length];

				System.out.println("" + k + "th SNP min=" + min[k] + ",max="
						+ max[k]);
			}
		}

		int k = 0;
		int goodRecords = 0;
		System.out.println("records.size=" + records.size());
		for (Enumeration<Record> rds = records.elements(); rds
				.hasMoreElements();) {
			Record rd = rds.nextElement();
			if (!rd.isFree())
				continue;
			rd.count = min[k];
			rd.min = min[k];
			rd.max = max[k];
			if (rd.min == rd.max) {
				goodRecords++;
			}
			System.out.println("min=" + rd.min + ",max=" + rd.max);

			k++;
		}
		int s = rectifySolutions(records, M, snps);
		System.out.println("" + goodRecords + " out of " + k
				+ " freeVariables are constant");
		printRecords(records);

		SOLUTIONS = s;
	}

	private static void printArray(Double[] a) {
		for (int i = 0; i < a.length; i++) {
			System.out.print("" + a[i] + ",");
		}
		System.out.println();
	}

	private static int rectifySolutions(Vector<Record> records,
			Vector<Double[]> M, SNP[] snps) throws Exception {
		Integer[] min = new Integer[records.size()];
		Integer[] max = new Integer[records.size()];
		long[] dim = new long[records.size()];
		for (int k = 0; k < min.length; k++) {
			min[k] = -1;
			max[k] = -1;
			dim[k] = 0;
		}

		long count = 1;
		int k = 0;
		for (Enumeration<Record> rds = records.elements(); rds
				.hasMoreElements();) {
			Record rd = rds.nextElement();
			if (rd.isFree()) {
				dim[k] = rd.max - rd.min + 1;
				count *= dim[k];
			}
			k++;
		}

		if (!isRealSolutionValid(records, M, snps)) {
			throw new Exception(
					"The real solution is not a valid solution of this problem");
		}

		int pass = 0;
		int inconsistent = 0;
		for (long n = 0; n < count; n++) {
			if (n % 10000 == 0)
				System.out.println("verifying " + n
						+ "th possible solution, total=" + count);
			long carry = count;
			long[] selects = new long[records.size()];
			long tmp = n;
			for (k = 0; k < dim.length; k++) {
				// System.out.println("k="+k+",len="+dim.length);
				if (dim[k] > 0) {
					carry = carry / dim[k];
					selects[k] = tmp / carry;
					tmp = tmp % carry;
				}
			}

			k = 0;
			boolean ok = true;
			for (Enumeration<Record> rds = records.elements(); rds
					.hasMoreElements();) {
				Record rd = rds.nextElement();
				if (!rd.isFree()) {
					int r = getRowByLeadIndex(M, k);
					Double[] row = M.get(r);
					double cc = 0;
					for (int i = k + 1; i < row.length - 1; i++) {
						if (row[i] != 0) {
							Record o = records.get(i);
							o.count = records.get(i).min + ((int) selects[i]);
							cc = cc + o.count * (-row[i]);
						}
					}
					cc = cc + row[row.length - 1];
					rd.count = check(cc);
					if (rd.count < 0) {
						ok = false;
						break;
					}
				}

				k++;
			}
			if (!ok)
				continue;

			// if (count == 1) {
			// printRecordsWithRealSolution(snps,records);
			// printMatrix(M);
			// }
			if (isConsistent(records, snps)) {
				pass++;
				for (int j = 0; j < records.size(); j++) {
					if (min[j] == -1) {
						min[j] = records.get(j).count;
					} else {
						int temp = records.get(j).count;
						if (min[j] > temp) {
							min[j] = temp;
						}
					}
					if (max[j] == -1) {
						max[j] = records.get(j).count;
					} else {
						int temp = records.get(j).count;
						if (max[j] < temp) {
							max[j] = temp;
						}
					}
				}
			} else {
				inconsistent++;
			}
		}
		System.out.println("" + pass + " out of " + count
				+ " solutions passed verification");
		System.out.println("" + inconsistent + " out of " + count
				+ " solutions are inconsistent");

		for (k = 0; k < min.length; k++) {
			Record rd = records.get(k);
			rd.min = min[k];
			rd.max = max[k];
			rd.count = min[k];
		}
		for (k = records.size() - 1; k >= 0; k--) {
			Record rd = records.get(k);
			if (rd.count == 0 && rd.min == rd.max) {
				records.remove(k);
			}
		}

		if (!hasRealSolution(snps, records)) {
			printRecordsWithRealSolution(snps, records);
			throw new Exception("Real solution is not in the solutions");
		}
		return pass;
	}

	private static boolean isRealSolutionValid(Vector<Record> records,
			Vector<Double[]> M, SNP[] snps) throws Exception {
		// test if free variables are correct
		int p = 0;
		for (int i = 0; i < records.size(); i++) {
			Record rd = records.get(i);
			if (rd.isFree()) {
				int c = calculateRecordCount(snps, rd);
				if (c < rd.min || c > rd.max) {
					System.out.println("error detected: " + i + "th varialble("
							+ p + "th free variable)=" + c + " out of range ["
							+ rd.min + "," + rd.max + "] record="
							+ rd.toString());
					return false;
				}
				p++;
			}
		}

		boolean ret = true;
		Integer[] solution = calculateASolution(records, M, snps);
		Integer[] realSolution = new Integer[records.size()];
		for (int i = 0; i < records.size(); i++) {
			Record rd = records.get(i);
			realSolution[i] = calculateRecordCount(snps, rd);
			if (realSolution[i].intValue() != solution[i].intValue()) {
				System.out.println("i=" + i + ",real=" + realSolution[i]
						+ ",fake=" + solution[i]);
				ret = false;
			}
		}
		if (!ret) {
			System.out.println("solution=");
			printArray(solution);
			System.out.println("real solution=");
			printArray(realSolution);
		}

		return ret;
	}

	private static void printArray(Integer[] a) {
		for (int i = 0; i < a.length; i++) {
			System.out.print("" + a[i] + ",");
		}
		System.out.println();
	}

	private static Integer[] calculateASolution(Vector<Record> records,
			Vector<Double[]> M, SNP[] snps) throws Exception {
		Integer[] solution = new Integer[records.size()];
		int k = 0;
		for (Enumeration<Record> rds = records.elements(); rds
				.hasMoreElements();) {
			Record rd = rds.nextElement();
			int r = getRowByLeadIndex(M, k);
			if (r == -1) {
				solution[k] = calculateRecordCount(snps, rd);
			}

			k++;
		}
		k = 0;
		for (Enumeration<Record> rds = records.elements(); rds
				.hasMoreElements();) {
			Record rd = rds.nextElement();
			int r = getRowByLeadIndex(M, k);
			if (r != -1) {
				Double[] row = M.get(r);
				double cc = 0;
				for (int i = k + 1; i < row.length - 1; i++) {
					if (row[i] != 0) {
						cc = cc + solution[i] * (-row[i]);
					}
				}
				cc = cc + row[row.length - 1];
				solution[k] = check(cc);
			}

			k++;
		}

		return solution;
	}

	private static void printRecordsWithRealSolution(SNP[] snps,
			Vector<Record> records) {
		int k = 0;
		for (Enumeration<Record> rds = records.elements(); rds
				.hasMoreElements();) {
			Record rd = rds.nextElement();
			System.out.println("" + k + ":r" + calculateRecordCount(snps, rd)
					+ "," + rd.toString());
			k++;
		}
	}

	private static boolean hasRealSolution(SNP[] snps, Vector<Record> records) {
		for (Enumeration<Record> rds = records.elements(); rds
				.hasMoreElements();) {
			Record rd = rds.nextElement();
			int count = calculateRecordCount(snps, rd);
			if (rd.min > count || rd.max < count) {
				return false;
			}
		}
		return true;
	}

	private static boolean isConsistent(Vector<Record> records, SNP[] snps) {
		int total = 0;
		for (Enumeration<Record> rds = records.elements(); rds
				.hasMoreElements();) {
			total += rds.nextElement().count;
		}
		if (total != CASE)
			return false;
		Record rd = records.get(0);
		for (int i = 0; i < rd.record.size() - 1; i++) {
			for (int k = i + 1; k < rd.record.size(); k++) {
				int ID1 = rd.getSNPIDByIndex(i);
				int ID2 = rd.getSNPIDByIndex(k);
				if (!isPairwiseFreqMatch(records, snps, ID1, ID2))
					return false;
			}
		}
		return true;
	}

	private static boolean isPairwiseFreqMatch(Vector<Record> records,
			SNP[] snps, int index1, int index2) {
		RSquare rs = snps[Math.max(index1, index2)].rsquares.get(Math.min(
				index1, index2));
		int cAC = 0;
		int cAD = 0;
		int cBC = 0;
		int cBD = 0;
		for (Enumeration<Record> e = records.elements(); e.hasMoreElements();) {
			Record record = e.nextElement();
			Allele s1 = record.record.get(record.getIndexBySNPID(rs.snp1
					.getID()));
			Allele s2 = record.record.get(record.getIndexBySNPID(rs.snp2
					.getID()));

			if (s1 == Allele.A && s2 == Allele.A) {
				cAC += record.count;
			} else if (s1 == Allele.A && s2 == Allele.T) {
				cAD += record.count;
			} else if (s1 == Allele.T && s2 == Allele.A) {
				cBC += record.count;
			} else if (s1 == Allele.T && s2 == Allele.T) {
				cBD += record.count;
			}
		}

		if (cAC != rs.iAC)
			return false;
		if (cAD != rs.iAD)
			return false;
		if (cBC != rs.iBC)
			return false;
		if (cBD != rs.iBD)
			return false;
		return true;
	}

	private static int getRowByLeadIndex(Vector<Double[]> M, int ldIndex) {
		int k = 0;
		for (Enumeration<Double[]> m = M.elements(); m.hasMoreElements();) {
			Double[] mm = m.nextElement();

			if (mm[ldIndex] == 1) {
				boolean yes = true;
				for (int i = 0; i < ldIndex; i++) {
					if (mm[i] != 0) {
						yes = false;
						break;
					}
				}
				if (yes)
					return k;
			}
			k++;
		}
		return -1;
	}

	private static Integer[] cpaMatlab(Double[] c, Double[][] A, Double[] b,
			Integer[] lb, Integer[] ub, Vector<Record> records) {
		MWNumericArray B = null; /* Stores input values a */
		Object[] result = null; /* Stores the result */
		Integer[] d = new Integer[2 * c.length];
		try {
			if (bnb == null) {
				bnb = new Matcher();
			}

			int[] dims = { b.length + 3, c.length + 1 };
			B = MWNumericArray.newInstance(dims, MWClassID.DOUBLE,
					MWComplexity.REAL);
			int[] index = { 1, 1 };
			for (int col = 0; col < dims[1] - 1; col++) {
				index[0] = 1;
				index[1] = col + 1;
				B.set(index, c[col]);
			}
			for (int col = 0; col < dims[1] - 1; col++) {
				index[0] = 2;
				index[1] = col + 1;
				B.set(index, lb[col]);
			}
			for (int col = 0; col < dims[1] - 1; col++) {
				index[0] = 3;
				index[1] = col + 1;
				B.set(index, ub[col]);
			}

			index[0] = 1;
			index[1] = dims[1];
			B.set(index, 100000); // MatIter
			index[0] = 2;
			B.set(index, 3600); // MaxTime
			index[0] = 3;
			B.set(index, CASE);

			for (int row = 0; row < dims[0] - 3; row++) {
				for (int col = 0; col < dims[1] - 1; col++) {
					index[0] = row + 4;
					index[1] = col + 1;
					B.set(index, A[row][col]);
				}
			}
			for (int row = 0; row < dims[0] - 3; row++) {
				index[0] = row + 4;
				index[1] = dims[1];
				B.set(index, b[row]);
			}

			// System.out.println(B);
			//result = bnb.linprog_bin_ex(1, B);
			// System.out.println("Result:");
			MWNumericArray re = (MWNumericArray) result[0];
			for (int k = 0; k < 2 * c.length; k++) {
				d[k] = (int) Math.round(re.getDouble(k + 1));
			}
		} catch (Exception e) {
			System.out.println("Exception: " + e.toString());
			System.exit(0);
		}

		finally {
			/* Free native resources */
			MWArray.disposeArray(B);
			MWArray.disposeArray(result);
		}
		return d;
	}

	protected void finalize() {
		if (bnb != null) {
			bnb.dispose();
		}
	}

	private static boolean isRealSolutionValidBeforeReduce(
			Vector<Record> records, Vector<Double[]> M, SNP[] snps) {
		Integer[] realSolution = new Integer[records.size()];
		for (int k = 0; k < records.size(); k++) {
			Record rd = records.get(k);
			realSolution[k] = calculateRecordCount(snps, rd);
			if (rd.getLB() > realSolution[k] || rd.getUB() < realSolution[k]) {
				System.out.println("error detected: count should be in ["
						+ rd.getLB() + "," + rd.getUB() + "], record="
						+ rd.toString() + ",count=" + realSolution[k]);
				return false;
			}
		}
		for (Enumeration<Double[]> rows = M.elements(); rows.hasMoreElements();) {
			Double[] row = rows.nextElement();
			double sum = 0;
			for (int k = 0; k < row.length - 1; k++) {
				sum += row[k] * realSolution[k];
			}
			if (Math.abs(sum - row[row.length - 1]) > MAXDEVIATION) {
				System.out.println("sum=" + sum);
				printArray(row);
				return false;
			}
		}
		return true;
	}

	private static void printMatrix(Vector<Double[]> M) {
		if (M.size() == 0)
			return;

		int rowCount = M.size();
		int colCount = M.get(0).length;
		System.out.println("Matrix(" + rowCount + "," + colCount + ")");
		for (int r = 0; r < rowCount; r++) {
			int c = 0;
			for (; c < colCount - 1; c++) {
				System.out.print("" + M.get(r)[c] + ",");
			}
			System.out.print("" + M.get(r)[c] + ";");
			if (r % 4 == 0)
				System.out.println();
		}
	}

	private static int check(double d) throws Exception {
		int i = (int) Math.round(d);
		if (Math.abs(i - d) > MAXDEVIATION) {
			throw new Exception("error detected: check failure, d=" + d);
		}
		return i;
	}

	private static void removeColFromMatrix(Vector<Double[]> M, int c) {
		int i = 0;
		for (Enumeration<Double[]> e = M.elements(); e.hasMoreElements();) {
			Double[] a = e.nextElement();
			Double[] b = new Double[a.length - 1];
			for (int k = 0; k < a.length; k++) {
				if (k < c) {
					b[k] = a[k];
				} else if (k > c) {
					b[k - 1] = a[k];
				}
			}
			M.set(i, b);

			i++;
		}
	}
	
	private static void computeSingleSNPFrequencies(SNP[] snps, int i) {
		Entries[] entries = new Entries[snps.length];
		for (int k = 0; k < snps.length; k++) {
			entries[k] = new Entries();
			snps[k].entries = entries[k];
		}
		Vector<Integer> bestIndexes = new Vector<Integer>();

		for (int k = 0; k < snps.length; k++) {
			if (k == i) {
				entries[k].appendWithoutDuplicate(new Entry(snps[k]
						.getCaseAFreq(), snps[k].getCaseTFreq()));
				continue;
			}
			if (snps[k].isSigSNP()) {
				entries[k].appendWithoutDuplicate(new Entry(snps[k]
						.getCaseAFreq(), snps[k].getCaseTFreq()));
				bestIndexes.add(k);
				continue;
			}
			boolean pass = false;

			long c1, c2, c3;
			long p1 = snps[i].getCaseAFreq();
			long p2 = snps[i].getCaseTFreq();
			for (c1 = 0; c1 <= p1; c1++) {
				// System.out.println("step 1");
				for (c2 = 0; c2 <= p2; c2++) {
					// System.out.println("step 2");
					long q1 = c1 + c2;
					long q2 = (p1 - c1) + (p2 - c2);
					// deal with rsquare
					int t = k > i ? k : i;
					int index = k > i ? i : k;
					double rsquare = snps[t].rsquares.get(index).rsquare;
					if (rsquare == -1)
						continue;
					double minRSQUARE = (int) (rsquare * RSQUARE_PRECISION)
							* 1.0 / RSQUARE_PRECISION;
					double maxRSQUARE = ((int) (rsquare * RSQUARE_PRECISION) + 1)
							* 1.0 / RSQUARE_PRECISION;

					double div = p1 * p2 * q1 * q2;
					if (div == 0)
						continue;
					rsquare = Math.pow(c1 * (p1 + p2) - p1 * q1, 2) / div;
					if (rsquare >= minRSQUARE && rsquare <= maxRSQUARE) {
						// System.out.println("step 3");
						double pvalue = snps[k].getPvalue();
						if (pvalue == -1)
							continue;
						double minPVALUE = (int) (pvalue * PVALUE_PRECISION)
								* 1.0 / PVALUE_PRECISION;
						double maxPVALUE = ((int) (pvalue * PVALUE_PRECISION) + 1)
								* 1.0 / PVALUE_PRECISION;

						long ctrlP = snps[k].getControlPopulation();
						for (c3 = 0; c3 <= ctrlP; c3++) {
							// System.out.println("step 4");
							long x11 = q1, x12 = c3, x21 = q2, x22 = ctrlP - c3;
							Allele maf = snps[k].getMAF();
							if (maf == Allele.A && (x11 > x21)) {
								continue;
							} else if (maf == Allele.T && (x11 < x21)) {
								continue;
							}
							long pp1 = x11 + x12, pp2 = x21 + x22, qq1 = x11
									+ x21, qq2 = x12 + x22;
							div = pp1 * qq1 * pp2 * qq2;
							if (div == 0)
								continue;
							// pvalue = Math.pow((x11*x22-x12*x21),
							// 2)*(pp1+pp2)/div;
							pvalue = Math.pow(x11 * (pp1 + pp2) - pp1 * qq1, 2)
									/ div;
							if (pvalue >= minPVALUE && pvalue <= maxPVALUE) {
								entries[k].appendWithoutDuplicate(new Entry(
										x11, x21));
								if (snps[k].getCaseAFreq() != x11
										|| snps[k].getCaseTFreq() != x21) {
									System.out
											.println("bad SNP" + k + " caseA="
													+ x11 + " caseT=" + x21
													+ " ctrlA=" + x12
													+ " ctrlT=" + x22);
									System.out.println("         caseA="
											+ snps[k].getCaseAFreq()
											+ " caseT="
											+ snps[k].getCaseTFreq()
											+ " ctrlA="
											+ snps[k].getControlAFreq()
											+ " ctrlT="
											+ snps[k].getControlTFreq());
									// System.exit(0);
								} else {
									pass = true;
								}
							}
						}
					}
				}
			}

			if (entries[k].getEntryNum() == 1) {
				bestIndexes.add(k);
			}

			if (pass == false) {
				System.out.println("SNP" + k
						+ " failed the frequency examination");
				System.out.println("         caseA=" + snps[k].getCaseAFreq()
						+ " caseT=" + snps[k].getCaseTFreq() + " ctrlA="
						+ snps[k].getControlAFreq() + " ctrlT="
						+ snps[k].getControlTFreq());
			}
		}

		int goodSNP = 0, lastGoodSNP = 0;
		while (goodSNP < snps.length) {
			goodSNP = 0;
			for (int k = 0; k < snps.length; k++) {
				// if (k == i) continue;
				int count = entries[k].getEntryNum();
				if (count <= 1) {
					goodSNP++;
					continue;
				}

				for (Enumeration<Integer> b = bestIndexes.elements(); b
						.hasMoreElements();) {
					Integer indexBest = b.nextElement();

					// System.out.println(entries[indexBest].getEntryNum());
					Entry bestEntry = entries[indexBest].entries.get(0);
					long x11 = 0, x21 = 0, x12 = 0, x22 = 0, p1 = 0, p2 = 0, q1 = 0, q2 = 0;
					p1 = bestEntry.getCaseA();
					p2 = bestEntry.getCaseT();
					boolean removable = true;
					for (Enumeration<Entry> e = entries[k].entries.elements(); e
							.hasMoreElements();) {
						Entry entry = e.nextElement();

						q1 = entry.getCaseA();
						q2 = entry.getCaseT();
						long min = p1 > q1 ? q1 : p1;
						for (x11 = 0; x11 <= min; x11++) {
							x12 = p1 - x11;
							x21 = q1 - x11;
							x22 = p2 - x21;

							if (x22 != q2 - x12)
								continue;
							if (x12 != q2 - x22)
								continue;
							if (x21 != p2 - x22)
								continue;
							if (x11 < 0 || x21 < 0 || x12 < 0 || x22 < 0) {
								continue;
							}
							Allele maf = snps[k].getMAF();
							if (maf == Allele.A && q1 > q2) {
								continue;
							} else if (maf == Allele.T && q1 < q2) {
								continue;
							}

							int t = k > indexBest ? k : indexBest;
							int index = k > indexBest ? indexBest : k;
							double rsquare = snps[t].rsquares.get(index).rsquare;
							double minRSQUARE = (int) (rsquare * RSQUARE_PRECISION)
									* 1.0 / RSQUARE_PRECISION;
							double maxRSQUARE = ((int) (rsquare * RSQUARE_PRECISION) + 1)
									* 1.0 / RSQUARE_PRECISION;
							double temp = RSquare.calculateRSQUARE(x11, x12,
									x21, x22);

							if (temp >= minRSQUARE && temp <= maxRSQUARE) {
								removable = false;
							}
						}

						if (removable) {
							count--;
							System.out.println("remove bad SNP" + k + " caseA="
									+ q1 + " caseT=" + q2);
							// if (q1 == snps[k].getCaseAFreq() && q2 ==
							// snps[k].getCaseTFreq()) {
							// System.out.println("victim deleted");
							// System.exit(0);
							// }
							entries[k].entries.remove(entry);
						}
					}
				}
				if (count > 1) {
					System.out.println("SNP" + k + " is really bad "
							+ entries[k].getEntryNum());
					for (int n = 0; n < entries[k].getEntryNum(); n++) {
						Entry ee = entries[k].entries.get(n);
						System.out.println(ee.toString());
					}
					System.out.println("real iCaseA=" + snps[k].getCaseAFreq()
							+ ",iCaseT=" + snps[k].getCaseTFreq());
					System.out.println("bestIndexes=" + bestIndexes.size());
				} else if (count == 1) {
					bestIndexes.add(k);
					// System.out.println("i="+i+",A="+snps[i].getCaseAFreq());
					// System.out.println("k="+k+",A="+snps[k].getCaseAFreq());
					// System.out.println("r2="+getRSquare(snps,i,k).rsquare);
					// System.out.println("SNP "+k+" is recovered");
					// System.exit(0); //test
				}
			}
			if (lastGoodSNP == goodSNP)
				break;
			lastGoodSNP = goodSNP;
		}

		boolean removed = true;
		while (removed) {
			removed = false;

			for (int k = 0; k < snps.length; k++) {
				if (entries[k].getEntryNum() > 1) {
					for (Enumeration<Entry> ee = entries[k].entries.elements(); ee
							.hasMoreElements();) {
						Entry e = ee.nextElement();
						for (int m = 0; m < snps.length; m++) {
							if (m == k)
								continue;
							if (entries[m].entries.size() == 0)
								continue;

							boolean pass = false;
							for (Enumeration<Entry> eee = entries[m].entries
									.elements(); eee.hasMoreElements();) {
								Entry f = eee.nextElement();
								RSquare rs = snps[Math.max(k, m)].rsquares
										.get(Math.min(k, m));
								double rsquare = rs.rsquare;
								int p = rs.snp1.getID();
								int q = rs.snp2.getID();
								long p1, p2, q1, q2;
								if (p == k) {
									p1 = e.getCaseA();
									p2 = e.getCaseT();
									q1 = f.getCaseA();
									q2 = f.getCaseT();
								} else {
									p1 = f.getCaseA();
									p2 = f.getCaseT();
									q1 = e.getCaseA();
									q2 = e.getCaseT();
								}
								long min = p1 > q1 ? q1 : p1;
								for (long x11 = 0; x11 <= min; x11++) {
									long x12 = p1 - x11;
									long x21 = q1 - x11;
									long x22 = p2 - x21;

									if (x22 != q2 - x12)
										continue;
									if (x12 != q2 - x22)
										continue;
									if (x21 != p2 - x22)
										continue;
									if (x11 < 0 || x21 < 0 || x12 < 0
											|| x22 < 0) {
										continue;
									}
									Allele maf = snps[p].getMAF();
									if (maf == Allele.A && p1 > p2) {
										continue;
									} else if (maf == Allele.T && p1 < p2) {
										continue;
									}
									maf = snps[q].getMAF();
									if (maf == Allele.A && q1 > q2) {
										continue;
									} else if (maf == Allele.T && q1 < q2) {
										continue;
									}

									double minRSQUARE = (int) (rsquare * RSQUARE_PRECISION)
											* 1.0 / RSQUARE_PRECISION;
									double maxRSQUARE = ((int) (rsquare * RSQUARE_PRECISION) + 1)
											* 1.0 / RSQUARE_PRECISION;
									double temp = RSquare.calculateRSQUARE(x11,
											x12, x21, x22);

									if (temp >= minRSQUARE
											&& temp <= maxRSQUARE) {
										pass = true;
										break;
									}
								}
								if (pass)
									break;
							}
							if (!pass) {
								System.out.println("remove bad SNP" + k
										+ " caseA=" + e.getCaseA() + " caseT="
										+ e.getCaseT());
								if (e.getCaseA() == snps[k].getCaseAFreq()
										&& e.getCaseT() == snps[k]
												.getCaseTFreq()) {
									System.out.println("victim deleted");
									System.exit(0);
								}
								entries[k].entries.remove(e);
								removed = true;
								break;
							}
						}
					}
				}
			}
		}

		int badSNP = 0;
		goodSNP = 0;
		for (int k = 0; k < snps.length; k++) {
			if (entries[k].getEntryNum() > 1) {
				badSNP++;
				snps[k].setBadSNP(true);
				System.out.println("SNP" + k + " is really bad "
						+ entries[k].getEntryNum());
				for (int n = 0; n < entries[k].getEntryNum(); n++) {
					Entry ee = entries[k].entries.get(n);
					System.out.println(ee.toString());
				}
				System.out.println("real iCaseA=" + snps[k].getCaseAFreq()
						+ ",iCaseT=" + snps[k].getCaseTFreq());
			} else if (entries[k].getEntryNum() == 1) {
				goodSNP++;
				System.out.println("SNP" + k + " is really good ");
			} else if (entries[k].getEntryNum() == 0) {
				System.out.println("SNP" + k + " is really disposable ");
				System.out.println("real iCaseA=" + snps[k].getCaseAFreq()
						+ ",iCaseT=" + snps[k].getCaseTFreq());
				// System.exit(0);
			}
		}

		System.out.println("" + goodSNP + " out of " + snps.length
				+ " are good");
		System.out.println("" + badSNP + " out of " + snps.length + " are bad");

		// System.exit(0);
	}

	private static SNP[] readSNPsFromFastaFile(String fileName, int cutSNP, int cutRec) {
		Vector<Allele[]> alleles = new Vector<Allele[]>();

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
					strLine = br.readLine();
					if (cutSNP > 0 && cutSNP < strLine.length())
						strLine = strLine.substring(0, cutSNP);

					Allele[] temp = getAlleleFromStr(strLine);
					Allele[] rewrite;
					if (firstLine == null) {
						firstLine = temp;
						secondLine = new Allele[temp.length];
						for (int m = 0; m < temp.length; m++)
							secondLine[m] = firstLine[m];

						firstLine_new = new Allele[firstLine.length];
						for (int i = 0; i < firstLine.length; i++)
							firstLine_new[i] = Allele.A;
						rewrite = firstLine_new;
					} else {
						// rewrite the string
						rewrite = new Allele[temp.length];
						for (int i = 0; i < temp.length; i++) {
							if (temp[i] == firstLine[i]) {
								rewrite[i] = Allele.A;
							} else {
								secondLine[i] = temp[i];
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
			all_snps[i] = new SNP("snp" + i, i, all.get(i), GWASAttack.CASE);
		}
		
		computeRSquares(all_snps);

		return all_snps;
	}

	private static Vector<RSquare> computeRSquares(SNP[] snps) {
		Vector<RSquare> all_rsquares = new Vector<RSquare>();
		for (int i = 0; i < snps.length - 1; i++) {
			for (int j = i + 1; j < snps.length; j++) {
				RSquare rsquare = new RSquare(snps[i], snps[j]);
				snps[i].appendRSquare(rsquare);
				snps[j].appendRSquare(rsquare);
				all_rsquares.add(rsquare);
			}
		}
		
		return all_rsquares;
	}

	private static Allele[] getAlleleFromStr(String str) {
		Allele[] record = new Allele[str.length()];
		for (int i = 0; i < record.length; i++) {
			record[i] = getAlleleFromChar(str.charAt(i));
		}
		return record;
	}

	private static Allele getAlleleFromChar(char c) {
		if (c == 'A' || c == 'a')
			return Allele.A;
		else if (c == 'T' || c == 't')
			return Allele.T;
		else if (c == 'G' || c == 'g')
			return Allele.G;
		else if (c == 'C' || c == 'c')
			return Allele.C;
		return Allele.U;
	}

	private static int calculateRecordCount(SNP[] snps, Record record) {
		Vector<Allele> alleles = record.record;

		int count = 0;
		for (int i = 0; i < CASE; i++) {
			boolean equal = true;
			int n = 0;
			for (Enumeration<Segment> segs = record.segments.elements(); segs
					.hasMoreElements();) {
				Segment seg = segs.nextElement();
				int tmp = n;
				for (int k = seg.begin; k < seg.begin + seg.length; k++) {
					if (snps[k].getAlleles().get(i) != alleles.get(tmp)) {
						equal = false;
						break;
					}
					tmp++;
				}
				if (!equal)
					break;
				n += seg.length;
			}
			if (equal)
				count++;
		}
		return count;
	}

	private static char getSignFromDouble(double r) {
		if (r > 0)
			return '+';
		else if (r == 0)
			return '0';
		else
			return '-';
	}

	public static void saveToFile(String fileName, Vector<Record> records)
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
		obj_out.writeObject(records);
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

	public static Vector<Record> readFromFile(String fileName)
			throws IOException, ClassNotFoundException {
		// Read from disk using FileInputStream.
		FileInputStream f_in = new FileInputStream(fileName);

		// Read object using ObjectInputStream.
		ObjectInputStream obj_in = new ObjectInputStream(f_in);

		// Read an object.
		Object obj = obj_in.readObject();

		// Is the object that you read in, say, an instance
		// of the Vector class?
		if (obj instanceof Vector) {
			// Cast object to a Vector
			Vector<Record> records = (Vector<Record>) obj;
			return records;
		} else { // ... the object is some other type ...
			return null;
		}
	}

	public static SNP[] readSNPsFromFile(String fileName) throws IOException,
			ClassNotFoundException {
		// Read from disk using FileInputStream.
		FileInputStream f_in = new FileInputStream(fileName);

		// Read object using ObjectInputStream.
		ObjectInputStream obj_in = new ObjectInputStream(f_in);

		// Read an object.
		Object obj = obj_in.readObject();

		// Is the object that you read in, say, an instance
		// of the Vector class?
		if (obj instanceof SNP[]) {
			// Cast object to a Vector
			SNP[] snps = (SNP[]) obj;
			return snps;
		} else { // ... the object is some other type ...
			return null;
		}
	}
}
