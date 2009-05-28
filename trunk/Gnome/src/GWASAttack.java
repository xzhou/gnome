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
	private Network net; 	//constraint net
	private Solver solver;	//constraint solver
	private Solution solution;	//one solution
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
		} 
		else if(approachChoice == 3)
		{
			//continue from 16%
			//please first run with method 1 and then continue
			SNP[] snps = readSNPsFromFile("snps.data");
			SignRecover.propagateSigns4(snps);
		}
		else {
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
					RSquare rs = UtilityFunctions.getRSquare(snps, index1, index2);

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

					double i1 = computeConfidenceValue(snps, rs, true);
					double i2 = computeConfidenceValue(snps, rs, false);
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
				UtilityFunctions.signRecoverRate(snps);
				System.out.println("undecidable=" + undecidable);
			}

			if (!propagating)
				break;
			propagating = false;
		}
		System.exit(0);
	}
	
	
	
	/*
	 * propagateSignOfRSquare
	 * given a sign recovered snps, propagate the signs further
	 */
	private static void propagateSignOfRSquare(SNP[] snps) throws Exception {
		boolean propagating = false;
		double c_i = 0.5;
		
		//repeatedly propagate
		while (true) {
			for (int index1 = 0; index1 < snps.length - 1; index1++) {
				for (int index2 = index1 + 1; index2 < snps.length; index2++) {
					RSquare rSquare = UtilityFunctions.getRSquare(snps, index1, index2);

					//if a r^2 is recovered, find another snp, 
					if (rSquare.isSignRecovered()) {
						//for any 3rd snp, and recover more signs
						for (int i = 0; i < snps.length; i++) {
							if (i == index1 || i == index2)
								continue;
							RSquare rs12 = UtilityFunctions.getRSquare(snps, index1, i);
							RSquare rs13 = UtilityFunctions.getRSquare(snps, index2, i);
							if (rs12.isSignRecovered() && rs13.isSignRecovered())
								continue;

							Vector<Double> pAB12 = new Vector<Double>();
							Vector<Double> pAB13 = new Vector<Double>();

							if (rs12.isSignRecovered()) {
								pAB12.add(rs12.pAB);
							} else {
								pAB12.add(rs12.v1);
								pAB12.add(rs12.v2);
							}

							if (rs13.isSignRecovered()) {
								pAB13.add(rs13.pAB);
							} else {
								pAB13.add(rs13.v1);
								pAB13.add(rs13.v2);
							}

							double deviation = 1e-2;
							//Validate each combination of p 
							for (int i1 = 0; i1 < pAB12.size(); i1++) {
								for (int i2 = 0; i2 < pAB13.size(); i2++) {
									
									boolean god = false;
									double p23 = rSquare.pAB, 
											p12 = pAB12.get(i1), 
											p13 = pAB13.get(i2);
									
									//?? why larger than one
									if (p12 > 1)
										p12 -= 1;
									if (p13 > 1)
										p13 -= 1;
									
									
									//god means the original sign
									if (Math.abs(p12 - rs12.pAB) < deviation
											&& Math.abs(p13 - rs13.pAB) < deviation) {
										god = true;
									}

									/*
									 * r sqauare sign is god
									 */
									
									/*
									 * angel is the computed sign
									 */
									
									boolean angel = UtilityFunctions.isMatchByMatlab(false,
											snps[i].getpA(), 
											snps[index1].getpA(), 
											snps[index2].getpA(), 
											p23, 
											p12, 
											p13);
									
									//this is consistent
									if (angel) {
										if (pAB12.get(i1) < 1)
											pAB12.set(i1, p12 + 1);
										if (pAB13.get(i2) < 1)
											pAB13.set(i2, p13 + 1);
									}
									//not consistent, print some information
									else if (god) {
										UtilityFunctions.isMatchByMatlab(true, 
												snps[i].getpA(),
												snps[index1].getpA(),
												snps[index2].getpA(), 
												p23, 
												p12,
												p13);
										// throw new
										// Exception("angel should be true because god is");
										System.err.println("angel should be true because god is");
									}
									// System.out.println("angel="+angel);
								}
							}

							
							//check recovered signs
							for (int i1 = 0; i1 < pAB12.size(); i1++) {
								double p12 = pAB12.get(i1);
								if (p12 < 1) {
									if (p12 == rs12.v1) {
										double interval = computeConfidenceValue(
												snps, rs12, false);
										if (interval >= c_i) {
											rs12.setSignRecovered(false, true);
											propagating = true;
										}
									} else if (p12 == rs12.v2) {
										double interval = computeConfidenceValue(
												snps, rs12, true);
										if (interval >= c_i) {
											rs12.setSignRecovered(true, true);
											propagating = true;
										}
									} else {
										throw new Exception(
												"error detected: p12=" + p12
														+ ",v1=" + rs12.v1
														+ ",v2=" + rs12.v2);
									}
								}
							}
							
							
							for (int i2 = 0; i2 < pAB13.size(); i2++) {
								double p13 = pAB13.get(i2);
								if (p13 < 1) {
									if (p13 == rs13.v1) {
										double interval = computeConfidenceValue(
												snps, rs13, false);
										if (interval >= c_i) {
											rs13.setSignRecovered(false, true);
											propagating = true;
										}
									} else if (p13 == rs13.v2) {
										double interval = computeConfidenceValue(
												snps, rs13, true);
										if (interval >= c_i) {
											rs13.setSignRecovered(true, true);
											propagating = true;
										}
									} else {
										throw new Exception(
												"error detected: p13=" + p13
														+ ",v1=" + rs13.v1
														+ ",v2=" + rs13.v2);
									}
								}
							}
						}
						UtilityFunctions.signRecoverRate(snps);
						UtilityFunctions.saveSNPsToFile("snps.dat", snps);
					}
				}
			}
			if (!propagating)
				break;
			UtilityFunctions.saveSNPsToFile("snps.dat", snps);
			System.out.println("propagating");
			propagating = false;
		}
	}
	


	//


	//compatible with how many snps
	public static double computeConfidenceValue(SNP[] snps, RSquare rs,
			boolean pos)
	{
		int index1 = rs.snp1.getID();
		int index2 = rs.snp2.getID();
		double p23 = rs.v2;
		
		//if it is positive
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
					
					angel = UtilityFunctions.isMatchByMatlab(false, snps[i].getpA(),
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
	//Simply use the signs of Hapmap as the signs for the real data
	private static void computeSignOfRSquareFromHapmap(SNP[] snps,
			SNP[] hapmap_snps) {
		int rec = 0, rev_rec = 0, true_rec = 0, true_rev_rec = 0;
		int hap_signs = 0, zeros = 0;
		double meanR = 0, meanNR = 0;
		System.out.println("SNPi SNPj r(recovered) r(real) r2");
		for (int index1 = 0; index1 < snps.length - 1; index1++) {
			for (int index2 = index1 + 1; index2 < snps.length; index2++) {
				RSquare rs = UtilityFunctions.getRSquare(snps, index1, index2); //
				RSquare hap_rs = UtilityFunctions.getRSquare(hapmap_snps, index1, index2);

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

	//6%
	private static void computeSignOfRSquare(SNP[] snps) throws Exception {
		int i = 0, count1 = 0, count2 = 0;
		for (int index1 = 0; index1 < snps.length - 1; index1++) {
			for (int index2 = index1 + 1; index2 < snps.length; index2++) {
				RSquare rs = UtilityFunctions.getRSquare(snps, index1, index2);
				double pA = rs.getpA();
				double pB = rs.getpB();
				double pa = 1 - pA;
				double pb = 1 - pB;
				double r = Math.sqrt(rs.rsquare * pA * pa * pB * pb);
				double v1 = r + pA * pB;
				double v2 = -r + pA * pB;
				rs.v1 = v1;
				rs.v2 = v2;

				double v = -1;
				double deviation = 0.01;
				
				//if v1 is invalid
				if (((v1 - pA) > deviation || (v1 - pB) > deviation)) {
					i++;
					v = v2;
				}
				
				//if v2 is in valid
				else if ((v2 - pA) > deviation || (v2 - pB) > deviation) {
					i++;
					v = v1;
				}
				
				else {
					if (pb - pA + v1 < -deviation) {
						i++;
						v = v2;
					} else if (pb - pA + v2 < -deviation) {
						i++;
						v = v1;
					}
				}
				
				//invalid
				if (v < 0)	continue;
				
				//invalid
				if ((v == v1 && rs.r < 0) || (v == v2 && rs.r > 0)) {
					System.err.println("incorect sign");
				}

				//sign recovered 
				rs.pAB = v;
				rs.setSignRecovered(true);
			}
		}
		UtilityFunctions.signRecoverRate(snps);
		System.out.println("i=" + i + ",count1=" + count1 + ",count2=" + count2);
	}

	//read us.freq.txt
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

			RSquare rsquare = UtilityFunctions.getRSquare(snps, i, j);

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



	private static void printRecords(Vector<Record> records) {
		int k = 0;
		for (Enumeration<Record> rds = records.elements(); rds
				.hasMoreElements();) {
			Record rd = rds.nextElement();
			System.out.println("" + k + ":" + rd.toString());
			k++;
		}
	}


	private static void printArray(Double[] a) {
		for (int i = 0; i < a.length; i++) {
			System.out.print("" + a[i] + ",");
		}
		System.out.println();
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
