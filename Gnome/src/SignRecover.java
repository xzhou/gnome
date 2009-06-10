import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.PrintStream;
import java.util.Vector;

import matcher.Matcher;

import com.mathworks.toolbox.javabuilder.MWArray;
import com.mathworks.toolbox.javabuilder.MWClassID;
import com.mathworks.toolbox.javabuilder.MWComplexity;
import com.mathworks.toolbox.javabuilder.MWNumericArray;
import lPSolver.LPSolver;


/**
 * 
 */

/**
 * @author xzhou
 * 
 */
public class SignRecover {
	
	//from the recovered signs, inferring other signs using 4 snps together
	/**
	 * using linear programming to propagate the sign
	 * 
	 * @param snps The snps of all data
	 * @return void
	 * 
	 */
	
	private static LPSolver lpSolver;
	
	public SignRecover(){
		try{
			lpSolver = new LPSolver();
		}catch(Exception e)
		{
			lpSolver = null;
		}
	}
	
	
	public static void propagateSigns4(SNP[] snps) throws Exception
	{
		boolean propagating = false;
		double confidence_level = 0.5;
		
		PrintStream matchResultOutStream = new PrintStream(
		        new BufferedOutputStream(new FileOutputStream(
		        		new File("LPresult.txt"))), true);
		
		PrintStream analysisFile = new PrintStream(
		        new BufferedOutputStream(new FileOutputStream(
				        new File("analysis.txt"))), true);
		
		PrintStream fpFile = new PrintStream(
		        new BufferedOutputStream(new FileOutputStream(
				        new File("fp.txt"))), true);
		
		PrintStream fnFile = new PrintStream(
		        new BufferedOutputStream(new FileOutputStream(
				        new File("fn.txt"))), true);
		
		long nTotalConbination = 0;
		long fp = 0;	//false positive
		long fn = 0;	//false negative
		
		
		
		long nRecovered = 0;
		long nFalseRecovered = 0;
		
		
		//repeatedly recover sign until no more sign can be recovered
		while(true)
		{
			//find 3 snps whose sign is recovered
			//partial order to reduce the complex
			for(int index2 = 0; index2 < snps.length - 1; index2++)
				for (int index3 = index2 + 1; index3 < snps.length - 1; index3++)
					for(int index4 = index3 + 1; index4 < snps.length - 1; index4++)
					{
						//remove different
						if(index2==index3 || index2 == index4 || index3 == index4)
							continue;
						//System.out.println("(" + index2 + "," + index3 + "," +index4 + ")");
						//get rsquare
						RSquare r2_23 = UtilityFunctions.getRSquare(snps, index2, index3);
						RSquare r2_24 = UtilityFunctions.getRSquare(snps, index2, index4);
						RSquare r2_34 = UtilityFunctions.getRSquare(snps, index3, index4);
						
						//we find 3 snps whose sign was recovered
						if(r2_23.isSignRecovered() && r2_24.isSignRecovered() && r2_34.isSignRecovered())
						{
							//find the fourth index
							for(int index1 = 0; index1 < snps.length-1; index1++)
							{
								if(index1 == index2 || index1 == index3 || index1 == index4)
									continue;

								RSquare r2_12 = UtilityFunctions.getRSquare(snps, index1, index2);
								RSquare r2_13 = UtilityFunctions.getRSquare(snps, index1, index3);
								RSquare r2_14 = UtilityFunctions.getRSquare(snps, index1, index4);

								if(r2_12.isSignRecovered() && r2_13.isSignRecovered() && r2_14.isSignRecovered())
								{
									continue;
								}
								matchResultOutStream.print("(" + index1 + "," +index2 + "," + index3 + "," +index4 + ")");	
								//output the correct sign
								int x = -1;	//correct sign of 12
								int y = -1; //correct sign of 13
								int z = -1; //correct sign of 14

								
								if(r2_12.sign == '+')	//record the correct v
									x = 0;
								else
									x = 1;
								
								if(r2_13.sign == '+')
									y = 0;
								else
									y = 1;
								
								if(r2_14.sign == '+')
									z = 0;
								else
									z = 1;
								
								Vector<Double> pAB12 = new Vector<Double>();
								Vector<Double> pAB13 = new Vector<Double>();
								Vector<Double> pAB14 = new Vector<Double>();
								
								if(r2_12.isSignRecovered())
									pAB12.add(r2_12.pAB);	//use the correct value
								else 
								{
									pAB12.add(r2_12.v1);
									pAB12.add(r2_12.v2);
								}
								
								if(r2_13.isSignRecovered())
									pAB13.add(r2_13.pAB);
								else
								{
									pAB13.add(r2_13.v1);
									pAB13.add(r2_13.v2);
								}
								
								if(r2_14.isSignRecovered())
									pAB14.add(r2_14.pAB);
								else
								{
									pAB14.add(r2_14.v1);
									pAB14.add(r2_14.v2);
								}
								
								matchResultOutStream.print("\t("+x+y+z+")\t\t");
								double deviation = 1e-2;
								
								//
								for(int i12 = 0; i12 < pAB12.size(); i12++)
								{
									for(int i13 = 0; i13 < pAB13.size(); i13++)
									{
										for(int i14 = 0; i14 < pAB14.size(); i14++)
										{
											boolean oracle = false;
											
											//for the three know pairwise frequency
											double p23 = r2_23.pAB;
											double p24 = r2_24.pAB;
											double p34 = r2_34.pAB;
											
											//for each new
											double p12 = pAB12.get(i12);
											double p13 = pAB13.get(i13);
											double p14 = pAB14.get(i14);
											
											//TODO remove the following code 
//											p23 = r2_23.getpAB();
//											p24 = r2_24.getpAB();
//											p34 = r2_34.getpAB();
//											p12 = r2_12.getpAB();
//											p13 = r2_13.getpAB();
//											p14 = r2_14.getpAB();
											
											//cancel the marker
											if(p12 > 1)
												p12 -= 1;
											if(p13 > 1)
												p13 -= 1;
											if(p14 > 1)
												p14 -= 1;
											
											//use the real data to validate
											if(Math.abs(p12 - r2_12.pAB) < deviation
													&& Math.abs(p13 - r2_13.pAB) < deviation
													&& Math.abs(p14 - r2_14.pAB) < deviation)
											{
												oracle = true;
											}								
											//check consistency
											boolean consistent = isConsistent4(false,
													snps[index1].getpA(),
													snps[index2].getpA(),
													snps[index3].getpA(),
													snps[index4].getpA(),
													p12, p13, p14,
													p23, p24, p34);
											
											//mark as consistent
											if(consistent)
											{	
												if(pAB12.get(i12) < 1)
													pAB12.set(i12, p12+1);
												if(pAB13.get(i13) < 1)
													pAB13.set(i13, p13 + 1);
												if(pAB14.get(i14) < 1)
													pAB14.set(i14, p14 + 1);
											}
											
											else if(oracle)
											{
												/*
												isConsistent4(true,
													snps[index1].getpA(),
													snps[index2].getpA(),
													snps[index3].getpA(),
													snps[index4].getpA(),
													p12, p13, p14,
													p23, p24, p34);
												*/
												//System.out.println();
											}
											matchResultOutStream.print("("+i12+i13+i14+":"+consistent+")\t");
											
											nTotalConbination ++;
											
											//false negative
											if(i12 == x && i13 == y && i14 == z && !consistent)
											{
												fn ++;
												fnFile.print("(" + index1 + "," +index2 + "," + index3 + "," +index4 + ")");
												fnFile.print("\t("+x+y+z+")\t");
												fnFile.print("("+i12+i13+i14+":"+consistent+")\t");
												fnFile.print("(" + snps[index1].getpA() + ", " + 
													snps[index2].getpA() + ", " +
													snps[index3].getpA() + ", " +
													snps[index4].getpA() + ", " +
													p12 + ", " + p13 + ", " +  p14 + ", " +
													p23 + ", " +  p24 + ", " + p34 + ")\n");
											}
											
											if( (i12 != x || i13 != y || i14 != z) && consistent)
											{
												fp ++;
												fpFile.print("(" + index1 + "," +index2 + "," + index3 + "," +index4 + ")");
												fpFile.print("\t("+x+y+z+")\t");
												fpFile.print("("+i12+i13+i14+":"+consistent+")\t");
												fpFile.print("(" + snps[index1].getpA() + ", " + 
													snps[index2].getpA() + ", " +
													snps[index3].getpA() + ", " +
													snps[index4].getpA() + ", " +
													p12 + ", " + p13 + ", " +  p14 + ", " +
													p23 + ", " +  p24 + ", " + p34 + ")\n");
											}
										} //end i14
									} //end i13
								} //end i12
								matchResultOutStream.println("\n");
								//for each snps, if not marked, it is eliminated
								for(int i12 = 0; i12 < pAB12.size(); i12++)
								{
									double p12 = pAB12.get(i12);
									if(p12<1)
									{
										if(p12 == r2_12.v1)
										{
											double myConfidenceLevel = UtilityFunctions.computeConfidenceValue(snps, r2_12, false);
											if(myConfidenceLevel >= confidence_level)
											{
												nRecovered ++;
												System.out.println("(" + index1 + "," + index2 + ") recovered\n");
												if (!r2_12.setSignRecovered(false, true))
													nFalseRecovered ++;
												//let the program continue
												propagating = true;
											}
										}
										else if(p12 == r2_12.v2)
										{
											double myConfidenceLevel = UtilityFunctions.computeConfidenceValue(snps, r2_12, true);
											if(myConfidenceLevel >= confidence_level)
											{
												nRecovered ++;
												System.out.println("(" + index1 + "," + index2 + ") recovered\n");
												if(!r2_12.setSignRecovered(true, true))
													nFalseRecovered ++;
												propagating = true;
											}
										}
										else
										{
											//error
											throw new Exception(
													"error detected: p12=" + p12
															+ ",v1=" + r2_12.v1
															+ ",v2=" + r2_12.v2);
										}
									}
								}
								
								//System.out.println("i13\n");
								for(int i13 = 0; i13 < pAB13.size(); i13++)
								{
									double p13 = pAB13.get(i13);
									if(p13<1)
									{
										if(p13 == r2_13.v1)
										{
											double myConfidenceLevel = UtilityFunctions.computeConfidenceValue(snps, r2_13, false);
											if(myConfidenceLevel >= confidence_level)
											{
												nRecovered ++;
												System.err.println("(" + index1 + "," + index3 + ") recovered\n");
												if (!r2_13.setSignRecovered(false, true))
													nFalseRecovered ++;
												propagating = true;
											}
										}
										else if(p13 == r2_13.v2)
										{
											double myConfidenceLevel = UtilityFunctions.computeConfidenceValue(snps, r2_13, true);
											if(myConfidenceLevel >= confidence_level)
											{
												nRecovered ++;
												System.err.println("(" + index1 + "," + index3 + ") recovered\n");
												if(!r2_13.setSignRecovered(true, true))
														nFalseRecovered ++;
												propagating = true;
											}
										}
										else
										{
											//error
											throw new Exception(
													"error detected: p12=" + p13
															+ ",v1=" + r2_13.v1
															+ ",v2=" + r2_13.v2);
										}
									}
								}	//for i13
								
								//System.out.println("i14\n");
								for(int i14 = 0; i14 < pAB14.size(); i14++)
								{
									double p14 = pAB14.get(i14);
									if(p14<1)
									{
										if(p14 == r2_14.v1)
										{
											//hog
											double myConfidenceLevel = UtilityFunctions.computeConfidenceValue(snps, r2_14, false);
											
											if(myConfidenceLevel >= confidence_level)
											{
												nRecovered ++;
												System.err.println("(" + index1 + "," + index4 + ") recovered\n");
												if(!r2_14.setSignRecovered(false, true))
													nFalseRecovered ++;
												propagating = true;
											}
										}
										else if(p14 == r2_14.v2)
										{
											double myConfidenceLevel = UtilityFunctions.computeConfidenceValue(snps, r2_14, true);
											if(myConfidenceLevel >= confidence_level)
											{
												nRecovered ++;
												System.err.println("(" + index1 + "," + index4 + ") recovered\n");
												if(!r2_14.setSignRecovered(true, true))
													nFalseRecovered++;
												propagating = true;
											}
										}
										else
										{
											//error
											throw new Exception(
													"error detected: p12=" + p14
															+ ",v1=" + r2_14.v1
															+ ",v2=" + r2_14.v2);
										}
									}
								} //end for i14	
							}							
							UtilityFunctions.signRecoverRate(snps);
							//UtilityFunctions.saveSNPsToFile("snps_2stage.dat", snps);
						} //end if recovered
						else
						{
							//do nothing, we need 3 snps with signs recovered
						}
						
						analysisFile.println("nRecovered = " + nRecovered + "nFalseRecovered= " + nFalseRecovered);
					} //end 3 for
			//if no more to propagate
			if(!propagating)
				break;
			//reset
			propagating = false;
		}//end while
		
		analysisFile.println("total = " + nTotalConbination);
		analysisFile.println("fp = " + fp);
		analysisFile.println("fn = " + fn);
		analysisFile.println("nRecovered = " + nRecovered + "nFalseRecovered= " + nFalseRecovered);
		analysisFile.close();
		
		matchResultOutStream.close();
	}
	
	//check if it is consistent using by solving linear constraints
	/*
	 * return value:	true:
	 * 						the constraint is satisfiable
	 * 					false
	 * 						not satisfiable, this is wrong
	 */
	public static boolean isConsistent4(boolean printInfo, Double pA1, Double pA2,
			Double pA3, Double pA4, double p12, double p13, double p14,
			double p23, double p24, double p34) {
		if(printInfo)
			System.out.println("isConsistent");
		
		MWNumericArray B = null;	/* Stores input values in a */
		Object[] result = null;		/*Stores the result*/
		
		boolean ret = false;
		try{
			if(lpSolver == null)
			{
				lpSolver = new LPSolver();
			}
			
			int[] dims = {1, 12};	//a column vector with 12 rows
			B = MWNumericArray.newInstance(dims, MWClassID.DOUBLE, MWComplexity.REAL);
			double deviation = -0.01;	// double devi = 0;
			insertArray(B, 1, 1, deviation);	//the first row is the deviation
			double pr = 0;
			if(printInfo)
				pr = 1;
			insertArray(B, 1, 2, pr);	//the second row is printable
			
			//insert the columns
			insertArray(B, 1, 3, pA1);	//the 3rd row is pA1 b1
			insertArray(B, 1, 4, pA2);	//b2
			insertArray(B, 1, 5, pA3);	//b3
			insertArray(B, 1, 6, pA4);	//b4
			insertArray(B, 1, 7, p12);	
			insertArray(B, 1, 8, p13);
			insertArray(B, 1, 9, p14);
			insertArray(B, 1, 10,p23);
			insertArray(B, 1, 11,p24);
			insertArray(B, 1, 12, p34);
			//System.out.println("-> call LPSolver");
			result = lpSolver.isConsistent16(1, B);
			
			MWNumericArray re = (MWNumericArray) result[0];
			if(re.getDouble(1) != 0)
				ret = true;
			else
				ret = false;
			
		}
		catch(Exception e)
		{
			e.printStackTrace();
		}
		finally{
			MWArray.disposeArray(B);
			MWArray.disposeArray(result);
		}
		//System.out.print("<- return LPSolver\n");
		return ret;
	}
	
	public static void insertArray(MWNumericArray B, int x, int y, double value)
	{
		int[] index = {x, y};
		B.set(index, value);
	}
	
	//from the recovered sings, inferring other pairwise signs using 3 snps as a block
	public static void propagateSigns3(SNP[] snps)
	{
		
	}
}
