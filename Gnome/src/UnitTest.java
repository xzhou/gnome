import java.io.*;

import lPSolver.LPSolver;

import com.mathworks.toolbox.javabuilder.MWArray;
import com.mathworks.toolbox.javabuilder.MWClassID;
import com.mathworks.toolbox.javabuilder.MWComplexity;
import com.mathworks.toolbox.javabuilder.MWNumericArray;


public class UnitTest {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		
		try {
			PrintStream ps = new PrintStream(
			        new BufferedOutputStream(new FileOutputStream(
			        new File("./javalog.txt"))), true);
			System.setOut(ps);
			testLPSolver();
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

	}

	
	public static void testLPSolver()
	{
		MWNumericArray B = null;	/* Stores input values in a */
		Object[] result = null;		/*Stores the result*/
		
		boolean ret = false;
		try{
			LPSolver lpSolver = new LPSolver();
			
			int[] dims = {1, 12};	//a column vector with 12 rows
			B = MWNumericArray.newInstance(dims, MWClassID.DOUBLE, MWComplexity.REAL);
			insertArray(B, 1, 1, 0.01);	//the first row is the deviation
			insertArray(B, 1, 2, 1);	//the second row is printable
			//insert the columns
			insertArray(B, 1, 3, 0.3);	//the 3rd row is pA1 b1
			insertArray(B, 1, 4, 0.4);	//b2
			insertArray(B, 1, 5, 0.1);	//b3
			insertArray(B, 1, 6, 0.2);	//b4
			insertArray(B, 1, 7, 0.2);	
			insertArray(B, 1, 8, 0.2);
			insertArray(B, 1, 9, 0.4);
			insertArray(B, 1, 10,0.5);
			insertArray(B, 1, 11,0.6);
			insertArray(B, 1, 12, 0.1);
			//System.out.println("-> call LPSolver");
			result = lpSolver.isConsistent4(1, B);
			
			MWNumericArray re = (MWNumericArray) result[0];
			if(re.getDouble(1) != 0)
			{
				System.out.println("Consistent");
				ret = true;
			}
			else
			{
				ret = false;
			}
			
		}
		catch(Exception e)
		{
			e.printStackTrace();
		}
		finally{
			MWArray.disposeArray(B);
			MWArray.disposeArray(result);
		}
		
	}
	
	public static void insertArray(MWNumericArray B, int x, int y, double value)
	{
		int[] index = {x, y};
		B.set(index, value);
	}
}
