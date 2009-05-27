import com.mathworks.toolbox.javabuilder.MWClassID;
import com.mathworks.toolbox.javabuilder.MWComplexity;
import com.mathworks.toolbox.javabuilder.MWNumericArray;


public class UnitTest {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		int[] dims = {1,8};
		MWNumericArray B = MWNumericArray.newInstance(dims, MWClassID.DOUBLE,
				MWComplexity.REAL);
		
		insertArray(B, 1, 1, 0.1);
		insertArray(B, 1, 2, 0.9);
		
		System.out.println(B.toString());
	}

	public static void insertArray(MWNumericArray B, int x, int y, double value)
	{
		int[] index = {x, y};
		B.set(index, value);
	}
}
