import jp.ac.kobe_u.cs.cream.DefaultSolver;
import jp.ac.kobe_u.cs.cream.IntVariable;
import jp.ac.kobe_u.cs.cream.Network;
import jp.ac.kobe_u.cs.cream.Solution;
import jp.ac.kobe_u.cs.cream.Solver;


public class SolverTest {
	public static void main(String[] args) {
		Test1();
	}

	private static void Test1() {
		Network net = new Network();

		IntVariable x00 = new IntVariable(net);
	    IntVariable x01 = new IntVariable(net);
	    IntVariable x02 = new IntVariable(net);
	    IntVariable x10 = new IntVariable(net);
	    IntVariable x11 = new IntVariable(net);
	    IntVariable x12 = new IntVariable(net);
	    IntVariable x20 = new IntVariable(net);
	    IntVariable x21 = new IntVariable(net);
	    IntVariable x22 = new IntVariable(net);

	    x00.ge(0); x01.ge(0); x02.ge(0);
	    x10.ge(0); x11.ge(0); x12.ge(0);
	    x20.ge(0); x21.ge(0); x22.ge(0);
	    
	    x00.add(x01).add(x02).add(x10).add(x11).add(x12).add(x20).add(x21).add(x22).equals(1064);
	    IntVariable lb = x00.add(x00).add(x01).add(x10);
	    lb.le(778);
	    lb.add(x11).ge(778);
	    x00.add(x00).add(x01).add(x01).add(x02).add(x02).add(x10).add(x11).add(x12).equals(1257);
	    x00.add(x00).add(x10).add(x10).add(x20).add(x20).add(x01).add(x11).add(x21).equals(1648);
	    
	    Solver solver = new DefaultSolver(net);

	    for (solver.start(); solver.waitNext(); solver.resume()) {
	        Solution solution = solver.getSolution();
	        int X00 = solution.getIntValue(x00);
	        int X01 = solution.getIntValue(x01);
	        int X02 = solution.getIntValue(x02);
	        int X10 = solution.getIntValue(x10);
	        int X11 = solution.getIntValue(x11);
	        int X12 = solution.getIntValue(x12);
	        int X20 = solution.getIntValue(x20);
	        int X21 = solution.getIntValue(x21);
	        int X22 = solution.getIntValue(x22);
	        System.out.println("x00 = " + X00 +
	        		"x01 = " + X01 +
	        		"x02 = " + X02 +
	        		"x10 = " + X10 +
	        		"x11 = " + X11 +
	        		"x12 = " + X12 +
	        		"x20 = " + X20 +
	        		"x21 = " + X21 +
	        		"x22 = " + X22);
	    }
	    solver.stop();
	}
}
