import java.util.Enumeration;
import java.util.Vector;

import jp.ac.kobe_u.cs.cream.IntVariable;
import jp.ac.kobe_u.cs.cream.Solution;
import jp.ac.kobe_u.cs.cream.SolutionHandler;
import jp.ac.kobe_u.cs.cream.Solver;

class FirstStepHandler implements SolutionHandler {
    public Vector<Record> records;
    public int solutions = 0;

    public FirstStepHandler(Vector<Record> records) {
    	this.records = records;
    }

    public synchronized void solved(Solver solver, Solution solution) {
    	if (solution != null) {
    		solutions++;
    		printHaplotypeSolution(solution,records);
    	}
    }
    
    private static void printHaplotypeSolution(Solution solution, Vector<Record> records) {
		int i = 0,recovered = 0,existed = 0;
		for (Enumeration<Record> e = records.elements();e.hasMoreElements();) {
			Record r = e.nextElement();
			int count = solution.getIntValue(r.constraint.count);
			
			if (r.count == -1) {
				r.count = count;
				r.min = count;
				r.max = count;
			} else {
				if (count < r.min) r.min = count;
				if (count > r.max) r.max = count;
			}
			
			if (r.count == r.min && r.count == r.max) {
				recovered++;
				if (r.count > 0) System.out.println(""+count+"-"+r.toString());
			} else if (r.min > 0) {
				existed++;
				System.out.println(""+count+"["+r.min+","+r.max+"]-"+r.toString());
			} else {
				System.out.println(""+count+"["+r.min+","+r.max+"]-"+r.toString());
			}
			                            
			i++;
		}
		System.out.println(""+recovered+" out of "+ i +" records recovered"); 
		System.out.println(""+existed+" out of "+ i +" records existed");
	}
}
