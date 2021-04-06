import DTPSolver.DTPSolver;

import java.io.IOException;
import java.util.ArrayList;


public class Main {

    public static void main(String[] args) throws IOException {
	// write your code here
        DTPSolver solver = new DTPSolver("instances/Range_100/ins_500_1.txt", 1);
//        DTPSolver solver = new DTPSolver("instances/test.txt", 0);
        solver.solve();
    }
}
