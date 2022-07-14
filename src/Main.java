import DTPSolver.DTPSolver;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;


public class Main {

    public static ArrayList<String> calibration(){
        var instances = new ArrayList<String>();
        instances.add("instances/Range_100/ins_050_1.txt");
        instances.add("instances/Range_100/ins_100_1.txt");
        instances.add("instances/Range_100/ins_200_1.txt");
        instances.add("instances/Range_100/ins_300_1.txt");
        instances.add("instances/Range_100/ins_400_1.txt");
        instances.add("instances/Range_100/ins_500_1.txt");

        instances.add("instances/Range_125/ins_50_1.txt");
        instances.add("instances/Range_125/ins_100_1.txt");
        instances.add("instances/Range_125/ins_200_1.txt");
        instances.add("instances/Range_125/ins_300_1.txt");
        instances.add("instances/Range_125/ins_400_1.txt");
        instances.add("instances/Range_125/ins_500_1.txt");

        instances.add("instances/Range_150/ins_50_1.txt");
        instances.add("instances/Range_150/ins_100_1.txt");
        instances.add("instances/Range_150/ins_200_1.txt");
        instances.add("instances/Range_150/ins_300_1.txt");
        instances.add("instances/Range_150/ins_400_1.txt");
        instances.add("instances/Range_150/ins_500_1.txt");
        return instances;
    }
    public static ArrayList<String> test(){
        var instances = new ArrayList<String>();
        instances.add("instances/Range_100/ins_050_1.txt");
//        instances.add("instances/Range_100/ins_050_2.txt");
//        instances.add("instances/Range_100/ins_050_3.txt");
        instances.add("instances/Range_100/ins_100_1.txt");
//        instances.add("instances/Range_100/ins_100_2.txt");
//        instances.add("instances/Range_100/ins_100_3.txt");
        instances.add("instances/Range_100/ins_200_1.txt");
//        instances.add("instances/Range_100/ins_200_2.txt");
//        instances.add("instances/Range_100/ins_200_3.txt");
        instances.add("instances/Range_100/ins_300_1.txt");
//        instances.add("instances/Range_100/ins_300_2.txt");
//        instances.add("instances/Range_100/ins_300_3.txt");
        instances.add("instances/Range_100/ins_400_1.txt");
//        instances.add("instances/Range_100/ins_400_2.txt");
//        instances.add("instances/Range_100/ins_400_3.txt");
        instances.add("instances/Range_100/ins_500_1.txt");
//        instances.add("instances/Range_100/ins_500_2.txt");
//        instances.add("instances/Range_100/ins_500_3.txt");
//
        instances.add("instances/Range_125/ins_50_1.txt");
//        instances.add("instances/Range_125/ins_50_2.txt");
//        instances.add("instances/Range_125/ins_50_3.txt");
        instances.add("instances/Range_125/ins_100_1.txt");
//        instances.add("instances/Range_125/ins_100_2.txt");
//        instances.add("instances/Range_125/ins_100_3.txt");
        instances.add("instances/Range_125/ins_200_1.txt");
//        instances.add("instances/Range_125/ins_200_2.txt");
//        instances.add("instances/Range_125/ins_200_3.txt");
        instances.add("instances/Range_125/ins_300_1.txt");
//        instances.add("instances/Range_125/ins_300_2.txt");
//        instances.add("instances/Range_125/ins_300_3.txt");
        instances.add("instances/Range_125/ins_400_1.txt");
//        instances.add("instances/Range_125/ins_400_2.txt");
//        instances.add("instances/Range_125/ins_400_3.txt");
//        instances.add("instances/Range_125/ins_500_1.txt");
//        instances.add("instances/Range_125/ins_500_2.txt");
//        instances.add("instances/Range_125/ins_500_3.txt");
//
//        instances.add("instances/Range_150/ins_50_1.txt");
//        instances.add("instances/Range_150/ins_50_2.txt");
//        instances.add("instances/Range_150/ins_50_3.txt");
//        instances.add("instances/Range_150/ins_100_1.txt");
//        instances.add("instances/Range_150/ins_100_2.txt");
//        instances.add("instances/Range_150/ins_100_3.txt");
//        instances.add("instances/Range_150/ins_200_1.txt");
//        instances.add("instances/Range_150/ins_200_2.txt");
//        instances.add("instances/Range_150/ins_200_3.txt");
//        instances.add("instances/Range_150/ins_300_1.txt");
//        instances.add("instances/Range_150/ins_300_2.txt");
//        instances.add("instances/Range_150/ins_300_3.txt");
//        instances.add("instances/Range_150/ins_400_1.txt");
//        instances.add("instances/Range_150/ins_400_2.txt");
//        instances.add("instances/Range_150/ins_400_3.txt");
//        instances.add("instances/Range_150/ins_500_1.txt");
//        instances.add("instances/Range_150/ins_500_2.txt");
//        instances.add("instances/Range_150/ins_500_3.txt");
        return instances;
    }


    public static void testAll() throws IOException {
        int run = 10;
        var resFile = "resTotal.csv";
        var instances = calibration();
        BufferedWriter bw = new BufferedWriter(new FileWriter(resFile, true));
        bw.write("\ninstance");
        for(var i=0; i<run; ++i){
            bw.write(","+ i + ",v_num,time");
        }
        bw.close();
        for(var instance : instances) {
            bw = new BufferedWriter(new FileWriter(resFile, true));
            bw.write("\n"+instance);
            bw.close();
            for(int i=0; i<run; ++i) {
                DTPSolver solver = new DTPSolver(instance, i);
                var sol = solver.solve();
                bw = new BufferedWriter(new FileWriter(resFile, true));
                bw.write("," + sol.tree.tree_weight + "," + sol.dominating_set.size() + "," + sol.time);
                bw.close();
            }
        }
    }
    public static void main(String[] args) throws IOException {
//        //parameters should be: instance_file seed time_limit tabu_tenure_range perturb_trigger_period perturb_strength
//        //For example instance.txt 1000 0 5-10 100 0.1-0.3
//        var instance_file = args[0];
//        var seed = Integer.parseInt(args[1]);
//        var time_limit = Double.parseDouble(args[2]);
//        var tt_range = args[3];
//        var ptp = Integer.parseInt(args[4]);
//        var ps = args[5];
//
//        var solver = new DTPSolver(instance_file, seed);
//        solver.setTime_limit(time_limit);
//        solver.set_tabu_tenure(tt_range);
//        solver.set_MAX_FAIL_COUNT(ptp);
//        solver.set_perturb_ratio(ps);
//        var sol = solver.solve();
//        System.out.println("FINAL BEST OBJ: " + sol.tree.tree_weight);
//    }
    /*The following is the old batch run main*/
	 //parameters should be: instance_file  time_limit run_count
        int run = Integer.parseInt(args[2]);
        var resFile = "resTotal.csv";
        var instance_file = args[0];
        var time_limit = Double.parseDouble(args[1]);

        BufferedWriter bw = new BufferedWriter(new FileWriter(resFile, true));
        bw.write("\ninstance");
        for(var i=0; i<run; ++i){
            bw.write(","+ i + ",v_num,time");
        }
        bw.write("\n"+instance_file);
        bw.close();

        for(int i=0; i<run; ++i) {
            DTPSolver solver = new DTPSolver(instance_file, i);
            solver.setTime_limit(time_limit);
            var sol = solver.solve();
            bw = new BufferedWriter(new FileWriter(resFile, true));
            bw.write("," + sol.tree.tree_weight + "," + sol.dominating_set.size() + "," + sol.time);
            bw.close();
        }
    }
}
