package DTPSolver;

import java.util.Set;
import java.util.TreeSet;

public class Solution implements Comparable{
    final Graph graph;
    final public Set<Vertex> dominating_set;
    final public int not_dominated_count;
    final public SpanningTree tree;
    final public double time;
    Solution(DTPSolver solver){
        solver.check_configuration();
        this.graph = solver.graph;
        this.dominating_set = new TreeSet<>(solver.X);
        this.tree = solver.tree;
        not_dominated_count = solver.X_minu.size();
        time = (System.currentTimeMillis() - solver.start_time)/1000.0;
    }

    @Override
    public int compareTo(Object o) {
        Solution sp = (Solution)o;
        int cmp = Integer.compare(dominating_set.size(), sp.dominating_set.size());

        if(cmp == 0) {
            var s1Iter = dominating_set.iterator();
            var s2Iter = sp.dominating_set.iterator();
            while (cmp == 0 && s1Iter.hasNext()) {
                Vertex s1v = s1Iter.next();
                Vertex s2v = s2Iter.next();
                cmp = Double.compare(s1v.index, s2v.index);
            }
        }
        return cmp;
    }

    @Override
    public boolean equals(Object o){
        Solution oo = (Solution) o;
        return dominating_set.equals(oo.dominating_set);
    }

    static public int compare_better(Solution s1, Solution s2){
        if(s1.tree == null && s2.tree == null){
            return Integer.compare(s1.not_dominated_count, s2.not_dominated_count);
        }else if(s1.tree == null){
            return 1;
        }else if(s2.tree == null){
            return -1;
        }else{
            return Double.compare(s1.tree.tree_weight, s2.tree.tree_weight);
        }

    }
}
