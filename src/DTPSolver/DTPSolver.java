package DTPSolver;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.*;


public class DTPSolver {
    final private Random random;
    private Graph graph;

    final Set<Vertex> X = new TreeSet<>();
    final Set<Vertex> X_plus = new TreeSet<>();
    final Set<Vertex> X_minu = new TreeSet<>();
    SpanningTree tree;

    int iter_count = 0;
    long start_time;

    Map<Integer, Set<Solution>> solutionPools = new TreeMap<>();
    Solution currBestSol;
    final int TABU_BASE = 10;
    final int TABU_VAR = 10;
    final int MAX_FAIL_COUNT = 100;

    public DTPSolver(String instance, int seed) throws IOException {
        random = new Random(seed);
        readInstance(instance);
    }

    public void readInstance(String instance) throws IOException {
        BufferedReader br = new BufferedReader(new FileReader(instance));
        String line = br.readLine();
        String[] data = line.split(" ");
        int nodeNum = Integer.parseInt(data[0]);
        int edgeNum = Integer.parseInt(data[1]);

        graph = new Graph(nodeNum);

        for (int i = 0; i < edgeNum; ++i) {
            line = br.readLine();
            data = line.split(" ");
            int sourceIndex = Integer.parseInt(data[0]);
            int sinkIndex = Integer.parseInt(data[1]);
            double weight = Double.parseDouble(data[2]);
            Vertex source = graph.get_vertex_by_index(sourceIndex);
            if (source == null) {
                source = new Vertex(sourceIndex);
                graph.add_vertex(source);
            }
            Vertex sink = graph.get_vertex_by_index(sinkIndex);
            if (sink == null) {
                sink = new Vertex(sinkIndex);
                graph.add_vertex(sink);
            }

            graph.add_edge(source, sink, weight);
        }

        br.close();
        graph.trim_to_size();
    }

    private int get_min_index_in_C(double[] C, boolean[] visited){
        var min_index = -1;
        for(int i=0; i<C.length; ++i){
            if(visited[i])continue;
            if(min_index == -1){
                min_index = i;
            }else if(C[i] < C[min_index]){
                min_index = i;
            }
        }
        return min_index;
    }

    private SpanningTree prims_algorithm(){
        double[] C = new double[graph.vertices.length];
        Edge[] E = new Edge[graph.vertices.length];
        boolean[] visited = new boolean[graph.vertices.length];
        Arrays.fill(C, Double.MAX_VALUE);
        Arrays.fill(E, null);
        Arrays.fill(visited, false);

        var tree_edges = new TreeSet<Edge>();
        var tree_weight = 0.0;
        var leaves = new TreeSet<Vertex>();

        while(tree_edges.size() < graph.vertices.length - 1){
            int min_index = get_min_index_in_C(C, visited);
            visited[min_index] = true;
            Vertex source = graph.vertices[min_index];
            if(E[min_index] != null) {
                tree_edges.add(E[min_index]);
                tree_weight += E[min_index].weight;
                leaves.add(source);
                leaves.remove(E[min_index].getOtherEdgeEnd(source));
            }
            for(Edge e : graph.edge_list[source.index]){
                Vertex sink = e.getOtherEdgeEnd(source);
                if(!visited[sink.index] && e.weight < C[sink.index]){
                    C[sink.index] = e.weight;
                    E[sink.index] = e;
                }
            }
        }

        var tree = new SpanningTree(tree_edges, tree_weight);
        tree.leaves = leaves;

        return tree;
    }

    private void init() {
        start_time = System.currentTimeMillis();

        SpanningTree init_tree = prims_algorithm();

        for (Vertex v : graph.vertices) {
            v.degree_to_X = graph.edge_list[v.index].size();
            v.is_in_X = true;
            X.add(v);
        }

        for(Vertex v : init_tree.leaves){
            remove_vertex_from_X(v);
        }
        tree = init_tree;

        System.out.println("Init tree: nodes_num="+X.size() + ", weight="+tree.tree_weight);
    }



    private void sampling(){
        for(;;){
            shrink_X();
            if(X_minu.size() == 0){
                tree = calc_spanning_tree(collect_edges_set_in_X());
            }else{
                tree = null;
            }
            System.out.println("Sampling K=" + X.size());
            boolean res = find_any_feasible_tree();
            var solSet = new TreeSet<Solution>();
            solSet.add(currBestSol);
            solutionPools.put(X.size(), solSet);
            System.out.print("\tFound solution for K=" + X.size()
                        + ", not_dominated=" + currBestSol.not_dominated_count);
            if(currBestSol.tree != null){
                System.out.println(", tree_weight=" + currBestSol.tree.tree_weight);
            }else{
                System.out.println();
            }
            if(!res)break;
        }
    }

    private Solution get_best_sol_from_pool(){
        double best_tree_weight = Double.MAX_VALUE;
        Set<Solution> best_solSet = null;
        for(var solSetEntry : solutionPools.entrySet()){
            int k = solSetEntry.getKey();
            Set<Solution> solSet = solSetEntry.getValue();
            Solution sol = solSet.iterator().next();
            if(sol.tree == null) continue;
            double weight = sol.tree.tree_weight;
            if( weight < best_tree_weight){
                best_tree_weight =weight;
                best_solSet = solSet;
            }
        }
        if(best_solSet == null){
            throw new Error("get best sol from pool error");
        }
        var iter = best_solSet.iterator();
        Solution bestSol = iter.next();
        int rand = random.nextInt(best_solSet.size());
        for(int i=0; i<rand; ++i){
            bestSol = iter.next();
        }
        return bestSol;
    }

    public void solve() {
        init();
        sampling();
        Solution bestSol = get_best_sol_from_pool();
        System.out.println("best tree weight:" + bestSol.tree.tree_weight
                + ", time:" + bestSol.time);

    }


    private void shrink_X() {
        find_cut_vertices();
        Vertex[] non_cut_points = X.stream().filter(v -> !v.is_cut).toArray(size -> new Vertex[size]);
        Vertex rand_remove_vertex = non_cut_points[random.nextInt(non_cut_points.length)];
        remove_vertex_from_X(rand_remove_vertex);
    }

    private void remove_vertex_from_X(Vertex v) {
        v.is_in_X = false;
        X.remove(v);
        X_plus.add(v);
        for (Edge e : graph.edge_list[v.index]) {
            Vertex u = e.getOtherEdgeEnd(v);
            u.degree_to_X--;
            if (!u.is_in_X && u.degree_to_X == 0) {
                X_plus.remove(u);
                X_minu.add(u);
            }
        }
    }

    private void add_vertex_to_X(Vertex v){
        v.is_in_X = true;
        X_plus.remove(v);
        X.add(v);
        for(Edge e : graph.edge_list[v.index]){
            Vertex u = e.getOtherEdgeEnd(v);
            u.degree_to_X++;
            if(!u.is_in_X && u.degree_to_X == 1){
                X_plus.add(u);
                X_minu.remove(u);
            }
        }
    }

    private TreeSet<Edge> collect_edges_set_in_X(){
        TreeSet<Edge> edges = new TreeSet<>();
        for (Vertex v : X) {
            for (Edge e : graph.edge_list[v.index]) {
                Vertex u = e.getOtherEdgeEnd(v);
                if (u.is_in_X) {
                    edges.add(e);
                }
            }
        }
        return edges;
    }


    private TreeSet<Edge> collect_edges_in_X_after_move(Vertex addInV, Vertex moveOutV){
        TreeSet<Edge> edges = collect_edges_set_in_X();
        for(Edge e : graph.edge_list[addInV.index]){
            Vertex u = e.getOtherEdgeEnd(addInV);
            if(u.is_in_X){
                edges.add(e);
            }
        }
        edges.removeAll(graph.edge_list[moveOutV.index]);
        return edges;
    }

    private ArrayList<Vertex> get_candidate_insert_v(){
        var candidates = new ArrayList<Vertex>();

        for(Vertex v : X_minu){
            for(Edge e : graph.edge_list[v.index]){
                Vertex u = e.getOtherEdgeEnd(v);
                if(u.degree_to_X > 0 && !u.is_in_X){
                    candidates.add(u);
                }
            }
        }
        return candidates;
    }

    private Move find_move(boolean is_calc_tree, ArrayList<Vertex> addInVs){
        Move bestMv = new Move(null, null, Integer.MAX_VALUE, null);
        Move bestMv_t = new Move(null, null, Integer.MAX_VALUE, null);
        int best_count = 0;
        int best_count_t = 0;
        find_cut_vertices();

        for(Vertex addInV : addInVs){
            for(Vertex moveOutV : X){
                if(moveOutV.is_cut)continue;
                if(addInV.degree_to_X == 1 && graph.get_edge_weight(addInV, moveOutV) < Double.MAX_VALUE){
                    continue;
                }
                Move mv = evaluate_move(addInV, moveOutV, is_calc_tree);
                if(addInV.tabu_iter <= iter_count) {
                    int cmp = compare_Move(mv, bestMv);
                    if (cmp < 0) {
                        bestMv = mv;
                        best_count = 1;
                    } else if (cmp == 0) {
                        best_count++;
                        if (random.nextInt(best_count) == 0) {
                            bestMv = mv;
                        }
                    }
                }else{
                    int cmp = compare_Move(mv, bestMv_t);
                    if(cmp < 0){
                        bestMv_t = mv;
                        best_count_t = 1;
                    }else if(cmp == 0){
                        best_count_t ++;
                        if(random.nextInt(best_count_t)==0){
                            bestMv_t = mv;
                        }
                    }
                }
            }
        }

        int cmp = compare_Move(bestMv, bestMv_t);
        if(bestMv.addInV == null || cmp > 0 && mv_improves_the_best(bestMv_t)){
            return bestMv_t;
        }else {
            return bestMv;
        }
    }

    private boolean mv_improves_the_best(Move mv){
        if(currBestSol == null || currBestSol.not_dominated_count > 0 && mv.delta < 0){
            return true;
        }else if(currBestSol.not_dominated_count == 0
                && mv.tree != null
                && mv.tree.tree_weight < currBestSol.tree.tree_weight){
            return true;
        }else{
            return false;
        }
    }

    private int compare_Move(Move a, Move b){
        if(a.delta == b.delta){
            return a.tree != null ? Double.compare(a.tree.tree_weight, b.tree.tree_weight) : 0;
        }else{
            return Integer.compare(a.delta, b.delta);
        }
    }

    private Move evaluate_move(Vertex addInV, Vertex moveOutV, boolean is_calc_tree){
        int delta = 0;
        for(Edge e : graph.edge_list[addInV.index]){
            Vertex v = e.getOtherEdgeEnd(addInV);
            if(!v.is_in_X && v.degree_to_X == 0){
                delta--;
            }
        }

        for(Edge e : graph.edge_list[moveOutV.index]){
            Vertex v = e.getOtherEdgeEnd(moveOutV);
            if(!v.is_in_X && v.degree_to_X == 1 && graph.get_edge_weight(v, addInV) == Double.MAX_VALUE){
                delta++;
            }
        }

        SpanningTree tree = null;

        if(is_calc_tree && X_minu.size() + delta == 0){
            tree = calc_spanning_tree(collect_edges_in_X_after_move(addInV, moveOutV));
        }

        return new Move(addInV, moveOutV, delta, tree);
    }


    private void make_move(Move mv){
        add_vertex_to_X(mv.addInV);
        remove_vertex_from_X(mv.moveOutV);
        tree = mv.tree;
        mv.moveOutV.tabu_iter = iter_count + TABU_BASE + random.nextInt(TABU_VAR);
    }

    private int current_configuration_compare_sol(Solution sol){
        if(sol == null) return -1;
        int cmp = Integer.compare(X_minu.size(), sol.not_dominated_count);
        if(cmp == 0 && tree != null){
            cmp = Double.compare(tree.tree_weight, sol.tree.tree_weight);
        }
        return cmp;
    }

    private int compare_sol(Solution s1, Solution s2){
        if(s1 == null && s2 != null){
            return 1;
        }else if(s1 != null && s2 == null){
            return -1;
        }
        int cmp = Integer.compare(s1.not_dominated_count, s2.not_dominated_count);
        if(cmp == 0 && s1.tree != null){
            cmp = Double.compare(s1.tree.tree_weight, s2.tree.tree_weight);
        }
        return cmp;
    }

    private boolean find_any_feasible_tree(){

        currBestSol = new Solution(this);
        for(;;++iter_count){
            Move mv = find_move(false, get_candidate_insert_v());
            if(mv.delta >= 0 || X_minu.isEmpty()){
                break;
            }
            make_move(mv);

        }

        if(!X_minu.isEmpty()){
            return false;
        }

        tree = calc_spanning_tree(collect_edges_set_in_X());
        currBestSol = new Solution(this);
        return true;
    }

//    private boolean local_search(boolean not_tabu){
//        boolean is_descending = true;
//        int fail_improve_count = 0;
//        long log_time = System.currentTimeMillis();
//        currBestSol = new Solution(this);
//        for(;;++iter_count){
//            Move mv = find_move();
//            if(mv.delta > 0
//                    || mv.delta == 0 && mv.tree == null
//                    || X_minu.size() == 0 && mv.tree.tree_weight >= tree.tree_weight){
//                is_descending = false;
//                if(not_tabu)break;
//            }
//            make_move(mv);
//
//            int cmp = current_configuration_compare_sol(currBestSol);
//            if(cmp < 0){
//                currBestSol = new Solution(this);
//                fail_improve_count = 0;
//            }else if(!is_descending){
//                fail_improve_count++;
//            }
//
//            long curr_time = System.currentTimeMillis();
//            if(curr_time - log_time > 1000) {
//                System.out.print("\titer=" + iter_count + ", X_m_size=" + X_minu.size());
//                if (tree != null) {
//                    System.out.println(", tree_weight=" + tree.tree_weight);
//                } else {
//                    System.out.println();
//                }
//                log_time = curr_time;
//            }
//
//            if(fail_improve_count >= MAX_FAIL_COUNT)break;
//
//            //check_configuration();
//        }
//
//        return X_minu.size() == 0;
//    }

    private SpanningTree calc_spanning_tree(TreeSet<Edge> candidate_edges) {
        int[] disjoint_set = new int[graph.vertices.length];
        Arrays.fill(disjoint_set, -1);

        var tree_edges = new TreeSet<Edge>();
        var tree_weight = 0.0;

        for (var edge : candidate_edges) {
            int i = edge.source.index;
            while (disjoint_set[i] >= 0) i = disjoint_set[i];
            int j = edge.sink.index;
            while (disjoint_set[j] >= 0) j = disjoint_set[j];
            if (i != j) {
                tree_edges.add(edge);
                tree_weight += edge.weight;
                if (i < j) {
                    disjoint_set[i] += disjoint_set[j];
                    disjoint_set[j] = i;
                } else {
                    disjoint_set[j] += disjoint_set[i];
                    disjoint_set[i] = j;
                }
            }
        }

        return new SpanningTree(tree_edges, tree_weight);
    }

    private int depth;
    private int num_root_child;
    private Vertex root;

    private void find_cut_vertices() {
        depth = 1;
        num_root_child = 0;
        for (Vertex v : X) {
            v.is_visited = false;
            v.dep = -1;
            v.low = -1;
            v.is_cut = false;
        }
        Vertex r = X.iterator().next();
        r.is_visited = true;
        root = r;
        cut_vertices_recur(r);
        if (num_root_child > 1) {
            r.is_cut = true;
        }
    }

    private void cut_vertices_recur(Vertex r) {
        r.is_visited = true;
        r.dep = depth;
        r.low = depth;
        ++depth;

        for (Edge e : graph.edge_list[r.index]) {
            Vertex tem = e.getOtherEdgeEnd(r);
            if (tem.is_in_X) {
                if (!tem.is_visited) {
                    cut_vertices_recur(tem);
                    r.low = Math.min(r.low, tem.low);
                    if (tem.low >= r.dep && r != root) {
                        r.is_cut = true;
                    } else if (r == root) {
                        ++num_root_child;
                    }
                } else {
                    r.low = Math.min(r.low, tem.dep);
                }
            }
        }
    }

    private void check_split(){
        var check_set = new TreeSet<Vertex>();
        check_set.addAll(X);
        check_set.addAll(X_plus);
        check_set.addAll(X_minu);
        if(check_set.size() != graph.vertices.length){
            throw new Error("check split error 1!");
        }

        var X_cpy = new TreeSet<>(X);
        var X_plus_cpy = new TreeSet<>(X_plus);
        var X_minu_cpy = new TreeSet<>(X_minu);

        X_cpy.removeAll(X_plus_cpy);
        if(X_cpy.size() != X.size()){
            throw new Error("check split error 2!");
        }
        X_plus_cpy.removeAll(X_minu_cpy);
        if(X_plus_cpy.size() != X_plus.size()){
            throw new Error("check split error 3!");
        }
        X_minu_cpy.removeAll(X_cpy);
        if(X_minu_cpy.size() != X_minu.size()){
            throw new Error("check split error 4!");
        }
    }

    private void check_consistency(){
        for(Vertex v : graph.vertices){
            if(v.is_in_X ^ X.contains(v)){
                throw new Error("check consistency error 1, " + v.index);
            }
            int calc_degree_in_X = 0;
            for(Edge e : graph.edge_list[v.index]){
                Vertex u = e.getOtherEdgeEnd(v);
                if(u.is_in_X){
                    calc_degree_in_X++;
                }
            }
            if(calc_degree_in_X != v.degree_to_X){
                throw new Error("check consistency error 2, " + v.index);
            }
            if((!v.is_in_X && v.degree_to_X > 0) ^ X_plus.contains(v)){
                throw new Error("check consistency error 3, " + v.index);
            }
            if((!v.is_in_X && v.degree_to_X == 0) ^ X_minu.contains(v)){
                throw new Error("check consistency error 4, " + v.index);
            }
        }
    }

    private void dfs_X_tree(Vertex parent, Vertex v, boolean[] visited){
        if(visited[v.index]){
            throw new Error("dfs_X_tree");
        }
        visited[v.index] = true;
        for(Edge e : graph.edge_list[v.index]){
            if(tree.tree_edges.contains(e)){
                Vertex u = e.getOtherEdgeEnd(v);
                if(u != parent) {
                    dfs_X_tree(v, u, visited);
                }
            }
        }
    }

    private void check_tree(){
        if(tree == null)return;
        for(Edge e : tree.tree_edges){
            if(!(X.contains(e.source) && X.contains(e.sink))){
                throw new Error("check tree error 1");
            }
        }
        if(tree.tree_edges.size() != X.size() - 1){
            throw new Error("check tree error 2");
        }
        boolean[] visited = new boolean[graph.vertices.length];
        dfs_X_tree(null, X.iterator().next(), visited);
    }

    private void check_configuration(){
        check_split();
        check_consistency();
        check_tree();
    }

    static class EdgeWeightComparator implements Comparator {

        @Override
        public int compare(Object o1, Object o2) {
            Edge e1 = (Edge) o1;
            Edge e2 = (Edge) o2;
            return Double.compare(e1.weight, e2.weight);
        }
    }

    public class SpanningTree {
        final public Set<Edge> tree_edges;
        final public double tree_weight;
        public Set<Vertex> leaves = null;

        SpanningTree(Set<Edge> edges, double weight) {
            tree_edges = edges;
            tree_weight = weight;
        }
    }

    private class Move {
        final public Vertex addInV;
        final public Vertex moveOutV;
        final public int delta;
        final public SpanningTree tree;

        Move(Vertex addInV, Vertex moveOutV, int delta, SpanningTree tree) {
            this.addInV = addInV;
            this.moveOutV = moveOutV;
            this.delta = delta;
            this.tree = tree;
        }
    }

    public class Solution implements Comparable{
        final Graph graph;
        final public Set<Vertex> dominating_set;
        final public int not_dominated_count;
        final public SpanningTree tree;
        final double time;
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
    }
}
