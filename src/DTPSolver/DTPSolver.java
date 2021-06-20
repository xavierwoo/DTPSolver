package DTPSolver;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.*;
import java.util.stream.Collectors;


public class DTPSolver {
    final private Random random;
    Graph graph;

    final Set<Vertex> X = new TreeSet<>();
    final Set<Vertex> X_plus = new TreeSet<>();
    final Set<Vertex> X_minu = new TreeSet<>();
    SpanningTree tree;

    int iter_count = 0;
    long start_time;
    double time_limit = 1000;

    Map<Integer, Set<Solution>> solutionPools = new TreeMap<>();
    Map<Integer, Integer> same_sol_again_count = new TreeMap<>();
    Solution currBestSol;
    final int TABU_BASE = 10;
    final int TABU_VAR = 50;
    final int MAX_FAIL_COUNT = 500;
    final int SOLVE_OFFSET = 2;
    final double PERTURB_STR_RATIO_MIN = 0.3;
    final double PERTURB_STR_RATIO_MAX = 0.5;
    final int PERTURB_TIMES_MAX = 3;
    final double CUT_RATIO = 0.05;

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
        System.out.println("Sampling...");
        for(;;){
            shrink_X();
            if(X_minu.size() == 0){
                tree = calc_spanning_tree(collect_edges_set_in_X());
            }else{
                tree = null;
            }
            boolean res = find_any_feasible_tree();
            if(!res)break;

            var solSet = new TreeSet<Solution>();
            solSet.add(currBestSol);
            solutionPools.put(X.size(), solSet);
            same_sol_again_count.put(X.size(), 0);
        }
        Solution sol = get_best_sol_from_pool();
        System.out.println("Best solution from sampling. tree_w=" + sol.tree.tree_weight +
                ", X_size=" + sol.dominating_set.size());
    }

    private Solution get_best_sol_from_pool(){
        double best_tree_weight = Double.MAX_VALUE;
        Set<Solution> best_solSet = null;
        for(var solSetEntry : solutionPools.entrySet()){
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
        return bestSol;
    }

    private int get_best_solution_X_size(){
        double best_weight = Double.MAX_VALUE;
        int best_X_size = 0;
        for(var solSetEntry : solutionPools.entrySet()){
            int X_size = solSetEntry.getKey();
            var sol = solSetEntry.getValue().iterator().next();
            if(sol.tree != null && sol.tree.tree_weight < best_weight){
                best_weight = sol.tree.tree_weight;
                best_X_size = X_size;
            }
        }
        return best_X_size;
    }

    private Solution get_random_solution_in_pool(int X_size){
        var solSet = solutionPools.get(X_size);
        if(solSet == null)return null;
        return get_random_item_in_set(solSet);
    }

    private <T> T get_random_item_in_set(Set<T> set){
        int randi = random.nextInt(set.size());
        var iter = set.iterator();
        T item = iter.next();
        for(int i=0; i < randi; ++i){
            item = iter.next();
        }
        return item;
    }

    private void calc_plus_minu_set(){
        for(Vertex v : graph.vertices){
            v.degree_to_X = 0;
            v.tabu_iter = 0;
            v.is_in_X = X.contains(v);
            for(Edge e : graph.edge_list[v.index]){
                Vertex u = e.getOtherEdgeEnd(v);
                if(X.contains(u)){
                    ++v.degree_to_X;
                }
            }
            if(!v.is_in_X) {
                if (v.degree_to_X > 0) {
                    X_plus.add(v);
                } else {
                    X_minu.add(v);
                }
            }
        }
    }

    private void rebuilt_curr_sol(Solution sol){
        X.clear();
        X.addAll(sol.dominating_set);
        X_plus.clear();
        X_minu.clear();
        calc_plus_minu_set();
        tree = sol.tree;
    }

    private void gen_random_configuration(int X_size){
        X.clear();
        X_plus.clear();
        X_minu.clear();
        int randi = random.nextInt(graph.vertices.length);
        X.add(graph.vertices[randi]);

        while(X.size() < X_size){
            var v_in_X = get_random_item_in_set(X);
            for(Edge e : graph.edge_list[v_in_X.index]){
                Vertex u = e.getOtherEdgeEnd(v_in_X);
                if(!X.contains(u)){
                    X.add(u);
                    break;
                }
            }
        }
        calc_plus_minu_set();
        if(X_minu.isEmpty()){
            tree = calc_spanning_tree(collect_edges_set_in_X());
        }else{
            tree = null;
        }
    }

    public Solution solve() {
        int[] try_sizes = new int[SOLVE_OFFSET * 2 + 1];
        init();
        sampling();
        Solution bestSol;

        while((System.currentTimeMillis() - start_time)/1000.0 < time_limit){
            int best_X_size = get_best_solution_X_size();
            try_sizes[0] = best_X_size;
            for(int i=1, pos=1; i<=SOLVE_OFFSET; i++, pos+=2){
                try_sizes[pos] = try_sizes[0] - i;
                try_sizes[pos+1] = try_sizes[0] + i;
            }
            //for(int X_size=best_X_size + SOLVE_OFFSET; X_size>=best_X_size-SOLVE_OFFSET; --X_size) {
            for(int X_size : try_sizes){
                if(!solutionPools.containsKey(X_size+1)){
                    System.out.println("Skip X_size="+X_size);
                    continue;
                }
                System.out.println("Try X_size="+X_size);
                Solution sol = get_random_solution_in_pool(X_size);

                if(sol == null){
                    gen_random_configuration(X_size);
                    currBestSol = new Solution(this);
                }else {
                    rebuilt_curr_sol(sol);
                    currBestSol = sol;
                }
                check_configuration();
                local_search();

                if(currBestSol.not_dominated_count > 0){
                    System.out.println("Failed for X_size="+X_size);
                    continue;
                }
                System.out.println("Best res for X_size="+X_size + ": tree_w=" + currBestSol.tree.tree_weight);
            }

        }

        bestSol = get_best_sol_from_pool();
        System.out.println("best tree weight:" + String.format("%.2f",bestSol.tree.tree_weight)
                + ", X_size:" + bestSol.dominating_set.size() + ", time:" + bestSol.time);
        return bestSol;
    }

    private Move get_random_move(){
        Vertex iv, rv;
        do{
            find_cut_vertices();
            iv = get_random_item_in_set(X_plus);
            var rv_candidate = X.stream().filter(v->!v.is_cut).collect(Collectors.toList());
            rv = rv_candidate.get(random.nextInt(rv_candidate.size()));
        }while(is_infeasible_move(iv, rv));
        return evaluate_move(iv, rv, true);
    }

    private void perturb(int strength){
        System.out.println("\tPerturb str="+strength);
        Solution sol = get_random_solution_in_pool(X.size());
        rebuilt_curr_sol(sol);

        for(int i=0; i<strength; ++i){
            Move mv = get_random_move();
            make_move(mv);
        }
        find_cut_vertices();
        check_configuration();
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
        graph.edge_list[moveOutV.index].forEach(edges::remove);
        return edges;
    }

    private ArrayList<Vertex> get_reduced_candidate_insert_v(){
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

    private boolean is_infeasible_move(Vertex addInV, Vertex moveOutV){
        return addInV.degree_to_X == 1 && graph.get_edge_weight(addInV, moveOutV) < Double.MAX_VALUE;
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
                if(is_infeasible_move(addInV,moveOutV)){
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

        SpanningTree new_tree = null;

        if(is_calc_tree && X_minu.size() + delta == 0){
            if(tree == null){
                new_tree = calc_spanning_tree(collect_edges_in_X_after_move(addInV, moveOutV));
            }else {
                //TODO:check if update improves
//                new_tree = update_current_tree(addInV, moveOutV);
                new_tree = calc_spanning_tree(collect_edges_in_X_after_move(addInV, moveOutV));
//
//                var test_tree = update_current_tree(addInV, moveOutV);
//                if(Math.abs(new_tree.tree_weight - test_tree.tree_weight)> 0.001){
//                    throw new Error("update calc spanning tree error");
//                }
            }
        }

        return new Move(addInV, moveOutV, delta, new_tree);
    }

    private SpanningTree update_current_tree(Vertex addInV, Vertex moveOutV){
        Graph g_tree = graph.gen_new_subgraph_from_tree(tree);

        update_calc_tree_add_v(X.size() + 1, g_tree, addInV);
        return update_calc_tree_remove_v(X.size(), g_tree, addInV, moveOutV);
    }

    private SpanningTree update_calc_tree_remove_v(int x_size, Graph g_tree, Vertex addInV,Vertex moveOutV){

        for(var e: graph.edge_list[moveOutV.index]){
            var other_v = e.getOtherEdgeEnd(moveOutV);
            if(other_v.is_in_X || other_v == addInV){
                g_tree.remove_edge(e);
            }
        }
        var curr_tree_edges = g_tree.get_all_edges_des();
        var candidate_edges = collect_edges_in_X_after_move(addInV, moveOutV);
        candidate_edges.removeAll(curr_tree_edges);

        return update_kruskal(g_tree, curr_tree_edges, candidate_edges);
    }

    private SpanningTree update_kruskal(Graph g_tree, Set<Edge> curr_tree_edges, Set<Edge> candidate_edges){
        int[] disjoint_set = new int[graph.vertices.length];
        Arrays.fill(disjoint_set, -1);
        prepare_disjoint_set(g_tree, disjoint_set);
        kruskal_subprocess(candidate_edges, disjoint_set, curr_tree_edges);

        var tree_weight = 0.0;
        for(Edge e : curr_tree_edges){
            tree_weight += e.weight;
        }
        return new SpanningTree(curr_tree_edges, tree_weight);
    }

    private void prepare_disjoint_set(Graph g_tree, int[] disjoint_set){
        boolean[] visited = new boolean[graph.vertices.length];
        for(Vertex v : X){
            if(visited[v.index])continue;
            visited[v.index] = true;
            dfs_prepare_disjoint_set(g_tree, v, disjoint_set, visited);
        }
    }

    private int dfs_prepare_disjoint_set(Graph g_tree, Vertex root, int[] disjoint_set, boolean[] visited){
        int vertices_count = 1;
        for(Edge e : g_tree.edge_list[root.index]){
            Vertex u = e.getOtherEdgeEnd(root);
            if(visited[u.index])continue;
            visited[u.index] = true;
            disjoint_set[u.index] = root.index;
            vertices_count += dfs_prepare_disjoint_set(g_tree, u, disjoint_set, visited);
        }
        if(disjoint_set[root.index] == -1){
            disjoint_set[root.index] = -vertices_count;
        }
        return vertices_count;
    }

    private void update_calc_tree_add_v(int x_size, Graph g_tree, Vertex addInV){
        for(var e : graph.edge_list[addInV.index]){
            var other_v = e.getOtherEdgeEnd(addInV);
            if(other_v.is_in_X) {
                g_tree.add_edge(e);
            }
        }
        var edges = g_tree.get_all_edges_des();
        int curr_edges_count = edges.size();

        boolean[] visited = new boolean[graph.vertices.length];
        for(Edge e : edges){
            Arrays.fill(visited, false);
            int count = DFS_subgraph_count(g_tree, e.source, e, visited);
            if(count != x_size) continue;
            g_tree.remove_edge(e);
            curr_edges_count--;
            if(curr_edges_count == x_size - 1)break;
        }
    }

    private int DFS_subgraph_count(Graph subgraph, Vertex start, Edge except_e, boolean[] visited){
        int count = 0;
        for(Edge e : subgraph.edge_list[start.index]){
            if(e == except_e) continue;
            Vertex r = e.getOtherEdgeEnd(start);
            if(visited[r.index]) continue;
            visited[r.index] = true;
            count ++;
            count += DFS_subgraph_count(subgraph, r, except_e, visited);
        }
        return count;
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
        if(cmp == 0 && tree != null && sol.tree != null){
            cmp = Double.compare(tree.tree_weight, sol.tree.tree_weight);
        }
        return cmp;
    }

    private boolean find_any_feasible_tree(){

        currBestSol = new Solution(this);
        for(;;++iter_count){
            Move mv = find_move(false, get_reduced_candidate_insert_v());
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

    private ArrayList<Vertex> prepare_candidate_addIns(){
        if(X_minu.isEmpty()){
            return new ArrayList<>(X_plus);
        }else{
            return get_reduced_candidate_insert_v();
        }
    }

    private boolean local_search(){
        int fail_improve_count = 0;
        long log_time = System.currentTimeMillis();
        if(X_minu.isEmpty() && tree == null){
            throw new Error("local search error1");
        }

        int per_str_min = (int)(X.size() * PERTURB_STR_RATIO_MIN);
        int per_str_max = (int)(X.size() * PERTURB_STR_RATIO_MAX);
        int per_times = 0;
        var solSet = solutionPools.get(X.size());
        if(solSet == null){
            solSet = new TreeSet<>();
            solSet.add(new Solution(this));
            solutionPools.put(X.size(), solSet);
            same_sol_again_count.put(X.size(), 0);
        }
        for(;(System.currentTimeMillis() - start_time)/1000.0 < time_limit;++iter_count){
            ArrayList<Vertex> candidate_addIns = prepare_candidate_addIns();
            Move mv = find_move(true, candidate_addIns);

            make_move(mv);
            fail_improve_count++;
            int cmp = current_configuration_compare_sol(currBestSol);
            if(cmp < 0){
                currBestSol = new Solution(this);
                fail_improve_count = 0;
                per_times = 0;
                if(solSet != null) {
                    solSet.clear();
                }
                solSet = new TreeSet<>();
                solutionPools.put(X.size(), solSet);
                solSet.add(currBestSol);
                same_sol_again_count.put(X.size(), 0);

            }else if(cmp == 0){
                var currSol = new Solution(this);
                if (solSet.contains(currSol)) {
                    same_sol_again_count.put(X.size(), same_sol_again_count.get(X.size()) + 1);
                } else {
                    solSet.add(currSol);
                }
            }

            long curr_time = System.currentTimeMillis();
            if(curr_time - log_time > 5000) {
                System.out.print("\titer=" + iter_count + ", X_size=" + X.size());
                if (tree != null) {
                    System.out.print(", tree_w=" + String.format("%.2f", tree.tree_weight) + ", best_tree_w="+String.format("%.2f", currBestSol.tree.tree_weight));
                } else {
                    System.out.print(", X_m_size=" + X_minu.size());
                }
                var bestSol = get_best_sol_from_pool();
                System.out.println(", ever_best_tree_w="+String.format("%.2f", bestSol.tree.tree_weight));
                log_time = curr_time;
            }

            if(per_times >= PERTURB_TIMES_MAX)break;
            if(fail_improve_count >= MAX_FAIL_COUNT){
                Solution bestSol = get_best_sol_from_pool();
                if(per_times >= PERTURB_TIMES_MAX/3
                        && (currBestSol.tree == null
                        ||
                        (currBestSol.tree.tree_weight - bestSol.tree.tree_weight)/bestSol.tree.tree_weight > CUT_RATIO)){
                    System.out.println("\tToo bad! dropped");
                    break;
                }
                Integer count = same_sol_again_count.get(X.size());
                int str = count== null? per_str_max : Math.min(per_str_min + same_sol_again_count.get(X.size()), per_str_max);
                perturb(str);
                per_times++;
                fail_improve_count =0;
            }

            //check_configuration();
        }

        return X_minu.size() == 0;
    }

    private SpanningTree calc_spanning_tree(TreeSet<Edge> candidate_edges) {
        int[] disjoint_set = new int[graph.vertices.length];
        Arrays.fill(disjoint_set, -1);

        var tree_edges = new TreeSet<Edge>();
        var tree_weight = 0.0;

        tree_weight = kruskal_subprocess(candidate_edges, disjoint_set, tree_edges);

        return new SpanningTree(tree_edges, tree_weight);
    }

    private double kruskal_subprocess(Set<Edge> candidate_edges, int[] disjoint_set, Set<Edge> tree_edges){
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
        return tree_weight;
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

    void check_configuration(){
        check_split();
        check_consistency();
        check_tree();
    }



    static private class Move {
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
}
