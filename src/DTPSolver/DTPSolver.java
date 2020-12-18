package DTPSolver;

import org.jgrapht.alg.interfaces.SpanningTreeAlgorithm.SpanningTree;
import org.jgrapht.alg.spanning.KruskalMinimumSpanningTree;
import org.jgrapht.graph.DefaultWeightedEdge;
import org.jgrapht.graph.SimpleWeightedGraph;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.*;
import java.util.stream.Collectors;

public class DTPSolver {
    final private Random random;
    final private SimpleWeightedGraph<Vertex, DefaultWeightedEdge> graph =
            new SimpleWeightedGraph<>(DefaultWeightedEdge.class);

    SimpleWeightedGraph<Vertex, DefaultWeightedEdge> X_graph;
    final Set<Vertex> X_plus = new TreeSet<>();
    final Set<Vertex> X_minu = new TreeSet<>();
    SpanningTree<DefaultWeightedEdge> currSpanningTree = null;

    Solution bestSolution = null;

    int iter_count = 0;
    final int tabu_tenure_min = 1;
    final int tabu_tenure_max = 10;
    Map<Vertex, Integer> tabuList = new TreeMap<>();

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

        TreeMap<Integer, Vertex> nodes = new TreeMap<>();

        while((line=br.readLine()) != null){
            data = line.split(" ");
            int sourceIndex = Integer.parseInt(data[0]);
            int targetIndex = Integer.parseInt(data[1]);
            double weight = Double.parseDouble(data[2]);
            Vertex source = nodes.get(sourceIndex);
            if(source == null){
                source = new Vertex(sourceIndex);
                nodes.put(sourceIndex, source);
                graph.addVertex(source);
            }
            Vertex target = nodes.get(targetIndex);
            if(target == null){
                target = new Vertex(targetIndex);
                nodes.put(targetIndex, target);
                graph.addVertex(target);
            }

            DefaultWeightedEdge e = graph.addEdge(source, target);
            graph.setEdgeWeight(e, weight);
        }

        if (nodeNum != graph.vertexSet().size()){
            throw new Error("instance vertices number error");
        }
        if (edgeNum != graph.edgeSet().size()){
            throw new Error("instance edges number error");
        }
        br.close();
    }

    private void init_X_graph(){
        X_graph = new SimpleWeightedGraph<>(DefaultWeightedEdge.class);
        for(Vertex v : graph.vertexSet()){
            X_graph.addVertex(v);
            v.degree_to_X_graph = graph.degreeOf(v);
            v.is_in_X_graph = true;
        }
        for(DefaultWeightedEdge e : graph.edgeSet()){
            Vertex source = graph.getEdgeSource(e);
            Vertex target = graph.getEdgeTarget(e);
            double weight = graph.getEdgeWeight(e);
            DefaultWeightedEdge ee = X_graph.addEdge(source, target);
            X_graph.setEdgeWeight(ee, weight);
        }
    }


    private int depth;
    private int num_root_child;
    private Vertex root;
    private void find_cut_vertices () {
        depth = 1;
        num_root_child = 0;
        for(Vertex v : X_graph.vertexSet()){
            v.is_visited = false;
            v.dep = -1;
            v.low = -1;
            v.is_cut = false;
        }
        Vertex r = X_graph.vertexSet().iterator().next();
        r.is_visited = true;
        root = r;
        cut_vertices_recur(r);
        if (num_root_child > 1){
            r.is_cut = true;
        }
    }

    private void cut_vertices_recur(Vertex r){
        r.is_visited = true;
        r.dep = depth;
        r.low = depth;
        ++depth;

        for(DefaultWeightedEdge e : graph.edgesOf(r)){
            Vertex tem = getOtherEdgeEnd(e, r);
            if(tem.is_in_X_graph){
                if(!tem.is_visited){
                    cut_vertices_recur(tem);
                    r.low = Math.min(r.low, tem.low);
                    if (tem.low >= r.dep && r != root){
                        r.is_cut = true;
                    }else if (r == root){
                        ++num_root_child;
                    }
                }else{
                    r.low = Math.min(r.low, tem.dep);
                }
            }
        }
    }

    private void init(){
        init_X_graph();
    }

    public void solve(){
        init();
        for(;;) {
            shrink_X_graph();
            local_search(true);
            if(!X_minu.isEmpty())break;
        }
    }

    private void updateBestSolution(){
        if(X_minu.size() < bestSolution.X_minu_size){
            bestSolution = new Solution(X_graph, X_minu.size(), currSpanningTree);
        }else if(X_minu.size() == bestSolution.X_minu_size){
            if(bestSolution.X_minu_size == 0){
                if(currSpanningTree.getWeight() < bestSolution.spanningTree.getWeight()){
                    bestSolution = new Solution(X_graph, X_minu.size(), currSpanningTree);
                }
            }
        }
    }

    private void local_search(boolean isSimpleLS){
        boolean isImprovingPhase = true;
        bestSolution = new Solution(X_graph, X_minu.size(), currSpanningTree);
        for(;;++iter_count){
            Move mv = findMove();
            int imprvType = isImprovingMv(mv);
            if(isSimpleLS && imprvType == 0)break;

            if(isImprovingPhase && imprvType == 0){
                isImprovingPhase = false;
            }else if(!isImprovingPhase && imprvType == 2){
                isImprovingPhase = true;
            }

            makeMove(mv);
            check_configuration();
            updateBestSolution();

            System.out.println("iter: " + iter_count
                    + ", Xs.Size= " + X_graph.vertexSet().size()
                    + ", Xm.Size = " + X_minu.size()
                    + (currSpanningTree == null ? "" : ", treeWeight = " + currSpanningTree.getWeight()));
        }
    }

    /***
     *
     * @param mv
     * @return 0 if it's not improving, 1 if it's improving, 2 if it improves the bestSolution
     */
    private int isImprovingMv(Move mv){
        if(mv.deltaX > 0){
            return 0;
        }else if (mv.deltaX < 0){
            if(mv.deltaX + X_minu.size() < bestSolution.X_minu_size){
                return 2;
            }else {
                return 1;
            }
        }else {
            if(mv.deltaX + X_minu.size() == 0) {
                double mvTreeWeight = mv.spanningTree.getWeight();
                if(mvTreeWeight < bestSolution.spanningTree.getWeight()){
                    return 2;
                }else if( mvTreeWeight < currSpanningTree.getWeight()){
                    return 1;
                }else{
                    return 0;
                }
            }else{
                return 0;
            }
        }
    }

    private Vertex getOtherEdgeEnd(DefaultWeightedEdge e, Vertex v){
        Vertex source = graph.getEdgeSource(e);
        if(source == v){
            return graph.getEdgeTarget(e);
        }else{
            return source;
        }
    }

    private void X_graph_remove_v(Vertex v, boolean update){
        if(X_graph.vertexSet().isEmpty()){
            throw new Error("X_graph vertex set is empty!");
        }
        X_graph.removeVertex(v);

        if(update) {
            v.is_in_X_graph = false;
            X_plus.add(v);
            for (DefaultWeightedEdge e : graph.edgesOf(v)) {
                Vertex u = getOtherEdgeEnd(e, v);
                --u.degree_to_X_graph;
                if (u.degree_to_X_graph == 0 && !u.is_in_X_graph) {
                    X_plus.remove(u);
                    X_minu.add(u);
                }
            }
        }
    }

    private void X_graph_insert_v(Vertex v, boolean update){
        X_graph.addVertex(v);
        List<Vertex> u_in_X_graph = null;
        List<Vertex> u_out_X_graph = null;
        if(update){
            u_in_X_graph = new ArrayList<>();
            u_out_X_graph = new ArrayList<>();
        }
        for(DefaultWeightedEdge e : graph.edgesOf(v)){
            Vertex u = getOtherEdgeEnd(e, v);
            if(u.is_in_X_graph){
                DefaultWeightedEdge edge = X_graph.addEdge(u, v);
                X_graph.setEdgeWeight(edge, graph.getEdgeWeight(e));
                if(update){
                    u_in_X_graph.add(u);
                }
            }else if(update){
                u_out_X_graph.add(u);
            }
        }

        if(update) {
            X_plus.remove(v);
            v.is_in_X_graph = true;

            for(Vertex u : u_in_X_graph){
                ++u.degree_to_X_graph;
            }
            for(Vertex u : u_out_X_graph){
                ++u.degree_to_X_graph;
                if(u.degree_to_X_graph == 1) {
                    X_minu.remove(u);
                    X_plus.add(u);
                }
            }
        }
    }

    private void shrink_X_graph(){
        find_cut_vertices();
        List<Vertex> nonCutPs = X_graph.vertexSet().stream().filter(a->!a.is_cut).collect(Collectors.toList());
        Vertex randomV = nonCutPs.get(random.nextInt(nonCutPs.size()));
        X_graph_remove_v(randomV, true);
        if(X_minu.isEmpty()){
            currSpanningTree = (new KruskalMinimumSpanningTree<>(X_graph)).getSpanningTree();
        }
    }

    private void check_set_split(){
        Set<Vertex> tmp = new TreeSet<>(X_graph.vertexSet());
        tmp.retainAll(X_plus);
        if(!tmp.isEmpty()){
            throw new Error("X^* and X^+ share elements");
        }

        tmp = new TreeSet<>(X_graph.vertexSet());
        tmp.retainAll(X_minu);
        if(!tmp.isEmpty()){
            throw new Error("X^* and X^- share elements");
        }

        tmp = new TreeSet<>(X_plus);
        tmp.retainAll(X_minu);
        if(!tmp.isEmpty()){
            throw new Error("X^+ and X^- share elements");
        }
    }

    private void check_consistency(){
        for (Vertex v: graph.vertexSet()){
            if (X_graph.vertexSet().contains(v)){
                if(!v.is_in_X_graph){
                    throw new Error("check_consistency 1");
                }
            }else if(X_plus.contains(v)){
                boolean is_dominated = false;
                for(DefaultWeightedEdge e : graph.edgesOf(v)){
                    Vertex u = getOtherEdgeEnd(e, v);
                    if(X_graph.vertexSet().contains(u)){
                        is_dominated = true;
                        break;
                    }
                }
                if(!is_dominated){
                    throw new Error("check_consistency 2");
                }
            }else if(X_minu.contains(v)){
                for(DefaultWeightedEdge e : graph.edgesOf(v)){
                    Vertex u = getOtherEdgeEnd(e, v);
                    if(X_graph.vertexSet().contains(u)){
                        throw new Error("check_consistency 3");
                    }
                }
            }else{
                throw new Error("check_consistency 4");
            }

            int d_x = 0;
            for(DefaultWeightedEdge e : graph.edgesOf(v)){
                Vertex u = getOtherEdgeEnd(e, v);
                if(X_graph.vertexSet().contains(u)){
                    ++d_x;
                }
            }
            if(d_x != v.degree_to_X_graph){
                throw new Error("check_consistency 5");
            }
        }
    }

    private void dfs_Xs(Vertex r){
        r.is_visited = true;
        for(DefaultWeightedEdge e : graph.edgesOf(r)){
            Vertex u = getOtherEdgeEnd(e, r);
            if(r.is_in_X_graph && !u.is_visited){
                dfs_Xs(u);
            }
        }
    }

    private void check_spanningTree(){
        for(Vertex v : X_graph.vertexSet()){
            v.is_visited = false;
        }

        Vertex v = X_graph.vertexSet().iterator().next();
        dfs_Xs(v);
        for(Vertex u : X_graph.vertexSet()){
            if(!u.is_visited){
                throw new Error("X_graph is not connected! ");
            }
        }

        if(!X_minu.isEmpty())return;

        if(currSpanningTree.getEdges().size() != X_graph.vertexSet().size() - 1){
            throw new Error("X_graph is not a tree! ");
        }
        for(Object o : currSpanningTree.getEdges()){
            DefaultWeightedEdge e = (DefaultWeightedEdge) o;
            Vertex source = X_graph.getEdgeSource(e);
            Vertex target = X_graph.getEdgeTarget(e);
            if(!graph.containsEdge(source, target)){
                throw new Error("spanning tree edge does not exist");
            }
        }
    }

    private void check_configuration(){
        check_set_split();
        check_consistency();
        check_spanningTree();
    }

    private int calcDeltaX(Vertex u, Vertex v){
        int deltaX = 0;
        for(DefaultWeightedEdge e : graph.edgesOf(u)){
            Vertex uu = getOtherEdgeEnd(e, u);
            if(uu.degree_to_X_graph == 0){
                --deltaX;
            }
        }

        for(DefaultWeightedEdge e : graph.edgesOf(v)){
            Vertex vv = getOtherEdgeEnd(e, v);
            if(vv.degree_to_X_graph == 1 && !graph.containsEdge(vv, u)){
                ++deltaX;
            }
        }
        return deltaX;
    }

    SpanningTree<DefaultWeightedEdge> calcNewWeight(Vertex iV, Vertex rV){

        X_graph_insert_v(iV, false);
        X_graph_remove_v(rV, false);

        var spanningTree = (new KruskalMinimumSpanningTree<>(X_graph)).getSpanningTree();

        X_graph_insert_v(rV, false);
        X_graph_remove_v(iV, false);
        check_configuration();
        return spanningTree;
    }

    private Move findMove(){
        Move bestMv = new Move(null, null, Integer.MAX_VALUE, null);
        find_cut_vertices();
        List<Vertex> iVs = X_graph.vertexSet().stream().filter(a->!a.is_cut).collect(Collectors.toList());
        for(Vertex u : X_plus){
            for(Vertex v : iVs){
                int deltaX = calcDeltaX(u,v);
                SpanningTree<DefaultWeightedEdge> spanningTree = null;
                if(deltaX + X_minu.size() == 0){
                    spanningTree = calcNewWeight(u, v);
                }
                Move mv = new Move(u, v, deltaX, spanningTree);
                if(bestMv.compareTo(mv) > 0){
                    bestMv = mv;
                }
            }
        }
        return bestMv;
    }

    private void makeMove(Move mv){
        X_graph_insert_v(mv.insertV, true);
        X_graph_remove_v(mv.removeV, true);
        currSpanningTree = mv.spanningTree;
        int tt = tabu_tenure_min + random.nextInt(tabu_tenure_max - tabu_tenure_min);
        tabuList.put(mv.removeV, tt);
    }


    class Solution implements Comparable{
        SimpleWeightedGraph<Vertex, DefaultWeightedEdge> X_graph;
        int X_minu_size;
        SpanningTree<DefaultWeightedEdge> spanningTree;

        Solution(SimpleWeightedGraph<Vertex, DefaultWeightedEdge> Xg, int X_minu_size, SpanningTree<DefaultWeightedEdge> sTree){
            X_graph = (SimpleWeightedGraph<Vertex, DefaultWeightedEdge>)Xg.clone();
            this.X_minu_size = X_minu_size;
            spanningTree = sTree;
        }

        @Override
        public int compareTo(Object o) {
//            if(o == null){
//                return -1;
//            }
//            if(o.getClass() != Solution.class){
//                throw new ClassCastException();
//            }
            Solution oo = (Solution)o;
            if(X_minu_size < oo.X_minu_size){
                return -1;
            }else if(X_minu_size > oo.X_minu_size){
                return 1;
            }else if(X_minu_size == 0){
                return Double.compare(spanningTree.getWeight(), oo.spanningTree.getWeight());
            }else{
                return 0;
            }
        }
    }

    class Move implements Comparable{
        Vertex insertV;
        Vertex removeV;
        int deltaX;
        SpanningTree<DefaultWeightedEdge> spanningTree;
        Move(Vertex iV, Vertex rV, int dX, SpanningTree<DefaultWeightedEdge> sT){
            insertV = iV;
            removeV = rV;
            deltaX = dX;
            spanningTree = sT;
        }

        @Override
        public int compareTo(Object o) {

            Move oo = (Move) o;

            if(deltaX < oo.deltaX){
                return -1;
            }else if(deltaX > oo.deltaX){
                return 1;
            }else if(spanningTree != null){
                return Double.compare(spanningTree.getWeight(), oo.spanningTree.getWeight());
            }else{
                return 0;
            }
        }
    }
}
