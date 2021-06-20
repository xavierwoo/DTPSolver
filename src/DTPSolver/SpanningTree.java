package DTPSolver;

import java.util.Set;

public class SpanningTree {
    final public Set<Edge> tree_edges;
    final public double tree_weight;
    public Set<Vertex> leaves = null;

    SpanningTree(Set<Edge> edges, double weight) {
        tree_edges = edges;
        tree_weight = weight;
    }
}
