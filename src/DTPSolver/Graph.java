package DTPSolver;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.Set;
import java.util.TreeSet;

public class Graph {
    public Vertex[] vertices;
    public ArrayList<Edge>[] edge_list;

    public Graph(int vertex_size){
        vertices = new Vertex[vertex_size];
        edge_list = new ArrayList[vertex_size];

        for(int i=0; i<edge_list.length; ++i){
            edge_list[i] = new ArrayList<>();
        }
    }

    public Vertex get_vertex_by_index(int index){
        return vertices[index];
    }

    public void add_vertex(Vertex vertex){
        vertices[vertex.index] = vertex;
    }

    public void add_edge(Vertex source, Vertex sink, double weight){

        Edge edge = new Edge(source, sink, weight);
        add_edge(edge);
    }

    public void add_edge(Edge e){
        if(get_edge(e.source, e.sink) != null){
            return;
        }
        edge_list[e.source.index].add(e);
        edge_list[e.sink.index].add(e);
    }

    public void remove_edge(Edge e){
        edge_list[e.source.index].remove(e);
        edge_list[e.sink.index].remove(e);
    }

    public Edge get_edge(Vertex source, Vertex sink){
        for(Edge e : edge_list[source.index]){
            Vertex u = e.getOtherEdgeEnd(source);
            if(u.index == sink.index){
                return e;
            }
        }
        return null;
    }

    public double get_edge_weight(Vertex source, Vertex sink){
        for(Edge e : edge_list[source.index]){
            Vertex u = e.getOtherEdgeEnd(source);
            if(u.index == sink.index){
                return e.weight;
            }
        }
        return Double.MAX_VALUE;
    }

    public void trim_to_size(){
        for(var list : edge_list){
            list.trimToSize();
        }
    }

    public Graph gen_new_subgraph_from_tree(SpanningTree tree){
        Graph new_graph = new Graph(vertices.length);
        System.arraycopy(vertices, 0, new_graph.vertices, 0, vertices.length);
        for(var e : tree.tree_edges){
            new_graph.add_edge(e);
        }
        return new_graph;
    }

    public Set<Edge> get_all_edges_des(){
        var edges = new TreeSet<>(new DescentCmp());
        for(var list : edge_list){
            edges.addAll(list);
        }
        return edges;
    }

    public Set<Edge> get_all_edges_inc(){
        var edges = new TreeSet<Edge>();
        for(var list : edge_list){
            edges.addAll(list);
        }
        return edges;
    }

    class DescentCmp implements Comparator<Edge>{
        @Override
        public int compare(Edge e1, Edge e2){
            int cmp = Double.compare(e2.weight, e1.weight);
            if(cmp != 0){
                return cmp;
            }else{
                cmp = Integer.compare(e1.source.index, e2.source.index);
                if (cmp != 0){
                    return cmp;
                }else {
                    return Integer.compare(e1.sink.index, e2.sink.index);
                }
            }
        }
    }
}
