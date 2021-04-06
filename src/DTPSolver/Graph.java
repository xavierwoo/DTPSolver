package DTPSolver;

import java.util.ArrayList;

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
        edge_list[source.index].add(edge);
        edge_list[sink.index].add(edge);
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
}
