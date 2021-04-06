package DTPSolver;

import java.util.Objects;

public class Vertex implements Comparable{
    final int index;
    int degree_to_X = 0;
    boolean is_in_X = false;
    boolean is_cut = false;

    int tabu_iter = 0;

    //used for determining cutting point
    boolean is_visited = false;
    int dep = 0;
    int low = 0;

    Vertex(int index){
        this.index = index;
    }

    @Override
    public String toString(){
        return String.valueOf(index);
    }

    @Override
    public int compareTo(Object o) {
//        if (o == null){
//            throw new NullPointerException();
//        }
//        if(o.getClass() != Vertex.class){
//            throw new ClassCastException();
//        }

        Vertex oo = (Vertex)o;
        return Integer.compare(index, oo.index);
    }

    @Override
    public int hashCode() {
        return Objects.hash(index);
    }

    @Override
    public boolean equals(Object o){
        if (o.getClass() != Vertex.class){
            return false;
        }

        return ((Vertex)o).index == index;
    }
}
