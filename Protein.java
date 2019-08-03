
import java.util.ArrayList;
import java.util.Iterator;
import java.util.LinkedList;



public class Protein {
    String pname;
    LinkedList<Node> neighbours=new LinkedList<>();
Protein()
{}
    
Protein(String name)
{pname=name;
 //neighbours=new LinkedList<>();
} 

///////// PROTEIN CLONING: CREATE AN IDENTICAL COPY OF PROTEIN /////////////
public Protein proteinCopy()
{Protein newProtein=new Protein();
 newProtein.pname=this.pname;
 for(Node tnode:this.neighbours)
 { Node newNode=new Node();
   newNode.nname=tnode.nname;
   newProtein.neighbours.add(newNode);
 }
   return newProtein;
}



  
}

class Complex
{ArrayList<String> cProtein=new ArrayList<>();
 ArrayList<String>core=new ArrayList<>();
 ArrayList<String>attachment=new ArrayList<>();
}
class Pair
{String pair1;
 String pair2;
 double weight;
}
class Node
{
String nname;
}
class Interactions
{ 
    String protein1;
    String protein2;
    String value;
}
class Gene{
String gname;
ArrayList<Double> expValue=new ArrayList<>();

}
class Graph {
ArrayList<Protein> ProteinChain=new ArrayList<>();

/////////// CREATES AN IDENTICAL COPY OF A GRAPH //////////
public Graph graphCopy()
{ Graph tgraph=new Graph();
  for(Protein protein:ProteinChain)
  { Protein tProtein;
    tProtein=protein.proteinCopy();
    tgraph.ProteinChain.add(protein);
  }
return tgraph;
}
///// TO FIND THE DENSITY OF A GRAPH  /////////////////
public double find_density()
{   double density=0.0;
    double tweight=0.0;
    double n=ProteinChain.size();
    for (Protein ProteinChain1 : ProteinChain) 
    { double pweight=0.0;
      for(Node tempNode:ProteinChain1.neighbours)
      {pweight=pweight++;
      }
      tweight=tweight+pweight;

    }
    density= tweight/(n*(n-1));
 return density;
}


/////////////////// ARRANGING GRAPH AGTER ADDING OR  REMOVING VERTICES//////// 
public void arrangeGraph()
{LinkedList<String> tlist=new LinkedList<>();
 for(Protein protein:ProteinChain)
 { tlist.add(protein.pname);
 }
for(Protein protein:ProteinChain)
 { LinkedList<Node> nodeList=new LinkedList<>();
   for(Node node:protein.neighbours)
   { if(tlist.contains(node.nname))
   { nodeList.add(node);
   }
   }
   protein.neighbours=nodeList;
 }
}
///////////// ARRANGING GRAPH SIMILAR TO arrangeGraph() /////////////
public Graph arrange()
{   ArrayList<String> graphProteins=new ArrayList<>();
    ArrayList<Pair> pairList=new ArrayList<>();
    Graph graph=new Graph();
    
   for(Protein protein:ProteinChain)
   { String nameP=protein.pname;
     graphProteins.add(nameP);
     for(Node node:protein.neighbours)
     { Pair tpair=new Pair();
       tpair.pair1=nameP;
       tpair.pair2=node.nname;
       pairList.add(tpair);
    } 
   }

     for(String nameP:graphProteins)
     { Protein protein=new Protein();
       protein.pname=nameP;
       ArrayList<String> tempNames=new ArrayList<>();
       for(Pair pair:pairList)
       { if(pair.pair1.contentEquals(nameP))
       { if(graphProteins.contains(pair.pair2))
           if(!tempNames.contains(pair.pair2))
           { Node node=new Node();
             node.nname=pair.pair2;
             protein.neighbours.add(node);
             tempNames.add(node.nname);
            }
       }
       else if(pair.pair2.contentEquals(nameP))
       { if(graphProteins.contains(pair.pair1))
        {if(!tempNames.contains(pair.pair1))
        { Node node=new Node();
          node.nname=pair.pair1;
          protein.neighbours.add(node);
          tempNames.add(node.nname);
        }}
       }
       }
       graph.ProteinChain.add(protein);
     }
return graph;
}
}