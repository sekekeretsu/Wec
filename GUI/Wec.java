
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.charset.Charset;
import java.nio.charset.StandardCharsets;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.LinkedList;
import java.util.Scanner;
import java.util.StringTokenizer;

/**
 *
 * @author seke14
 */
public class Wec {
      static String weighted_data;   //="C:\\\\coach\\\\collins2007.txt";
      static String gene_data;       //="C:\\\\coach\\\\gene.txt";
      static String benchmark_data;  //="C:\\\\coach\\\\sgd.txt";
      static double balanceT;
      static double weightT;
      static double enrichT;
      static double filterT;
      
    public static void main(String[] args) throws IOException {
         weighted_data=args[0];
         gene_data=args[1];
         benchmark_data=args[2];
         balanceT=Double.parseDouble(args[3]);
       weightT=Double.parseDouble(args[4]);
       enrichT=Double.parseDouble(args[5]);
      filterT=Double.parseDouble(args[6]);
      //   balanceT=0.7;
      //   weightT=0.3;
      //  enrichT=0.8;
      //  filterT=0.8;
        ProteinComplex PC=new ProteinComplex();
        PC.distinctProteins(weighted_data);
        PC.readBenchmarkComplex(benchmark_data);
        PC.gepreader(gene_data);
        PC.matchProtienToGene();
        PC.removeInvalidInteraction();
        PC.neighborhoodProteins();
        PC.complexIdentification(balanceT,weightT);
        PC.redundencyFilter(filterT);
        PC.complexEnrichmentProcess(enrichT);
        PC.displayPredictedComplex();
        PC.Compare();
        PC.CompareToFindRealComplex();
       
        PC.writeDetectedComplex();
        PC.displayComplexCollections();
        PC.evaluation();   
          System.out.println();
                    System.out.println("*****************************************************************************************************************************");
                    System.out.println("******************     PLEASE CHECK WORKING DIRECTORY FOR OUTPUT FILES       **************************");
                    System.out.println("*****************************************************************************************************************************");
    }
}


class ProteinComplex{ 
   ArrayList<Interactions> InteractionList=new ArrayList<>();
   ArrayList<Interactions> InteractionList_filtered=new ArrayList<>();
   ArrayList<String> dist_element=new ArrayList<>();
   ArrayList<Protein> neighbourList=new ArrayList<>();
   ArrayList<Gene>   geneList=new ArrayList<>();
   ArrayList<String> dist_genename=new ArrayList<>();
   ArrayList<String> dist_protein=new ArrayList<>();
   ArrayList<String> nomatch=new ArrayList<>();
   ArrayList<Graph> similarNeighbourList=new ArrayList<>();
   ArrayList<Graph> coreList=new ArrayList<>();
   ArrayList<Graph> coreListFiltered=new ArrayList<>();
   ArrayList<Complex> predictedComplexList=new ArrayList<>();
   ArrayList<Complex> benchmarkComplexList=new ArrayList<>();
   ArrayList<Complex> matchComplexList=new ArrayList<>(); 
    ArrayList<Complex> matchRealComplexList=new ArrayList<>();
    ArrayList<Complex> coveredRealComplex=new ArrayList<>();
    ArrayList<Complex> coveredRealComplex1=new ArrayList<>();
    ArrayList<Complex> coreComplexList=new ArrayList<>();
   
    
    ////////// PROCESS TO IDENTIFY PROTEIN CLUSTERS OR POTEINTIAL COMPLEXES   ///////////////////////////////
    public void complexIdentification(double alpha,double threshold)
   {  // double alpha=0.7 ;
      // double threshold=0.3;
       for(Protein protein:neighbourList)
       { Graph newGraph=new Graph();
        Protein oprotein=protein.proteinCopy();
        for(Node node:oprotein.neighbours)
        { Protein nprotein=findProtein(node.nname);
         double eCC=findEcc(nprotein,oprotein);
         if(eCC==0)
         {continue;}
         else
         { Gene gene1=findGene(node.nname);
           Gene gene2=findGene(oprotein.pname);
           if(gene1==null||gene2==null)
           {continue;}
           else
           {
               double pCC=evalSimilarity(gene1,gene2);
               double tweight;
               //if(pCC>=0)
               //{
               tweight=(eCC*alpha) + ((1-alpha)*pCC);    
               if(tweight>=threshold)
               { newGraph.ProteinChain.add(nprotein);
               }
               //}
               /*else
               {tweight=eCC;
                  if(tweight>=0.65)
               { newGraph.ProteinChain.add(nprotein);
               }
               }
               */
           }
         }
        }
        if(newGraph.ProteinChain.size()>=2)
        {  newGraph.arrangeGraph();
            coreList.add(newGraph);
       }
       }
   }
    
    //////////////////// FILTERING REDUNDANT COMPLEXES //////////////////////
    public void redundencyFilter(double filterT)
   { coreListFiltered.add(coreList.get(0));
      for(Graph tcore1:coreList)
      { double NA_max=0.0;
        Graph max=new Graph();
        for(Graph tcore2:coreListFiltered)
        { double na;
        na=neighbourAffinity(tcore1,tcore2);
        if(na>NA_max)
        {
            NA_max=na;
            max=tcore2;
        }
        else
        {continue;
        }
        
        }
        if(NA_max < filterT)
        { coreListFiltered.add(tcore1);
        }
        else
        {
            double den1=findGraphEcc(tcore1);
            double den2=findGraphEcc(max);
           // double size1=tcore1.ProteinChain.size();
            //double size2=max.ProteinChain.size();
            if((den1)>=(den2))
            { coreListFiltered.add(tcore1);
              coreListFiltered.remove(max);
            }
           }
          }
       System.out.println("size of filtered list:"+coreListFiltered.size());
   }
    
    ////// EDGE CLUSTERING COEFFICIENT(ECC) : PROCESS TO EVALUATE THE ECC OF TWO PROTEINS   ///////
   double findEcc(Protein nprotein,Protein oprotein)
   { double ecc=0.0;
     double countCommon=0.0;
     ArrayList<String> tempList=new ArrayList<>();
     for(Node tnode:oprotein.neighbours)
     { tempList.add(tnode.nname);
     }
     for(Node tnode:nprotein.neighbours)
     {if(tempList.contains(tnode.nname))
     {countCommon++;}
     }
     if(countCommon>=0)
     { double min=0.0;
       double sizenprotein=nprotein.neighbours.size();
       double sizeoprotein=oprotein.neighbours.size();
       if(sizenprotein<=sizeoprotein)
       {min=sizenprotein;}
       else
       { min=sizeoprotein;
       }
        ecc=countCommon/(min-1);
       }
   return ecc;
   }
   //////// TO FIND TOTAL EDGE CLUSTERING COFFICIENT OF EDGES IN A GRAPH ////////
   double findGraphEcc(Graph graph)
   { ArrayList<String> pnameList=new ArrayList<>();
     double total=0.0;
     for(Protein protein:graph.ProteinChain)
     { pnameList.add(protein.pname);
     }
     for(String tname:pnameList)
     { Protein tprotein=findProtein(tname);
      for(Node tnode:tprotein.neighbours)
      { if(pnameList.contains(tnode.nname))
      {   Protein nprotein=findProtein(tnode.nname);
          double tweight=findEcc(nprotein,tprotein);
          total=total+tweight;
      }
      }
     }
     total=total/2;
   return total;
   }
   
   ////////////// WRITING INTERACTION WITH VALID EXPRESSION DATA TO FILE ///////////////
    public void writeInteractions()
    { try{
    File f=new File("output_filtered.txt");
    f.createNewFile();
    FileWriter Fw=new FileWriter(f);
    for(Interactions interaction:InteractionList_filtered)
    {Fw.write(interaction.protein1+"\t"+interaction.protein2);
     Fw.write("\n");
    }
    Fw.flush();
    Fw.close();
    }
   catch(Exception e){
            System.out.println("EXCEPTION"+e.getMessage());
        }
    }
    
    /////////// DISPLAYING PREDICTED COMPLEXES  ////////////////////////
    public void displayPredictedComplex() throws IOException
    { System.out.println(".................displaying predicted complex:...................  "); 
    File f=new File("output_PredictedComplex.txt");
       f.createNewFile();
       FileWriter fw=new FileWriter(f); 
    for(Complex tComplex:predictedComplexList)
    { for(String tprotein:tComplex.cProtein)
      {System.out.print(tprotein+"\t"); 
      fw.write(tprotein);
      fw.write("\t");
      }
      System.out.println();
      fw.write("\n");
      }
       fw.flush();
       fw.close();
      } 
    
    
    ////////////  EVLUATION MATRICS: EVALUATING PRECISION RECALL AND F-MEASURE ///////////////
    public void evaluation()
    {
              double fmeasure;
        double detected=0;
        double predicted=0;
        double benchmark=0;
        double realMatched=0;
        double pre=0;
        double recall=0;
        double no_interaction;
        double no_protein;
        double no_gene;
    
         no_interaction=InteractionList_filtered.size();
         no_protein=dist_protein.size();
         no_gene=geneList.size();
        benchmark=benchmarkComplexList.size();
        predicted=predictedComplexList.size();
        detected=matchComplexList.size();
        realMatched=matchRealComplexList.size();
        pre=detected/predicted;
        recall=realMatched/benchmark;
        fmeasure=2*(pre*recall)/(pre+recall);
        System.out.println();
        System.out.println("............    FINAL OUTPUT       ...........");
        System.out.println();
        System.out.println("No of cores before filtering :"+coreList.size());
        System.out.println("No of cores after  filtering :"+coreListFiltered.size());
        System.out.println("no of complexes in reference data : "+benchmarkComplexList.size());
        System.out.println("no of predicted complex   : "+predictedComplexList.size());
        System.out.println("no of detected complex    : "+matchComplexList.size());
        System.out.println("no of realmatch complex   : "+matchRealComplexList.size()+"\n");
        System.out.println("                precision : "+pre);
        System.out.println("                    recall: "+recall);
        System.out.println("                f-measure : "+fmeasure);
        System.out.println();
        System.out.println();
        
       
        
    }
  
  /////////  DISPLAY SIZE OF REAL, PREDICTED ,MATCHED COMPLEXES INFORMATION///////////
      void displayComplexCollections()
     { System.out.println();
     System.out.println("Total number of real complexes in the reference data : "+benchmarkComplexList.size());
     System.out.println("Total number of predicted complexes                  : "+predictedComplexList.size());
     System.out.println("No of predicted complexes having a match in reference data  :"+matchComplexList.size());
     System.out.println("No of real complex that has match with predicted complexes  :"+coveredRealComplex.size());
    
    
     }
      
      /////// WRITE DETECTED COMPLEXES TO FILE ///////////////
         public void writeDetectedComplex() throws IOException
   {
       File f=new File("output_DetectedComplexS.txt");
       f.createNewFile();
       FileWriter fw=new FileWriter(f);
       
       System.out.println();
       System.out.println("ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo");
       System.out.println("oooooooooooooooooooooooooo  WRITEIN AND DISPLAYING DETECTED COMPLEXES    oooooooooooooooooooooooooooooooo");
       System.out.println("ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo");
       for(Complex complex:matchComplexList)
      { for(String protein:complex.cProtein)
        { System.out.print(protein);
          fw.write(protein);
        System.out.print("\t");
        fw.write("\t");}
      System.out.println();
      fw.write("\n");
      }
       fw.flush();
       fw.close();
        System.out.println("covered real complex list:"+coveredRealComplex.size());
        System.out.println("detected complex list: "+ matchComplexList.size());
        System.out.print("the no of proteins complexes:"+ predictedComplexList.size());
   }
   
         /////////////////////  TO DISPLAY DETECTED COMPLEXES OR MATCHED COMPLEXES /////////////////////
   public void displayDetectedComplex()
   {
       System.out.println();
       System.out.println("oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo");
       System.out.println("ooooooooooooooo  DISPLAYING DETECTED COMPLEXES ooooooooooooooooooooooooooooooooooooooooooooooooo");
       System.out.println("oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo");
       for(Complex complex:matchComplexList)
      { for(String protein:complex.cProtein)
        { System.out.print(protein);
        System.out.print("\t");}
      System.out.println();
      }
        System.out.println("covered real complex list:"+coveredRealComplex.size());
         System.out.println("detected complex list: "+ matchComplexList.size());
       System.out.print("the no of proteins complexes:"+ predictedComplexList.size());
   }
   
   //////// TO ADD HIGHLY CONNECTED PROTEINS A POTENTIAL COMPLEX TO ENRICH A COMPLEX ///////////////////
  public void complexEnrichmentProcess(double enrichT)
  {
   for(Graph tgraph:coreListFiltered)
   {Graph graph1=new Graph();
    Complex newComplex=new Complex();
   ///////////////////////////...linkedlist or array...?????...../////////////
    LinkedList<String> coreProteinList=new LinkedList<>();
    LinkedList<String> attachment=new LinkedList<>();
    double coreListSize=tgraph.ProteinChain.size();
    for(Protein protein1:tgraph.ProteinChain)
    { coreProteinList.add(protein1.pname);
      newComplex.cProtein.add(protein1.pname);
      Protein protein2=findProtein(protein1.pname);
      protein2=protein2.proteinCopy();  
      graph1.ProteinChain.add(protein2);
    }
    LinkedList<String> candidateProteins=new LinkedList<>();
    for(Protein protein1:graph1.ProteinChain)
    { for(Node node1:protein1.neighbours)
     { if(candidateProteins.contains(node1.nname)||coreProteinList.contains(node1.nname))
        { continue;
        }
     else
     { candidateProteins.add(node1.nname);
     }
     }
    }
    for(String pname:candidateProteins)
    { double intersection=0.0;
      for(Protein protein:graph1.ProteinChain)
      { for(Node node1:protein.neighbours)
        { if(pname.contentEquals(node1.nname))
        { intersection++;
          break;
        }
        }
      }
      if((intersection*intersection)>=(enrichT*coreListSize))
      {attachment.add(pname);
       newComplex.cProtein.add(pname);
      }
    }
    predictedComplexList.add(newComplex);
   }
  }
   
  
   /////////// TO FIND THE NEIGHBORHOOD AFFINITY OF COMPLEXES TO BE USED IN FILTERING REDUNDANCY /////////////////////
   double neighbourAffinity(Graph graph1,Graph graph2)
   {double result=0.0;
   double size1=graph1.ProteinChain.size();
   double size2=graph2.ProteinChain.size();
   double match=0.0;
   if(size1>=size2)
   {  LinkedList<String> ProteinList=new LinkedList<>();
      for(Protein protein1:graph1.ProteinChain)
      {ProteinList.add(protein1.pname);
      }
   for(Protein protein2:graph2.ProteinChain)
      { if(ProteinList.contains(protein2.pname))
      { match++;
      }
      }
   }
   else
   {
       LinkedList<String> ProteinList2=new LinkedList<>();
       for(Protein protein2:graph2.ProteinChain)
       {
           ProteinList2.add(protein2.pname);
       }
       for(Protein protein1:graph1.ProteinChain)
       { if(ProteinList2.contains(protein1.pname))
       { match++;
       }
       }
   }
   double denom=size1*size2;
   double NA=(match*match)/denom;
   result=NA;
    return result;
   }
   
   
   ///////////// TO REMOVE INTERACTIONS CONTAINING PROTEINS WITHOUR GENE EXPRESSION DATA //////////////
   public void removeInvalidInteraction()
   {   int count=0;
       System.out.println(" the size of the old [[before]] filtering Interaction LIst is:"+ InteractionList.size());
       for(int x=0;x<InteractionList.size();x++)
        {
            if(nomatch.contains(InteractionList.get(x).protein1)||nomatch.contains(InteractionList.get(x).protein2))
                {count++;continue;}
            else
                {InteractionList_filtered.add(InteractionList.get(x));
                }
        }
   System.out.println(" the size of the old [[after]] filtering Interaction LIst is:"+ InteractionList.size());
   System.out.println(" the size of the new Interaction LIst is:"+ InteractionList_filtered.size());
   System.out.println(" count :b "+count);
   
   for(Interactions i1:InteractionList_filtered)
   {  
       String name1=i1.protein1;
       String name2=i1.protein2;
      
          if(!(dist_protein.contains(name1)))
          dist_protein.add(name1);
          if(!(dist_protein.contains(name2)))
          dist_protein.add(name2);
   }
   }
    ////// FINDING PROTEIN : RETURNS A PROTEIN WITH NAME AND ITS CORRESPONDING NEIGHBORING PROTEINS.///////////
   Protein findProtein(String pname)
   {Protein protein=new Protein();
    for(int x=0;x<neighbourList.size();x++)
     {   
        if(neighbourList.get(x).pname.contentEquals(pname))
        {
         protein= neighbourList.get(x);
         break;
        }
        if(x==(neighbourList.size()-1)&&protein.pname==null)
        {System.out.println("the protein is not found:"+ pname);
           }
        }
    Protein proteinR;
    proteinR=protein.proteinCopy();
    return proteinR;
    }
   
   //////////// RETURNING A GENE AND ITS CORRESPONDING EXPRESSION VALUE //////////
   Gene findGene(String gName)
   { Gene returnGene=null;
     for(int x=0;x<geneList.size();x++)
     { Gene tgene=geneList.get(x);
       if(tgene.gname.contentEquals(gName))
       { returnGene=tgene;break;} //////////// break???
     }
   return returnGene;
   }
   
   ////// TO EVALUATE THE SIMILARITY IN EXPRESSION OF TWO PROTEINS/GENES /////////// 
  Double evalSimilarity(Gene gene1,Gene gene2)
   {Double result;
    int size=36;
    Double x,y;
    Double xSum=0.0;
    Double ySum=0.0;
    Double xSumSquare=0.0;
    Double ySumSquare=0.0;
    Double num,Denom,Denom1;
    Double xy=0.0;
    for(int i=0;i<size;i++)
    {
        x=gene1.expValue.get(i);
        y=gene2.expValue.get(i);
        xSum=xSum+x;
        ySum=ySum+y;
        xSumSquare=xSumSquare+(x*x);
        ySumSquare=ySumSquare+(y*y);
        xy=xy+(x*y);
    }
    num=(size*xy)-(xSum*ySum);
    Denom1=((size*xSumSquare)-(xSum*xSum))*((size*(ySumSquare))-(ySum*ySum));
    Denom=Math.sqrt(Denom1);
    result=num/Denom;
    
   // System.out.println("num:"+num);
    
    //System.out.println("denom:"+Denom);
    System.out.println("similarity : "+result);
    return result;
   }
   //////// FINDING PROTEINS HAVING NO GENE INFORMATION ////////
   public void matchProtienToGene()
   { for(int x=0;x<dist_element.size();x++)
   { String temp=dist_element.get(x);
     if(!(dist_genename.contains(temp)))
      { nomatch.add(temp);
        }
   }
   }
   //////////// DISPLAYING GENE INFORMATIONS ///////////////
     void displayGenes()
     { for(Gene a:geneList)
     { 
        System.out.print(a.gname);
        for(Double b:a.expValue)
        {System.out.print(" "+b);}
                System.out.println();
     }
          System.out.println("the sizeof geneList or total genes from GEP       : "+geneList.size());
          System.out.println("the sizeof distinct elements from weightd PPI     : "+dist_element.size());
          System.out.println("the sizeof  dist gene from GEP                    : "+dist_genename.size());
          System.out.println("the sizeof proteins form PPI mathing genes from GEP: "+dist_protein.size());
          System.out.println("the sizeof Proteins with no match or expression   : "+nomatch.size());
     }
     
     ////////// INPUT : READING GENE EXPRESSION DATA FROM FILE //////////////////
     public void gepreader(String fFileName) throws IOException
    {
       fFilePath_GeneExp=Paths.get(fFileName);
       try (Scanner scanner =  new Scanner(fFilePath_GeneExp, ENCODING.name())){
         while (scanner.hasNextLine()){
                Gene tempGene=new Gene();
                 String aLine=scanner.nextLine();
                 StringTokenizer st = new StringTokenizer(aLine);
                 String firstToken;
                 firstToken=st.nextToken();
                 String nameToken=st.nextToken();
                 dist_genename.add(nameToken);
                 tempGene.gname=nameToken;
                   while (st.hasMoreTokens()) {
                   tempGene.expValue.add(Double.parseDouble(st.nextToken()));
           }  
                  geneList.add(tempGene);
      }  
    }
      
    }
     /////// DISPLAYING PROEIN PROTEIN INTERACTIONS /////////
     void displayIntractions()
     {   System.out.println("size of the interactions:"+InteractionList.size());
         System.out.println("size of the interactions:"+dist_element.size());
         for(int i=0;i<InteractionList.size();i++)
         { 
          System.out.println(i+" "+InteractionList.get(i).protein1+" "+InteractionList.get(i).protein2+" "+InteractionList.get(i).value);
         }
     }
     
   //// FORMS NEIGHBORHOOD GRAPH : FIND THE NGEIGHBORING PROTEINS THAT HAS A DIRECT INTERACTION WITH THE SELECTED PROTEIN //////
    public void neighborhoodProteins() throws IOException
    {
   for(int i=0;i<dist_protein.size();i++)
    {
         String tdist_protein=dist_protein.get(i);
         Protein protein=new Protein(tdist_protein);
         protein.neighbours=findNeighbors(tdist_protein);    
         neighbourList.add(protein);
    
    }
    }
   
    ///////////// FINDING THE DIRECT NEIGHBORS OF A PROTEIN  //////////////////
    LinkedList<Node> findNeighbors(String passedProteinname) throws IOException
    {  LinkedList<Node> templl=new LinkedList<>();
       LinkedList<String> tempNeighbours=new LinkedList<>();
        for(int x=0;x<InteractionList_filtered.size();x++)
        { String name1=InteractionList_filtered.get(x).protein1;
          String name2=InteractionList_filtered.get(x).protein2;
          String value=InteractionList_filtered.get(x).value;
           
          if(name1.contentEquals(name2))
          { continue;
          }
         
          if(passedProteinname.contentEquals(name1))
          {
              if(tempNeighbours.contains(name2))
                {continue;}
              else    
                 {   Node tnode=new Node();
                     tnode.nname=name2;
                    // tnode.weight=Double.parseDouble(value);
                     templl.add(tnode);
                     tempNeighbours.add(name2);
                 }
          }
          else if(passedProteinname.contentEquals(name2))
          {
              if(tempNeighbours.contains(name1))
              {continue;}
              else
              {
              Node tnode=new Node();
              tnode.nname=name1;
             // tnode.weight=Double.parseDouble(value);
              templl.add(tnode);
              tempNeighbours.add(name1);
             }
          }
        }
        return templl;
    }
    ////// DISPLAY THE NEIGHBORHOOD PROTEINS //////////////
    void displayProteinNeighbour()
    {  for(int i=0;i<neighbourList.size();i++)
    {System.out.println( "PROTEIN ::::"+neighbourList.get(i).pname);
         for(int j=0;j<neighbourList.get(i).neighbours.size();j++)
         {System.out.print( neighbourList.get(i).neighbours.get(j).nname+" ");}
    System.out.println( "\n");
    }
    }
    
    ////////INPUT : READING THE PPI NETWORK DATA FROM FILE ////////
    public  void distinctProteins(String aFileName) throws IOException {
       
         fFilePath=Paths.get(aFileName);
         try (Scanner scanner =  new Scanner(fFilePath, ENCODING.name())){
         while (scanner.hasNextLine())
         {        Interactions tInteraction=new Interactions();
                
             String aLine=scanner.nextLine();
                         String name1=new String();
                         String name2=new String();
                         String value=new String();
                         StringTokenizer st=new StringTokenizer(aLine);
                     if (st.hasMoreTokens())
                         name1 = st.nextToken();
                      if(st.hasMoreTokens())
                         name2 = st.nextToken();
                      if(st.hasMoreTokens())
                         value=st.nextToken();
                      tInteraction.protein1=name1;
                      tInteraction.protein2=name2;
                      tInteraction.value=value;
                      InteractionList.add(tInteraction);
                  if(!dist_element.contains(name1.trim()))
                      dist_element.add(name1.trim());
                  if(!dist_element.contains(name2.trim()))
                      dist_element.add(name2.trim());
         
         }System.out.print(" end line ");
       }  
}
 
    /////////// INPUT : READING THE BENCHMARK COMPLEX DATA /////////
    public void readBenchmarkComplex(String bname) throws IOException
    { fFilePath_benchmark=Paths.get(bname);
     try(Scanner scanner=new Scanner(fFilePath_benchmark,ENCODING.name()))
     { 
     while(scanner.hasNextLine())
     { String aLine;aLine=scanner.nextLine();
       Complex tComplex=new Complex();   
     StringTokenizer stzer=new StringTokenizer(aLine);
        while(stzer.hasMoreTokens())
        {String thisString=stzer.nextToken();
         tComplex.cProtein.add(thisString);
        } benchmarkComplexList.add(tComplex);
     }
     }
    }
    /////////// FINDING EVERY REAL PROTEINS THAT HAS A MATCH WITH A PREDICTED COMPLEX //////////////////////
     public void CompareToFindRealComplex()
    {   System.out.println();
        System.out.println("...............comparing real complexes with predicted complexes...........");
        System.out.println();
        for(Complex pComplex:benchmarkComplexList)
        { LinkedList<String> pComplexElements=new LinkedList<>();
           double size1=pComplex.cProtein.size();
                for(String string:pComplex.cProtein)
                 { pComplexElements.add(string);
                 }
             
        Complex maxComplex=new Complex();
        double maxCloseness=0.0;
    for(Complex bComplex:predictedComplexList)
        {    double match=0.0;
             for(String string1:bComplex.cProtein)
              {    
                 if(pComplexElements.contains(string1))
                 {     match++;
                 }
             }
             double size2=bComplex.cProtein.size();
             double prod1=match*match;
             double prod2=size1*size2;
             double closeness=prod1/prod2;
             if(closeness>maxCloseness)
             { maxCloseness=closeness;
               maxComplex=bComplex;
             }
        }
    if(maxCloseness>0.255)
    { matchRealComplexList.add(pComplex);
      coveredRealComplex1.add(maxComplex);
    }
    }
    }
     /////// FINDING EVERY PREDICTED COMPLEXES THAT HAS A MATCH IN THE REAL COMPLEX BENCHMARK DATA OR REFERENCE DATA ////////
    public void Compare()
    { System.out.println("............comparing predicted complexes with real complexes............");
   
        for(Complex pComplex:predictedComplexList)
        { LinkedList<String> pComplexElements=new LinkedList<>();
           double size1=pComplex.cProtein.size();
                for(String string:pComplex.cProtein)
                 { pComplexElements.add(string);
                 }
             
        Complex maxComplex=new Complex();
        double maxCloseness=0.0;
    for(Complex bComplex:benchmarkComplexList)
        {    double match=0.0;
             for(String string1:bComplex.cProtein)
              {    
                 if(pComplexElements.contains(string1))
                 {     match++;
                 }
             }
             double size2=bComplex.cProtein.size();
             double prod1=match*match;
             double prod2=size1*size2;
             double closeness=prod1/prod2;
             if(closeness>maxCloseness)
             { maxCloseness=closeness;
               maxComplex=bComplex;
             }
        }
    if(maxCloseness>0.255)
    { matchComplexList.add(pComplex);
      coveredRealComplex.add(maxComplex);
    }
    }
    }
 private Path fFilePath_benchmark;
    private Path fFilePath_GeneExp;
    private static Path fFilePath;
    private final static Charset ENCODING = StandardCharsets.UTF_8;
}////..............CLASS CLOSED HERE.............//////
