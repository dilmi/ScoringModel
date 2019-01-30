/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package lncRNAMutations;

import static java.lang.Double.NaN;
import java.util.Random;
import java.util.Vector;

/**
 *
 * @author dilmiperera
 */
public class GeneMapBootstrap {
    String lncRNAToGene="";
    String geneExpression="";
    int distanceForMap=10000;
    
   
    public static void main(String[] args) {
        getCorrelation();
        //System.wait();
        bootstrap(1000);
    }
    
    public static void checkName() {
        String lncRNA[][]=DataSetReader.readDataSet("C:\\Users\\Dilmi\\Desktop\\ncRNA\\lncRNADisease-red.txt", 7, "\t");
        String genes[][]=DataSetReader.readDataSet("C:\\Users\\Dilmi\\Desktop\\ncRNA\\Gene.txt", 1, "\t");
         boolean found=false;
        for (int i = 0; i < lncRNA.length; i++) {
            found=false;
            
            for (int j = 0; j < genes.length; j++) {
                //System.out.println(""+ lncRNA[i][4]);
                if(genes[j][0].equalsIgnoreCase(lncRNA[i][4])){
                    lncRNA[i][6]=genes[j][0];
                    found=true;
                    break;
                }
                
            }
            if(!found && lncRNA[i][5]!=null) {
                System.out.println(""+lncRNA[i][5]);
                String alias[]=lncRNA[i][5].split(";");
                for (int j = 0; j < genes.length; j++) {
                    for (int k = 0; k < alias.length; k++) {
                        if(genes[j][0].equalsIgnoreCase(alias[k].trim())){
                            lncRNA[i][6]=genes[j][0];
                            break;
                        }
                    }
                }
                    
            }
        }
        DataSetWriter.writeToFile2D("C:\\Users\\Dilmi\\Desktop\\ncRNA\\lncRNADisease-red-gene.txt", lncRNA);
 
    }
    
    public static void geneMap(String[] args) {
        String lncRNA[][]=DataSetReader.readDataSet("", 4, "\t");
        String genes[][]=DataSetReader.readDataSet("", 4, "\t");
        
        
        for (int i = 0; i < 10; i++) {
            
        }
        
        for (int i = 0; i < lncRNA.length; i++) {
            
        }
               
    }
    
    public static void getCorrelation() {
        String lncRNAtoGene[][]=DataSetReader.readDataSet("C:\\Users\\Dilmi\\Desktop\\ncRNA\\lncRNADisease-red-gene-edited-100kb-mapped-to-GENCODE-reduced.txt", 3, "\t");
        String geneExp[][]=DataSetReader.readDataSet("C:\\Users\\Dilmi\\Desktop\\ncRNA\\expressionData_allSamples2.txt", 376, "\t");
        
        for (int i = 0; i < lncRNAtoGene.length; i++) {
            Vector<Double> lncRNA=new Vector<Double>();
            Vector<Double> gene=new Vector<Double>();
            String lncRNAname=lncRNAtoGene[i][0];
            String geneName=lncRNAtoGene[i][1];
            for (int j = 1; j < geneExp.length; j++) {
                if(lncRNAname.equalsIgnoreCase(geneExp[j][0])){
                    //lncRNA=new String[geneExp[0].length-1];
                    for (int k = 1; k < geneExp[0].length; k++) {
                        lncRNA.add(Double.parseDouble(geneExp[j][k]));
                    }
                }
                if(geneName.equalsIgnoreCase(geneExp[j][0]) || geneName.replace("-", "").equalsIgnoreCase(geneExp[j][0])){
                    //gene=new String[geneExp[0].length-1];
                    for (int k = 1; k < geneExp[0].length; k++) {
                        gene.add(Double.parseDouble(geneExp[j][k]));
                    }
                }
            }
            if(!lncRNA.isEmpty() && !gene.isEmpty() ){
                double corr = calcCorrelation(lncRNA, gene);
                lncRNAtoGene[i][2]=String.valueOf(corr);
            }
        }
        DataSetWriter.writeToFile2D("C:\\Users\\Dilmi\\Desktop\\ncRNA\\lncRNADisease-100kb-GENCODE-correlation.txt", lncRNAtoGene);
 
    }

    
    public static void bootstrap(int iterations) {
        
        String lncRNAtoGene[][]=DataSetReader.readDataSet("C:\\Users\\Dilmi\\Desktop\\ncRNA\\lncRNADisease-100kb-GENCODE-correlation.txt", 3+iterations, "\t");
        String geneExp[][]=DataSetReader.readDataSet("C:\\Users\\Dilmi\\Desktop\\ncRNA\\expressionData_allSamples2.txt", 376, "\t");
        String geneList[][]=DataSetReader.readDataSet("C:\\Users\\Dilmi\\Desktop\\ncRNA\\Gene.txt", 1, "\t");
        for (int n = 0; n < iterations; n++) {
        
        for (int i = 0; i < lncRNAtoGene.length; i++) {
            
            if(!lncRNAtoGene[i][2].equals("null")){
                
            Vector<Double> lncRNA=new Vector<Double>();
            Vector<Double> gene=new Vector<Double>();
            String lncRNAname=lncRNAtoGene[i][0];
            int random=new Random().nextInt(geneList.length - 1) + 1;
            String geneName=geneList[random][0];
            for (int j = 1; j < geneExp.length; j++) {
                if(lncRNAname.equalsIgnoreCase(geneExp[j][0])){
                    //lncRNA=new String[geneExp[0].length-1];
                    for (int k = 1; k < geneExp[0].length; k++) {
                        try{
                            lncRNA.add(Double.parseDouble(geneExp[j][k]));
                        }
                        catch(Exception e){
                            lncRNA.add(NaN);
                        }
                        
                    }
                }
                if(geneName.equalsIgnoreCase(geneExp[j][0]) || geneName.replace("-", "").equalsIgnoreCase(geneExp[j][0])){
                    //gene=new String[geneExp[0].length-1];
                    for (int k = 1; k < geneExp[0].length; k++) {
                        
                        if(geneExp[j][k]!=null){
                            gene.add(Double.parseDouble(geneExp[j][k]));
                        }
                        
                    }
                }
            }
            if(!lncRNA.isEmpty() && !gene.isEmpty() && lncRNA.size()==gene.size()){
                double corr = calcCorrelation(lncRNA, gene);
                lncRNAtoGene[i][n+3]=String.valueOf(corr);
            }}
        }
        
        }
        
        DataSetWriter.writeToFile2D("C:\\Users\\Dilmi\\Desktop\\ncRNA\\lncRNADisease-100kb-GENCODE-correlation-bootstrap.txt", lncRNAtoGene);
        
        //String lncRNAtoGene[][]=DataSetReader.readDataSet("C:\\Users\\Dilmi\\Desktop\\ncRNA\\lncRNADisease-100kb-GENCODE-correlation-bootstrap.txt", 3+iterations, "\t");
        
        String lncRNAtoGenepval[][]=DataSetReader.readDataSet("C:\\Users\\Dilmi\\Desktop\\ncRNA\\lncRNADisease-100kb-GENCODE-correlation.txt", 4, "\t");
        
        for (int i = 0; i < lncRNAtoGene.length; i++) {
            
            if(!lncRNAtoGene[i][2].equals("null")){
                double value=Double.parseDouble(lncRNAtoGene[i][2]);
                Vector<Double> dist=new Vector<Double>();
                for (int j = 3; j < iterations; j++) {
                    if(lncRNAtoGene[i][j]!=null && !lncRNAtoGene[i][j].equalsIgnoreCase("NaN") &&!lncRNAtoGene[i][j].equalsIgnoreCase("null") ){
                        dist.add(Double.parseDouble(lncRNAtoGene[i][j]));
                    }
                    
                   
                }
                lncRNAtoGenepval[i][3]=String.valueOf(pValCalc(dist,value));
            }
            }
        DataSetWriter.writeToFile2D("C:\\Users\\Dilmi\\Desktop\\ncRNA\\lncRNADisease-100kb-GENCODE-correlation-bootstrap-with-pval.txt", lncRNAtoGenepval);
        
    }
    
    public static double pValCalc(Vector<Double> distV,double value) {
        
        
        double pVal=0.0;
        double sum=0;
        double sumsq=0;
        int n=distV.size();
        double t;
        int df;
        double x;
        for (int i = 0; i < distV.size(); i++) {
            x=distV.get(i);
            sum=sum+x;
            sumsq=sumsq+x*x;
        }
        
        double mean=sum/n;
        double std=Math.sqrt((sumsq-(sum * sum)/n)/(n-1));
        
        t = ((value - mean)/std)*Math.sqrt(n);
        System.out.println(""+t);
        df = n - 1;
        pVal=calcPval(t, df, 2);
        return pVal;
    }
    
        public static double calcPval(double tva, int df, int side) {
        double ttp = TtoP(tva, df);
        if (side == 1) {
            return ttp / 2;
        } else {
            return ttp;
        }
    }

    public static double TtoP(double t, int df) {
        double tsq = t * t;
        double p = 0;
        double abst = Math.abs(t);
        if (df == 1) {
            p = 1 - 2 * Math.atan(abst) / Math.PI;
        } else if (df == 2) {
            p = 1 - abst / Math.sqrt(tsq + 2);
        } else if (df == 3) {
            p = 1 - 2 * (Math.atan(abst / Math.sqrt(3)) + abst * Math.sqrt(3) / (tsq + 3)) / Math.PI;
        } else if (df == 4) {
            p = 1 - abst * (1 + 2 / (tsq + 4)) / Math.sqrt(tsq + 4);
        } else {
            double z = TtoZ(abst, df);
            if (df > 4) {
                p = Norm_p(z);
            } else {
                p = Norm_p(z);
            }
        }
        return p;
    }

    public static double TtoZ(double t, int df) {
        double A9 = df - 0.5;
        double B9 = 48 * A9 * A9;
        double T9 = t * t / df, Z8, P7, B7, z;

        if (T9 >= 0.04) {
            Z8 = A9 * Math.log(1 + T9);
        } else {
            Z8 = A9 * (((1 - T9 * 0.75) * T9 / 3 - 0.5) * T9 + 1) * T9;
        }
        P7 = ((0.4 * Z8 + 3.3) * Z8 + 24) * Z8 + 85.5;
        B7 = 0.8 * Math.pow(Z8, 2) + 100 + B9;
        z = (1 + (-P7 / B7 + Z8 + 3) / B9) * Math.sqrt(Z8);
        return z;
    }

    public static double Norm_p(double z) {
        double absz = Math.abs(z);
        double a1 = 0.0000053830;
        double a2 = 0.0000488906;
        double a3 = 0.0000380036;
        double a4 = 0.0032776263;
        double a5 = 0.0211410061;
        double a6 = 0.0498673470;
        double p = (((((a1 * absz + a2) * absz + a3) * absz + a4) * absz + a5) * absz + a6) * absz + 1;
        p = Math.pow(p, -16);
        return p;
    }
    
    public static double calcCorrelation(Vector<Double> xVect, Vector<Double> yVect) {
    double meanX = 0.0, meanY = 0.0;
    for(int i = 0; i < xVect.size(); i++)
    {
        meanX += xVect.elementAt(i);
        meanY += yVect.elementAt(i);
    }

    meanX /= xVect.size();
    meanY /= yVect.size();

    double sumXY = 0.0, sumX2 = 0.0, sumY2 = 0.0;
    for(int i = 0; i < xVect.size(); i++)
    {
      sumXY += ((xVect.elementAt(i) - meanX) * (yVect.elementAt(i) - meanY));
      sumX2 += Math.pow(xVect.elementAt(i) - meanX, 2.0);
      sumY2 += Math.pow(yVect.elementAt(i) - meanY, 2.0);
    }

    return (sumXY / (Math.sqrt(sumX2) * Math.sqrt(sumY2)));
  }
    
}
