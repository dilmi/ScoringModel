/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package lncRNAMutations;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Vector;
import org.apache.commons.math3.stat.inference.MannWhitneyUTest;

/**
 *
 * @author Dilmi
 */
public class GeneMaplncaNet {
    //Lung - 570 , 107
    //Breast - 1216 , 214
    //Colon - 506 , 88
    public static void main(String[] args) {
//        mapExpression();
//        mapMutations();
        compareGroups();
//        mapMiTlncRNAtolncaNet();
    }
    
//    public static void mapMiTlncRNAtolncaNet() {
//        String lncRNAtoGene[][]=DataSetReader.readDataSet("C:\\Users\\Dilmi\\Desktop\\Lowy\\Non-codingRNA\\Lung\\LUAD.cancercensus.lnc.txt", 7, "\t");
//        String lncRNAID[][]=DataSetReader.readDataSet("C:\\Users\\Dilmi\\Desktop\\Lowy\\Non-codingRNA\\Lung\\combined-Lung-TCGA-in-miT-UCSC-all-sample-lncRNA.txt",5 , "\t");
//        String lncRNAtoGenetoSample[][]=new String[lncRNAID.length][3];
//                int count=0;
//        for (int i = 0; i < lncRNAID.length; i++) {
//            if(i%1000==0){System.out.println(""+i);}
//            for (int j = 0; j < lncRNAtoGene.length; j++) {
//                
//                if(lncRNAID[i][2].equals(lncRNAtoGene[j][2])){
//                    //System.out.println("True");
//                    lncRNAtoGenetoSample[count][0]=lncRNAID[i][0];
//                    lncRNAtoGenetoSample[count][1]=lncRNAtoGene[j][1];
//                    lncRNAtoGenetoSample[count][2]=lncRNAID[i][2];
//                    count++;
//                    break;
//                }
//                
//            }
//        }
//        DataSetWriter.writeToFile2D("C:\\Users\\Dilmi\\Desktop\\Lowy\\Non-codingRNA\\Lung\\combined-Colon-TCGA-in-miT-UCSC-all-sample-lncRNA-gene.txt", lncRNAtoGenetoSample);
//    }
//    
    public static void mapExpression(){
        String lncRNAtoGene[][]=DataSetReader.readDataSet("C:\\Users\\Dilmi\\Desktop\\Lowy\\Non-codingRNA\\Colon\\Colon-miT-lncaNet-to-fill.txt", 88, "\t");
        String geneExp[][]=DataSetReader.readDataSet("C:\\Users\\Dilmi\\Desktop\\Lowy\\Non-codingRNA\\Colon\\TMMnorm.txt",506 , "\t");
        for (int i = 1; i < geneExp.length; i++) {
            for (int j = 1; j < lncRNAtoGene.length; j++) {
                if(geneExp[i][0].equals(lncRNAtoGene[j][0])){
                    for (int k = 1; k < geneExp[0].length; k++) {
                        for (int l = 1; l < lncRNAtoGene[0].length; l++) {
                            if(geneExp[0][k].equals(lncRNAtoGene[0][l])){
                                lncRNAtoGene[j][l]=geneExp[i][k];
                                break;
                            }
                        }
                    }
                }
                
            }
        }
        DataSetWriter.writeToFile2D("C:\\Users\\Dilmi\\Desktop\\Lowy\\Non-codingRNA\\Colon\\Colon-miT-lncaNet-to-filled.txt", lncRNAtoGene);
       
    }
    
    public static void mapMutations() {
        String lncRNAtoSample[][]=DataSetReader.readDataSet("C:\\Users\\Dilmi\\Desktop\\Lowy\\Non-codingRNA\\Colon\\combined-Colon-TCGA-in-miT-UCSC-all-sample-lncRNA-to-gene.txt", 3, "\t");
        String mutationMap[][]=DataSetReader.readDataSet("C:\\Users\\Dilmi\\Desktop\\Lowy\\Non-codingRNA\\Colon\\Colon-miT-to-mutations.txt", 88, "\t");
        for (int i = 0; i < lncRNAtoSample.length; i++) {
            for (int j = 1; j < mutationMap.length; j++) {
                if(lncRNAtoSample[i][1].equals(mutationMap[j][0])){
                    for (int k = 1; k < mutationMap[0].length; k++) {
                        if(lncRNAtoSample[i][0].equals(mutationMap[0][k])){
                            mutationMap[j][k]="1";
                        }
                    }
                }
            }
        }
        DataSetWriter.writeToFile2D("C:\\Users\\Dilmi\\Desktop\\Lowy\\Non-codingRNA\\Colon\\Colon-to-lncRNA-mutaitonMap-miT.txt", mutationMap);
    }
    
    public static void compareGroups() {
        String lncRNAtoGene[][]=DataSetReader.readDataSet("C:\\Users\\Dilmi\\Desktop\\Lowy\\Non-codingRNA\\Colon\\Colon-miT-lncaNet-to-filled.txt", 88, "\t");
        String mutationMap[][]=DataSetReader.readDataSet("C:\\Users\\Dilmi\\Desktop\\Lowy\\Non-codingRNA\\Colon\\Colon-to-lncRNA-mutaitonMap-miT.txt", 88, "\t");
        String output[][]=DataSetReader.readDataSet("C:\\Users\\Dilmi\\Desktop\\Lowy\\Non-codingRNA\\Colon\\Colon-lncaNet-miT-map-p-value.txt", 7, "\t");
        
        for (int i = 1; i < output.length; i++) {
            String lncRNA=output[i][0];
            String gene=output[i][1];
            String mutationMatrix[]=new String [mutationMap[0].length-1];
            String geneExpressionMatrix[]=new String [mutationMap[0].length-1];
            
            for (int j = 1; j < lncRNAtoGene.length; j++) {
                if(gene.equals(lncRNAtoGene[j][0]) ){
                    //System.out.print(""+lncRNA);
                    //System.out.println("\t"+gene);
                    for (int k = 1; k < lncRNAtoGene[0].length; k++) {
                        geneExpressionMatrix[k-1]=lncRNAtoGene[j][k];
                    }
                    break;
                }
            }
            for (int j = 1; j < mutationMap.length; j++) {
                if(lncRNA.equals(mutationMap[j][0])  ){
                    
                    for (int k = 1; k < mutationMap[0].length; k++) {
                        mutationMatrix[k-1]=mutationMap[j][k];
                        //System.out.println(""+mutationMatrix[k-1]);
                    }
                    break;
                }
            }
            
            ArrayList<Double> mutants=new ArrayList<Double>();
            ArrayList<Double> wildtype=new ArrayList<Double>();
            ArrayList<Double> normal=new ArrayList<Double>();
            
            for (int j = 0; j < mutationMatrix.length; j++) {
                if(mutationMatrix[j]==null){
                    System.out.println(""+j);
                    System.out.println(""+lncRNA);
                
                    break;
                }
                   // System.out.println(""+j);
                   // System.out.println(""+lncRNA);
                
                else{
                    //System.out.println("true");
                    if(mutationMatrix[j].equals("1")){
                        mutants.add(Double.parseDouble(geneExpressionMatrix[j]));
                    }

                    else if(mutationMatrix[j].equals("0")){
                        wildtype.add(Double.parseDouble(geneExpressionMatrix[j]));
                    }

                    else if(mutationMatrix[j].equals("-1")){
                        normal.add(Double.parseDouble(geneExpressionMatrix[j]));
                    }
                }
            }
            
            if(mutants.size()>0){
                ///output[i][2]=String.valueOf(pValCalc(wildtype, mutants));
                double mutantsA[]=new double [mutants.size()];
                for (int j = 0; j < mutantsA.length; j++) {
                    mutantsA[j]=(double)mutants.get(j);
                }
                double wildtypeA[]=new double [wildtype.size()];
                for (int j = 0; j < wildtypeA.length; j++) {
                    wildtypeA[j]=(double)wildtype.get(j);
                }
                //output[i][2]=String.valueOf(new MannWhitneyUTest().mannWhitneyUTest(wildtypeA, mutantsA));
                output[i][2]=String.valueOf(pValCalc(wildtype, mutants));
                output[i][3]=String.valueOf(mutants.size());
                
                Collections.sort(mutants);
                if (!mutants.isEmpty()) {
                        if ((mutants.size() % 2) == 1) {
                            output[i][4] = String.valueOf(mutants.get(((mutants.size() + 1) / 2) - 1).doubleValue());
                        } else if ((mutants.size() % 2) == 0) {
                            output[i][4] = String.valueOf((mutants.get((mutants.size() / 2) - 1).doubleValue() + mutants.get(((mutants.size() / 2) + 1) - 1).doubleValue()) / 2);
                        }
                }
                
                Collections.sort(wildtype);
                if (!wildtype.isEmpty()) {
                        if ((wildtype.size() % 2) == 1) {
                            output[i][5] = String.valueOf(wildtype.get(((wildtype.size() + 1) / 2) - 1).doubleValue());
                        } else if ((wildtype.size() % 2) == 0) {
                            output[i][5] = String.valueOf((wildtype.get((wildtype.size() / 2) - 1).doubleValue() + wildtype.get(((wildtype.size() / 2) + 1) - 1).doubleValue()) / 2);
                        }
                    }
                
                Collections.sort(normal);
                if (!normal.isEmpty()) {
                        if ((normal.size() % 2) == 1) {
                            output[i][6] = String.valueOf(normal.get(((normal.size() + 1) / 2) - 1).doubleValue());
                        } else if ((normal.size() % 2) == 0) {
                            output[i][6] = String.valueOf((normal.get((normal.size() / 2) - 1).doubleValue() + normal.get(((normal.size() / 2) + 1) - 1).doubleValue()) / 2);
                        }
                    }
            }
                
            
            
        }
        
        DataSetWriter.writeToFile2D("C:\\Users\\Dilmi\\Desktop\\Lowy\\Non-codingRNA\\Colon\\Colon-lncRNADB-miT-map-p-value-foldChange-filled2.txt", output);
    }
    
    public static double pValCalc(ArrayList<Double> sample1,ArrayList<Double> sample2) {
        double sum1=0;
        double sumsq1=0;
        int n1=sample1.size();
        
        double x1;
        for (int i = 0; i < sample1.size(); i++) {
            x1=sample1.get(i);
            sum1=sum1+x1;
            sumsq1=sumsq1+x1*x1;
        }
        
        double mean1=sum1/n1;
        double var1=(sumsq1-(sum1 * sum1)/n1)/(n1-1);
        //System.out.print(""+mean1);
        
        double sum2=0;
        double sumsq2=0;
        int n2=sample2.size();
        
        double x2;
        for (int i = 0; i < sample2.size(); i++) {
            x2=sample2.get(i);
            sum2=sum2+x2;
            sumsq2=sumsq2+x2*x2;
        }
        
        double mean2=sum2/n2;
        double var2=(sumsq2-(sum2 * sum2)/n2)/(n2-1);
        //System.out.println("\t"+mean2);
        
        double pVal=0.0;
        double t;
        int df;
        double var=(var1/n1+var2/n2);
        t = (mean2 - mean1)/
                Math.sqrt(var);
        System.out.println(""+Math.sqrt(var));
        df = n1+n2-2;
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
