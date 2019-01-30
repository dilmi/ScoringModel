/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package lncRNAMutations;

import java.io.File;

/**
 *
 * @author Dilmi
 */
public class retrieveFiles {
    public static void getGeneExp() {
        String fileList[][]=DataSetReader.readDataSet("/home/dilmi/data/Combined/temp/Expression/BRCA/ListOfFiles.txt", 2, "\t");
        String geneExp[][]=DataSetReader.readDataSet("/home/dilmi/data/Combined/temp/Expression/BRCA/template-geneExp.txt", 1219, "\t");
        File folder = new File("/home/dilmi/data/Combined/temp/Expression/BRCA/RNASeqV2/UNC__IlluminaHiSeq_RNASeqV2/Level_3");
        
        File[] listOfFiles = folder.listFiles();

        for (int i = 1; i < geneExp[0].length; i++) {
            if(geneExp[0][i].equals(fileList[i][1])){
                //System.out.println("True");
                String file=fileList[i][0];
                for (int j = 0; j < listOfFiles.length; j++) {
                    if(listOfFiles[j].getName().contains(file) && listOfFiles[j].getName().contains(".rsem.genes.normalized_results")){
                        //System.out.println("True1");
                        String geneExpFile[][]=DataSetReader.readDataSet(listOfFiles[j].getAbsolutePath(), 2, "\t"); 
                        for (int k = 1; k < geneExp.length; k++) {
                            if(geneExp[k][0].equals(geneExpFile[k][0])){
                                //System.out.println("True2");
                                geneExp[k][i]=geneExpFile[k][1];
                            }
                            
                        }
                        break;

                    }
                }
            }
        }
        DataSetWriter.writeToFile2D("/home/dilmi/data/Combined/temp/Expression/BRCA/GeneExpression.txt", geneExp);
        
    }
    

        
    
}
