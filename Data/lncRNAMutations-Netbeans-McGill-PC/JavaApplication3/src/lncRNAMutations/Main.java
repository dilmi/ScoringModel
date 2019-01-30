/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package lncRNAMutations;



/**
 *
 * @author dilmiperera
 */
public class Main {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        //GeneMapBootstrap.getCorrelation();
        //System.wait();
        //GeneMapBootstrap.bootstrap(1000);
//        getFinalList();
//        MutSig.mut();
        //retrieveFiles.getGeneExp();
        /*String mutations[][]=DataSetReader.readDataSet("/Users/dilmiperera/Desktop/Lowy-Dilmi-PC/sampleTolncRNA.csv", 3, ",");
        String mutationToSample[][]=DataSetReader.readDataSet("/Users/dilmiperera/Desktop/Lowy-Dilmi-PC/mutationToSample.csv", 48, ",");
        
        
        for (int i = 1; i < mutationToSample.length; i++) {
            for (int j = 1; j < mutationToSample[0].length; j++) {
                mutationToSample[i][j]="0";
            }
        }
        
        for (int i = 1; i < mutations.length; i++) {
//            System.out.println(""+i);
            for (int j = 1; j < mutationToSample.length; j++) {
                
                if(mutations[i][1].equals(mutationToSample[j][0])){
                    for (int k = 1; k < mutationToSample[1].length; k++) {
                        if(mutationToSample[0][k].equals(mutations[i][0])){
                            mutationToSample[j][k]="1";
                            
                            break;
                        }
                    }
                    
                    break;
                }
            }
        }
        
        DataSetWriter.writeToFile2D("/Users/dilmiperera/Desktop/Lowy-Dilmi-PC/mutationToSample-edited.tsv", mutationToSample);*/
    }
    
    
    public static void getFinalList(){
        String finalList[][]=DataSetReader.readDataSet("/Users/dilmiperera/Desktop/Lowy-Dilmi-PC/Jason-sampleList-summary.txt", 2, "\t");
        String currentList[][]=DataSetReader.readDataSet("/Users/dilmiperera/Desktop/Lowy-Dilmi-PC/currentList.txt", 3, "\t");
        for (int i = 0; i < currentList.length; i++) {
            boolean tumourFound=false;
            String samplesNotFound="-";
            String currentSamples[]=currentList[i][1].split(",");
            System.out.println(""+i);
            for (int j = 0; j < finalList.length; j++) {
                if(finalList[j][0].equals(currentList[i][0])){
                    tumourFound= true;
                    System.out.println("True");
                    
                    String finalSamples[]=finalList[j][1].split(",");
                    for (int k = 0; k < currentSamples.length; k++) {
                        boolean sampleFound=false;
                        for (int l = 0; l < finalSamples.length; l++) {
                            if(currentSamples[k].equals(finalSamples[l])){
                                sampleFound=true;
                                break;
                            }
                        }
                        if(!sampleFound){
                            if(!currentList[j][0].equals("Lung") && !currentSamples[k].contains("TCGA")){ 
                                if(samplesNotFound.equals("-")){
                                    samplesNotFound=currentSamples[k];
                                }
                                else{
                                    String temp=samplesNotFound;
                                    samplesNotFound=temp+","+currentSamples[k];
                                }
                            }
                        }
                    }
                    currentList[i][2]=samplesNotFound;
                    break;
                }
            }
            if(!tumourFound){
                currentList[i][2]=samplesNotFound;
            }
        }
        DataSetWriter.writeToFile2D("/Users/dilmiperera/Desktop/Lowy-Dilmi-PC/finalList.txt", currentList);
    }
    
}
