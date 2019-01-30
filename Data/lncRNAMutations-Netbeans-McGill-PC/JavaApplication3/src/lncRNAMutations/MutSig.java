/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package lncRNAMutations;


/**
 *
 * @author Dilmi
 */
public class MutSig {
    public static void mut() {
        String mutations[][]=DataSetReader.readDataSet("/home/dilmi/data/Combined/Promoter_mutations_paper_mutations_final_for_fasta_edited.fa", 4, "\t");
        String template[][]=DataSetReader.readDataSet("/home/dilmi/data/Combined/MutSig.txt", 1166, "\t");
        
        for (int j = 3; j < template.length ; j++) {
                for (int k = 3; k < template[5].length; k++) {
                    template[j][k]="0";
                }
            }
        
        
        for (int i = 0; i < 1000; i++) {
            String sample=mutations[i][0];
            String baseChange=mutations[i][2];
            String trinucleotide=mutations[i][1];
            String count="1";
            
            
            
            for (int j = 3; j < template[5].length; j++) {
                if(template[2][j].equals(sample)){
                    
                    for (int k = 3; k < template.length; k++) {
                        if(template[k][0].equals(baseChange)&& template[k][1].equals(trinucleotide)){
                            
                            try {
                                //System.out.println("g");
                                String temp=template[k][j];
                                template[k][j]=String.valueOf(Double.parseDouble(count)+Double.parseDouble(temp));
                            }
                            catch(Exception e){
                                //System.out.println(""+e.getMessage());
                            }
                            
                        }
                    }
                    
                }
            }
        }
        DataSetWriter.writeToFile2D("/home/dilmi/data/Combined/MutSig-filled3.txt", template);
    }
    
}
