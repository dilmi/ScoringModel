/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package promoter;

/**
 *
 * @author Dilmi
 */
public class Promoter {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        //mut();
        String template[][]=DataSetReader.readDataSet("C:\\Users\\Dilmi\\Desktop\\Lowy\\Promoter Mutations paper\\MutSig-filled3.txt", 1166, "\t");
        for (int i = 3; i < template.length; i++) {
            for (int j = 3; j < template[4].length; j++) {
                String temp=template[i][j];
                template[i][j]=String.valueOf(Double.parseDouble(temp)/(Double.parseDouble(template[1][j])));
            }
        }
        DataSetWriter.writeToFile2D("C:\\Users\\Dilmi\\Desktop\\Lowy\\Promoter Mutations paper\\MutSig-filled-normalizeed.txt", template);
        
    }
    public static void mut() {
        //C:\Users\Dilmi\Desktop\Lowy\Promoter Mutations paper\
        String mutations[][]=DataSetReader.readDataSet("/home/dilmi/data/Combined/Promoter_mutations_paper_mutations_final_for_fasta_edited.fa", 4, "\t");
        String template[][]=DataSetReader.readDataSet("/home/dilmi/data/Combined/MutSig.txt", 1166, "\t");
        //String mutations[][]=DataSetReader.readDataSet("C:\\Users\\Dilmi\\Desktop\\Lowy\\Promoter Mutations paper\\temp125.txt", 4, "\t");
        //String template[][]=DataSetReader.readDataSet("C:\\Users\\Dilmi\\Desktop\\Lowy\\Promoter Mutations paper\\MutSig.txt", 1166, "\t");
        for (int j = 3; j < template.length ; j++) {
                for (int k = 3; k < template[5].length; k++) {
                    template[j][k]="0";
                }
            }
        
        int temcount=0;
        
        for (int i = 0; i < mutations.length; i++) {
            boolean found=false;
            String sample=mutations[i][0];
            String baseChange=mutations[i][2];
            String count="1";
            String trinucleotide=mutations[i][1];
            String inverseTri="";
            for (int j = 0; j < 3; j++) {
                inverseTri=inverseTri.concat(getInverse(mutations[i][1].substring(j, j+1)));
            }
            String inverseBaseChange=getInverse(baseChange.substring(0, 1))+">"+getInverse(baseChange.substring(2));
            
            
           
            
            
            
            
            for (int j = 3; j < template[4].length; j++) {
                if(template[2][j].equals(sample)){
                    
                    for (int k = 3; k < template.length; k++) {
                        if((template[k][0].equals(baseChange)&& template[k][1].equals(trinucleotide)) || (template[k][0].equals(inverseBaseChange)&& template[k][1].equals(inverseTri))){
                            
                            try {
                                //System.out.println("g");
                                String temp=template[k][j];
                                template[k][j]=String.valueOf(Double.parseDouble(count)+(Double.parseDouble(temp)));
                                temcount++;
                                found=true;
                                //System.out.println(""+temcount);
                            }
                            catch(Exception e){
                                //System.out.println(""+e.getMessage());
                            }
                            
                        }
                    }
                    
                }
            }
            if(!found){
                
            }
            
        }
        DataSetWriter.writeToFile2D("/home/dilmi/data/Combined/MutSig-filled3.txt", template);
        //DataSetWriter.writeToFile2D("C:\\Users\\Dilmi\\Desktop\\Lowy\\Promoter Mutations paper\\MutSig-filled4.txt", template);
    }
    
    public static String getInverse(String base){
        if(base.equals("A")){
            return "T";
        }
        else if(base.equals("C")){
            return "G";
        }
        else if(base.equals("G")){
            return "C";
        }
        else if(base.equals("T")){
            return "A";
        }
        
        else{
            return "error";
        }
    }
}
