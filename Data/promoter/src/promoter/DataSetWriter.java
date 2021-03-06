/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package promoter;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.PrintWriter;

/**
 *
 * @author Dilmi
 */
public class DataSetWriter {
    public DataSetWriter(){
        
    }

    public static void writeToFile1D(String link, double[]data){
        try{
        // Create file
            FileWriter fstream = new FileWriter(link);
            PrintWriter out = new PrintWriter(fstream);
            for(int i = 0 ; i < data.length ; i++){
                out.println(Double.toString(data[i]));

            }
        //Close the output stream
        out.close();
        }catch (Exception e){//Catch exception if any
          System.err.println("Error: " + e.getMessage());
        }
    }
    public static void writeToFile1D(String link, String[]data){
            try{
            // Create file
                FileWriter fstream = new FileWriter(link);
                PrintWriter out = new PrintWriter(fstream);
                for(int i = 0 ; i < data.length ; i++){
                    out.println(data[i]);

                }
            //Close the output stream
            out.close();
            }catch (Exception e){//Catch exception if any
              System.err.println("Error: " + e.getMessage());
            }
        }
    public static void writeToFile2D(String link, double[][]data){
        try{
        // Create file
            FileWriter fstream = new FileWriter(link);
            PrintWriter out = new PrintWriter(fstream);
            for(int i = 0 ; i < data.length ; i++){
                for(int j = 0 ; j< data[0].length ; j++){
                    out.println(Double.toString(data[i][j]));
                }
            }
        //Close the output stream
        out.close();
        }catch (Exception e){//Catch exception if any
          System.err.println("Error: " + e.getMessage());
        }
    }
    
     public static void writeToFile2D(String link, String[][]data){
        try{
        // Create file
            FileWriter fstream = new FileWriter(link);
            PrintWriter out = new PrintWriter(fstream);
            for(int i = 0 ; i < data.length ; i++){
                for(int j = 0 ; j< data[i].length ; j++){
                    out.print(data[i][j]+"\t");
                }
                out.println();
            }
        //Close the output stream
        out.close();
        }catch (Exception e){//Catch exception if any
          System.err.println("Error: " + e.getMessage());
        }
    }
    
    public static void writeToFile(String link, double data){
        try{
        // Create file
            FileWriter fstream = new FileWriter(link);
            PrintWriter out = new PrintWriter(fstream);
            out.println(data);
            
            //Close the output stream
        out.close();
        }catch (Exception e){//Catch exception if any
          System.err.println("Error: " + e.getMessage());
        }
    }

    public static void writeObjectToFile(String link, Object object){
        try {
            FileOutputStream fout = new FileOutputStream(link);
            ObjectOutputStream oos = new ObjectOutputStream(fout);
            oos.writeObject(object);
            oos.close();
            fout.close();
        }
        catch(Exception e){
            e.printStackTrace();
        }
    }

  
    
}
