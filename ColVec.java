
package hsap.csr;

public class ColVec {
    public int[] J;
    public double[] D;
    private int numRows;
    
    public ColVec(int[] Js, double[] Ds){
        J = Js;
        D = Ds;
        numRows = 0;
        for(int k=0; k<J.length; k++){
            if(numRows<J[k]+1) numRows = J[k]+1;
        }
    }
    
    public ColVec(int[] Ds){
        int maxIndex=0, count=0;
        for(int i=0; i<Ds.length; i++){
            if(Ds[i]!=0){
                maxIndex = i;
                count++;
            }
        }
        J = new int[count];
        D = new double[count];
        int index=0;
        for(int i=0; i<Ds.length; i++){
            if(Ds[i]!=0){
                J[index]=i;
                D[index++]=Ds[i];
            }
        }
        numRows = maxIndex+1;
    }
    
    public int numRows(){
        return numRows;
    }
    
    public int numCols(){
        return 1;
    }
    
    public String toString(){
        String s = "";
        s+="[";
        for(int i=0; i<J.length; i++){
            s+=J[i]+" ";
        }
        s=s.substring(0,s.length()-1);
        s+="]\n";
        s+="[";
        for(int i=0; i<D.length; i++){
            s+=D[i]+" ";
        }
        s=s.substring(0,s.length()-1);
        return s+"]";
    }
}
