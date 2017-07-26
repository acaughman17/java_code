
package hsap.csr;

public class TradMatrix{
    private double[][] mat;
    
    public TradMatrix(double[][] m){mat = m;}
    
    public TradMatrix(CSRMatrix m){
        mat = new double[m.numRows()][m.numCols()];
        for(int i=0; i<m.I.length-1; i++){
            for(int k=m.I[i]; k<m.I[i+1]; k++){
                mat[i][m.J[k]] = m.D[k];
            }
        }
    }
    
    public double at(int row, int col){return mat[row][col];}
    
    public static TradMatrix randomSparseMatrix(int numEntries, int numCols, int numRows, double[] weights){
        if(numEntries > numRows*numCols){numEntries = numRows*numCols;}
        double[][] data = new double[numRows][numCols];
        int column = (int)(Math.random()*numCols);
        int offset = 1;
        data[numRows-1][column] = weights[(int)(Math.random()*weights.length)];
        if(column!=numCols-1){
            int coord = (int)(Math.random()*numRows);
            while(coord!=numRows-1){
                coord = (int)(Math.random()*numRows);
            }
            data[coord][numCols-1] = weights[(int)(Math.random()*weights.length)];
            offset++;
        }
        for(int i=0; i<numEntries-offset; i++){
            int row = (int)(Math.random()*numRows), col = (int)(Math.random()*numCols);
            while(data[row][col]!=0){
                row = (int)(Math.random()*numRows); col = (int)(Math.random()*numCols);
            }
            data[row][col] = weights[(int)(Math.random()*weights.length)];
        }
        return new TradMatrix(data);
    }
    
    public static TradMatrix randomSparseGraph(int numEntries, int n, double[] weights){
        if(numEntries > n*n){numEntries = n*n;}
        double[][] data = new double[n][n];
        for(int i=0; i<numEntries; i++){
            int row = (int)(Math.random()*n), col = (int)(Math.random()*n);
            while(data[row][col]!=0 || row>col || row==col){
                row = (int)(Math.random()*n); col = (int)(Math.random()*n);
            }
            int index = (int)(Math.random()*weights.length);
            data[row][col] = weights[index]; data[col][row] = weights[index];
        }
        return new TradMatrix(data);
    }

    public int numRows(){return mat.length;}

    public int numCols(){return mat[0].length;}
    
    public String toString(){
        String s = "[";
        for(int i=0; i<numRows(); i++){
            for(int j=0; j<numCols(); j++){s+=(int)mat[i][j]+"\t";}
            s=s.substring(0,s.length()-1);
            s+="]\n[";
        }
        return s.substring(0,s.length()-2);
    }
}
