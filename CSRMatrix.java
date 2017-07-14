
package hsap.csr;

public class CSRMatrix{
    public int[] I;
    public int[] J;
    public double[] D;
    private int numCols;
    private int numRows;

    public CSRMatrix(int[] Is, int[] Js, double[] Ds, int numR, int numC){
        I = Is;
        J = Js;
        D = Ds;
        numRows = numR;
        numCols = numC;
    }
    
    public CSRMatrix(TradMatrix m){
        I = new int[m.numRows()+1];
        java.util.List<Integer> js = new java.util.ArrayList<>();
        java.util.List<Double> ds = new java.util.ArrayList<>();
        for(int i=0; i<m.numRows(); i++){
            for(int j=0; j<m.numCols(); j++){
                if(m.at(i,j)!=0){
                    I[i+1]++;
                    js.add(j);
                    ds.add(m.at(i,j));
                }
            }
        }
        for(int i=2; i<I.length; i++){
            I[i] += I[i-1];
        }
        J = new int[js.size()];
        D = new double[ds.size()];
        for(int k=0; k<js.size(); k++){
            J[k] = js.get(k);
            D[k] = ds.get(k);
        }
        
        numRows = m.numRows();
        numCols = m.numCols();
    }
    
    public CSRMatrix copy(){
        int[] nI = new int[I.length], nJ = new int[J.length]; double[] nD = new double[D.length];
        System.arraycopy(I, 0, nI, 0, I.length);
        System.arraycopy(J, 0, nJ, 0, J.length);
        System.arraycopy(D, 0, nD, 0, D.length);
        return new CSRMatrix(nI, nJ, nD, numRows(), numCols());
    }
    
    public CSRMatrix(ColVec v){
        I = new int[v.numRows()+1];
        for(int i=0; i<v.J.length; i++){
            I[v.J[i]+1]++;
        }
        for(int i=2; i<I.length; i++){
            I[i] += I[i-1];
        }
        J = new int[v.J.length];
        D = v.D;
        
        numRows = v.numRows();
        numCols = 1;
    }
    
    public CSRMatrix transpose(){
        //**********************************************
        /*int[] colCount = new int[numCols];
        for(int i=0; i<I.length-1; i++){
            for(int j=I[i]; j<I[i+1]; j++){
                colCount[J[j]]++;
            }
        }
        
        int[] newI = new int[numCols+1];
        newI[0] = 0;
        for(int i=1; i<newI.length; i++){
            newI[i] = colCount[i-1]+newI[i-1];
        }*/
        //**********************************************
        int[] nI = new int[numCols+1];
        for(int k=0; k<J.length; k++){
            nI[J[k]+1]++;
        }
        for(int i=2; i<nI.length; i++){
            nI[i] += nI[i-1];
        }
        //**********************************************
        int[] nJ = new int[nI[nI.length-1]];
        double[] nD = new double[nI[nI.length-1]];
        
        int[] counter = new int[numCols];
        for(int i=0; i<I.length-1; i++){
            for(int k=I[i]; k<I[i+1]; k++){
                nJ[nI[J[k]]+counter[J[k]]]=i;
                nD[nI[J[k]]+counter[J[k]]]=D[k];
                counter[J[k]]++;
            }
        }
        
        return new CSRMatrix(nI, nJ, nD, numCols, numRows);
    }

    public int numRows() {
        return numRows;
    }

    public int numCols() {
        return numCols;
    }
    
    public CSRMatrix neg(CSRMatrix m) {
        int[] nI = new int[I.length], nJ = new int[J.length]; double[] nD = new double[D.length];
        System.arraycopy(I, 0, nI, 0, I.length);
        System.arraycopy(J, 0, nJ, 0, J.length);
        for(int i=0; i<D.length; i++){nD[i] = -D[i];}
        return new CSRMatrix(nI, nJ, nD, numRows, numCols);
    }
    
    public CSRMatrix add(CSRMatrix m) {
        if(m.numRows() != numRows() || m.numCols() != numCols()){
            throw new Error("dimension mismatch");
        }
        
        int[] I0=I, J0=J, I1=m.I, J1=m.J;
        double[] D0=D, D1=m.D;
        
        int[] nI = new int[I0.length];
        for(int i=0; i<I0.length-1; i++){
            int k0 = I0[i], k1 = I1[i];
            while(k0<I0[i+1] && k1<I1[i+1]){
                if(J0[k0]<J1[k1]){k0++;} else if(J0[k0]>J1[k1]){k1++;} else {k0++; k1++;}
                nI[i+1]++;
            }
            nI[i+1] += I0[i+1]-k0 + I1[i+1]-k1;
        }
        
        for(int i=2; i<nI.length; i++){nI[i] += nI[i-1];}
        
        int[] nJ = new int[nI[nI.length-1]]; double[] nD = new double[nI[nI.length-1]];
        
        for(int i=0; i<I0.length-1; i++){
            int count=nI[i], k0 = I0[i], k1 = I1[i];
            while(k0<I0[i+1] && k1<I1[i+1]){
                if(J0[k0]<J1[k1]){nJ[count] = J0[k0]; nD[count] = D0[k0]; k0++;} else 
                if(J0[k0]>J1[k1]){nJ[count] = J1[k1]; nD[count] = D1[k1]; k1++;} else 
                {nJ[count] = J0[k0]; nD[count] = D0[k0]+D1[k1]; k0++; k1++;}
                count++;
            }
            while(k0<I0[i+1]){nJ[count] = J0[k0]; nD[count] = D0[k0]; k0++; count++;}
            while(k1<I1[i+1]){nJ[count] = J1[k1]; nD[count] = D1[k1]; k1++; count++;}
        }
        return new CSRMatrix(nI, nJ, nD, numRows, numCols);
    }
    
    public CSRMatrix sub(CSRMatrix m) {
        if(m.numRows() != numRows() || m.numCols() != numCols()){
            throw new Error("dimension mismatch");
        }
        
        int[] I0=I, J0=J, I1=m.I, J1=m.J;
        double[] D0=D, D1=m.D;
        
        int[] nI = new int[I0.length];
        for(int i=0; i<I0.length-1; i++){
            int k0 = I0[i], k1 = I1[i];
            while(k0<I0[i+1] && k1<I1[i+1]){
                if(J0[k0]<J1[k1]){k0++;} else if(J0[k0]>J1[k1]){k1++;} else {k0++; k1++;}
                nI[i+1]++;
            }
            nI[i+1] += I0[i+1]-k0 + I1[i+1]-k1;
        }
        
        for(int i=2; i<nI.length; i++){nI[i] += nI[i-1];}
        
        int[] nJ = new int[nI[nI.length-1]]; double[] nD = new double[nI[nI.length-1]];
        
        for(int i=0; i<I0.length-1; i++){
            int count=nI[i], k0 = I0[i], k1 = I1[i];
            while(k0<I0[i+1] && k1<I1[i+1]){
                if(J0[k0]<J1[k1]){nJ[count] = J0[k0]; nD[count] = D0[k0]; k0++;} else 
                if(J0[k0]>J1[k1]){nJ[count] = J1[k1]; nD[count] = -D1[k1]; k1++;} else 
                {nJ[count] = J0[k0]; nD[count] = D0[k0]-D1[k1]; k0++; k1++;}
                count++;
            }
            while(k0<I0[i+1]){nJ[count] = J0[k0]; nD[count] = D0[k0]; k0++; count++;}
            while(k1<I1[i+1]){nJ[count] = J1[k1]; nD[count] = -D1[k1]; k1++; count++;}
            if(count>nI[i+1])System.out.println(count+" "+nI[i+1]);
        }
        return new CSRMatrix(nI, nJ, nD, numRows, numCols);
    }
    
    public CSRMatrix mult(CSRMatrix m) {
        if(m.numCols() != numRows()){
            throw new Error("dimension mismatch");
        }
        int[] I0=I, J0=J;
        double[] D0=D;
        int[] I1=m.I, J1=m.J;
        double[] D1=m.D;
        int[] nI=new int[numRows()+1];
        
        for(int i0=0;i0<I0.length-1;i0++){
            boolean[] flag=new boolean[m.numCols()];
            for(int k0=I0[i0];k0<I0[i0+1];k0++){
                for(int k1=I1[J0[k0]];k1<I1[J0[k0]+1];k1++){
                    flag[J1[k1]]=true;
                }
            }
            for(int i=0;i<flag.length;i++){
                if(flag[i]){
                    nI[i0+1]++;
                }
            }
        }
        
        
        for(int i=2;i<nI.length;i++){
            nI[i]+=nI[i-1];
        }
        int[] nJ = new int[nI[nI.length-1]];
        double[] nD = new double[nI[nI.length-1]];
        int count=0;
        for(int i0=0; i0<I0.length-1; i0++){
            int[] temp = new int[m.numCols()];
            for(int k0=I0[i0]; k0<I0[i0+1]; k0++){
                for(int k1=I1[J0[k0]]; k1<I1[J0[k0]+1]; k1++){
                    temp[J1[k1]] += D0[k0]*D1[k1]; //wtf... int += double?
                }
            }
            for(int i=0;i<m.numCols();i++){
                if(temp[i]!=0){
                    nJ[count]=i;
                    nD[count]=temp[i];
                    count++;
                }
            }
        }        
        return new CSRMatrix(nI,nJ,nD,numRows,m.numCols());
    }
    
    public ColVec mult(ColVec v){
        class PairWeight implements Comparable<PairWeight>{int x, y, w; PairWeight(int xx, int yy, int ww){x = xx; y = yy; w = ww;}
            @Override public int compareTo(PairWeight p) {if(p.x>x) return -1; if(p.x<x) return 1; if(p.y>y) return -1; if(p.y<y) return 1; return 0;}}
        
        java.util.List<PairWeight> pairs = new java.util.ArrayList<PairWeight>();
        
        int[] I0 = I, J0 = J; double[] D0 = D;
        int[] J1 = v.J; double[] D1 = v.D;
        
        int maxRow = 0;
        for(int i0=0; i0<I0.length-1; i0++){
            int k0 = I0[i0], k1 = 0, sum = 0;
            while(k0<I0[i0+1] && k1<J1.length)
                if(J0[k0]>J1[k1]) k1++; 
                else if(J0[k0]<J1[k1]) k0++; 
                else sum += D0[k0++]*D1[k1];
            if(sum!=0){
                pairs.add(new PairWeight(i0, 0, sum));
                maxRow = i0>maxRow ? i0 : maxRow;
            }
        }
        
        int[] nJ = new int[pairs.size()]; double[] nD = new double[pairs.size()];
        for(int i=0; i<pairs.size(); i++){
            PairWeight p = pairs.get(i);
            nJ[i] = p.x;
            nD[i] = p.w;
        }
        return new ColVec(nJ, nD);
    }
    
    public CSRMatrix[] decompose(){
        //0:lower tri, 1:dia, 2:upper tri
        int[] nI0 = new int[I.length], nI1 = new int[I.length], nI2 = new int[I.length];
        
        int pos;
        for(int i=0; i<I.length-1; i++){
            pos = 0;
            for(int k=I[i]; k<I[i+1]; k++){
                if(pos==0)if(i<J[k]){pos = 2;} else if(i==J[k]){pos = 1;}
                switch (pos){
                    case 0: nI0[i+1]++;break;
                    case 1: nI1[i+1]++;break;
                    case 2: nI2[i+1]++;
                }
                if(pos==1)pos++;         
            }
        }
        for(int i=2; i<I.length; i++){
            nI0[i] += nI0[i-1];
            nI1[i] += nI1[i-1];
            nI2[i] += nI2[i-1];
        }
        int c0=0, c1=0, c2=0;
        int[] nJ0 = new int[nI0[nI0.length-1]], nJ1 = new int[nI1[nI1.length-1]], nJ2 = new int[nI2[nI2.length-1]];
        double[] nD0 = new double[nI0[nI0.length-1]], nD1 = new double[nI1[nI1.length-1]], nD2 = new double[nI2[nI2.length-1]];
        for(int i=0; i<I.length-1; i++){
            pos = 0;
            for(int k=I[i]; k<I[i+1]; k++){
                if(pos==0)if(i<J[k]){pos = 2;} else if(i==J[k]){pos = 1;}
                switch (pos){
                    case 0: nJ0[c0] = J[k];nD0[c0] = D[k];c0++;break;
                    case 1: nJ1[c1] = J[k];nD1[c1] = D[k];c1++;break;
                    case 2: nJ2[c2] = J[k];nD2[c2] = D[k];c2++;
                }
                if(pos==1)pos++;     
            }
        }
        
        return new CSRMatrix[] {new CSRMatrix(nI0, nJ0, nD0, numRows, numCols),new CSRMatrix(nI1, nJ1, nD1, numRows, numCols),new CSRMatrix(nI2, nJ2, nD2, numRows, numCols)}; 
    }
    
    public static CSRMatrix randomSparseMatrix(int numEntries, int numCols, int numRows, double[] weights){
        /*class PairWeight implements Comparable<PairWeight>{
            int x, y; double w;
            PairWeight(int xx, int yy, double ww){x = xx; y = yy; w = ww;}
            @Override public int compareTo(PairWeight p) {if(p.x>x) return -1; if(p.x<x) return 1; if(p.y>y) return -1; if(p.y<y) return 1; return 0;}
        }
        if(numEntries > numRows*numCols){
            numEntries = numRows*numCols;
        }
        PairWeight[] pairs = new PairWeight[numEntries];
        pairs[0] = new PairWeight(numRows-1,(int)(Math.random()*numCols),weights[(int)(Math.random()*weights.length)]);
        int coord = 0;
        while(coord!=numRows-1){
            coord = (int)(Math.random()*numRows);
        }
        pairs[1] = new PairWeight(coord,numCols-1,Math.random()*weights[(int)(Math.random()*weights.length)]);
        
        for(int i=0; i<numEntries-2; i++){
            boolean contains = false; 
            int row = (int)(Math.random()*numRows);
            int col = (int)(Math.random()*numCols);
            for(int j=0; j<i+2; j++)if(pairs[j].x == row && pairs[j].y==col)contains=true;
            while(!contains){
                row = (int)(Math.random()*numRows);
                col = (int)(Math.random()*numCols);
                for(int j=0; j<i+2; j++)if(pairs[j].x == row && pairs[j].y==col)contains=true;
            }
            pairs[i+2] = new PairWeight(row,col,Math.random()*weights[(int)(Math.random()*weights.length)]);
        }
        java.util.Arrays.sort(pairs);
        int[] I = new int[numRows];
        I[numRows-1]=pairs.length;
        int[] J = new int[pairs.length];
        double[] D = new double[pairs.length];
        
        for(int i=0; i<pairs.length; i++){
            PairWeight p = pairs[i];
            I[p.x]++;
            J[i] = p.y;
            D[i] = p.w;
        }
        
        return new CSRMatrix(I,J,D);*/
        return new CSRMatrix(TradMatrix.randomSparseMatrix(numEntries, numRows, numCols, weights));
    }

    public void scale(double d){
        for(int i=0; i<D.length; i++){
            D[i] *= d;
        }
    }
    
    public static CSRMatrix importFromTxtSorted(java.io.File f){
        class PairWeight implements Comparable<PairWeight>{int x, y; double w; PairWeight(int xx, int yy, double ww){x = xx; y = yy; w = ww;}
            @Override public int compareTo(PairWeight p) {if(p.x>x) return -1; if(p.x<x) return 1; if(p.y>y) return -1; if(p.y<y) return 1; return 0;}}
        
        java.util.List<PairWeight> pairs = new java.util.ArrayList<PairWeight>();
        
        String line;
        String delim = " ";
        try (java.io.BufferedReader br = new java.io.BufferedReader(new java.io.FileReader(f))) {

            while ((line = br.readLine()) != null && line.charAt(0)!='\t') {
                int x = Integer.parseInt(line.substring(0,line.indexOf(delim)));
                line = line.substring(line.indexOf(delim));
                int y = Integer.parseInt(line.substring(0,line.indexOf(delim)));
                line = line.substring(line.indexOf(delim));
                double w = Double.parseDouble(line);
                pairs.add(new PairWeight(x,y,w));
            }
            br.close();
        } catch (java.io.IOException e) {
            e.printStackTrace();
        }
    
        int maxRow = pairs.get(pairs.size()-1).x;
        int[] nI = new int[maxRow+2], nJ = new int[pairs.size()]; double[] nD = new double[pairs.size()];
        for(int i=0; i<pairs.size(); i++){
            PairWeight p = pairs.get(i);
            nI[p.x+1]++;
            nJ[i] = p.y;
            nD[i] = p.w;
        }
        for(int i=2; i<maxRow+1; i++){
            nI[i] += nI[i-1];
        }
        return new CSRMatrix(nI, nJ, nD, 0, 0);
    }

    public String toString(){
        String s = "";
        s+="[";
        for(int i=0; i<I.length; i++){
            s+=I[i]+" ";
        }
        s=s.substring(0,s.length()-1);
        s+="]\n";
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
