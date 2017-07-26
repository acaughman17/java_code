
package hsap.csr;

public class CSRMatrix{
    public int[] I, J; public double[] D;
    public int numCols, numRows;

    public CSRMatrix(int[] Is, int[] Js, double[] Ds, int numR, int numC){I = Is; J = Js; D = Ds; numRows = numR; numCols = numC;}
    
    public CSRMatrix(TradMatrix m){
        I = new int[m.numRows()+1];
        java.util.List<Integer> js = new java.util.ArrayList<>(2*m.numRows());
        java.util.List<Double> ds = new java.util.ArrayList<>(2*m.numRows());
        for(int i=0; i<m.numRows(); i++){
            for(int j=0; j<m.numCols(); j++){
                if(m.at(i,j)!=0){I[i+1]++; js.add(j); ds.add(m.at(i,j));}
            }
        }
        for(int i=2; i<I.length; i++){I[i] += I[i-1];}
        J = new int[js.size()]; D = new double[ds.size()];
        for(int k=0; k<js.size(); k++){J[k] = js.get(k); D[k] = ds.get(k);}
        numRows = m.numRows(); numCols = m.numCols();
    }
    
    public CSRMatrix copy(){
        int[] nI = new int[I.length], nJ = new int[J.length]; double[] nD = new double[D.length];
        System.arraycopy(I, 0, nI, 0, I.length); System.arraycopy(J, 0, nJ, 0, J.length); System.arraycopy(D, 0, nD, 0, D.length);
        return new CSRMatrix(nI, nJ, nD, numRows(), numCols());
    }
    
    public static CSRMatrix identity(int n){
        int[] nI = new int[n+1], nJ = new int[n]; double[] nD = new double[n];
        for(int i=0; i<n; i++){nI[i+1] = i+1; nJ[i] = i; nD[i] = 1;}
        return new CSRMatrix(nI, nJ, nD, n, n);
    }
    
    public CSRMatrix(ColVec v){
        I = new int[v.numRows()+1]; J = new int[v.J.length]; D = v.D;
        for(int i=0; i<v.J.length; i++){I[v.J[i]+1]++;}
        for(int i=2; i<I.length; i++){I[i] += I[i-1];}
        numRows = v.numRows(); numCols = 1;
    }
    
    public CSRMatrix transpose(){
        int[] nI = new int[numCols+1];
        for(int k=0; k<J.length; k++){nI[J[k]+1]++;}
        for(int i=2; i<nI.length; i++){nI[i] += nI[i-1];}
        int[] nJ = new int[nI[nI.length-1]]; double[] nD = new double[nI[nI.length-1]];
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

    public int numRows(){return numRows;}

    public int numCols(){return numCols;}
    
    public CSRMatrix neg(CSRMatrix m){
        int[] nI = new int[I.length], nJ = new int[J.length]; double[] nD = new double[D.length];
        System.arraycopy(I, 0, nI, 0, I.length); System.arraycopy(J, 0, nJ, 0, J.length); for(int i=0; i<D.length; i++){nD[i] = -D[i];}
        return new CSRMatrix(nI, nJ, nD, numRows, numCols);
    }
    
    public CSRMatrix add(CSRMatrix m){
        if(m.numRows() != numRows() || m.numCols() != numCols()){throw new Error("dimension mismatch");}       
        int[] I0=I, J0=J, I1=m.I, J1=m.J; double[] D0=D, D1=m.D;
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
    
    public CSRMatrix sub(CSRMatrix m){
        if(m.numRows() != numRows() || m.numCols() != numCols()){throw new Error("dimension mismatch");}
        int[] I0=I, J0=J, I1=m.I, J1=m.J; double[] D0=D, D1=m.D;
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
        }
        return new CSRMatrix(nI, nJ, nD, numRows, numCols);
    }
    
    public CSRMatrix mult(CSRMatrix m) {
        if(m.numRows != numCols){throw new Error("dimension mismatch: "+numRows+"x"+numCols+" times "+m.numRows+"x"+m.numCols);}
        int[] I0=I, I1=m.I, J0=J, J1=m.J; double[] D0=D, D1=m.D;
        int[] nI=new int[numRows()+1];
        for(int i0=0;i0<I0.length-1;i0++){
            boolean[] flag=new boolean[m.numCols()];
            for(int k0=I0[i0];k0<I0[i0+1];k0++){
                for(int k1=I1[J0[k0]];k1<I1[J0[k0]+1];k1++){
                    flag[J1[k1]]=true;
                }
            }
            for(int i=0;i<flag.length;i++){
                if(flag[i]){nI[i0+1]++;}
            }
        }
        for(int i=2;i<nI.length;i++){nI[i]+=nI[i-1];} 
        int[] nJ = new int[nI[nI.length-1]]; double[] nD = new double[nI[nI.length-1]];
        int count=0;
        for(int i0=0; i0<I0.length-1; i0++){
            double[] temp = new double[m.numCols()];
            boolean[] flag=new boolean[m.numCols()];
            for(int k0=I0[i0]; k0<I0[i0+1]; k0++){
                for(int k1=I1[J0[k0]]; k1<I1[J0[k0]+1]; k1++){
                    temp[J1[k1]] += D0[k0]*D1[k1];
                    flag[J1[k1]] = true;
                }
            }
            for(int i=0;i<m.numCols();i++){
                if(flag[i]){
                    nJ[count]=i;
                    nD[count]=temp[i];
                    count++;
                }
            }
        }        
        return new CSRMatrix(nI,nJ,nD,numRows,m.numCols);
    }
    
    //needs rewrite
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
    
    //0:lower tri, 1:dia, 2:upper tri
    public CSRMatrix[] decompose(){
        int[] nI0 = new int[I.length], nI1 = new int[I.length], nI2 = new int[I.length];
        for(int i=0; i<I.length-1; i++){
            int pos = 0;
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
            nI0[i] += nI0[i-1]; nI1[i] += nI1[i-1]; nI2[i] += nI2[i-1];
        }
        int c0=0, c1=0, c2=0;
        int[] nJ0 = new int[nI0[nI0.length-1]], nJ1 = new int[nI1[nI1.length-1]], nJ2 = new int[nI2[nI2.length-1]];
        double[] nD0 = new double[nI0[nI0.length-1]], nD1 = new double[nI1[nI1.length-1]], nD2 = new double[nI2[nI2.length-1]];
        for(int i=0; i<I.length-1; i++){
            int pos = 0;
            for(int k=I[i]; k<I[i+1]; k++){
                if(pos==0)if(i<J[k]){pos = 2;} else if(i==J[k]){pos = 1;}
                switch (pos){
                    case 0: nJ0[c0] = J[k];nD0[c0] = D[k];c0++; break;
                    case 1: nJ1[c1] = J[k];nD1[c1] = D[k];c1++; break;
                    case 2: nJ2[c2] = J[k];nD2[c2] = D[k];c2++;
                }
                if(pos==1)pos++;     
            }
        }
        return new CSRMatrix[] {new CSRMatrix(nI0, nJ0, nD0, numRows, numCols),new CSRMatrix(nI1, nJ1, nD1, numRows, numCols),new CSRMatrix(nI2, nJ2, nD2, numRows, numCols)}; 
    }
    
    public int rowSum(int index){return I[index+1]-I[index];}
    
    public double weightedRowSum(int index){double sum=0; for(int i=I[index]; i<I[index+1]; i++){sum+=D[i];} return sum;}
    
    //Luby's algorithm
    public int[][] matchingSet(){
        int[][] M = new int[J.length][2];
        int count = 0;
        float[] W = new float[J.length];
        for(int i=0;i<I.length-1;i++){
            for(int k=I[i];k<I[i+1];k++){
                /*double wi = weightedRowSum(i), wj = weightedRowSum(J[k]);
                W[k]=1/((float)(Math.random()+rowSum(i)*wi/(wi+wj)+rowSum(J[k])*wj/(wi+wj))); //10% time dif*/
                W[k]=1/((float)(Math.random()+rowSum(i)+rowSum(J[k])));
            }
        }
        boolean found = true;
        while(found){
            found=false;
            for(int i=0;i<I.length-1;i++){
                for(int k=I[i]; k<I[i+1]; k++){
                    if(W[k]!=0){
                        found=true;
                        int j=J[k];
                        boolean isMax=true;
                        for(int k0=I[i];k0<I[i+1] && isMax;k0++){
                            if(W[k0]>W[k]) isMax=false;
                        }
                        for(int k0=I[j];k0<I[j+1] && isMax;k0++){
                            if(W[k0]>W[k]) isMax=false;
                        }
                        if(isMax){
                            int min = (i<j) ? i : j;
                            M[count][0]=min; M[count][1]=i+j-min; count++;
                            //erase
                            for(int i0=0;i0<I.length-1;i0++){
                                for(int k0=I[i0];k0<I[i0+1];k0++){
                                      if(i0==i || i0==j || J[k0]==i || J[k0]==j) W[k0]=0;
                                }
                            }
                        }
                    }
                }
            }
        }
        int[][] set = new int[count][2];
        for(int i=0; i<count; i++){set[i] = M[i];}
        return set;
      }
    
    public static CSRMatrix partition(CSRMatrix A, int numParts){
        int[][] m;
        CSRMatrix p = identity(A.numRows);
        do {
            m =  A.matchingSet();
            CSRMatrix np = interpolationMatrix(A,m);
            A = mult(np,A.toLaplacian(),np.transpose()).fromLaplacian();//P * L * P^T = nL
            p = np.mult(p);
        } while(m.length!=0 && p.numRows > numParts);
        return p;
    }
    
    public CSRMatrix toLaplacian(){
        int[] nI = new int[I.length], nJ = new int[J.length+I.length-1]; double[] nD = new double[D.length+I.length-1];
        for(int i=1; i<I.length; i++){nI[i] = I[i]+i;}
        int count=0;
        for(int i=0;i<I.length-1; i++){
            boolean passed=false;
            for(int k=I[i]; k<I[i+1]; k++){
                if(!passed && J[k]>i){
                    double sum=0; for(int k0=I[i]; k0<I[i+1]; k0++){sum+=D[k0];}
                    nJ[k+count] = i; nD[k+count] = sum;
                    count++;
                    passed=true;
                }
                nJ[k+count] = J[k]; nD[k+count] = -D[k];
            }
            if(!passed){
                double sum=0; for(int k0=I[i]; k0<I[i+1]; k0++){sum+=D[k0];}
                nJ[I[i+1]+count] = i; nD[I[i+1]+count] = sum;
                count++;
            }
        }
        return new CSRMatrix(nI, nJ, nD, numRows, numCols);
    }
    
    public CSRMatrix fromLaplacian(){
        int[] nI = new int[I.length], nJ = new int[J.length-(I.length-1)]; double[] nD = new double[D.length-(I.length-1)];
        for(int i=1; i<I.length; i++){nI[i] = I[i]-i;}
        int count=0;
        for(int i=0;i<I.length-1; i++){
            boolean passed=false;
            for(int k=nI[i]; k<nI[i+1]; k++){
                if(!passed && J[k+count]>=i){
                    count++;
                    passed=true;
                }
                nJ[k] = J[k+count]; nD[k] = -D[k+count];
            }
            if(!passed){count++;}
        }
        return new CSRMatrix(nI, nJ, nD, numRows, numCols);
    }
    
    public static CSRMatrix interpolationMatrix(CSRMatrix A, int[][] set){
        int n = A.numRows;
        boolean[] accounted= new boolean[n];
        for(int i=0; i<set.length; i++){
            accounted[set[i][0]]=true;
            accounted[set[i][1]]=true;
        } 
        int[] leftover=new int[n-2*set.length];
        int count=0;
        for(int i=0;i<accounted.length;i++){
            if(!accounted[i]){
                leftover[count]=i;
                count++;
            }
        } 
        int[] nI=new int[n-set.length+1], nJ=new int[n]; double[] nD=new double[n];
        for(int i=0;i<set.length;i++){
            nI[i+1]=2*(i+1);
            nJ[2*i]=set[i][0];
            nJ[2*i+1]=set[i][1];
            nD[2*i]=1;
            nD[2*i+1]=1;
        }
        count=2*set.length;
        for(int i=set.length+1;i<n-set.length+1;i++){
            nI[i]=nI[i-1]+1;
            nJ[count]=leftover[i-(set.length+1)];
            nD[count]=1;
            count++;
        }
        return new CSRMatrix(nI,nJ,nD,nI.length-1,A.I.length-1);
    }
    
    public static CSRMatrix mult(CSRMatrix... Ms){
        CSRMatrix ret = Ms[0];
        for(int i=1; i<Ms.length; i++){ret = ret.mult(Ms[i]);}
        return ret;
    }
    
    public static CSRMatrix randomSparseMatrix(int numEntries, int numCols, int numRows, double[] weights){return new CSRMatrix(TradMatrix.randomSparseMatrix(numEntries, numRows, numCols, weights));}

    public void scale(double d){for(int i=0; i<D.length; i++){D[i] *= d;}}
    
    //rewrite maybe
    public static CSRMatrix importFromTxtSorted(java.io.File f){
        class PairWeight implements Comparable<PairWeight>{int x, y; double w; PairWeight(int xx, int yy, double ww){x = xx; y = yy; w = ww;}
            @Override public int compareTo(PairWeight p) {if(p.x>x) return -1; if(p.x<x) return 1; if(p.y>y) return -1; if(p.y<y) return 1; return 0;}}
        java.util.List<PairWeight> pairs = new java.util.ArrayList<PairWeight>();
        String line, delim = " ";
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
        } catch (java.io.IOException e) {e.printStackTrace();}
        int maxRow = pairs.get(pairs.size()-1).x;
        int[] nI = new int[maxRow+2], nJ = new int[pairs.size()]; double[] nD = new double[pairs.size()];
        for(int i=0; i<pairs.size(); i++){
            PairWeight p = pairs.get(i);
            nI[p.x+1]++; nJ[i] = p.y; nD[i] = p.w;
        }
        for(int i=2; i<maxRow+1; i++){nI[i] += nI[i-1];}
        return new CSRMatrix(nI, nJ, nD, 0, 0);
    }

    @Override
    public String toString(){
        String s = "[";
        for(int i=0; i<I.length; i++){s+=I[i]+"\t";}
        s=s.substring(0,s.length()-1)+"]\n[";
        for(int i=0; i<J.length; i++){s+=J[i]+"\t";}
        s=s.substring(0,s.length()-1)+"]\n[";
        for(int i=0; i<D.length; i++){s+=D[i]+"\t";}
        return s.substring(0,s.length()-1)+"]";
    }
    
    public static CSRMatrix connectedComponents(CSRMatrix A){return partition(A,1);}
}
