public static CSRMatrix add(CSRMatrix A, CSRMatrix B){
        int[] ai=A.getI();
        int[] aj=A.getJ();
        int[] ad=A.getD();
        int[] bi=B.getI();
        int[] bj=B.getJ();
        int[] bd=B.getD();
        int count=-1;
        for(int i=0;i<ai.length-1;i++){
            int ak=ai[i];
            int bk=bi[i];
            while(ak<ai[i+1]||bk<bi[i+1]){
                count++;
                if(aj[ak]>bj[bk]){
                    bk++;
                }
                else if(aj[ak]<bj[bk]){
                    ak++;
                }
                else{
                    ak++;
                    bk++;
                }
            }
            while(ak<ai[i+1]){
                count++;
                ak++;
            }
            while(bk<bi[i+1]){
                count++;
                bk++;
            }
        }
        int[] ci=new int[ai.length];
        int[] cj=new int[count];
        int[] cd=new int[count];
        count=-1;
        for(int i=0;i<ai.length-1;i++){
            int ak=ai[i];
            int bk=bi[i];
            while(ak<ai[i+1]&&bk<bi[i+1]){
                count++;
                if(aj[ak]>bj[bk]){
                    ci[i+1]++;
                    cj[count]=bj[bk];
                    cd[count]=bd[bk];
                    bk++;
                }
                else if(aj[ak]<bj[bk]){
                    ci[i+1]++;
                    cj[count]=aj[ak];
                    cd[count]=ad[ak];
                    ak++;
                }
                else{
                    ci[i+1]++;
                    cj[count]=aj[ak];
                    cd[count]=ad[ak]+bd[bk];
                    ak++;
                    bk++;
                }
            }
            while(ak<ai[i+1]){                
                count++;
                ci[i+1]++;
                cj[count]=aj[ak];
                cd[count]=ad[ak];
                ak++;
            }
            while(bk<bi[i+1]){
                count++;
                ci[i+1]++;
                cj[count]=bj[bk];
                cd[count]=bd[bk];
                bk++;
            }
        }
        
        for(int i=0;i<ci.length-1;i++){
            ci[i+1]+=ci[i];
        }
        CSRMatrix C = new CSRMatrix(ci,cj,cd);
        return C;
    }
