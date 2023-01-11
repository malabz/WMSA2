package cn.malab.lab.WMSA2.psa;

public class profileAlign {
    public static String[] Align(String[] A, String[] B, char[] alphabet, int kk) {
        PSA psa;
        if (A.length == 1 && B.length == 1) {
            psa = new PSA(A[0], B[0]);
        } else {
            psa = new PSA(A, B, alphabet, kk);
        }
        return psa.getAlign();
    }
}
