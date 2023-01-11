package cn.malab.lab.WMSA2.psa;

public class PSA {

    private String[] strings;

    public PSA(String A, String B) {
        SeqAlign(A, B);
    }

    public PSA(String[] A, String[] B, char[] alphabet, int kk) {
        PrfAlign(A, B, alphabet, kk);
    }

    public String[] getAlign() {
        return strings;
    }

    private void SeqAlign(String A, String B) {
        double diff = getDiff(A.length(), B.length());
        if (diff > 0.4) {
            DP dp = new DP(A, B);
            strings = dp.getStrAlign();
        } else if (diff > 0.1) {
            Kband kd = new Kband(A, B);
            strings =  kd.getStrAlign();
        } else {
            FMAlign fmAlign = new FMAlign(A, B);
            strings = fmAlign.getStrAlign();
        }
    }

    private void PrfAlign(String[] A, String[] B, char[] alphabet, int kk) {
        double diff = getDiff(A[0].length(), B[0].length());
        if (diff > 0.4) {
            multiDP mdp = new multiDP(A, B, alphabet);
            strings = mdp.getStrsAlign();
        } else if (diff > 0.1) {
            multiKband mkd = new multiKband(A, B, alphabet, kk);
            strings =  mkd.getStrsAlign();
        } else {
            FastMultiAlign fma = new FastMultiAlign(A, B, alphabet, kk, 0, 0);
            strings = fma.getStrsAlign();
        }
    }

    private double getDiff(int l1, int l2) {
        return (double) Math.abs(l1 - l2) / Math.max(l1, l2);
    }

}
