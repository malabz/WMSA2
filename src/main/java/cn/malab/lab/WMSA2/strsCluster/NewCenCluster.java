package cn.malab.lab.WMSA2.strsCluster;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import cn.malab.lab.WMSA2.io.string;
import cn.malab.lab.WMSA2.measure.kmer;

public class NewCenCluster {
    private final double sim;
    private final boolean silent;
    private final String[] strs;
    private Map<Integer, int[]> clusters;

    public NewCenCluster(String[] strs, double sim, boolean silent) {
        this.sim = sim;
        this.strs = strs;
        this.silent = silent;
        this.genClusters();
    }

    public Map<Integer, int[]> getClusters() {
        return this.clusters;
    }

    private void genClusters() {
        int remainder = strs.length;
        clusters = new HashMap<>();
        boolean[] clued = new boolean[strs.length];
        int[][] kmers = kmer.Counterk(strs, 6);

        while (remainder >= 1) {
            String output = "  clu : " + clusters.size() + "  remain : " + remainder + "      ";
            if (!silent) {
                System.out.print(output);
            }
            int idxc = pickLongestID(clued);
            int[] simIds = pickSimilarIds(kmers, clued, idxc);
            clusters.put(idxc, simIds);
            remainder -= (simIds.length + 1);
            if (!silent) {
                System.out.print(string.repeat("\b", output.length()));
            }
        }
    }

    private int[] pickSimilarIds(int[][] kmers, boolean[] clued, int idxc) {
        clued[idxc] = true;
        List<Integer> res = new ArrayList<>();
        for (int i = 0; i < strs.length; i++) {
            if (!clued[i] && kmer.getDistance(kmers[idxc], kmers[i]) >= sim) {
                res.add(i);
                clued[i] = true;
            }
        }
        int[] resInt = new int[res.size()];
        int i = 0;
        for (Integer r : res) {
            resInt[i++] = r;
        }
        return resInt;
    }

    private int pickLongestID(boolean[] clued) {
        int res = -1, maxlen = -1;
        for (int i = 0; i < strs.length; i++) {
            if (!clued[i] && strs[i].length() > maxlen) {
                res = i;
                maxlen = strs[i].length();
            }
        }
        return res;
    }
}
