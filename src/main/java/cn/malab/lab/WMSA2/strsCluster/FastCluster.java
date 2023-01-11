package cn.malab.lab.WMSA2.strsCluster;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class FastCluster {
    private final double sim;
    private final String[] strs;
    private final int[] lens;
    private final boolean silent;
    private Map<Integer, int[]> clusters;

    public FastCluster(String[] strs, double sim, boolean silent) {
        this.sim = sim;
        this.strs = strs;
        this.silent = silent;
        this.lens = getLens(strs);
        this.genClusters();
    }

    /**
     * return the idxs of each cluster
     * 
     * @return clusters
     */
    public Map<Integer, int[]> getClusters() {
        return clusters;
    }

    private int[] getLens(String[] strs) {
        int[] res = new int[strs.length];
        for (int i = 0; i < strs.length; i++) {
            res[i] = strs[i].length();
        }
        return res;
    }

    /**
     * find the longest one
     */
    private int pickLongest() {
        int res = 0;
        for (int i = 1; i < lens.length; i++) {
            res = lens[i] > lens[res] ? i : res;
        }
        return res;
    }

    private int[] LensToFind(int idxc) {
        int length = (int) (lens[idxc] * sim);
        lens[idxc] = -1;
        List<Integer> res = new ArrayList<>();
        for (int i = 0; i < lens.length; i++) {
            if (lens[i] >= length) {
                res.add(i);
                lens[i] = -1;
            }
        }
        int[] resInt = new int[res.size() + 1];
        int i = 1;
        resInt[0] = idxc;
        for (int r : res)
            resInt[i++] = r;
        return resInt;
    }

    private void genClusters() {
        int remainder = this.strs.length;
        clusters = new HashMap<>();
        // 先构建长度类
        Map<Integer, int[]> LensClusters = new HashMap<>();
        while (remainder > 0) {
            int idxc = pickLongest();
            LensClusters.put(idxc, LensToFind(idxc));
            remainder -= (LensClusters.get(idxc).length);
        }
        // 对每个长度类再进行聚类
        for (int key : LensClusters.keySet()) {
            int[] idxs = LensClusters.get(key);
            String[] strings = new String[idxs.length];
            for (int i = 0; i < strings.length; i++) {
                strings[i] = strs[idxs[i]];
            }
            // CenCluster ccluster = new CenCluster(strings, sim, silent);
            NewCenCluster ccluster = new NewCenCluster(strings, sim, silent);
            Map<Integer, int[]> temp = ccluster.getClusters();
            for (int k : temp.keySet()) {
                int[] clus = temp.get(k);
                for (int i = 0; i < clus.length; i++) {
                    clus[i] = idxs[clus[i]];
                }
                clusters.put(idxs[k], clus);
            }
        }
        if (!silent) {
            System.out.println("clusters: " + clusters.size() + "                    ");
        }
    }
}
