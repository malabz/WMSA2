package cn.malab.lab.WMSA2.strsCluster;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import cn.malab.lab.WMSA2.io.string;
import cn.malab.lab.WMSA2.measure.kmer;

public class CenCluster {
    private final double sim;
    private final boolean silent;
    private String[] strs;
    private int[] lens;
    private Map<Integer, int[]> clusters;

    /**
     * gen clusters
     *
     * @param strs
     * @param sim
     * @param silent
     */
    public CenCluster(String[] strs, double sim, boolean silent) {
        this.sim = sim;
        this.strs = strs;
        this.lens = getLens(strs);
        this.silent = silent;
        genClusters();
    }

    /**
     * gen clusters (aligned)
     */
    public CenCluster(String[] strs, int[] lens, double sim, boolean silent) {
        this.sim = sim;
        this.strs = strs;
        this.lens = lens;
        this.silent = silent;
        genClusters();
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
     * gen clusters
     */
    private void genClusters() {
        clusters = new HashMap<>();
        // a map of change idx --> origin idx
        int[] idxMap = new int[strs.length];
        for (int i = 0; i < strs.length; i++) {
            idxMap[i] = i;
        }

        while (strs.length > 1) {
            String output = "  clu : " + clusters.size() + "  remain : " + strs.length + "      ";
            if (!silent) {
                System.out.print(output);
            }

            int idxc = pickLongest();
            // starDist sDist = new starDist(strs, !iter);
            kmer sDist = new kmer(strs, 6);
            int[] simOnes = pickSimilars(sDist.getDismatrix1D(idxc), idxc);
            simOnes = delOutliers(idxMap, idxc, simOnes);
            clusters.put(idxMap[idxc], simOnes);
            idxMap[idxc] = -1;

            // 剩余数量
            int rnums = idxMap.length - simOnes.length - 1;
            if (rnums == 0) {
                if (!silent) {
                    System.out.print(string.repeat("\b", output.length()));
                }
                break;
            }
            String[] newStrs = new String[rnums];
            int[] newlens = new int[rnums], newIdx = new int[rnums];
            for (int i = 0, j = 0; i < idxMap.length; i++) {
                if (idxMap[i] != -1) {
                    newStrs[j] = strs[i];
                    newlens[j] = lens[i];
                    newIdx[j++] = idxMap[i];
                }
            }
            this.strs = newStrs;
            this.lens = newlens;
            idxMap = newIdx;
            if (!silent) {
                System.out.print(string.repeat("\b", output.length()));
            }
        }
        if (idxMap.length == 1) {
            clusters.put(idxMap[0], new int[0]);
        }
    }

    private int[] delOutliers(int[] idxMap, int idxc, int[] clusters) {
        int[] resInt = new int[clusters.length];
        for (int i = 0; i < clusters.length; i++) {
            resInt[i] = idxMap[clusters[i]];
            idxMap[clusters[i]] = -1;
        }
        return resInt;
        // if (clusters.length <= 2) {
        //     int[] resInt = new int[clusters.length];
        //     for (int i = 0; i < clusters.length; i++) {
        //         resInt[i] = idxMap[clusters[i]];
        //         idxMap[clusters[i]] = -1;
        //     }
        //     return resInt;
        // }
        // FMAlign fmAlign = new FMAlign(strs[idxc].replaceAll("-", ""));
        // double[] scores = new double[clusters.length];
        // double ave = 0.0;
        // for (int i = 0; i < scores.length; i++) {
        //     fmAlign.AlignStrB(strs[clusters[i]].replaceAll("-", ""));
        //     scores[i] = score.sp(fmAlign.getStrAlign()[0], fmAlign.getStrAlign()[1]);
        //     ave += scores[i];
        // }
        // ave /= scores.length;
        // List<Integer> res = new ArrayList<>();
        // for (int i = 0; i < scores.length; i++) {
        //     if (scores[i] >= ave || scores[i] >= 0.9) {
        //         res.add(clusters[i]);
        //     }
        // }
        // int[] resInt = new int[res.size()];
        // int i = 0;
        // for (int value : res) {
        //     resInt[i++] = idxMap[value];
        //     idxMap[value] = -1;
        // }
        // return resInt;
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

    /**
     * pick similar ones
     *
     * @param dismatrix
     * @return idxs
     */
    private int[] pickSimilars(double[] dismatrix, int idxc) {
        List<Integer> res = new ArrayList<>();
        for (int i = 0; i < dismatrix.length; i++) {
            if (dismatrix[i] >= sim) {
                res.add(i);
            }
        }

        int[] resInt = new int[res.size()];
        int i = 0;
        for (Integer r : res) {
            resInt[i++] = r >= idxc ? r + 1 : r;
        }

        return resInt;
    }
}
