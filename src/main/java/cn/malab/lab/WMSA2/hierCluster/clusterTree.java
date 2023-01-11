package cn.malab.lab.WMSA2.hierCluster;

import java.util.List;
import java.util.Map;
import java.util.ArrayList;

import cn.malab.lab.WMSA2.strsCluster.CenCluster;
import cn.malab.lab.WMSA2.strsCluster.CenCluster;
import cn.malab.lab.WMSA2.strsCluster.FastCluster;
import cn.malab.lab.WMSA2.measure.strsdist;

public class clusterTree {

    private final String[] strs;
    private boolean aligned = false;
    private int[] lens;
    public List<int[]> TreeList;

    public clusterTree(String[] strs) {
        this.strs = strs;
        this.TreeList = new ArrayList<>();
        genTree();
    }

    public clusterTree(String[] strs, String[] strsed, int[] lens) {
        this.strs = strs;
        this.lens = lens;
        this.aligned = true;
        this.TreeList = new ArrayList<>();
        // genTreeAligned(strsed);
    }

    private void genTree() {
        // 第一步先完成初步聚类
        Map<Integer, int[]> cluIdx = Cluster();
        int[][] cenNames = new int[cluIdx.size()][];
        int globaln = strs.length;
        int icN = 0;
        // 第二步完成构建子树
        for (int k : cluIdx.keySet()) {
            if (cluIdx.get(k).length < 1) {
                cenNames[icN++] = new int[] { k, k };
                continue;
            }
            String[] cstrs = new String[cluIdx.get(k).length + 1];
            cstrs[0] = strs[k];
            for (int i = 1; i < cstrs.length; i++)
                cstrs[i] = strs[cluIdx.get(k)[i - 1]];
            int[] idxs = new int[cluIdx.get(k).length + 1];
            idxs[0] = k;
            System.arraycopy(cluIdx.get(k), 0, idxs, 1, idxs.length - 1);
            // 1. compute distance
            strsdist sdist = new strsdist(cstrs, "kmer");
            // 2. genTree
            effupgma htree = new effupgma(sdist.getDismatrix2D(), idxs, globaln);
            TreeList.addAll(htree.TreeList);
            globaln += cluIdx.get(k).length;
            cenNames[icN++] = new int[] { globaln - 1, k };
        }
        if (cenNames.length < 2) {
            return;
        }
        // 2. build father tree
        String[] fstrs = new String[cenNames.length];
        int[] idxs = new int[cenNames.length];
        for (int i = 0; i < fstrs.length; i++) {
            fstrs[i] = strs[cenNames[i][1]];
            idxs[i] = cenNames[i][0];
        }
        strsdist sdist = new strsdist(fstrs, "kmer");
        effupgma htree = new effupgma(sdist.getDismatrix2D(), idxs, globaln);
        TreeList.addAll(htree.TreeList);
    }

    // TODO:
    // private void genTreeAligned(String[] strsed) {
    //     // 第一步先完成初步聚类
    //     Map<Integer, int[]> cluIdx = Cluster();
    //     int[][] cenNames = new int[cluIdx.size()][];
    //     int globaln = strs.length;
    //     int icN = 0;
    //     // 第二步完成构建子树
    //     for (int k : cluIdx.keySet()) {
    //         if (cluIdx.get(k).length < 1) {
    //             cenNames[icN++] = new int[] { k, k };
    //             continue;
    //         }
    //         String[] cstrs = new String[cluIdx.get(k).length + 1];
    //         cstrs[0] = strs[k];
    //         for (int i = 1; i < cstrs.length; i++)
    //             cstrs[i] = strs[cluIdx.get(k)[i - 1]];
    //         int[] idxs = new int[cluIdx.get(k).length + 1];
    //         idxs[0] = k;
    //         System.arraycopy(cluIdx.get(k), 0, idxs, 1, idxs.length - 1);
    //         // 1. compute distance
    //         strsdist sdist;
    //         if (aligned) {
    //             sdist = new strsdist(cstrs, "aligned");
    //         } else {
    //             sdist = new strsdist(cstrs, "kmer");
    //         }
    //         // 2. genTree
    //         effupgma htree = new effupgma(sdist.getDismatrix2D(), idxs, globaln);
    //         TreeList.addAll(htree.TreeList);
    //         globaln += cluIdx.get(k).length;
    //         cenNames[icN++] = new int[] { globaln - 1, k };
    //     }
    //     if (cenNames.length < 2) {
    //         return;
    //     }
    //     // 2. build father tree
    //     String[] fstrs = new String[cenNames.length];
    //     int[] idxs = new int[cenNames.length];
    //     for (int i = 0; i < fstrs.length; i++) {
    //         fstrs[i] = strs[cenNames[i][1]];
    //         idxs[i] = cenNames[i][0];
    //     }
    //     strsdist sdist;
    //     if (aligned) {
    //         sdist = new strsdist(fstrs, "aligned");
    //     } else {
    //         sdist = new strsdist(fstrs, "kmer");
    //     }
    //     effupgma htree = new effupgma(sdist.getDismatrix2D(), idxs, globaln);
    //     TreeList.addAll(htree.TreeList);
    // }
    
    private Map<Integer, int[]> Cluster() {
        if (!aligned) {
            FastCluster fastCluster = new FastCluster(strs, 0.85, true);
            // c not in int[]
            return fastCluster.getClusters();
        } else {
            // 比对过的来建树
            CenCluster cenCluster = new CenCluster(strs, lens, 0.9, true);
            return cenCluster.getClusters();
        }
    }
}
