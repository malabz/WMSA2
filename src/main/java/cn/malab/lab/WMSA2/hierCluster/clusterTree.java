package cn.malab.lab.WMSA2.hierCluster;

import java.util.List;
import java.util.Map;
import java.util.ArrayList;

import cn.malab.lab.WMSA2.strsCluster.FastCluster;
import cn.malab.lab.WMSA2.measure.strsdist;

public class clusterTree {

    private final String[] strs;
    private node rootNode; 
    private List<int[]> TreeList;

    public node getRootNode() {
        return rootNode;
    }

    public List<int[]> getTreeList() {
        return TreeList;
    }

    public clusterTree(String[] strs) {
        this.strs = strs;
        this.TreeList = new ArrayList<>();
        genTree();
    }

    public clusterTree(String[] strs, String[] labels) {
        this.strs = strs;
        genNode(labels);
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
            TreeList.addAll(htree.getTreeList());
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
        TreeList.addAll(htree.getTreeList());
    }

    private void genNode(String[] labels) {
        // 第一步先完成初步聚类
        Map<Integer, int[]> cluIdx = Cluster();
        int[][] cenNames = new int[cluIdx.size()][];
        node[] nodes = new node[cluIdx.size()];
        int globaln = strs.length;
        int icN = 0;
        // 第二步完成构建子树
        for (int k : cluIdx.keySet()) {
            // 类中序列为1
            if (cluIdx.get(k).length < 1) {
                nodes[icN] = new leafnode(labels[k], icN);
                cenNames[icN++] = new int[] { k, k };
                continue;
            }
            // 类中序列大于1
            // 取出序列和标签
            String[] cstrs = new String[cluIdx.get(k).length + 1];
            String[] lbs = new String[cluIdx.get(k).length + 1];
            cstrs[0] = strs[k];
            lbs[0] = labels[k];
            for (int i = 1; i < cstrs.length; i++) {
                cstrs[i] = strs[cluIdx.get(k)[i - 1]];
                lbs[i] = labels[cluIdx.get(k)[i - 1]];
            }
            // 取出索引值
            int[] idxs = new int[cluIdx.get(k).length + 1];
            idxs[0] = k;
            System.arraycopy(cluIdx.get(k), 0, idxs, 1, idxs.length - 1);
            // 建树
            // 1. compute distance
            strsdist sdist = new strsdist(cstrs, "kmer");
            // 2. genTree
            effupgma htree = new effupgma(sdist.getDismatrix2D(), lbs, idxs, globaln);
            
            // 更新全局索引、树节点、中心序列的索引
            globaln += cluIdx.get(k).length;
            nodes[icN] = htree.getRootNode();
            cenNames[icN++] = new int[] { globaln - 1, k };
        }
        if (cenNames.length < 2) {
            if (nodes.length > 0) rootNode = nodes[0];
            return;
        }
        // 2. build father tree
        String[] fstrs = new String[cenNames.length];
        String[] flbs = new String[cenNames.length];
        int[] idxs = new int[cenNames.length];
        for (int i = 0; i < fstrs.length; i++) {
            fstrs[i] = strs[cenNames[i][1]];
            flbs[i] = nodes[i].toString();
            idxs[i] = i;
        }
        strsdist sdist = new strsdist(fstrs, "kmer");
        effupgma htree = new effupgma(sdist.getDismatrix2D(), flbs, idxs, cenNames.length);
        
        rootNode = htree.getRootNode();
    }

    
    private Map<Integer, int[]> Cluster() {
        FastCluster fastCluster = new FastCluster(strs, 0.85, true);
        // c not in int[]
        return fastCluster.getClusters();
    }
}
