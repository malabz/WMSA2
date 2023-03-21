package cn.malab.lab.WMSA2.hierCluster;

import cn.malab.lab.WMSA2.measure.strsdist;

public class guidetree {
    private final String[] strs,strsed;
    private final String mode;

    /**
     * mode: "nj" or "upgma" or "cluster"
     */
    public guidetree(String[] strs, String mode) {
        this.strs = strs;
        this.strsed = null;
        this.mode = mode;
    }

    /**
     * mode: "nj" or "upgma" or "cluster"
     */
    public guidetree(String[] strs, String[] strsed, String mode) {
        this.strs = strs;
        this.strsed = strsed;
        this.mode = mode;
    }

    /**
     * gen an order list of the treealign
     */
    public int[][] genTreeList(boolean silent) {
        // compute similarity matrix
        if (mode.equals("cluster")) {
            clusterTree cTree = new clusterTree(strs);
            return cTree.getTreeList().toArray(new int[0][]);
        }
        strsdist sdist;
        if (strsed != null) {
            sdist = new strsdist(strsed, "aligned");
        } else {
            sdist = new strsdist(strs, "kmer");
        }
        double[][] simMatrix = sdist.getDismatrix2D();
        if (mode.equals("upgma")) {
            effupgma htree = new effupgma(simMatrix);
            return htree.getTreeList().toArray(new int[0][]);
        } else if (mode.equals("nj")) {
            NeighborJoining nJtree = new NeighborJoining(simMatrix, silent);
            return nJtree.TreeList.toArray(new int[0][]);
        } else {
            throw new IllegalArgumentException("unknown mode: " + mode);
        }
    }
}
