package cn.malab.lab.WMSA2.measure;

import cn.malab.lab.WMSA2.psa.STAlign;

/**
 * 计算序列间的最长公共子序列的距离
 */
public class lcs {
    String[] strs;

    public lcs(String[] strs) {
        this.strs = strs;
    }

    /**
     * 计算序列间两两之间的不相似度
     */
    public double[][] getDismatrix() {
        double[][] dismatrix = new double[strs.length][strs.length];
        for (int i = 0; i < strs.length; i++) {
            STAlign sAlign = new STAlign(strs[i]);
            for (int j = i + 1; j < strs.length; j++) {
                dismatrix[i][j] = dismatrix[j][i] = 1 - sAlign.getSimstrB(strs[j]);
            }
        }
        return dismatrix;
    }
}