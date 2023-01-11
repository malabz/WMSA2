package cn.malab.lab.WMSA2.msa;

import java.util.Date;
import java.text.SimpleDateFormat;
import java.util.HashMap;

import cn.malab.lab.WMSA2.hierCluster.guidetree;
import cn.malab.lab.WMSA2.io.string;
import cn.malab.lab.WMSA2.measure.*;
import cn.malab.lab.WMSA2.psa.profileAlign;

public class treeAlign {
    private String[] straligned;
    private final String Treemode;
    private int[] orders;
    private int[][] treelist;
    private final int num, kk;
    private final boolean silent, aligned;

    /**
     * Choose one Gen tree mode "nj" or "upgma" or "cluster" and choose silent or
     * output
     */
    public treeAlign(String[] strs, String treemode, boolean silent) {
        // 序列数目
        this.num = strs.length;
        // 是否输出比对过程
        this.silent = silent;
        // 数据是否已经被对齐
        this.aligned = false;
        // 构建指导树的模式
        this.Treemode = treemode;
        // 得到一个合适的k值初始值
        this.kk = score.getK(strs, false);
        // 对齐序列
        Align(strs);
        reOrder();
    }


    public treeAlign(String[] strs, int[][] treelist, boolean silent) {
        // 序列数目
        this.num = strs.length;
        // 是否输出比对过程
        this.silent = silent;
        // 数据是否已经被对齐
        this.aligned = false;
        // 构建指导树的模式
        this.Treemode = "";
        this.treelist = treelist;
        // 得到一个合适的k值初始值
        this.kk = score.getK(strs, false);
        // 对齐序列
        Align(strs);
        reOrder();
    }

    /**
     * Choose one Gen tree mode "nj" or "upgma" or "cluster" and choose silent or
     * output
     */
    public treeAlign(String[] strs, String treemode, boolean silent, int kk) {
        // 序列数目
        this.num = strs.length;
        // 是否输出比对过程
        this.silent = silent;
        // 数据是否已经被对齐
        this.aligned = false;
        // 构建指导树的模式
        this.Treemode = treemode;
        // 得到一个合适的k值初始值
        this.kk = kk;
        // 对齐序列
        Align(strs);
        reOrder();
    }

    /**
     * Choose one Gen tree mode "nj" or "upgma" or "cluster" and choose silent or
     * output
     */
    public treeAlign(String[] strs, String[] strsed, boolean aligned) {
        this.num = strs.length;
        this.silent = false;
        this.aligned = aligned;
        this.Treemode = "upgma";
        this.straligned = strsed;
        // 得到一个合适的k值初始值
        this.kk = score.getK(strs, false);
        Align(strs);
        reOrder();
    }

    /**
     * To get the alignment results.
     */
    public String[] getStrsAlign() {
        return this.straligned;
    }

    public int getK() {
        return kk;
    }

    private void Align(String[] strs) {
        // 1.build tree
        SimpleDateFormat sdf = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss");
        if (!silent && !this.Treemode.equals("")) {
            System.out.println("[" + sdf.format(new Date()) + "] building the " + Treemode + " tree");
        }
        int[][] treeList = genTreeList(strs);
        HashMap<Integer, String[]> strsList = new HashMap<>();
        HashMap<Integer, int[]> labelsList = new HashMap<>();
        if (!silent && !this.Treemode.equals("")) {
            System.out.println("[" + sdf.format(new Date()) + "] Done.");
        }
        // 2.align
        if (!silent) {
            System.out.println("[" + sdf.format(new Date()) + "] aligning...");
        }
        char[] alphabet = kmer.Counter(strs);
        int len = treeList.length, i = 0;
        for (int[] readyAlign : treeList) {
            String outToScreen = "  " + (i + 1) + " / " + len;
            if (!silent) {
                System.out.print(outToScreen);
            }
            String[] strsA = getStrsList(strs, strsList, readyAlign[0]);
            String[] strsB = getStrsList(strs, strsList, readyAlign[1]);
            strsList.put(readyAlign[2], profileAlign.Align(strsA, strsB, alphabet, kk));
            labelsList.put(readyAlign[2], combineLabels(labelsList, readyAlign[0], readyAlign[1]));
            strsList.remove(readyAlign[0]);
            strsList.remove(readyAlign[1]);
            labelsList.remove(readyAlign[0]);
            labelsList.remove(readyAlign[1]);

            i++;
            if (!silent) {
                System.out.print(string.repeat("\b", outToScreen.length()));
            }
        }
        this.orders = labelsList.get(treeList[treeList.length - 1][2]);
        this.straligned = strsList.get(treeList[treeList.length - 1][2]);
        if (!silent) {
            System.out.println("[" + sdf.format(new Date()) + "] Done.");
        }
    }

    private int[] combineLabels(HashMap<Integer, int[]> labelsList, int l1, int l2) {
        int[] listL1 = l1 < num ? new int[] { l1 } : labelsList.get(l1);
        int[] listL2 = l2 < num ? new int[] { l2 } : labelsList.get(l2);
        int[] res = new int[listL1.length + listL2.length];
        System.arraycopy(listL1, 0, res, 0, listL1.length);
        System.arraycopy(listL2, 0, res, listL1.length, listL2.length);
        return res;
    }

    private void reOrder() {
        int i = 0;
        String[] temp = new String[num];
        for (int j : orders) {
            temp[j] = straligned[i++];
        }
        straligned = temp;
    }

    // private int[] getLens(String[] strs) {
    //     int[] res = new int[strs.length];
    //     for (int i = 0; i < strs.length; i++) {
    //         res[i] = strs[i].length();
    //     }
    //     return res;
    // }

    private int[][] genTreeList(String[] strs) {
        if (this.treelist != null) {
            return this.treelist;
        }
        // 未对齐的和对齐的两种生成树的方式
        if (!aligned) {
            guidetree gTree = new guidetree(strs, this.Treemode);
            return gTree.genTreeList(silent);
        } else {
            // clusterTree cTree = new clusterTree(straligned, getLens(strs));
            // return cTree.TreeList.toArray(new int[0][]);
            guidetree gTree = new guidetree(strs, straligned, this.Treemode);
            return gTree.genTreeList(silent);
        }
    }

    private String[] getStrsList(String[] strs, HashMap<Integer, String[]> strsList, int key) {
        return key < this.num ? new String[] { strs[key] } : strsList.remove(key);
    }
}