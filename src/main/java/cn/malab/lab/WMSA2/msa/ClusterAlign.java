package cn.malab.lab.WMSA2.msa;

import java.io.IOException;
import java.text.SimpleDateFormat;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.Date;
import java.util.HashMap;
import java.util.Map;

import cn.malab.lab.WMSA2.hierCluster.guidetree;
import cn.malab.lab.WMSA2.io.string;
import cn.malab.lab.WMSA2.psa.PSA;
import cn.malab.lab.WMSA2.psa.profileAlign;
import cn.malab.lab.WMSA2.strsCluster.*;
import cn.malab.lab.WMSA2.measure.*;

public class ClusterAlign {
    private final int kk, numsStrs;
    private int[] order;
    private final char[] alphabet;
    private String[] labels, strsAligned;
    private Map<Integer, int[]> cluIdx;
    // mode1 "c" or "t" mode2 "c" or "t"
    private String mode1 = "t", mode2 = "t1";
    private double sim;
    public String[] getStrsAlign() {
        reOrder();
        return strsAligned;
    }

    public ClusterAlign(String[] strs, double sim, String modein, String modeout) {
        mode1 = modein;
        mode2 = modeout;
        this.sim = sim;
        numsStrs = strs.length;
        kk = score.getK(strs, false);
        alphabet = kmer.Counter(strs);
        strsAligned = new String[strs.length];
        multiAlign(cluAlign(getCluster(strs)));
        // new reAlign(strsAligned);
    }

    /**
     * get cluster from cdhitfilew
     */
    public ClusterAlign(String[] strs, String[] labels, String cdhitfile) throws IOException {
        numsStrs = strs.length;
        kk = score.getK(strs, false);
        alphabet = kmer.Counter(strs);
        strsAligned = new String[strs.length];
        this.labels = labels;
        multiAlign(cluAlign(getCluCdhit(strs, cdhitfile)));
    }

    /**
     * reorder the strings
     */
    private void reOrder() {
        int j = 0;
        String[] res = new String[strsAligned.length];
        for (int i : order) {
            res[i] = strsAligned[j++];
        }
        strsAligned = res;
    }

    /**
     * gen the cluster by own
     */
    private Map<Integer, String[]> getCluster(String[] strs) {
        SimpleDateFormat sdf = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss");
        System.out.print("[" + sdf.format(new Date()) + "] ");
        System.out.print("clustering ... ");
        FastCluster fastCluster = new FastCluster(strs, this.sim, false);
        cluIdx = fastCluster.getClusters();
        HashMap<Integer, String[]> multiStrs = new HashMap<>();
        for (int c : cluIdx.keySet()) {
            int[] newOne = new int[cluIdx.get(c).length + 1];
            newOne[0] = c;
            System.arraycopy(cluIdx.get(c), 0, newOne, 1, newOne.length - 1);
            cluIdx.put(c, newOne);

            String[] newStrs = new String[cluIdx.get(c).length];
            for (int i = 0; i < newStrs.length; i++) {
                newStrs[i] = strs[cluIdx.get(c)[i]];
            }
            multiStrs.put(c, newStrs);
        }
        System.out.println("[" + sdf.format(new Date()) + "] Done.");
        return multiStrs;
    }

    /**
     * get cluster strings [][] by cdhitdfile
     */
    private Map<Integer, String[]> getCluCdhit(String[] strs, String cdhitfile) throws IOException {
        exCDhit cdhit = new exCDhit();
        String[][] names = cdhit.readClstr(cdhitfile);
        Map<String, Integer> idxmap = new HashMap<>();
        Map<Integer, String[]> multiStrs = new HashMap<>();
        cluIdx = new HashMap<>();
        Pattern p = Pattern.compile("(>\\s*\\w+).*");
        for (int i = 0; i < labels.length; i++) {
            Matcher m = p.matcher(labels[i]);
            if (m.matches())
                idxmap.put(m.group(1), i);
            else
                throw new NumberFormatException();
        }
        for (String[] name : names) {
            String[] temp = new String[name.length];
            int[] tempidx = new int[name.length];
            for (int j = 0; j < name.length; j++) {
                tempidx[j] = idxmap.get(name[j]);
                temp[j] = strs[tempidx[j]];
            }
            multiStrs.put(tempidx[0], temp);
            cluIdx.put(tempidx[0], tempidx);
        }
        return multiStrs;
    }

    private Map<Integer, String[]> cluAlign(Map<Integer, String[]> multiStrs) {
        mode1 = "c";
        // if (multiStrs.size() <= 1) {
        //     int length = 0;
        //     for (int c : multiStrs.keySet())
        //         length = multiStrs.get(c).length;
        //     mode1 = length <= 3000 ? "t" : "c";
        // }
        int i = 1;
        SimpleDateFormat sdf = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss");
        System.out.print("[" + sdf.format(new Date()) + "] ");
        System.out.print("aligning each cluster\n");
        for (int c : multiStrs.keySet()) {
            String outToScreen = "  " + (i++) + "/" + multiStrs.size();
            System.out.print(outToScreen);
            if (multiStrs.get(c).length == 2) {
                PSA psa = new PSA(multiStrs.get(c)[0], multiStrs.get(c)[1]);
                multiStrs.put(c, psa.getAlign());
            } else if (multiStrs.get(c).length < 2) {
                System.out.print(string.repeat("\b", outToScreen.length()));
                continue;
            } else if (mode1.equals("c") || multiStrs.get(c).length >= 3000) {
                centerAlign cAlign = new centerAlign(multiStrs.get(c), true);
                multiStrs.put(c, cAlign.getStrsAlign());
            } 
            // else if (mode1.equals("t")) {
            //     treeAlign tAlign = new treeAlign(multiStrs.get(c), "upgma", true, kk);
            //     multiStrs.put(c, tAlign.getStrsAlign());
            // } 
            else {
                throw new IllegalArgumentException("Unknown mode " + mode1);
            }
            System.out.print(string.repeat("\b", outToScreen.length()));
        }
        System.out.println("[" + sdf.format(new Date()) + "] Done.");
        return multiStrs;
    }

    private void multiAlign(Map<Integer, String[]> multiStrsed) {
        if (numsStrs > 10000) {
            mode2 = "t1";
        }
        switch (mode2) {
            case "t1":
                TreeAlign1(multiStrsed);
                break;
            case "t2":
                TreeAlign2(multiStrsed);
                break;
            default:
                throw new IllegalArgumentException("Unknown mode " + mode2);
        }
    }

    /**
     * 类中取1，再来比对
     */
    public void TreeAlign1(Map<Integer, String[]> multiStrs) {
        SimpleDateFormat sdf = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss");
        System.out.print("[" + sdf.format(new Date()) + "] ");
        System.out.println("combining all cluster");
        if (multiStrs.size() < 2) {
            for (int i : multiStrs.keySet()) {
                strsAligned = multiStrs.get(i);
                order = cluIdx.get(i);
            }
        } else {
            // cluIdx 索引映射 cluIdx [centerIdxc, ...] N + 1
            // multiStrs 类的字符串映射
            String[] centerStrs = new String[multiStrs.size()];
            int[] keys = new int[multiStrs.size()];
            int i = 0;
            for (int key : multiStrs.keySet()) {
                centerStrs[i] = PickRealOne(multiStrs.get(key));
                keys[i++] = key;
            }
            treeAlign tAlign = new treeAlign(centerStrs, "cluster", true, kk);
            centerStrs = tAlign.getStrsAlign();
            for (i = 0; i < centerStrs.length; i++) {
                insertGapStrings(centerStrs[i], multiStrs.get(keys[i]));
            }
            order = new int[numsStrs];
            i = 0;
            for (int key : keys) {
                int len = cluIdx.get(key).length;
                System.arraycopy(cluIdx.remove(key), 0, order, i, len);
                System.arraycopy(multiStrs.remove(key), 0, strsAligned, i, len);
                i += len;
            }
        }
        System.out.println("[" + sdf.format(new Date()) + "] Done.");
    }

    /**
     * 生成一颗中心序列的树，再来树比对
     */
    private void TreeAlign2(Map<Integer, String[]> multiStrs) {
        SimpleDateFormat sdf = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss");
        System.out.print("[" + sdf.format(new Date()) + "] ");
        System.out.println("combining all clusters");
        if (multiStrs.size() < 2) {
            for (int i : multiStrs.keySet()) {
                strsAligned = multiStrs.get(i);
                order = cluIdx.get(i);
            }
        } else {
            Map<Integer, int[]> newcluIdx = new HashMap<>();
            Map<Integer, Integer> idxmap = new HashMap<>();
            Map<Integer, String[]> NewStrsed = new HashMap<>();
            int[][] treeList = GenList(multiStrs, idxmap);
            int i = 0;
            for (int[] readyAlign : treeList) {
                String outToScreen = "  " + (i++) + "/" + treeList.length;
                System.out.print(outToScreen);
                String[] strsA, strsB;
                int[] IdxA, IdxB, IdxC;
                if (readyAlign[0] < idxmap.size()) {
                    strsA = multiStrs.remove(idxmap.get(readyAlign[0]));
                    IdxA = cluIdx.remove(idxmap.get(readyAlign[0]));
                } else {
                    strsA = NewStrsed.remove(readyAlign[0]);
                    IdxA = newcluIdx.remove(readyAlign[0]);
                }
                if (readyAlign[1] < idxmap.size()) {
                    strsB = multiStrs.remove(idxmap.get(readyAlign[1]));
                    IdxB = cluIdx.remove(idxmap.get(readyAlign[1]));
                } else {
                    strsB = NewStrsed.remove(readyAlign[1]);
                    IdxB = newcluIdx.remove(readyAlign[1]);
                }
                NewStrsed.put(readyAlign[2], profileAlign.Align(strsA, strsB, alphabet, kk));
                IdxC = new int[IdxA.length + IdxB.length];
                System.arraycopy(IdxA, 0, IdxC, 0, IdxA.length);
                System.arraycopy(IdxB, 0, IdxC, IdxA.length, IdxB.length);
                newcluIdx.put(readyAlign[2], IdxC);

                System.out.print(string.repeat("\b", outToScreen.length()));
            }
            strsAligned = NewStrsed.get(treeList[treeList.length - 1][2]);
            order = newcluIdx.get(treeList[treeList.length - 1][2]);
        }
        System.out.println("[" + sdf.format(new Date()) + "] Done.");
    }

    private int[][] GenList(Map<Integer, String[]> multiStrs, Map<Integer, Integer> idxMap) {
        String[] strs = new String[multiStrs.size()];
        int tempIdx = 0;
        for (int key : multiStrs.keySet()) {
            strs[tempIdx] = multiStrs.get(key)[0].replaceAll("-", "");
            idxMap.put(tempIdx++, key);
        }
        guidetree gTree = new guidetree(strs, "cluster");
        return gTree.genTreeList(true);
    }

    private void insertGapStrings(String str, String[] strings) {
        int[] mark = new int[strings[0].length() + 1];
        int i = 0, nums = 0;
        for (char c : str.toCharArray()) {
            if (c == '-')
                nums++;
            else {
                mark[i++] = nums;
                nums = 0;
            }
        }
        mark[i] = nums;
        for (i = 0; i < strings.length; i++) {
            strings[i] = insertGap(mark, strings[i]);
        }
    }

    private String insertGap(int[] mark, String seq) {
        assert mark.length == seq.length() + 1;
        StringBuilder seqGap = new StringBuilder();
        int len = mark.length;
        for (int i = 0; i < len; i++) {
            seqGap.append(string.repeat("-", mark[i]));
            if (i < len - 1)
                seqGap.append(seq.charAt(i));
        }
        return seqGap.toString();
    }

    private String PickRealOne(String[] strings) {
        char[] res = new char[strings[0].length()];
        for (int i = 0; i < res.length; i++) {
            HashMap<Character, Integer> charCounter = new HashMap<>();
            for (String string : strings) {
                char tempc = string.charAt(i);
                charCounter.put(tempc, charCounter.containsKey(tempc) ? charCounter.get(tempc) + 1 : 1);
            }
            charCounter.remove('-');
            char resc = 'n';
            int num = 0;
            for (char c : charCounter.keySet()) {
                if (charCounter.get(c) > num) {
                    num = charCounter.get(c);
                    resc = c;
                }
            }
            res[i] = resc;
        }
        return new String(res);
    }
}
