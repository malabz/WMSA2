package cn.malab.lab.WMSA2.msa;

import cn.malab.lab.WMSA2.io.string;
import cn.malab.lab.WMSA2.psa.FMAlign;
import cn.malab.lab.WMSA2.psa.Kband;
import cn.malab.lab.WMSA2.psa.STAlign;

public class centerAlign {
    private final String[] strsaligned;
    private final String mode;
    private final int num;
    private boolean silent = false;

    /**
     *
     * @param strs
     * @param mode "kband" "suffix" "fmindex"
     */
    public centerAlign(String[] strs, String mode) {
        this.num = strs.length;
        this.mode = mode;
        this.strsaligned = Align(strs);
    }

    /**
     * mode "fmindex"
     * 
     * @param strs
     */
    public centerAlign(String[] strs) {
        this.num = strs.length;
        this.mode = "fmindex";
        this.strsaligned = Align(strs);
    }

    /**
     * no output in silence mode "fmindex"
     * 
     * @param strs
     * @param silent
     */
    public centerAlign(String[] strs, boolean silent) {
        this.num = strs.length;
        this.mode = "fmindex";
        this.silent = silent;
        this.strsaligned = Align(strs);
    }

    public String[] getStrsAlign() {
        return this.strsaligned;
    }

    private String[][] alignToAll(String[] strs, int maxcol) {
        String[][] resultsAlign = new String[num - 1][2];
        switch (this.mode) {
            case "suffix":
                STAlign stAlign = new STAlign(strs[maxcol]);
                for (int i = 0; i < num - 1; i++) {
                    String outToScreen = "    " + (i + 1) + " / " + (num - 1);
                    if (!silent)
                        System.out.print(outToScreen);
                    int j = i >= maxcol ? i + 1 : i;
                    stAlign.AlignStrB(strs[j]);
                    resultsAlign[i] = stAlign.getStrAlign();
                    if (!silent)
                        System.out.print(string.repeat("\b", outToScreen.length()));
                }
                if (!silent)
                    System.out.println("    " + (num - 1) + " / " + (num - 1));
                break;
            case "kband":
                for (int i = 0; i < num - 1; i++) {
                    String outToScreen = "    " + (i + 1) + " / " + (num - 1);
                    if (!silent)
                        System.out.print(outToScreen);
                    int j = i >= maxcol ? i + 1 : i;
                    Kband kbAlign = new Kband(strs[maxcol], strs[j]);
                    resultsAlign[i] = kbAlign.getStrAlign();
                    if (!silent)
                        System.out.print(string.repeat("\b", outToScreen.length()));
                }
                System.out.println("    " + (num - 1) + " / " + (num - 1));
                break;
            case "fmindex":
                FMAlign fmAlign = new FMAlign(strs[maxcol]);
                for (int i = 0; i < num - 1; i++) {
                    String outToScreen = "    " + (i + 1) + " / " + (num - 1);
                    if (!silent)
                        System.out.print(outToScreen);
                    int j = i >= maxcol ? i + 1 : i;
                    fmAlign.AlignStrB(strs[j]);
                    resultsAlign[i] = fmAlign.getStrAlign();
                    if (!silent)
                        System.out.print(string.repeat("\b", outToScreen.length()));
                }
                if (!silent)
                    System.out.println("    " + (num - 1) + " / " + (num - 1));
                break;
            default:
                throw new IllegalArgumentException("unkown mode: " + mode);
        }
        return resultsAlign;
    }

    private int[] markGapsCenSeq(String[][] resultsAlign, String centerSeq) {
        int[] markInsertion = new int[centerSeq.length() + 1];
        for (String[] str2 : resultsAlign) {
            int i = 0, counter = 0;
            for (char c : str2[0].toCharArray()) {
                if (c == '-')
                    counter++;
                else {
                    markInsertion[i] = Math.max(markInsertion[i], counter);
                    counter = 0;
                    i++;
                }
            }
            markInsertion[i] = Math.max(markInsertion[i], counter);
        }
        return markInsertion;
    }

    private String[] insertGapTostrs(String[][] resultsAlign, String centerSeq, int[] markInsertion, int maxcol) {
        String[] newStrsaligned = new String[num];
        int idxinsert = 0;
        newStrsaligned[maxcol] = insertGap(markInsertion, centerSeq);

        for (String[] str2 : resultsAlign) {
            String outToScreen = "    " + (idxinsert + 1) + " / " + (num - 1);
            if (!silent)
                System.out.print(outToScreen);
            char[] tempA = str2[0].toCharArray();
            int[] mark = new int[tempA.length + 1];
            int pi = 0, pj = 0, total = 0;
            for (char c : tempA) {
                if (c == '-')
                    total++;
                else {
                    mark[pi++] = markInsertion[pj++] - total;
                    while (total != 0) {
                        pi++;
                        total--;
                    }
                }
            }
            mark[pi] = markInsertion[pj] - total;
            if (idxinsert >= maxcol) {
                newStrsaligned[++idxinsert] = insertGap(mark, str2[1]);
            } else {
                newStrsaligned[idxinsert++] = insertGap(mark, str2[1]);
            }
            if (!silent)
                System.out.print(string.repeat("\b", outToScreen.length()));
        }
        return newStrsaligned;
    }

    private String[] Align(String[] strs) {
        // 1. find the maxcol
        if (!silent)
            System.out.println("\nfind the center seq");
        int maxcol = findMaxcol(strs);
        if (!silent)
            System.out.println("\n    No." + maxcol);
        if (!silent)
            System.out.println("\nAlign the center seq\n");
        // 2. alignOneToAll
        String[][] resultsAlign = alignToAll(strs, maxcol);
        String centerSeq = strs[maxcol];
        // 3. mark the gaps in center seq
        int[] markInsertion = markGapsCenSeq(resultsAlign, centerSeq);
        // 4. insert the gap
        if (!silent)
            System.out.println("\nAlign all seqs\n");
        return insertGapTostrs(resultsAlign, centerSeq, markInsertion, maxcol);
    }

    private String insertGap(int[] mark, String seq) {
        StringBuilder seqGap = new StringBuilder();
        int len = mark.length;
        for (int i = 0; i < len; i++) {
            seqGap.append(string.repeat("-", mark[i]));
            if (i < len - 1)
                seqGap.append(seq.charAt(i));
        }
        return seqGap.toString();
    }

    private int findMaxcol(String[] strs) {
        int max = 0;
        for (int i = 0; i < num; i++) {
            if (strs[i].length() > strs[max].length())
                max = i;
        }
        return max;
    }
}
