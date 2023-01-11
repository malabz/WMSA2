package cn.malab.lab.WMSA2.measure;

import java.util.Arrays;
import java.util.HashMap;

/**
 * 统计kmer，计算字母表的大小，计算字符之间的距离
 */
public class kmer {
    private final String[] strs;
    private final int num, k;

    /**
     * @param k kmer的步长
     */
    public kmer(String[] strs, int k) {
        this.k = k;
        this.strs = strs;
        this.num = strs.length;
    }

    /**
     * the number of alphabet is smaller than 5
     *
     * @return alphabet
     */
    public static char[] Counter(String[] strs) {
        HashMap<Character, Integer> alphabet = new HashMap<>();
        for (String str : strs) {
            char[] chars = str.toCharArray();
            for (char c : chars) {
                alphabet.put(c, alphabet.containsKey(c) ? alphabet.get(c) + 1 : 1);
            }
        }

        char[] alpbet = new char[alphabet.keySet().size()];
        int i = 0;
        for (char c : alphabet.keySet()) {
            alpbet[i++] = c;
        }
        if (alpbet.length <= 4) {
            return alpbet;
        }
        HashMap<Integer, Character> newAl = new HashMap<>();
        char[] keys = new char[alpbet.length];
        int[] values = new int[alpbet.length];
        i = 0;
        for (char key : alphabet.keySet()) {
            keys[i] = key;
            values[i++] = alphabet.get(key);
        }
        i = 0;
        for (int value : values) {
            while (newAl.containsKey(value)) {
                value++;
            }
            values[i] = value;
            newAl.put(value, keys[i++]);
        }
        Arrays.sort(values);
        char[] newalphabet = new char[4];
        for (i = values.length - 1; i > values.length - 5; i--) {
            newalphabet[values.length - 1 - i] = newAl.get(values[i]);
        }
        return newalphabet;
    }

    /**
     * to get the idx matrix of k mer
     *
     * @return kmer[][]
     */
    public static int[][] Counterk(String[] strs, int k) {
        int num = strs.length;
        HashMap<String, Integer> wordsIdx = new HashMap<>();
        int temp = 0;
        for (String str : strs) {
            // for (int i = 0; i <= str.length() - this.k; i += this.k) {
            for (int i = 0; i <= str.length() - k; i++) {
                String kstr = str.substring(i, i + k);
                if (!wordsIdx.containsKey(kstr)) {
                    wordsIdx.put(kstr, temp++);
                }
            }
        }
        int[][] kmerAll = new int[num][wordsIdx.keySet().size()];
        for (int i = 0; i < num; i++) {
            String str = strs[i];
            // for (int j = 0; j <= str.length() - this.k; j += this.k) {
            for (int j = 0; j <= str.length() - k; j++) {
                kmerAll[i][wordsIdx.get(str.substring(j, j + k))]++;
            }
        }
        return kmerAll;
    }

    /**
     *  得到中心序列与其他序列之间的距离
     *  @param center
     *  @return distance
     */
    public double[] getDismatrix1D(int center) {
        int[][] kmerAll = Counterk(strs, k);
        double[] dists = new double[num-1];
        for (int i = 0; i < center; i++) {
            dists[i] = getDistance(kmerAll[center], kmerAll[i]);
        }
        for (int i = center + 1; i < num; i++) {
            dists[i-1] = getDistance(kmerAll[center], kmerAll[i]);
        }
        return dists;
    }

    /**
     * compute the dis-similarity between the strs
     * 
     * @return distance[][]
     */
    public double[][] getDismatrix() {
        int[][] kmerAll = Counterk(strs, k);
        double[][] dismatrix = new double[num][num];
        for (int i = 0; i < num; i++) {
            for (int j = i + 1; j < num; j++) {
                // 不相似度
                dismatrix[i][j] = dismatrix[j][i] = 1 - getDistance(kmerAll[i], kmerAll[j]);
            }
        }
        return dismatrix;
    }

    public static double getDistance(int[] kmerA, int[] kmerB) {
        double sim;
        long respoint = 0;
        long lenA = 0, lenB = 0;
        for (int i = 0; i < kmerA.length; i++) {
            respoint += (long) kmerA[i] * kmerB[i];
            lenA += (long) kmerA[i] * kmerA[i];
            lenB += (long) kmerB[i] * kmerB[i];
        }
        sim = (double) respoint / Math.sqrt(lenA * lenB);
        return sim;
    }

}
