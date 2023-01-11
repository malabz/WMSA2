package cn.malab.lab.WMSA2.psa;

import cn.malab.lab.WMSA2.io.string;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;

public class multiDP {

    private final int ms = 7, mis = -3;
    private final int d = 15, e = 2;
    private final String[] A, B;
    private String[] alignA, alignB;
    private final int numA, numB;
    // private int lenA, lenB;
    private int[][] alphA, alphB;
    private boolean state = false;
    HashMap<Character, Integer> mapAlph;

    /**
     * 
     * @param A
     * @param B
     * @param alphabet
     */
    public multiDP(String[] A, String[] B, char[] alphabet) {
        if (A[0].length() > B[0].length() && B[0].length() > 0) {
            this.A = B;
            this.B = A;
            this.numA = B.length;
            this.numB = A.length;
            this.state = true;
        } else {
            this.A = A;
            this.B = B;
            this.numA = A.length;
            this.numB = B.length;
        }
        this.InitalAB();
        this.CountAlphabet(alphabet);
        this.Align();
    }

    /**
     * To get the aligned results.
     */
    public String[][] getStrAlign() {
        if (state)
            return new String[][] { alignB, alignA };
        return new String[][] { alignA, alignB };
    }

    public String[] getStrsAlign() {
        String[] res = new String[alignA.length + alignB.length];
        if (state) {
            System.arraycopy(alignB, 0, res, 0, alignB.length);
            System.arraycopy(alignA, 0, res, alignB.length, alignA.length);
        } else {
            System.arraycopy(alignA, 0, res, 0, alignA.length);
            System.arraycopy(alignB, 0, res, alignA.length, alignB.length);
        }
        return res;
    }

    /**
     * init the result arrays
     */
    private void InitalAB() {
        this.alignA = new String[numA];
        this.alignB = new String[numB];
        for (int i = 0; i < numA; i++) {
            this.alignA[i] = "";
        }
        for (int i = 0; i < numB; i++) {
            this.alignB[i] = "";
        }
    }

    private void CountAlphabet(char[] alphabet) {
        mapAlph = new HashMap<>();
        for (int i = 0; i < alphabet.length; i++)
            mapAlph.put(alphabet[i], i);
        mapAlph.put('-', alphabet.length);
        mapAlph.put('n', alphabet.length + 1);
        alphA = new int[A[0].length()][alphabet.length + 2];
        alphB = new int[B[0].length()][alphabet.length + 2];
        for (int i = 0; i < A[0].length(); i++) {
            for (String s : A) {
                char c = mapAlph.containsKey(s.charAt(i)) ? s.charAt(i) : 'n';
                alphA[i][mapAlph.get(c)] += 1;
            }
        }
        for (int i = 0; i < B[0].length(); i++) {
            for (String s : B) {
                char c = mapAlph.containsKey(s.charAt(i)) ? s.charAt(i) : 'n';
                alphB[i][mapAlph.get(c)] += 1;
            }
        }
    }

    private int Match(int idxm, int idxn) {
        int[] tempA = alphA[idxm], tempB = alphB[idxn];
        int results = 0, len = tempA.length;
        for (int i = 0; i < len - 2; i++)
            results += (tempA[i] * tempB[i]);
        results *= (this.ms - this.mis);
        results += (numA - tempA[len - 1] - tempA[len - 2]) * (numB - tempB[len - 1] - tempB[len - 2]) * this.mis;
        results += (numA - tempA[len - 2]) * (tempB[len - 1] * this.ms - tempB[len - 2] * this.e);
        results -= tempA[len - 2] * (numB - tempB[len - 2]) * this.e;
        results += tempA[len - 1] * (numB - tempB[len - 1] - tempB[len - 2]) * this.ms;
        return results / (numA * numB);
    }

    private float[][][] Init(int m, int n) {
        float[][][] p = new float[3][m + 1][n + 1];
        for (float[][] l : p) {
            for (int i = 0; i < m + 1; i++) {
                Arrays.fill(l[i], Float.NEGATIVE_INFINITY);
            }
        }
        p[0][0][0] = 0;
        for (int j = 1; j < n + 1; j++) {
            p[1][0][j] = -e * j;
        }
        for (int i = 1; i < m + 1; i++) {
            p[2][i][0] = -e * i;
        }
        return p;
    }

    private float Max3(float p0, float p1, float p2) {
        return Math.max(p0, Math.max(p1, p2));
    }

    private void Align() {
        // n >= m
        int m = A[0].length();
        int n = B[0].length();
        // special instances
        if (m == 0 && n == 0) {
            return;
        } else if (m == 0) {
            Arrays.fill(this.alignA, string.repeat("-", n));
            this.alignB = this.B;
            return;
        } else if (n == 0) {
            Arrays.fill(this.alignB, string.repeat("-", m));
            this.alignA = this.A;
            return;
        }
        float[][][] p = Init(m, n);

        for (int i = 1; i < m + 1; i++) {
            for (int j = 1; j < n + 1; j++) {
                // p[0] : A[i] ~ B[j]
                p[0][i][j] = Max3(p[0][i - 1][j - 1], p[1][i - 1][j - 1], p[2][i - 1][j - 1]) + Match(i - 1, j - 1);
                // p[1] : B[j] ~ -
                p[1][i][j] = Math.max(p[0][i][j - 1] - d, p[1][i][j - 1] - e);
                // p[2] : A[j] ~ -
                p[2][i][j] = Math.max(p[0][i - 1][j] - d, p[2][i - 1][j] - e);
            }
        }
        TraceBack(p, m, n);
    }

    private void TraceBack(float[][][] p, int m, int n) {
        int channel = p[1][m][n] >= p[0][m][n] ? (p[1][m][n] >= p[2][m][n] ? 1 : 2)
                : (p[0][m][n] >= p[2][m][n] ? 0 : 2);
        int i = m, j = n;
        List<Integer> trace = new ArrayList<>();
        while (i > 0 || j > 0) {
            if (channel == 1 && j > 0) {
                channel = -1;
                if (p[1][i][j] == p[1][i][j - 1] - e)
                    channel = 1;
                else if (i > 0 && j > 1 && p[1][i][j] == p[0][i][j - 1] - d)
                    channel = 0;
                trace.add(1);
                j--;
            } else if (channel == 0 && i > 0 && j > 0) {
                channel = -1;
                if (i > 1 && j > 1 && p[0][i][j] == p[0][i - 1][j - 1] + Match(i - 1, j - 1)) {
                    channel = 0;
                } else if (j > 1 && p[0][i][j] == p[1][i - 1][j - 1] + Match(i - 1, j - 1)) {
                    channel = 1;
                } else if (i > 1 && p[0][i][j] == p[2][i - 1][j - 1] + Match(i - 1, j - 1)) {
                    channel = 2;
                }
                trace.add(0);
                i--;
                j--;
            } else if (channel == 2 && i > 0) {
                channel = -1;
                if (p[2][i][j] == p[2][i - 1][j] - e)
                    channel = 2;
                else if (i > 1 && j > 0 && p[2][i][j] == p[0][i - 1][j] - d)
                    channel = 0;
                trace.add(2);
                i--;
            } else {
                throw new IllegalStateException("channel = " + channel);
            }
        }
        this.Bach(trace.toArray(new Integer[0]));
    }
    
    private void Bach(Integer[] traces) {
        char[][] alignA = new char[numA][traces.length];
        char[][] alignB = new char[numB][traces.length];
        int idxi = this.A[0].length(), idxj = this.B[0].length(), i = traces.length - 1;
        for (int trace : traces) {
            switch (trace) {
            case 0:
                for (int idxA = 0; idxA < numA; idxA++) {
                    alignA[idxA][i] = this.A[idxA].charAt(idxi - 1);
                }
                for (int idxB = 0; idxB < numB; idxB++) {
                    alignB[idxB][i] = this.B[idxB].charAt(idxj - 1);
                }
                idxi--;
                idxj--;
                break;
            case 1:
                for (int idxA = 0; idxA < numA; idxA++) {
                    alignA[idxA][i] = '-';
                }
                for (int idxB = 0; idxB < numB; idxB++) {
                    alignB[idxB][i] = this.B[idxB].charAt(idxj - 1);
                }
                idxj--;
                break;
            case 2:
                for (int idxA = 0; idxA < numA; idxA++) {
                    alignA[idxA][i] = this.A[idxA].charAt(idxi - 1);
                }
                for (int idxB = 0; idxB < numB; idxB++) {
                    alignB[idxB][i] = '-';
                }
                idxi--;
                break;
            }
            i--;
        }
        for (i = 0; i < numA; i++) {
            this.alignA[i] = new String(alignA[i]);
        }
        for (int j = 0; j < numB; j++) {
            this.alignB[j] = new String(alignB[j]);
        }
    }
}
