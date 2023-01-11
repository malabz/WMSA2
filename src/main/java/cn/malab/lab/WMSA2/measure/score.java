package cn.malab.lab.WMSA2.measure;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;

import cn.malab.lab.WMSA2.msa.centerAlign;
import cn.malab.lab.WMSA2.sample.sampleStrings;

public class score {
    /**
     * compute the sps score
     * 
     * @return score
     * @throws IOException
     */
    public static double avgSps(String path) throws IOException {
        int length = 0;
        try (BufferedReader br = new BufferedReader(new FileReader(path))) {
            String temp;
            StringBuilder line = new StringBuilder();
            
            while ((temp = br.readLine()) != null) {
                if (line.length() > 0 && temp.charAt(0) == '>') {
                    length = line.length();
                    break;
                } else if (temp.length() > 0 && temp.charAt(0) != '>') {
                    line.append(temp);
                }
            }
        }
        // "A" "G" "C" "T/U" "-" "N"
        System.out.println("     length: " + length);
        int size = 0;
        int[][] alnums = new int[length][6];

        try (BufferedReader br = new BufferedReader(new FileReader(path))) {
            String temp;
            StringBuilder line = new StringBuilder();
            
            while ((temp = br.readLine()) != null) {
                if (temp.length() > 0 && temp.charAt(0) == '>') {
                    size++;
                    int idx = 0;
                    for (char c : line.toString().toCharArray()) {
                        if (c == 'A' || c == 'a') {
                            alnums[idx][0]++;
                        } else if (c == 'G' || c == 'g') {
                            alnums[idx][1]++;
                        } else if (c == 'C' || c == 'c') {
                            alnums[idx][2]++;
                        } else if (c == 'U' || c == 'u' || c == 'T' || c == 't') {
                            alnums[idx][3]++;
                        } else if (c == '-') {
                            alnums[idx][4]++;
                        } else {
                            alnums[idx][5]++;
                        }
                        idx++;
                    }

                    line = new StringBuilder();
                } else if (temp.length() > 0 && temp.charAt(0) != '>') {
                    line.append(temp);
                }
            }
            int idx = 0;
            for (char c : line.toString().toCharArray()) {
                if (c == 'A' || c == 'a') {
                    alnums[idx][0]++;
                } else if (c == 'G' || c == 'g') {
                    alnums[idx][1]++;
                } else if (c == 'C' || c == 'c') {
                    alnums[idx][2]++;
                } else if (c == 'U' || c == 'u' || c == 'T' || c == 't') {
                    alnums[idx][3]++;
                } else if (c == '-') {
                    alnums[idx][4]++;
                } else {
                    alnums[idx][5]++;
                }
                idx++;
            }
        }
        System.out.println("       nums: " + size);
        double avgscore = 0;
        int ms = 1, mis = -1, gap = -2, odd = size, even = size - 1;
        if ((size & 1) == 0) {
            odd--;
            even++;
        }
        even >>= 1;
        for (int i = 0; i < length; i++) {
            long match = 0, mismatch, gapmis;
            for (int j = 0; j < 4; j++) {
                match += (long)alnums[i][j] * (alnums[i][j] - 1);
            }
            match >>= 1;
            match += (long) alnums[i][5] * (alnums[i][0] + alnums[i][1] + alnums[i][2] + alnums[i][3]);
            match += (long) (alnums[i][5] >> 1) * (alnums[i][5] - 1);
            mismatch = (alnums[i][0] + alnums[i][1]) * ((long)alnums[i][2] + alnums[i][3]);
            mismatch += (long) alnums[i][0] * alnums[i][1] + (long)alnums[i][2] * alnums[i][3];
            gapmis = (alnums[i][0] + alnums[i][1] + alnums[i][2] + alnums[i][3] + alnums[i][5]) * (long)alnums[i][4];

            avgscore += ((match / (double)odd) / (double) even) * (double)ms;
            avgscore += ((mismatch / (double)odd) / (double) even) * (double)mis;
            avgscore += ((gapmis / (double)odd) / (double) even) * (double)gap;
        }

        return avgscore;
    }


    /**
     * try to get an appropriate k
     * 
     * @param strs
     * @param sampled
     * @return k
     */
    public static int getK(String[] strs, boolean sampled) {
        if (!sampled) {
            strs = sampleStrings.getSampleStrs(strs);
        }
        centerAlign cAlign = new centerAlign(strs, true);
        strs = cAlign.getStrsAlign();
        int nums = strs.length, gap = 0;
        for (int i = 0; i < nums; i++) {
            for (int j = i + 1; j < nums; j++) {
                gap += countGap(strs[i], strs[j]);
            }
        }
        return Math.max(gap / (nums * nums / 2), 1);
    }

    private static int countGap(String A, String B) {
        int nums = 0, len = A.length();
        int lenA = 0, lenB = 0;
        for (int i = 0; i < len; i++) {
            if (A.charAt(i) == B.charAt(i) && A.charAt(i) == '-')
                nums++;
            if (A.charAt(i) != '-')
                lenA++;
            if (B.charAt(i) != '-')
                lenB++;
        }
        return Math.min(len - lenA - nums, len - lenB - nums);
    }
}
