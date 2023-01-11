package cn.malab.lab.WMSA2.io;

import java.io.*;
import java.text.SimpleDateFormat;
import java.util.*;

/**
 * read or write the fasta file
 */
public class Fasta {
    /**
     * 
     * @param path
     * @return String[lables[], Strings[]]
     * @throws IOException
     */
    public static String[][] readFasta(String path) throws IOException {
        List<String> strs = new ArrayList<>();
        List<String> labels = new ArrayList<>();

        try (BufferedReader br = new BufferedReader(new FileReader(path))) {
            String temp;
            StringBuilder line = new StringBuilder();
            while ((temp = br.readLine()) != null) {
                if (temp.length() > 0 && temp.charAt(0) == '>') {
                    if (line.length() != 0) {
                        strs.add(line.toString().toLowerCase());
                        line = new StringBuilder();
                    }
                    labels.add(temp.toLowerCase());
                } else if (temp.length() > 0) {
                    line.append(temp);
                }
            }
            strs.add(line.toString().toLowerCase());
        }
        return new String[][] { labels.toArray(new String[0]), strs.toArray(new String[0]) };
    }

    public static void countInfo(String[] strs) {
        Map<Character, Integer> map = new HashMap<Character, Integer>() {
            {
                put('a', 0);
                put('g', 0);
                put('c', 0);
                put('t', 0);
                put('u', 0);
                put('-', 0);
                put('n', 0);
            }
        };
        Map<Character, Integer> Amchar = new HashMap<>();
        long length = 0, minl = Long.MAX_VALUE, maxl = 0;
        int n = 0, gap = 0;
        for (int i = 0; i < strs.length; i++) {
            length += strs[i].length();
            minl = Math.min(minl, strs[i].length());
            maxl = Math.max(maxl, strs[i].length());
            int tempn = 0;
            for (char c : strs[i].toCharArray()) {
                if (!map.containsKey(c)) {
                    tempn++;
                    n++;
                    if (Amchar.containsKey(c))
                        Amchar.put(c, Amchar.get(c) + 1);
                    else
                        Amchar.put(c, 1);
                } else if (c == '-') {
                    tempn++;
                    map.put(c, map.get(c) + 1);
                }
            }
            if (tempn > 0) {
                strs[i] = strs[i].replaceAll("[-]", "");
                strs[i] = strs[i].replaceAll("[^AGCTUNagctun]", "n");
                gap += map.get('-');
                map.put('-', 0);
            }
        }
        SimpleDateFormat sdf = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss");
        System.out.println("[" + sdf.format(new Date()) + "] decected " + strs.length + " seqs");
        System.out.print("[" + sdf.format(new Date()) + "] max/min/ave len: " + maxl + "/" + minl + "/" + length / strs.length);
        if (gap > 0) {
            System.out.print("  gaps(deleted):" + gap);
        }
        if (n > 0) {
            System.out.print("  ambiguous characters:" + n);
            for (char c : Amchar.keySet()) {
                System.out.print(" " + c + ":" + Amchar.get(c));
            }
        }
        System.out.println();
    }

    /**
     * 
     * @param strings
     * @param labels
     * @param path
     */
    public static void writeFasta(String[] strings, String[] labels, String path, boolean width) throws IOException {
        try (Writer write = new FileWriter(path); BufferedWriter bw = new BufferedWriter(write)) {
            for (int i = 0; i < strings.length; i++) {
                bw.write(labels[i] + "\n");
                int j = 1;
                if (width) {
                    for (; j * 60 < strings[i].length(); j++) {
                        bw.write(strings[i].substring(j * 60 - 60, j * 60) + "\n");
                    }
                    bw.write(strings[i].substring(j * 60 - 60) + "\n");
                } else {
                    bw.write(strings[i].toUpperCase() + "\n");
                }
            }
        }
    }
}