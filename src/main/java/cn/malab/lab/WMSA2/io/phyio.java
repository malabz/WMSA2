package cn.malab.lab.WMSA2.io;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.Writer;
import java.util.Deque;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.Map;

public class phyio {
    /**
     * read Newick file and generate TreeList 
     */
    public static int[][] readAndGenTreeList(String phypath, String[] labels, int N) throws Exception {
        StringBuilder line = new StringBuilder();
        try (BufferedReader br = new BufferedReader(new FileReader(phypath))) {
            String temp;
            while ((temp = br.readLine()) != null) {
                line.append(temp);
            }
        }
        int[][] treeList = new int[N-1][3];
        if (labels != null && labels.length > 0) {
            N = labels.length;
            int global_N = N, idx = 0;

            Map<String, Integer> map = new HashMap<>();
            for (int i = 0; i < N; i++) {
                map.put(labels[i].substring(1), i);
            }
            
            Deque<Integer> stack = new LinkedList<>();
            for (int i = 0; i < line.length(); i++) {
                if (line.charAt(i) == ':' && line.charAt(i-1) != ')') {
                    int tmp = i - 1;
                    while (line.charAt(tmp) != '(' && line.charAt(tmp) != ',') {
                        tmp--;
                    }
                    stack.push(map.get(line.substring(tmp + 1, i)));
                } else if (line.charAt(i) == ')') {
                    treeList[idx++] = new int[]{stack.pop(), stack.pop(), global_N};
                    stack.push(global_N++);
                }
            }
            if (idx != N - 1) {
                throw new Exception("idx = " + idx);
            }
        } else {
            int global_N = N, idx = 0;
            Deque<Integer> stack = new LinkedList<>();
            int numsToAgger = 0;
            for (int i = 0; i < line.length(); i++) {
                if (Character.isDigit(line.charAt(i))) {
                    int tmp = i + 1;
                    while (Character.isDigit(line.charAt(tmp))) {
                        tmp++;
                    }
                    stack.push(Integer.parseInt(line.substring(i, tmp)) -  1);
                    i = tmp - 1;
                } else if (line.charAt(i) == ')') {
                    while (numsToAgger-- > 1) {
                        treeList[idx++] = new int[]{stack.pop(), stack.pop(), global_N};
                        stack.push(global_N++);
                    }
                    treeList[idx++] = new int[]{stack.pop(), stack.pop(), global_N};
                    stack.push(global_N++);
                } else if (line.charAt(i) == '(') {
                    numsToAgger = 0;
                } else if (line.charAt(i) == ',') {
                    numsToAgger++;
                }
            }
            if (idx != N - 1) {
                throw new Exception("idx = " + idx);
            }
        }
        return treeList;
    }

    public static void writeNewick(String path, String tree) throws IOException {
        try (Writer write = new FileWriter(path); BufferedWriter bw = new BufferedWriter(write)) {
            bw.write(tree + "\n");
        }
    }
}
