package cn.malab.lab.WMSA2.hierCluster;

import java.util.ArrayList;
import java.util.Deque;
import java.util.LinkedList;
import java.util.List;


public class effupgma {
    private final double[][] dmatrix;
    private final node[] nodes;
    private final int[] nums;
    private final int len;
    private node rootNode;
    private int remainder;
    private int global_n;
    private double headValue, secondValue;
    
    private List<int[]> TreeList;

    public List<int[]> getTreeList() {
        return TreeList;
    }

    public effupgma(double[][] matrix) {
        this.dmatrix = matrix;
        this.len = matrix.length;
        this.remainder = this.len;
        this.global_n = this.len;
        this.nums = new int[this.len];
        this.nodes = new node[this.len];
        this.TreeList = new ArrayList<>();

        for (int i = 0; i < len; i++) {
            nums[i] = 1;
            this.nodes[i] = new leafnode("" + i, i);
        }

        this.genTree();
    }

    public effupgma(double[][] matrix, int[] idxs, int global_n) {
        this.dmatrix = matrix;
        this.len = matrix.length;
        this.remainder = this.len;
        this.global_n = global_n;
        this.nums = new int[this.len];
        this.nodes = new node[this.len];
        this.TreeList = new ArrayList<>();

        for (int i = 0; i < len; i++) {
            this.nums[i] = 1;
            this.nodes[i] = new leafnode("" + idxs[i], i, idxs[i]);
        }

        this.genTree();
    }

    /**
     * 输出指导树用
     * @param matrix
     * @param names
     * @param idxs
     * @param global_n
     */
    public effupgma(double[][] matrix, String[] names, int[] idxs, int global_n) {
        this.dmatrix = matrix;
        this.len = matrix.length;
        this.remainder = this.len;
        this.global_n = global_n;
        this.nums = new int[this.len];
        this.nodes = new node[this.len];
        this.TreeList = new ArrayList<>();

        for (int i = 0; i < len; i++) {
            this.nums[i] = 1;
            this.nodes[i] = new leafnode(names[i], i, idxs[i]);
        }

        this.genTree();
    }

    public node getRootNode() {
        return rootNode;
    }

    private void genTree() {
        while (remainder >= 2) {
            Deque<node> queue = new LinkedList<>();
            // 先加入两个元素
            for (int i = 0; i < len; i++) {
                if (nodes[i] != null) {
                    // 添加第一个元素
                    queue.push(nodes[i]);
                    GetHeadValue(queue);
                    // 添加第二个元素
                    assert queue.peek() != null;
                    queue.push(nodes[queue.peek().getMinimumIdx()]);
                    GetHeadValue(queue);
                    break;
                }
            }
            while (headValue != secondValue) {
                // 只需要留下最后一个元素即可
                queue.removeLast();
                assert queue.peek() != null;
                queue.push(nodes[queue.peek().getMinimumIdx()]);
                GetHeadValue(queue);
            }
            // 合并两个节点
            combineNodes(queue.pop(), queue.pop());
            remainder--;
        }
    }


    private void combineNodes(node n1, node n2) {
        int i1 = n1.getIdx(), i2 = n2.getIdx();
        double min = this.dmatrix[i1][i2];
        n1.setLen(min / 2 - n1.getDistance());
        n2.setLen(min / 2 - n2.getDistance());        
        nodes[i1] = new midnode(n1, n2, i1, global_n);
        nodes[i2] = null;
        this.rootNode = nodes[i1];
        // 更新list
        TreeList.add(new int[] {n1.getNum(), n2.getNum(), global_n++});
        // 更新距离
        renewDist(i1, i2);
        // 更新nums值
        nums[i1] += nums[i2];
        nums[i2] = 0;
    }

    /**
     * renew whole area of the distance matrix
     */
    private void renewDist(int idxi, int idxj) {
        for (int i = 0; i < idxi; i++) {
            if (nodes[i] == null) {
                continue;
            }
            dmatrix[i][idxi] = (dmatrix[i][idxi] * nums[idxi] + dmatrix[i][idxj] * nums[idxj]) / (nums[idxi] + nums[idxj]);
            dmatrix[idxi][i] = dmatrix[i][idxi];
        }
        for (int j = idxi + 1; j < len; j++) {
            if (nodes[j] == null) {
                continue;
            }
            dmatrix[idxi][j] = (dmatrix[idxi][j] * nums[idxi] + dmatrix[idxj][j] * nums[idxj]) / (nums[idxi] + nums[idxj]);
            dmatrix[j][idxi] = dmatrix[idxi][j];
        }
        for (int i = 0; i < len; i++) {
            dmatrix[i][idxj] = dmatrix[idxj][i] = Double.POSITIVE_INFINITY;
        }
    }

    private void GetHeadValue(Deque<node> queue) {
        secondValue = headValue;
        node cur = queue.peek();
        assert cur != null;
        int idx = cur.getIdx(), target = -1;
        double minValue = 2.0d;
        for (int j = 0; j < len; j++) {
            if (j != idx && dmatrix[idx][j] < minValue) {
                minValue = dmatrix[idx][j];
                target = j;
            }
        }
        cur.setMinimumIdx(target);
        cur.setMinimumLength(minValue);
        headValue = minValue;
    }
}
