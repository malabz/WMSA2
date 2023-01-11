package cn.malab.lab.WMSA2.hierCluster;

public class leafnode extends node {
    private final String name;
    private final int idx, num;

    public leafnode(String name, int idx) {
        this.name = name;
        this.idx = idx;
        this.num = idx;
    }

    public leafnode(String name, int idx, int num) {
        this.name = name;
        this.idx = idx;
        this.num = num;
    }

    @Override
    public String toString() {
        return name;
    }

    @Override
    public double getDistance() {
        return 0.0;
    }

    @Override
    public int getIdx() {
        return idx;
    }

    @Override
    public int getNum() {
        return num;
    }
}
