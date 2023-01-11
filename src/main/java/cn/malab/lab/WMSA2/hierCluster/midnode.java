package cn.malab.lab.WMSA2.hierCluster;

public class midnode extends node {
    private final node a, b;
    private final int idx, num;

    /**
     * for efficient upgma
     * @param a
     * @param b
     * @param idx
     * @param num
     */
    public midnode(node a, node b, int idx, int num) {
        this.a = a;
        this.b = b;
        this.idx = idx;
        this.num = num;
    }

    /**
     * for normal upgma
     * @param a
     * @param b
     * @param num
     */
    public midnode(node a, node b, int num) {
        this.a = a;
        this.b = b;
        this.idx = num;
        this.num = num;
    }

    public String toString() {
        return "(" + a + ":" + a.getLen() + ", " + b + ":" + b.getLen() + ")";
    }

    @Override
    public double getDistance() {
        return a.getLen() + a.getDistance();
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
