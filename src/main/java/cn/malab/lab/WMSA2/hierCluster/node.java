package cn.malab.lab.WMSA2.hierCluster;

public abstract class node {
    private int minimumIdx;
    private double branchlength, minimumLength;

    public abstract int getNum();
    public abstract int getIdx();
    public abstract double getDistance();

    public double getLen() {
        return branchlength;
    }

    public void setLen(double length) {
        branchlength = length;
    }

    public double getMinimumLength() {
        return minimumLength;
    }

    public void setMinimumLength(double minimumLength) {
        this.minimumLength = minimumLength;
    }

    public int getMinimumIdx() {
        return minimumIdx;
    }

    public void setMinimumIdx(int minimumIdx) {
        this.minimumIdx = minimumIdx;
    }
}
