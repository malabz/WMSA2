package cn.malab.lab.WMSA2.Main;

import java.io.FileWriter;
import java.io.IOException;
import java.io.Writer;
import java.text.SimpleDateFormat;
import java.util.Date;

import cn.malab.lab.WMSA2.hierCluster.guidetree;
import cn.malab.lab.WMSA2.io.Fasta;
import cn.malab.lab.WMSA2.io.phyio;
import cn.malab.lab.WMSA2.measure.score;
import cn.malab.lab.WMSA2.msa.ClusterAlign;
import cn.malab.lab.WMSA2.msa.centerAlign;
import cn.malab.lab.WMSA2.msa.treeAlign;

public class Main {
    private static String mode;
    private static String infile;
    private static String outfile;
    private static String newickfile;
    private static String newickmode;

    public static void main(String[] args) throws Exception {
        // 1. 解析输入参数
        parse(args);
        // 2. 打印输入参数
        print_args();

        // 参数
        SimpleDateFormat sdf = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss");
        String[][] res;
        String[] labels, strs;

        // 5. 匹配比对模式
        switch (mode) {
            // 树比对
            case "tree":
                res = readData(sdf);
                labels = res[0];
                strs = res[1];
                // cluster模式：StarTree模式构建指导树
                treeAlign talign = new treeAlign(strs, "cluster", false);
                // 得到比对结果
                String[] strsTal = talign.getStrsAlign();
                // 将比对结果写入到文件中
                Fasta.writeFasta(strsTal, labels, outfile, true);
                break;
            // 星比对
            case "center":
                res = readData(sdf);
                labels = res[0];
                strs = res[1];
                centerAlign calign = new centerAlign(strs, "fmindex");
                String[] strsCal = calign.getStrsAlign();
                Fasta.writeFasta(strsCal, labels, outfile, true);
                break;
            // 混合策略比对
            case "mix":
                res = readData(sdf);
                labels = res[0];
                strs = res[1];
                ClusterAlign clalign = new ClusterAlign(strs, "t", "t2");
                String[] strsClal = clalign.getStrsAlign();
                Fasta.writeFasta(strsClal, labels, outfile, true);
                break;
            case "startree":
                res = readData(sdf);
                strs = res[1];
                // clusterTree cTree = new clusterTree(strs);
                guidetree gtree = new guidetree(strs, "upgma");
                gtree.genTreeList(true);
                break;
            case "withguidetree":
                res = readData(sdf);
                labels = res[0];
                strs = res[1];
                int[][] treelist;
                if (newickmode.equalsIgnoreCase("mbed")) {
                    treelist = phyio.readAndGenTreeList(newickfile, labels, labels.length);
                } else if (newickmode.equalsIgnoreCase("parttree")) {
                    treelist = phyio.readAndGenTreeList(newickfile, null, labels.length);
                } else {
                    throw new IllegalArgumentException("unkown mode: " + newickmode);
                }
                treeAlign NAlign = new treeAlign(strs, treelist, false);
                // 得到比对结果
                String[] strsaln = NAlign.getStrsAlign();
                // 将比对结果写入到文件中
                Fasta.writeFasta(strsaln, labels, outfile, true);
                break;
            case "spscore":
                System.out.println("Avg SpScore: " + score.avgSps(infile));
                break;
            default:
                args_help();
                throw new IllegalArgumentException("unkown mode: " + mode);
        }
    }

    private static String[][] readData(SimpleDateFormat sdf) throws IOException {
        System.out.println("[" + sdf.format(new Date()) + "] Reading data");
        String[][] res = Fasta.readFasta(infile);
        System.out.println("[" + sdf.format(new Date()) + "] Done.");
        // 打印序列的长度、数目信息
        if (!mode.equals("startree")) {
            Fasta.countInfo(res[1]);
        }
        return res;
    }

    private static void parse(String[] args) throws IOException {
        // 比对方式选择 -m
        // 输入文件位置 -i
        // 输出文件位置 -o
        // newick文件位置 -ni
        // newick文件模式 -nm
        if (args.length == 0 || args.length > 10) {
            args_help();
            System.exit(0);
        }
        for (int i = 0; i < args.length; i++) {
            if (args[i].equals("-m") && args.length > i + 1) {
                if (   args[i + 1].equalsIgnoreCase("mix")
                    || args[i + 1].equalsIgnoreCase("tree")
                    || args[i + 1].equalsIgnoreCase("center")
                    || args[i + 1].equalsIgnoreCase("spscore") 
                    || args[i + 1].equalsIgnoreCase("startree") 
                    || args[i + 1].equalsIgnoreCase("withGuideTree")) {
                    mode = args[++i].toLowerCase();
                } else {
                    args_help();
                    throw new IllegalArgumentException("unknown mode: " + args[i + 1]);
                }
            } else if (args[i].equals("-i") && args.length > i + 1) {
                infile = args[++i];
            } else if (args[i].equals("-o") && args.length > i + 1) {
                outfile = args[++i];
            } else if (args[i].equals("-ni") && args.length > i + 1) {
                newickfile = args[++i];
            } else if (args[i].equals("-nm") && args.length > i + 1) {
                newickmode = args[++i];
            } else {
                args_help();
                System.exit(0);
            }
        }

        if (infile == null) {
            args_help();
            throw new IllegalArgumentException("infile path is null");
        } else if (outfile == null && !(mode.equals("spscore") || mode.equals("startree"))) {
            args_help();
            throw new IllegalArgumentException("outfile path is null");
        } else if (mode == null) {
            mode = "mix";
        }
        if (!(mode.equals("spscore") || mode.equals("startree"))) {
            try (Writer ignored = new FileWriter(outfile)) {}
        }
    }

    private static void args_help() {
        System.out.println("\nusage: java -jar " + " [-m] mode [-i] path [-o] path");
        System.out.println();
        System.out.println("  necessary arguments: ");
        System.out.println("    -i  Input file path (nucleotide sequences in fasta format)");
        System.out.println("    -o  Output file path");
        System.out.println();
        System.out.println("  optional arguments: ");
        System.out.println("    -m  three align option (default mode: Mix)");
        System.out.println("         1. Tree   more accurate but slower");
        System.out.println("         2. Mix    less accurate but faster");
        System.out.println();
    }

    private static void print_args() {
        System.out.println("\n**");
        System.out.println("** mode: " + mode);
        System.out.println("** infile: " + infile);
        System.out.println("** outfile: " + outfile);
        System.out.println("**");
        System.out.println();
    }
}
