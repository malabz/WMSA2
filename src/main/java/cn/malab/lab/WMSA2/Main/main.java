package cn.malab.lab.WMSA2.Main;

import cn.malab.lab.WMSA2.io.Fasta;
import cn.malab.lab.WMSA2.msa.ClusterAlign;
import cn.malab.lab.WMSA2.msa.treeAlign;

import java.io.FileWriter;
import java.io.IOException;
import java.io.Writer;
import java.text.SimpleDateFormat;
import java.util.Date;

public class Main {
    private static String mode;
    private static String infile;
    private static String outfile;
    private static double sim = 0.9;

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
            case "pro":
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
            // 混合策略比对
            case "win":
                res = readData(sdf);
                labels = res[0];
                strs = res[1];
                ClusterAlign clalign = new ClusterAlign(strs, sim, "t", "t2");
                String[] strsClal = clalign.getStrsAlign();
                Fasta.writeFasta(strsClal, labels, outfile, true);
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
        Fasta.countInfo(res[1]);
        return res;
    }

    private static void parse(String[] args) throws IOException {
        // 比对方式选择 -m
        // 输入文件位置 -i
        // 输出文件位置 -o
        if (args.length == 0 || args.length > 8) {
            args_help();
            System.exit(0);
        }
        for (int i = 0; i < args.length; i++) {
            if (args[i].equals("-m") && args.length > i + 1) {
                if (args[i + 1].equalsIgnoreCase("pro") || args[i + 1].equalsIgnoreCase("win")) {
                    mode = args[++i].toLowerCase();
                } else {
                    args_help();
                    throw new IllegalArgumentException("unknown mode: " + args[i + 1]);
                }
            } else if (args[i].equals("-i") && args.length > i + 1) {
                infile = args[++i];
            } else if (args[i].equals("-o") && args.length > i + 1) {
                outfile = args[++i];
            } else if (args[i].equals("-s") && args.length > i + 1) {
                sim = Double.parseDouble(args[++i].trim());
            } else {
                args_help();
                throw new IllegalArgumentException("unknown: " + args[i]);
            }
        }

        if (infile == null) {
            args_help();
            throw new IllegalArgumentException("infile path is null");
        } else if (outfile == null) {
            args_help();
            throw new IllegalArgumentException("outfile path is null");
        } else if (mode == null) {
            mode = "pro";
        }

        // 判断输出是否可写入
        try (Writer ignored = new FileWriter(outfile)) {}

    }

    private static void args_help() {
        System.out.println("\nusage: java -jar WMSA2.jar " + " [-m] mode [-s] sim [-i] path [-o] path");
        System.out.println();
        System.out.println("  necessary arguments: ");
        System.out.println("    -i  Input file path (nucleotide sequences in fasta format)");
        System.out.println("    -o  Output file path");
        System.out.println();
        System.out.println("  optional arguments: ");
        System.out.println("    -m  three align option (default mode: Pro)");
        System.out.println("         1. Pro   more accurate but slower");
        System.out.println("         2. Win   less accurate but faster");        
        System.out.println("    -s  the similarity of the cluster (used in Win mode. default: 0.9)");
        System.out.println();
    }

    private static void print_args() {
        System.out.println("\n**");
        System.out.println("** mode: " + mode);
        System.out.println("** infile: " + infile);
        System.out.println("** outfile: " + outfile);
        if (mode.equals("win")) {
            System.out.println("** similarity: " + sim);
        }
        System.out.println("**");
        System.out.println();
    }
}

