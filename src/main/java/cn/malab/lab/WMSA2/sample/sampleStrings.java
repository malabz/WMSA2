package cn.malab.lab.WMSA2.sample;


public class sampleStrings {

    public static String[] getSampleStrs(String[] strs) {
        int len = strs.length / 10;
        if (len <= 1) {
            return strs;
        }
        String[] sampledStrs = new String[20];
        for (int i = 0, j = 0; i < 10; i++) {
            int[] temp = LongShort(strs, i * len, (i + 1) * len);
            assert temp != null;
            sampledStrs[j++] = strs[temp[0]];
            sampledStrs[j++] = strs[temp[1]];
        }
        return sampledStrs;
    }

    private static int[] LongShort(String[] strs, int start, int end) {
        if (start >= strs.length) {
            return null;
        }
        end = Math.min(end, strs.length);
        int longOne = start, shortOne = end - 1;
        for (int j = start; j < end; j++) {
            if (strs[j].length() > strs[longOne].length()) {
                longOne = j;
            } else if (strs[j].length() <= strs[shortOne].length()) {
                shortOne = j;
            }
        }
        return new int[] { longOne, shortOne };
    }

}
