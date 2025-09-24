package raj.cbm.util;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Arrays;

/**
 * Utility to load deterministic RNG seeds from a plain-text file into a
 * 2-D array indexed as {@code seeds[scenarioIndex][runIndex]}.
 *
 * <h3>Expected file format</h3>
 * The file should contain exactly {@code scenarios} lines. Each line must have
 * at least {@code runs} whitespace-separated integers (parsable as {@code long}),
 * one per run. A common convention is 30 integers per line.
 *
 * <pre>
 * # Example seeds.txt (3 scenarios × 5 runs shown)
 *  101  202  303  404  505
 *  111  222  333  444  555
 *  123  234  345  456  567
 * </pre>
 *
 * <p><b>Notes & pitfalls</b></p>
 * <ul>
 *   <li>This loader does no padding or truncation beyond reading the first {@code runs}
 *       tokens per line; if a line has fewer than {@code runs} tokens or the file has
 *       fewer than {@code scenarios} lines, you'll get an {@link IndexOutOfBoundsException}.</li>
 *   <li>Invalid numeric tokens will raise a {@link NumberFormatException}.</li>
 *   <li>IO failures (missing file, permissions, etc.) propagate as {@link IOException}.</li>
 * </ul>
 *
 * @author rchoudhary
 * @since 1.0
 */
public class Seeds {

    /**
     * Loads a seeds matrix from {@code seeds.txt}-style file.
     *
     * @param path       path to the seeds file
     * @param scenarios  number of scenario rows to read (lines)
     * @param runs       number of seeds per scenario (columns) to read from each line
     * @return a {@code long[scenarios][runs]} matrix where {@code [i][r]} is the seed for
     *         scenario {@code i} and run {@code r}
     * @throws IOException if the file cannot be read
     * @throws IndexOutOfBoundsException if the file has fewer lines than {@code scenarios}
     *                                   or any line has fewer than {@code runs} tokens
     * @throws NumberFormatException if a token cannot be parsed as a {@code long}
     */
    public static long[][] load(String path, int scenarios, int runs) throws IOException {
        var lines = Files.readAllLines(Paths.get(path));
        long[][] s = new long[scenarios][runs];
        if (lines.size() < scenarios) {
            throw new IllegalArgumentException("seeds.txt has fewer lines (" + lines.size() + ") than scenarios (" + scenarios + ")");
        }
        for (int i = 0; i < scenarios; i++) {
            String[] toks = lines.get(i).trim().split("\\s+");
            if (toks.length < runs) {
                throw new IllegalArgumentException("seeds.txt line " + (i+1) + " has only " + toks.length + " integers; expected at least " + runs);
            }
            for (int r = 0; r < runs; r++) {
                s[i][r] = Long.parseLong(toks[r]);
            }
        }
        return s;
    }
}
