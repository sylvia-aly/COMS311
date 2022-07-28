package edu;

import java.util.ArrayList;

/**
 * @author: Sylvia Nguyen
 * 
 *          A class to implement a dynamic programming algorithm for computing
 *          the edit distance between two sequences and their alignment with
 *          this mimimum distance.
 *
 *          The algorithm computes the cost matrix in a reverse order so that
 *          the optimal alignment can be recovered in a forward order, with its
 *          first pair reported and saved at the beginning of each array.
 *
 *          Acknowledgements
 *
 *          This template is created by Xiaoqiu Huang. This project was
 *          conceived during discussions with the other instructors (David
 *          Fernandez-Baca and Seyedehzahra Hosseini) for Com S 311 in Spring
 *          2022.
 */

// import java.lang.reflect.Array;
// import java.util.Arrays;
// import java.lang.Math;
// import java.util.Random;

public class Alignment {
	private int editDistance; // The edit distance
	private int editnum = 0; // The number of differences in the alignment
	private int alignlen; // The length of the alignment
	private char[] top; // The top row of the alignment
	private char[] mid; // The middle row of the alignment
	private char[] bot; // The bottom row of the alignment
	private int[][] Smat; // Cost matrix S
	private char[] A; // Sequence A with its first letter in A[1]
	private int M; // Its length (== m)
	private char[] B; // Sequence B with its first letter in B[1]
	private int N; // Its length (== n)

	/*
	 * Input parameter seqone is a reference to a Sequence object containing the
	 * name, length and character array of one sequence. Input parameter seqtwo is
	 * for another sequence. Input parameter param is a reference to a
	 * SubstitutionGap object containing a match score, a mismatch penalty and a gap
	 * penalty. If seqone is null, throws an RuntimeException exception with the
	 * message "Sequence object one is null." If seqtwo is null, throws an
	 * RuntimeException exception with the message "Sequence object two is null." If
	 * param is null, throws an RuntimeException exception with the message
	 * "SubstitutionGap object is null."
	 * 
	 * This method computes the cost matrix in a reverse order and produces the
	 * alignment by doing traceback in a forward order, with the first alignment
	 * pair saved at the beginning of the three character arrays (top, mid and bot).
	 * 
	 * If at most one of the two sequences is empty (of length 0), the method is
	 * expected to work correctly by producing an alignment made of only deletions
	 * or insertions. However, if both sequences are of length 0, then the method
	 * throws an RuntimeException exception with the message
	 * "Both sequences are empty."
	 * 
	 * The method obtains one sequence and its length from the first parameter
	 * seqone, and saves them in the char[] array A and the int scalar variable M,
	 * with the first letter of the sequence stored at index 1 of A. Similarly, the
	 * other sequence and its length from the second parameter seqtwo are saved in
	 * the char[] array B and the int scalar variable N, with the first letter of
	 * the sequence stored at index 1 of B.
	 * 
	 * Next the method obtains a match score, a mismatch penalty, and an
	 * insertion/deletion penalty from the third parameter param. If letters a and b
	 * are identical, then sigma(a, b) is the match score; otherwise, sigma(a, b) is
	 * the mismatch penalty. The insertion/deletion penalty is denoted by r in the
	 * assignment specification.
	 * 
	 * The method computes the cost matrix S(i, j) for 0 <= i <= M and 0 <= j <= N,
	 * according to the algorithm in the assignment specification, where M and N are
	 * the sequence lengths m and n in the assignment specification, and A[i]
	 * contains the letter a_i and B[j] contains the letter b_j in the
	 * specification. The matrix S is implemented by using the int [][] array
	 * variable Smat of M+1 rows and N+1 columns. After the matrix S is computed,
	 * the method saves S(0, 0), the edit distance, in the int scalar variable
	 * editDistance.
	 * 
	 * Next the method generates an optimal sequence alignment by using the
	 * traceback algorithm in the assignment specification. To save the alignment,
	 * the three char[] arrays top, mid and bot of sufficient length, say, M+N+2,
	 * are created. Let (x_1, y_1) (x_2, y_2) ... (x_k, y_k) be all pairs of the
	 * alignment of length k in left-to-right order. Then, for 1 <= j <= k, top[j-1]
	 * is set to x_j, bot[j-1] is set to y_j, and mid[j-1] is set to '|' if x_j ==
	 * y_j and ' ' (a single space) otherwise. Each of the remaining slots in each
	 * array is set to ' '. The alignment length variable alignlen is set to k. The
	 * variable editnum, which denotes the number of differences in the alignment,
	 * is set to the number of pairs (x_j, y_j) such that x_j != y_j.
	 */
	public Alignment(Sequence seqone, Sequence seqtwo, SubstitutionGap param) {

		// char array
		
		if(seqone == null ) {
			throw new RuntimeException("Sequence object one is null.");
		}
		else if(seqtwo == null) {
			throw new RuntimeException("Sequence object two is null.");
		}
		else if(param == null) {
			throw new RuntimeException("SubstitutionGap object is null.");
		}
		
		
		
		if(seqone.getLength() == 0 && seqtwo.getLength() == 0) {
			throw new RuntimeException("Both sequences are empty.");
		}
		
		
		
		A = seqone.getSeq();

		// length a
		M = seqone.getLength();

		// char array
		B = seqtwo.getSeq();

		// length b
		N = seqtwo.getLength();

		int r = param.getgap();

		creatingMatrix(r, param);

		editDistance = Smat[0][0];

		creatingOptimalAlignment(r);

	}
// creating the matrix to help calculate edit distance and applying traceback procedures to
	public void creatingMatrix(int r, SubstitutionGap param) {

		Smat = new int[M + 1][N + 1];

		Smat[M][N] = 0;

		for (int j = N - 1; j >= 0; j--) {
			Smat[M][j] = Smat[M][j + 1] + r;
		}
		for (int i = M - 1; i >= 0; i--) {
			Smat[i][N] = Smat[i + 1][N] + r;

			for (int j = N - 1; j >= 0; j--) {

				int sigma = param.getmatch();

				if (A[i + 1] != B[j + 1]) {
					sigma = param.getmismat();
				}

				Smat[i][j] = Math.min(Math.min(Smat[i + 1][j + 1] + sigma, Smat[i + 1][j] + r), Smat[i][j + 1] + r);

			}
		}

	}
//creating the optimal alignment - helps to calculate the number of differences and length of alignment 
	public void creatingOptimalAlignment(int r) {

		int i = 0;
		int j = 0;

		// OptimalAlignment
		ArrayList<Character> topList = new ArrayList<Character>();
		ArrayList<Character> midList = new ArrayList<Character>();
		ArrayList<Character> botList = new ArrayList<Character>();

		while (i <= M && j <= N) {

			if (i == M && j == N) {
				break;
			}

			if (j == N || i < M && Smat[i][j] == Smat[i + 1][j] + r) {
				topList.add(A[i + 1]);
				midList.add('-');
				botList.add(' ');
				editnum++;
				i++;
				continue;

			}
			if (i == M || Smat[i][j] == Smat[i][j + 1] + r) {
				topList.add(' ');
				midList.add('-');
				botList.add(B[j + 1]);
				editnum++;
				j++;
				continue;

			}

			if (A[i + 1] == B[j + 1]) {
				topList.add(A[i + 1]);
				midList.add('|');
				botList.add(B[j + 1]);
			}

			else {
				topList.add(A[i + 1]);
				midList.add(' ');
				botList.add(B[j + 1]);
				editnum++;
			}
			i++;
			j++;

		}
		alignlen = topList.size();

		top = new char[topList.size()];
		mid = new char[midList.size()];
		bot = new char[botList.size()];

		for (int x = 0; x < top.length; x++) {

			top[x] = topList.get(x);
			mid[x] = midList.get(x);
			bot[x] = botList.get(x);

		}

	}

	// TODO: This is the only method that you need to implement in this assignment.
	// Please follow the instructions given above this method.
	// You are required to implement the algorithm in the specification,
	// where S(m, n) is set to 0 and S(0, 0) is the edit distance;
	// You should not implement other algorithm variants where S(0, 0) is set
	// to 0 and S(m, n) is the edit distance. This requirement is enforced in this
	// template.

// You should not modify the methods below except the main() method.

// Returns a String object containing one of the two sequences and its length:
// sequence A followed by its length if kind == 'A', and
// sequence B followed by its length otherwise.
	public String getSequence(char kind) {
		int j, len;
		char[] seq;

		if (kind == 'A') {
			seq = A;
			len = M;
		} else {
			seq = B;
			len = N;
		}
		StringBuilder out = new StringBuilder();
		for (j = 1; j <= len; j++) {
			out.append(String.format("%c", seq[j]));
		}
		out.append(String.format(" of length %d", len));
		return out.toString();
	}

// Returns the edit distance between the two sequences.
	public int getEditdistance() {
		return editDistance;
	}

// Returns the length of the alignment of the two sequences.
	public int getAlignmentLength() {
		return alignlen;
	}

// Returns the number of differences between the two sequences.
	public int getNumberOfDifference() {
		return editnum;
	}

// Returns a String object containing one row of the optimal alignment:
// top row if kind == 'T', middle row if kind == 'M', and bottom row if kind == 'B'.
// Otherwise, it throws a RuntimeException exception.
	public String getAlignment(char kind) {
		if (kind == 'T')
			return new String(top, 0, alignlen);
		if (kind == 'M')
			return new String(mid, 0, alignlen);
		if (kind == 'B')
			return new String(bot, 0, alignlen);
		throw new RuntimeException("Alignment: Unexpected type character");
	}

// Returns the cost matrix.
	public int[][] getMatrix() {
		return Smat;
	}

// The function returns a String object containing
// each value in the cost matrix. The values are reported first in order of row index i 
// and then in order of column index j in the form of (i, j, Smat[i][j])
// so that the row and column coordinates of a value are easily recognized.
// This function is for your own use to check on each matrix,
// so any format is OK.
	public String matrixtoString() {
		int i, j;
		int[] row;

		if (Smat[M][N] != 0)
			throw new RuntimeException("matrixtoString: did not set S(m, n) to 0");

		StringBuilder res = new StringBuilder();

		res.append(String.format("Matrix S \n\n"));
		res.append(String.format(" Seq2     "));
		for (j = 0; j <= N; j++)
			if (j > 0) {
				res.append(String.format(" %c   ", B[j]));
			} else {
				res.append(String.format("       "));
			}

		res.append(String.format("\n"));

		{
			res.append(String.format("  column  "));
		}
		for (j = 0; j <= N; j++) {
			res.append(String.format(" %3d ", j));
		}

		res.append(String.format("\n"));
		res.append(String.format("Seq1\n"));
		for (i = 0; i <= M; i++) {
			row = Smat[i];
			if (i > 0) {
				res.append(String.format("%c row %2d  ", A[i], i));
			} else {
				res.append(String.format("  row %2d  ", i));
			}
			for (j = 0; j <= N; j++) {
				res.append(String.format("%5d", row[j]));
			}

			res.append(String.format("\n\n"));
		}
		return res.toString();
	} // matrixtoString()

// A main() method for illustrating a small example.
	public static void main(String[] args) {

		SubstitutionGap param = new SubstitutionGap();
		Sequence seqone = new Sequence("A", "ACGTCGAGCTA");
		Sequence seqtwo = new Sequence("B", "ACCTCGACTA");
		Alignment align = new Alignment(seqone, seqtwo, param);
		System.out.println("Match score: " + param.getmatch());
		System.out.println("Mismatch penalty: " + param.getmismat());
		System.out.println("Gap penalty: " + param.getgap());
		System.out.println("Sequence A: " + align.getSequence('A'));
		System.out.println("Sequence B: " + align.getSequence('B'));

		// how do you figure out the edit distance (today)
		//
		System.out.println("Edit distance: " + align.getEditdistance());

		// how do you figure out the number of differences (today)
		// figure out how to calculate the number of differneces
		System.out.println("Number of differences: " + align.getNumberOfDifference());

		// The optimal alignment - traceback (wednesday)
		System.out.println("Length of alignment: " + align.getAlignmentLength());
		System.out.println("Top row: " + align.getAlignment('T'));
		System.out.println("Mid row: " + align.getAlignment('M'));
		System.out.println("Bot row: " + align.getAlignment('B'));
		System.out.println();

		// The matrix (today)
		// There's the pseudocode, use that as a reference
		// I still am not sure how the matrix is formed
		System.out.println(align.matrixtoString());

	} // main

	/*
	 * Below is the output from main() on the two short DNA sequences.
	 * 
	 * Match score: 0 Mismatch penalty: 20 Gap penalty: 15 Sequence A: ACGTCGAGCTA
	 * of length 11 Sequence B: ACCTCGACTA of length 10 Edit distance: 35 Number of
	 * differences: 2 Length of alignment: 11 Top row: ACGTCGAGCTA Mid row: ||
	 * ||||-||| Bot row: ACCTCGA CTA
	 * 
	 * Matrix S
	 * 
	 * Seq2 A C C T C G A C T A column 0 1 2 3 4 5 6 7 8 9 10 Seq1 row 0 35 50 45 60
	 * 75 90 105 120 135 150 165
	 * 
	 * A row 1 50 35 30 45 60 75 90 105 120 135 150
	 * 
	 * C row 2 65 50 35 30 45 60 75 90 105 120 135
	 * 
	 * G row 3 60 45 30 15 30 45 60 75 90 105 120
	 * 
	 * T row 4 75 60 45 30 15 30 45 60 75 90 105
	 * 
	 * C row 5 90 75 60 45 30 15 30 45 60 75 90
	 * 
	 * G row 6 75 80 65 50 35 30 15 30 45 60 75
	 * 
	 * A row 7 90 75 60 45 30 15 20 15 30 45 60
	 * 
	 * G row 8 105 90 75 60 45 30 15 0 15 30 45
	 * 
	 * C row 9 120 105 90 75 60 45 30 15 0 15 30
	 * 
	 * T row 10 135 120 105 90 75 60 45 30 15 0 15
	 * 
	 * A row 11 150 135 120 105 90 75 60 45 30 15 0
	 * 
	 */

} // Alignment

// A class for alignment parameters
class SubstitutionGap {
	int match = 0; // a bonus score for a pair of identical DNA letters
	int mismat = 20; // a penalty score for a pair of different DNA letters
	int gapext = 15; // a penalty score for each letter in a gap

	public void setmatch(int mat) {
		match = mat;
	}

	public void setmismat(int mis) {
		mismat = mis;
	}

	public void setgap(int gap) {
		gapext = gap;
	}

	public int getmatch() {
		return match;
	}

	public int getmismat() {
		return mismat;
	}

	public int getgap() {
		return gapext;
	}
} // SubstitutionGap

// A class for an object holding a sequence
class Sequence {
	private String name; // the name of a sequence
	private int length; // the length of the sequence
	private char[] seq; // an array holding the sequence from index 1 to index length

	public Sequence(String aname, String aseq) {
		if (aname == null)
			throw new NullPointerException("No sequence name");
		if (aseq == null)
			throw new NullPointerException("No sequence");
		name = new String(aname);
		length = aseq.length();
		seq = new char[length + 1];
		seq[0] = ' ';
		for (int j = 1; j <= length; j++)
			seq[j] = aseq.charAt(j - 1);
	}

	public String getName() {
		return name;
	}

	public int getLength() // returns length.
	{
		return length;
	}

	public char[] getSeq() // returns seq.
	{
		return seq;
	}
} // Sequence