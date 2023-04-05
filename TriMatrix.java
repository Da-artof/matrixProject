/**
 * PROJECT III: TriMatrix.java
 *
 * This file contains a template for the class TriMatrix. Not all methods are
 * implemented. Make sure you have carefully read the project formulation
 * before starting to work on this file. You will also need to have completed
 * the Matrix class.
 *
 * Remember not to change the names, parameters or return types of any
 * variables in this file!
 *
 * The function of the methods and instance variables are outlined in the
 * comments directly above them.
 */
import java.util.Random;

public class TriMatrix extends Matrix {
    /**
     * An array holding the diagonal elements of the matrix.
     */
    private double[] diag;

    /**
     * An array holding the upper-diagonal elements of the matrix.
     */
    private double[] upper;

    /**
     * An array holding the lower-diagonal elements of the matrix.
     */
    private double[] lower;
    
    /**
     * Constructor function: should initialise m and n through the Matrix
     * constructor and set up the data array.
     *
     * @param N  The dimension of the array.
     */
    public TriMatrix(int N) {
	super(N,N);
	if (N > 0) {
		diag = new double[N];
		lower = new double[N-1];
		upper = new double[N-1];
	} else {
		throw new MatrixException("Matrix dimension must be greater than 0");
	}
    }
    
    /**
     * Getter function: return the (i,j)'th entry of the matrix.
     *
     * @param i  The location in the first co-ordinate.
     * @param j  The location in the second co-ordinate.
     * @return   The (i,j)'th entry of the matrix.
     */
    public double getIJ(int i, int j) {
        if (i < 0 || j < 0) {
		throw new MatrixException("Matrix indices must be non-negative.");
	} else if (i >= this.m || j >= this.n) {
                throw new MatrixException("Matrix indices must be less than matrix dimensions.");
        } else if (i == j) {
		return diag[i];
        } else if (i == j - 1) {
		return upper[i];
        } else if (i == j + 1) {
        	return lower[j];
	} else {
		return 0;
	}
    }
    
    /**
     * Setter function: set the (i,j)'th entry of the data array.
     *
     * @param i    The location in the first co-ordinate.
     * @param j    The location in the second co-ordinate.
     * @param val  The value to set the (i,j)'th entry to.
     */
    public void setIJ(int i, int j, double val) {
        if (i == j) {
                diag[i] = val;
        } else if (i == j - 1) {
                upper[i] = val;
        } else if (i == j + 1) {
                lower[j] = val;
        } else if (val == 0) {
        } else {
                throw new MatrixException("Tridiagonal Matrix entries must be 0 if not on any of the 3 central diagonals");
        }
    }
    
    /**
     * Return the determinant of this matrix.
     *
     * @return The determinant of the matrix.
     */
    public double determinant() {
	TriMatrix A = this.decomp();
	double det = A.diag[0];
	for (int i = 1; i < this.m; i++){
		det = det * A.getIJ(i,i);
	}
	return det;
    }
    
    /**
     * Returns the LU decomposition of this matrix. See the formulation for a
     * more detailed description.
     * 
     * @return The LU decomposition of this matrix.
     */
    public TriMatrix decomp() {
	TriMatrix A = new TriMatrix(this.m);
	for (int i = 0; i < this.m - 1; i++) {
		A.upper[i] = this.upper[i];
	}	
	A.diag[0] = this.diag[0];
	A.lower[0] = this.lower[0] / this.diag[0];
	for (int i = 1; i < this.m; i++) {
		if (i < this.m - 1){
			A.diag[i] = this.diag[i] - (A.lower[i-1]*A.upper[i-1]);
			A.lower[i] = this.lower[i] / A.diag[i];
		} else {
			A.diag[i] = this.diag[i] - (A.lower[i-1]*A.upper[i-1]);
		}
	}
	return A;
    }

    /**
     * Add the matrix to another matrix A.
     *
     * @param A  The Matrix to add to this matrix.
     * @return   The sum of this matrix with the matrix A.
     */
    public Matrix add(Matrix A){
	double sum;
	if (A.m == A.n && A.m == this.m) {	
		if (A instanceof TriMatrix) {
			TriMatrix B = new TriMatrix(this.m);
			for (int i = 0; i < this.m; i++) {
				for (int j = 0; j < this.m; j++) {
					if (i == j){
						sum = A.getIJ(i,j) + this.getIJ(i,j);
						B.setIJ(i,j,sum);
					} else if (i == j - 1){
                                                sum = A.getIJ(i,j) + this.getIJ(i,j);
                                                B.setIJ(i,j,sum);
                                        } else if (i == j + 1){
                                                sum = A.getIJ(i,j) + this.getIJ(i,j);
                                                B.setIJ(i,j,sum);
                                        }
	        	        }
			}
			return B;
		} else {
	                return A.add(this);
		}
	} else {
		throw new MatrixException("Matrix dimensions must be equal for addition");
	}	
    }
    
    /**
     * Multiply the matrix by another matrix A. This is a _left_ product,
     * i.e. if this matrix is called B then it calculates the product BA.
     *
     * @param A  The Matrix to multiply by.
     * @return   The product of this matrix with the matrix A.
     */
    public Matrix multiply(Matrix A) {
        if (this.n == A.m) {
                GeneralMatrix C = new GeneralMatrix(this.m,A.n);
                for (int i = 0; i < A.m; i++) {
                        for (int j = 0; j < A.n; j++) {
	               	        double sum = 0;
             	            	for (int k = 0; k < this.n; k++) {
                      	                double prod = this.getIJ(i,k) * A.getIJ(k,j);
                       		        sum = sum + prod;
                       	        }
                                C.setIJ(i,j,sum);
                	}
		}
                return C;
        } else {
                throw new MatrixException("No. of columns of left matrix must equal no. of rows of right matrix.");
        }
    }
    
    /**
     * Multiply the matrix by a scalar.
     *
     * @param a  The scalar to multiply the matrix by.
     * @return   The product of this matrix with the scalar a.
     */
    public Matrix multiply(double a) {
	TriMatrix A = new TriMatrix(this.n);
	for (int i = 0; i < this.n; i++) {
		for (int j = 0; j < this.n; j++) {
			A.setIJ(i,j,this.getIJ(i,j)*a);
        	}
	}
	return A;
    }

    /**
     * Populates the matrix with random numbers which are uniformly
     * distributed between 0 and 1.
     */
    public void random() {
        for (int i = 0; i < this.m; i++) {
                for (int j = 0; j < this.n; j++) {
			if (i == j) {	
                                Random r = new Random();
                                double R = r.nextDouble();
                                this.setIJ(i,j,R);
			} else if (j == i + 1) {
                                Random r = new Random();
                                double R = r.nextDouble();
                                this.setIJ(i,j,R);
                        } else if (j == i - 1) {
                                Random r = new Random();
                                double R = r.nextDouble();
                                this.setIJ(i,j,R);
                        } else {
				this.setIJ(i,j,0);
			}
                }
        }
    }
    
    /*
     * Your tester function should go here.
     */
    public static void main(String[] args) {
	TriMatrix A = new TriMatrix(3);
	TriMatrix B = new TriMatrix(3);
	GeneralMatrix C = new GeneralMatrix(3,3);
	A.random();
	B.random();
	C.random();
	System.out.println("TriMat A = " + A.toString());
	System.out.println("TriMat B = " + B.toString());
	System.out.println("GenMat C = " + C.toString());
	System.out.println("Mat A (0,0) = " + A.getIJ(0,0));	
	A.setIJ(0,0,1);	
	System.out.println(A.getIJ(0,0));
	System.out.println("Det(A) = " + A.determinant());
	System.out.println("A + B = " + A.add(B).toString());
	System.out.println("A + C = " + A.add(C).toString());
	System.out.println("A * B = " + A.multiply(B).toString());
	System.out.println("A * C = " + A.multiply(C).toString());
	System.out.println("2 * A = " + A.multiply(2).toString());	
    }
}
