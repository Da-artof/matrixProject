/*
 * PROJECT III: GeneralMatrix.java
 *
 * This file contains a template for the class GeneralMatrix. Not all methods
 * are implemented. Make sure you have carefully read the project formulation
 * before starting to work on this file. You will also need to have completed
 * the Matrix class.
 *
 * Remember not to change the names, parameters or return types of any
 * variables in this file!
 *
 * The function of the methods and instance variables are outlined in the
 * comments directly above them.
 */

import java.util.Arrays;
import java.util.Random;
public class GeneralMatrix extends Matrix {
    /**
     * This instance variable stores the elements of the matrix.
     */
    private double[][] data;

    /**
     * Constructor function: should initialise m and n through the Matrix
     * constructor and set up the data array.
     *
     * @param m  The first dimension of the array.
     * @param n  The second dimension of the array.
     */
    public GeneralMatrix(int m, int n) throws MatrixException {
	super(m,n);
	if (m < 1 || n < 1) {
		throw new MatrixException("Matrix dimension must be greater than 0.");
	}else {
		data = new double[m][n];
	}
    }

    /**
     * Constructor function. This is a copy constructor; it should create a
     * copy of the matrix A.
     *
     * @param A  The matrix to create a copy of.
     */
    public GeneralMatrix(GeneralMatrix A) {
	super(A.m,A.n);
	data = new double[A.m][A.n];
	int i = 0, j = 0;
	while(i < A.m) {
		while (j < A.n) {
			data[i][j] = A.getIJ(i,j);
			j++;
		}
		j = 0;
		i++;
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
	if (i < 0 || j < 0){
		throw new MatrixException("Matrix indices must be 0 or greater");
	} else if (i >= m || j >= n) {
		throw new MatrixException("Matrix indices must be less than no. of rows/columns");
	} else {
		return this.data[i][j];
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
	if (i < 0 || j < 0){
                throw new MatrixException("Matrix indices must be 0 or greater");
	} else if (i >= m || j >= n) {
                throw new MatrixException("Matrix indices must be less than no. of rows/columns");
        } else {
		this.data[i][j] = val;
        }
    }
    
    /**
     * Return the determinant of this matrix.
     *
     * @return The determinant of the matrix.
     */
    public double determinant() {
	double[] d = new double[1];
	GeneralMatrix LU = this.decomp(d);
	double det = LU.getIJ(0,0);
	for(int i = 1; i < this.m; i++) {
		det = det*LU.getIJ(i,i);
	}
	det = det * d[0];
	return det;
    }

    /**
     * Add the matrix to another matrix A.
     *
     * @param A  The Matrix to add to this matrix.
     * @return   The sum of this matrix with the matrix A.
     */
    public Matrix add(Matrix A) {
	if ((A.m == this.m) && (A.n == this.n)) {
		GeneralMatrix B = new GeneralMatrix(A.m,A.n);
		for (int i = 0; i < m; i++){
			for (int j = 0; j < n; j++){
				double sum = this.getIJ(i,j) + A.getIJ(i,j);
				B.setIJ(i,j,sum);
			}
		}
	return B;
	}else {
		throw new MatrixException("Dimensions of matrices must be equal for addition.");
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
	}else   {
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
        GeneralMatrix C = new GeneralMatrix(this.m,this.n);
	for (int i = 0; i < this.m; i++) {
                for (int j = 0; j < this.n; j++) {
                 	double prod = this.getIJ(i,j) * a;
      	             	C.setIJ(i,j,prod);
		}
	}
	return C;
    }


    /**
     * Populates the matrix with random numbers which are uniformly
     * distributed between 0 and 1.
     */
    public void random() {
        for (int i = 0; i < this.m; i++) {
                for (int j = 0; j < this.n; j++) {
			Random r = new Random();
			double R = r.nextDouble();
			this.setIJ(i,j,R);
		}
	}
    }

    /**
     * Returns the LU decomposition of this matrix; i.e. two matrices L and U
     * so that A = LU, where L is lower-diagonal and U is upper-diagonal.
     * 
     * On exit, decomp returns the two matrices in a single matrix by packing
     * both matrices as follows:
     *
     * [ u_11 u_12 u_13 u_14 ]
     * [ l_21 u_22 u_23 u_24 ]
     * [ l_31 l_32 u_33 u_34 ]
     * [ l_41 l_42 l_43 l_44 ]
     *
     * where u_ij are the elements of U and l_ij are the elements of l. When
     * calculating the determinant you will need to multiply by the value of
     * d[0] calculated by the function.
     * 
     * If the matrix is singular, then the routine throws a MatrixException.
     *
     * This method is an adaptation of the one found in the book "Numerical
     * Recipies in C" (see online for more details).
     * 
     * @param d  An array of length 1. On exit, the value contained in here
     *           will either be 1 or -1, which you can use to calculate the
     *           correct sign on the determinant.
     * @return   The LU decomposition of the matrix.
     */
    public GeneralMatrix decomp(double[] d) {
        // This method is complete. You should not even attempt to change it!!
        if (n != m)
            throw new MatrixException("Matrix is not square");
        if (d.length != 1)
            throw new MatrixException("d should be of length 1");
        
        int           i, imax = -10, j, k; 
        double        big, dum, sum, temp;
        double[]      vv   = new double[n];
        GeneralMatrix a    = new GeneralMatrix(this);
        
        d[0] = 1.0;
        
        for (i = 1; i <= n; i++) {
            big = 0.0;
            for (j = 1; j <= n; j++)
                if ((temp = Math.abs(a.data[i-1][j-1])) > big)
                    big = temp;
            if (big == 0.0)
                throw new MatrixException("Matrix is singular");
            vv[i-1] = 1.0/big;
        }
        
        for (j = 1; j <= n; j++) {
            for (i = 1; i < j; i++) {
                sum = a.data[i-1][j-1];
                for (k = 1; k < i; k++)
                    sum -= a.data[i-1][k-1]*a.data[k-1][j-1];
                a.data[i-1][j-1] = sum;
            }
            big = 0.0;
            for (i = j; i <= n; i++) {
                sum = a.data[i-1][j-1];
                for (k = 1; k < j; k++)
                    sum -= a.data[i-1][k-1]*a.data[k-1][j-1];
                a.data[i-1][j-1] = sum;
                if ((dum = vv[i-1]*Math.abs(sum)) >= big) {
                    big  = dum;
                    imax = i;
                }
            }
            if (j != imax) {
                for (k = 1; k <= n; k++) {
                    dum = a.data[imax-1][k-1];
                    a.data[imax-1][k-1] = a.data[j-1][k-1];
                    a.data[j-1][k-1] = dum;
                }
                d[0] = -d[0];
                vv[imax-1] = vv[j-1];
            }
            if (a.data[j-1][j-1] == 0.0)
                a.data[j-1][j-1] = 1.0e-20;
            if (j != n) {
                dum = 1.0/a.data[j-1][j-1];
                for (i = j+1; i <= n; i++)
                    a.data[i-1][j-1] *= dum;
            }
        }
        
        return a;
    }

    /*
     * Your tester function should go here.
     */
    public static void main(String[] args) {
	System.out.println("Expected: same 2*2 matrix twice, different 2*2 matrix, 2*3 matrix");
	GeneralMatrix B = new GeneralMatrix(2,2);
	GeneralMatrix D = new GeneralMatrix(2,3);
	GeneralMatrix E = new GeneralMatrix(2,2);	
	B.random();
	D.random();
	E.random();
	GeneralMatrix C = new GeneralMatrix(B);
	System.out.println(B.toString());
	System.out.println(C.toString());
	System.out.println(E.toString());
	System.out.println(D.toString());
	System.out.println("Expected: top left entry of first matrix");
	System.out.println(B.getIJ(0,0));
	System.out.println("Expected:  2*2 matrix with 2 in bottom right corner");
	B.setIJ(1,1,2);
	System.out.println(B.toString());
	System.out.println("Expected: det of first matrix");
	System.out.println(B.determinant());
	System.out.println("Expected: sum of 2*2 matrices, 2*3 matrix, 2*2 matrix");
	System.out.println(B.add(E).toString());
	System.out.println(B.multiply(D).toString());
	System.out.println(B.multiply(E).toString());
	System.out.println("Expected: 2* first 2*2 matrix");
	System.out.println(B.multiply(2).toString());		
    }
}
