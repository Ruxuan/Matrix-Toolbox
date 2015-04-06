// Matrix Toolbox by Ru Li

import java.lang.Math;

public class matrix {
	double[][] data;
	double[] b;
	int rows = 0;
	int cols = 0;
	double[][] rref;
	double[] rref_b;
	int rank;
	double det = 1;
	double det_d = 1;
	double[][] transpose;
	double[][] inverse;
	double[] soln;
	
	// Constructor
	public matrix(double[][] array, double[] b, int r, int c){
		this.rows = r;
		this.cols = c;
		this.data = array;
		this.rref = new double[r][c];
		for (int i = 0; i < this.rows ; ++i) {
			for (int j = 0; j < this.cols ; ++j) {
				this.rref[i][j] = array[i][j];
			}
		}
		this.b    = b;
		this.rref_b = new double[rows];
		for (int i = 0; i < this.rows ; ++i) {
			this.rref_b[i] = b[i];
		}
	}
	
	// Display methods **************************************
	// prints matrix
	public void matrix_display() {			
		System.out.println("---Data------------");
		for (int i = 0; i < this.rows ; ++i) {
			System.out.print("| ");
			for (int j = 0; j < this.cols ; ++j) {
				System.out.printf("%4.2f ", this.data[i][j]);
			}
			System.out.printf(" | %6.2f |", this.b[i]);
			System.out.println();
		}
		System.out.println("-------------------");			
	}
	
	//printf this.rref
	public void rref_display() {
		//this.reduced_row_echelon_form()
		System.out.println("--RREF-------------");
		for (int i = 0; i < this.rows ; ++i) {
			System.out.print("| ");
			for (int j = 0; j < this.cols ; ++j) {
				System.out.printf("%4.2f ", this.rref[i][j]);
			}
			System.out.printf(" | %6.2f |", this.rref_b[i]);
			System.out.println();
		}
		System.out.println("-------------------");	
	}
	
	// prints this.soln
	public void soln_display() {
		System.out.println("--SOLN-------------");
		for (int i = 0; i < this.cols; ++i) {
			System.out.printf("x%d = %6.4f\n", i+1, this.soln[i]);
		}
		System.out.println("-------------------");
	}
	
	// prints the determinant
	public void det_display() {
		System.out.println("--DET--------------");
		System.out.printf("%6.4f\n", this.det);
		System.out.println("-------------------");
	}
	
	// prints the transpose
	public void transpose_display() {
		if (this.rows != this.cols) return;
		System.out.println("--TRANSPOSE--------");
		for (int i = 0; i < this.rows ; ++i) {
			System.out.print("| ");
			for (int j = 0; j < this.cols ; ++j) {
				System.out.printf("%4.2f ", this.transpose[i][j]);
			}
			System.out.println("|");
		}
		System.out.println("-------------------");	
	}
	
	// prints the inverse
	public void inverse_display() {
		if (this.rows != this.cols) return;
		System.out.println("--INV--------------");
		for (int i = 0; i < this.rows ; ++i) {
			System.out.print("| ");
			for (int j = 0; j < this.cols ; ++j) {
				System.out.printf("%4.2f ", this.inverse[i][j]);
			}
			System.out.println("|");
		}
		System.out.println("-------------------");
	}
	
	// prints the rank
	public void rank_display() {
		System.out.println("--RANK-------------");
		System.out.printf("Rank = %d\n", this.rank);
		System.out.println("-------------------");
	}
	
	//***********************************************************
	// Matrix Operations ****************************************
	// this * other
	public double[][] matrix_multiplication(matrix other) {
		double[][] product = new double[this.rows][other.cols];
		if (this.cols != other.rows) {
			System.out.println("Unable to perform multiplication of matrices");
			return product;
		}
		int n_r = this.rows;
		int n_c = other.cols;
		
		for (int i = 0; i < n_r; ++i) {
			for (int j = 0; j < n_c; ++j) {
				int dot = 0;
				for (int k = 0; k < this.cols; ++k) {
					dot += other.data[k][j] * this.data[i][k];
				}
				product[i][j] = dot;
			}
		}
		return product;
	}
	
	// this + other
	public double[][] matrix_addition(matrix other) {
		double[][] sum = new double[this.rows][this.cols];
		if (this.cols != other.cols || this.rows != other.cols) {
			System.out.println("Unable to add matrices");
			return sum;
		}
		
		int n_r = this.rows;
		int n_c = this.cols;
		
		for (int i = 0; i < n_r; ++i) {
			for (int j = 0; j < n_c; ++j) {
				sum[i][j] = this.data[i][j] + other.data[i][j];
			}
		}
		return sum;
	}
	
	// scales this by c
	public void ro_scale(double[][] m, double[] b, int r, double c) {
		for (int i = 0; i < this.cols ; ++i) {
			m[r][i] *= c;
		}
		b[r] *= c;
	}
			
	// switches the positions of r1 and r2
	public void ro_switch(double[][] m, double[] b, int r1, int r2) {
		double temp[] = m[r1];
		m[r1] = m[r2];
		m[r2] = temp;
		
		double tmp = b[r1];
		b[r1] = this.b[r2];
		b[r2] = tmp;
	}
	
	// Adds/subtracts r1 by r2 scaled by c
	public void ro_add(double[][] m, double[] b, int r1, int r2, double c) {
		for (int i = 0; i < this.cols ; ++i) {
			m[r1][i] += c * m[r2][i];
		}
		b[r1] += b[r2]*c;
	}
	
	// Row Reducing Algorithms ****************************************
	// row echelon form (Gauss Jordan Elimination)
	public void gauss_jordan_elimination() {
		int min = Math.min(this.rows, this.cols);
		
		for (int k = 0; k < min; ++k) {
			// Find the k-th pivot
			int i_max = k;
			for (int i = k+1; i < this.rows; ++i) {
	            if (Math.abs(this.rref[i][k]) > Math.abs(this.rref[i_max][k])) {
	                i_max = i;
	            }
	        }
			
			// Error detection
			if (this.rref[i_max][k] == 0) {
				System.out.println("Error: Matrix is Singular or Inconsistent!");
				break;
			}
			
			if (k != i_max) this.det_d *= -1;
			this.ro_switch(this.rref, this.rref_b, k, i_max);
			
			// Do for all rows below pivot
			for (int i = k + 1; i < this.rows ; ++i) {
				double temp = this.rref[i][k] / this.rref[k][k];
				// Do for a all remaining elements in current row
				for (int j = k + 1; j < this.cols ; ++j) {
					this.rref[i][j] -= this.rref[k][j] * temp;
				}
				// Fill lower triangular matrix with zeros
				this.rref[i][k] = 0;
				this.rref_b[i] -= temp * rref_b[k];
			}
		}			
	}
	
	// Row Reduced Echelon form
	public void reduced_row_echelon_form() {
		int lead = 0;
		int rowCount = this.rows;
		int colCount = this.cols;
		for (int r = 0; r < rowCount; ++r) {
			if (lead >= colCount) return;
			int i = r;
			while (this.rref[i][lead] == 0){
				++i;
				if (i == rowCount) {
					i = r;
					++lead;
					if (colCount == lead) {
						return;
					}
				}
			}
			this.ro_switch(this.rref, this.rref_b, i,r);	
			if (this.rref[r][lead] != 0) {
				this.ro_scale(this.rref, this.rref_b, r, 1/this.rref[r][lead]);
			}
			
			for (int j = 0; j < rowCount; ++j) {
				if (j == r) continue;
				double val = -this.rref[j][lead];
				this.ro_add(this.rref, this.rref_b, j, r, val);
			}
			++lead;
		}
	}
	// Data Computation ***************************************
	// Computes all data
	public void find_soln() {
		this.gauss_jordan_elimination();
		this.back_substitution();
		this.determinant();
		this.rank();
		//the functions above requires row-echelon-form first
		this.inverse();
		this.transpose();
		this.reduced_row_echelon_form();
	}		
	
	// Use back substitution after gauss-jordan elimination
	public void back_substitution() {
		// Back substitution to find solutions
		int min = Math.min(this.rows, this.cols);
		this.soln = new double[min];
		
		for (int a = min - 1; a >= 0; --a) {
			double temp = this.rref_b[a];
			for (int c = min - 1; c > a; --c) {
				temp -= this.soln[c] * this.rref[a][c];
			}
			this.soln[a] = temp / this.rref[a][a];
		}
	}
	
	// Finds determinant of matrix
	public void determinant() {
		if (this.rows != this.cols) {
			this.det = 0;
		} else {
			for (int i = 0; i < this.rows; ++i) {
				this.det *= this.rref[i][i];
			}
		}
		this.det *= this.det_d;
	}
	
	// Finds transpose of matrix
	public void transpose() {
		if (this.rows != this.cols) {
			this.transpose = new double[0][0];
		} else {
			this.transpose = new double[rows][cols];
			for (int i = 0; i < this.rows ; ++i) {
				for (int j = 0; j < this.cols ; ++j) {
					this.transpose[i][j] = this.data[j][i];
				}
			}
		}
	}
	
	// Finds rank of a matrix
	public void rank() {
		this.rank = 0;
		for (int i = 0; i < this.rows; ++i) {
			int rec = 0;
			for (int j = 0; j < this.cols; ++j) {
				if (this.rref[i][j] != 0) {
					rec = 1;
				}
			}
			if (rec == 1) ++this.rank;
		}
	}
	
	// Finds inverse of a matrix
	public void inverse() {
		if (this.det == 0) return;
		else this.inverse = new double[rows][cols];
		
		int rowCount = this.rows;
		int colCount = this.cols;
		
		for (int inv = 0 ; inv < rowCount; ++inv) {
			double[][] data_copy = new double[rowCount][colCount];
			for (int i = 0; i < rowCount ; ++i) {
				for (int j = 0; j < colCount ; ++j) {
					data_copy[i][j] = this.data[i][j];
				}
			}
			
			double[] colInv = new double [rowCount];
			colInv[inv] = 1;
			
			int lead = 0;
			for (int r = 0; r < rowCount; ++r) {
				if (lead >= colCount) return;
				int i = r;
				while (data_copy[i][lead] == 0){
					++i;
					if (i == rowCount) {
						i = r;
						++lead;
						if (colCount == lead) {
							return;
						}
					}
				}
				this.ro_switch(data_copy, colInv, i, r);
				if (data_copy[r][lead] != 0) {
					this.ro_scale(data_copy, colInv, r, 1/data_copy[r][lead]);
				}
				for (int j = 0; j < rowCount; ++j) {
					if (j == r) continue;
					double val = -data_copy[j][lead];
					this.ro_add(data_copy, colInv, j, r, val);
				}
				++lead;
			}
			for (int i = 0; i < rowCount; ++i) {
				this.inverse[i][inv] = colInv[i];
			}
		}
	}
	//__________________________________________________________
	
	public static void main(String[] args) {
		int input_rows = 3;
		int input_cols = 3;
		double [][] input_data = new double[][] {{1,2,3},
                					             {0,1,4},
                					             {5,6,0}};
		double [] input_aug = new double[] {14, 14, 17};
		
		matrix standard = new matrix(input_data, input_aug, input_rows, input_cols);
		standard.matrix_display();
		standard.find_soln();
		standard.rref_display();
		standard.soln_display();
		standard.det_display();
		standard.transpose_display();
		standard.rank_display();
		standard.inverse_display();
	}
}