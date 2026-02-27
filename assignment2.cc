/**
 * @file assignment1.cc
 * @brief Skeleton of code for 5614 Assignment 2.
 * 
 * 5614 Assignment 2. This file contains the bones of the class
 * definition plus the main function. The main function that you
 * finish with should be the same as provided.
 * Please add the rest of the code as instructed in the document. 
 * Don't forget to use comments.
 * 	
 * @author R. Morrin
 * @version 7.0
 * @date 2026-02-15
 */
#include <iostream> 			// Needed for I/O
#include <cmath> 			// Needed for maths functions
#include <iomanip> 			// Output formatting
#include <print>
#include <fstream>


constexpr double MYPI {4 * std::atan(1.0)};
//Use the below if not using g++ and you get errors with line above
//constexpr double MYPI = M_PI;

/**
 * @class Gaussian
 * @brief A class to represent a normal distribution
 *
 * This class represents a Gaussian distribution.
 */
class Gaussian
{
    public:
	Gaussian (); 				// Default Constructor    // WRITE DEFINITION INLINE
	Gaussian (double const mean, double const stdev); 	// Overloaded constructor // WRITE DEFINITION INLINE
	double get_mu() const; 						// WRITE DEFINITON INLINE !!
	double get_sigma() const; 					// WRITE DEFINITON INLINE !!
	double normalised(double const x) const; 				// return normalized z value  // WRITE DEFINITION INLINE

	// Write remaining member function definitions outside of class body
	Gaussian(const Gaussian& rhs); 					// Copy constructor
	Gaussian & operator=(const Gaussian &rhs); 			// Copy assignment operator
	double pdf(double const x) const; 					// return pdf at x
	double cdf(double const x) const; 				// return cdf at x using std::erfc
    double inverse_cdf(double z) const;             // inverse cdf value at x
	void print_parameters() const; 					// Print distribution parameters

	// You don't need to modify below here in class definition.
	~Gaussian() { 							// Destructor
	    std::println("Destroying object with mu = {} stdev = {}", mu, sigma);
	}

    private:
	double mu; 							///< Holds mean of Gaussian
	double sigma; 							///< Holds standard deviation of Gaussian
};



/// TODO
// Put member function defintions below here
//


/// TODO
// Non-Member Functions
// These are just another way to get same functionality.
// void print_parameters(const Gaussian & dist);
// double pdf(const Gaussian & dist, double x);
// double cdf(const Gaussian & dist, const double x);


/// TODO Write the function defintion for this function below main
void generate_file_of_rngs(const Gaussian& dist, std::string filename, int seed, int N);

// Main function for 5614 Assigment 2

int main()
{
    Gaussian A;  			// Create object of class Gaussian via default constructor
    Gaussian B {1, 2}; 			// Create object of class Gaussian via overloaded constructor

    // Create a list. Is actually ok to use auto here
    auto list = {-8., -3.0, -2.0, -1.0, -0.5, 0.0, 0.5, 1.0, 2.0, 3.0};

    // Calculate and print cumulative density function values for A for a number of input values.
    // Most of these lines are for formatting output. We will cover this later.
    A.print_parameters();
    std::println("{:-<60}", "");
    std::println("{0: >2}x{0:>11}Phi(x)", "");

    for (auto  i : list) {
        std::println("{0: >3.1f}\t{1: >18.16f}", i, A.cdf(i));
    }
    std::println("{:-<60}\n", "");


    // Calculate and print cumulative density function values for B for a number of input values.
    B.print_parameters();
    std::println("{:-<60}", "");
    std::println("{0: >2}x{0:>11}Phi(x)", "");
    for (auto  i : list) {
        std::println("{0: >3.1f}\t{1: >18.16f}", i, B.cdf(i));
    }
    std::println("{:-<60}", "");

    //
    // Use the non-member functions
    //

    Gaussian C {2,5}; 			// Create another object of class Gaussian
    print_parameters(C); 		// Print the parameters of D using non-member function

    // Print parameters of A, followed by an example CDF calculation using the free functions.
    std::println("Checking Free functions");
    std::println("PDF of A at x=1 is {:.12f}", pdf(A,1.0));
    std::println("CDF of A at x=2 is {:.12f}\n", cdf(A,2.0));

    // Copy-and-swap
    Gaussian temp {A};
    A = B; 				// Should call copy assignment operator you wrote
    B = temp;

    print_parameters(A);
    std::println("Inverse CDF for A with Phi=0.4  is {:.6f}"  , A.inverse_cdf(0.4));
    std::println("Inverse CDF for A with Phi=0.5  is {:.6f}\n", inverse_cdf(A, 0.5));
    print_parameters(B);
    std::println("Inverse CDF for B with Phi=0.25 is {:.6f}"  , B.inverse_cdf(0.25));
    std::println("Inverse CDF for B with Phi=0.68 is {:.6f}\n", inverse_cdf(B, 0.68));


    // Generate file of Gaussian distributed RNGs
    std::println("Generating file of 1'000'000 normally distributed numbers\n");
    generate_file_of_rngs(A, "output.txt", 12345, 1'000'000);

    return 0;
}

