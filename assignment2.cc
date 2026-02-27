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
	Gaussian (): mu {0}, sigma{1} { 				// Default Constructor    // WRITE DEFINITION INLINE
        std::print("Constructing default with mean 0.0 and stdev 1.0\n" );
    }
	Gaussian (double const mean, double const stdev): mu{mean}, sigma{stdev}{ 	// Overloaded constructor // WRITE DEFINITION INLINE
        std::print("Constructing with mean {} and stdev {}", mu, sigma);
    }
	double get_mu() const { 						// WRITE DEFINITON INLINE !!
        return mu;
    }
	double get_sigma() const {					// WRITE DEFINITON INLINE !!
        return sigma;
    }
	double normalised(double const x) const{				// return normalized z value  // WRITE DEFINITION INLINE
        return (x - mu) / sigma;
    }

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

Gaussian::Gaussian(const Gaussian& rhs) : mu{rhs.mu}, sigma{rhs.sigma} {
    std::print("Gaussian Copy Constructor called");
}

Gaussian & Gaussian::operator=(const Gaussian &rhs) {
    if (this == &rhs) { //Check to make sure not trying to copy the itself
        return *this; 
    }

    mu = rhs.mu;
    sigma = rhs.sigma;

    return *this; //returns reference to current gaussian object.

}

void Gaussian::print_parameters() const {
    std::print("Normal distribution with mean {} and standard deviation {}\n", this.mu, this.sigma);
}

double Gaussian::pdf(double const x) const {
    double z = normalised(x); 
    double coefficient = 1.0 / (sigma * std::sqrt(2.0 * MYPI));
    return coefficient * std::exp(-0.5 * z * z);
}

double Gaussian::cdf(double const x) const {
    double z = normalised(x);
    // Implements Eq. 2: 0.5 * erfc(-z / sqrt(2))
    return 0.5 * std::erfc(-z / std::sqrt(2.0));
}

double Gaussian::inverse_cdf(double z) const {
    // Constants from Abramowitz & Stegun
    double const c0 = 2.515517;
    double const c1 = 0.802853;
    double const c2 = 0.010328;
    
    double const d1 = 1.432788;
    double const d2 = 0.189269;
    double const d3 = 0.001308;
    
    // Determine which half of the distribution we are in
    bool is_lower_half = (z < 0.5);
    
    // Uif z >= 0.5 use 1 - z
    double p;

    if (is_lower_half) {
        p = z;
    } else {
        p = 1.0 - z;
    }
    
    // Calculate t = sqrt(ln(1 / p^2)) which simplifies to sqrt(-2 * ln(p))
    double t = std::sqrt(-2.0 * std::log(p));
    
    // Horner's Method for the numerator and denominator polynomials
    double numerator = c0 + t * (c1 + t * c2);
    double denominator = 1.0 + t * (d1 + t * (d2 + t * d3));
    
    // Calculate standard normal value
    double z_std = t - (numerator / denominator);
    
    // Apply symmetry for the lower half of the distribution
    if (is_lower_half) {
        z_std = -z_std;
    }
    
    // Transform standard z_std back to x in N(mu, sigma^2)
    return z_std * sigma + mu;
}


/// TODO
// Non-Member Functions
// These are just another way to get same functionality.
void print_parameters(const Gaussian & dist){
    std::print("Normal distribution with mean {} and standard deviation {}\n", dist.get_mu(), dist.get_sigma());
}
double pdf(const Gaussian & dist, double x){
    double z = dist.normalised(x); 
    double coefficient = 1.0 / (dist.get_sigma() * std::sqrt(2.0 * MYPI));
    return coefficient * std::exp(-0.5 * z * z);
}


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
    //std::println("CDF of A at x=2 is {:.12f}\n", cdf(A,2.0)); Commented out this line as not implementing free function version of inverse_cdf

    // Copy-and-swap
    Gaussian temp {A};
    A = B; 				// Should call copy assignment operator you wrote
    B = temp;

    print_parameters(A);
    std::println("Inverse CDF for A with Phi=0.4  is {:.6f}"  , A.inverse_cdf(0.4));
    //std::println("Inverse CDF for A with Phi=0.5  is {:.6f}\n", inverse_cdf(A, 0.5)); Commented out this and line below as not implementing free function version of inverse_cdf
    print_parameters(B);
    std::println("Inverse CDF for B with Phi=0.25 is {:.6f}"  , B.inverse_cdf(0.25));
    //std::println("Inverse CDF for B with Phi=0.68 is {:.6f}\n", inverse_cdf(B, 0.68));


    // Generate file of Gaussian distributed RNGs
    std::println("Generating file of 1'000'000 normally distributed numbers\n");
    generate_file_of_rngs(A, "output.txt", 12345, 1'000'000);

    return 0;
}



