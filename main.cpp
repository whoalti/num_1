#include <iostream> 
#include <cmath>    
#include <iomanip>  
#include <string>

double f(double x) {
    return std::sinh(x) - 12.0 * std::tanh(x) - 0.311;
}

double f_prime(double x) {
    double cosh_x = std::cosh(x);
    return cosh_x - 12.0 / (cosh_x * cosh_x);
}

double phi(double x) {
    return std::asinh(12.0 * std::tanh(x) + 0.311);
}

double f_second_prime(double x) {
    double cosh_x = std::cosh(x);
    double tanh_x = std::tanh(x);
    return std::sinh(x) + 24.0 * tanh_x / (cosh_x * cosh_x);
}

void calculateAprioriRelaxation(double a, double b, double x0, double x_star, double epsilon) {
    double m1 = std::abs(f_prime(a));
    double M1 = std::abs(f_prime(b));
    
    if (m1 > M1) std::swap(m1, M1);

    if (m1 <= 1e-9) {
        std::cout << "\n--- A Priori Relaxation Estimate ---" << std::endl;
        std::cout << "Error: First derivative is zero or too small on the interval." << std::endl;
        return;
    }

    double q_o = (M1 - m1) / (M1 + m1);
    double numerator = log(std::abs(x0 - x_star) / epsilon);
    double denominator = log(1.0 / q_o);
    int steps = static_cast<int>(floor(numerator / denominator)) + 1;

    std::cout << "\n--- A Priori Relaxation Estimate ---" << std::endl;
    std::cout << "Interval: [" << a << ", " << b << "]" << std::endl;
    std::cout << "m1 = " << m1 << ", M1 = " << M1 << std::endl;
    std::cout << "Convergence factor q_o = " << q_o << std::endl;
    std::cout << "Predicted iterations: " << steps << std::endl;
}



void calculateAprioriNewton(double a, double b, double x0, double x_star, double epsilon) {
    double m1 = std::abs(f_prime(a));
    double M2 = std::abs(f_second_prime(b));

    if (m1 <= 1e-9) {
        std::cout << "\n--- A Priori Newton's Estimate ---" << std::endl;
        std::cout << "Error: First derivative is zero or too small on the interval." << std::endl;
        return;
    }

    double q_newton = (M2 * std::abs(x0 - x_star)) / (2.0 * m1);
    
    std::cout << "\n--- A Priori Newton's Estimate ---" << std::endl;
    std::cout << "Interval: [" << a << ", " << b << "]" << std::endl;
    std::cout << "m1 = " << m1 << ", M2 = " << M2 << std::endl;
    std::cout << "Convergence factor q = " << q_newton << std::endl;

    if (q_newton >= 1.0) {
        std::cout << "Convergence condition q < 1 is NOT met. Prediction is not reliable." << std::endl;
        return;
    }

    double numerator = log(std::abs(x0 - x_star) / epsilon);
    double denominator = log(1.0 / q_newton);
    double inner_log = numerator / denominator + 1.0;
    int steps = static_cast<int>(floor(log2(inner_log))) + 1;
    
    std::cout << "Predicted iterations: " << steps << std::endl;
}

void newtonsMethod(double x0, double epsilon, int max_iter) {
    std::cout << "--- Newton's Method ---" << std::endl;
    double x = x0;
    for (int i = 0; i < max_iter; ++i) {
        double fx = f(x);
        double fpx = f_prime(x);

        if (std::abs(fpx) < 1e-12) {
            std::cout << "Error: Derivative is zero. Cannot proceed." << std::endl;
            return;
        }

        double x_next = x - fx / fpx;

        std::cout << "Iteration " << i + 1 << ": x = " << x_next
                  << ", |x_next - x| = " << std::abs(x_next - x) << std::endl;

        if (std::abs(x_next - x) < epsilon) {
            std::cout << "\nConverged successfully!" << std::endl;
            std::cout << "The root is approximately: " << x_next << std::endl;
            std::cout << "Number of iterations: " << i + 1 << std::endl;
            return;
        }
        x = x_next;
    }
    std::cout << "\nFailed to converge within " << max_iter << " iterations." << std::endl;
}

void relaxationMethod(double x0, double epsilon, int max_iter) {
    std::cout << "\n--- Relaxation Method ---" << std::endl;
    double x = x0;
    for (int i = 0; i < max_iter; ++i) {
        double x_next = phi(x);

        std::cout << "Iteration " << i + 1 << ": x = " << x_next
                  << ", |x_next - x| = " << std::abs(x_next - x) << std::endl;

        if (std::abs(x_next - x) < epsilon) {
            std::cout << "\nConverged successfully!" << std::endl;
            std::cout << "The root is approximately: " << x_next << std::endl;
            std::cout << "Number of iterations: " << i + 1 << std::endl;
            return;
        }
        x = x_next;
    }
     std::cout << "\nFailed to converge within " << max_iter << " iterations." << std::endl;
}


// The equation has 3 roots: -3.14986, 0.02828, 3.20207
int main() { 
    double tolerance = 1e-4;

    std::cout << "Enter tolerance (default is 1e-4): ";
    std::string line;
    std::getline(std::cin, line);

    if (!line.empty()) {
        try {
            tolerance = std::stod(line);
        } catch (const std::invalid_argument&) {
            std::cout << "Invalid input. Using default tolerance 1e-4." << std::endl;
            tolerance = 1e-4;
        } catch (const std::out_of_range&) {
            std::cout << "Input out of range. Using default tolerance 1e-4." << std::endl;
            tolerance = 1e-4;
        }
    }

    int max_iterations = 20;

    double a = 2.5;
    double b = 3.5;
    double initial_guess_apriori = 3.5;
    double root_approx = 3.202071; 
    
    std::cout << std::fixed << std::setprecision(6);
    
    calculateAprioriRelaxation(a, b, initial_guess_apriori, root_approx, tolerance);
    calculateAprioriNewton(a, b, initial_guess_apriori, root_approx, tolerance);

    std::cout << "\n========================================\n" << std::endl;
    std::cout << "--- Practical Calculation ---" << std::endl;

    newtonsMethod(initial_guess_apriori, tolerance, max_iterations);
    relaxationMethod(initial_guess_apriori, tolerance, max_iterations);

    return 0;
}