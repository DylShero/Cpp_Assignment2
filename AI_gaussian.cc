#include <cmath>
#include <stdexcept>
#include <random>
#include <string>

/**
 * @brief Represents a Gaussian (Normal) distribution N(μ, σ²)
 */
class Gaussian {
public:
    // ── Constructors ─────────────────────────────────────────────────────────

    /// Default: standard normal N(0, 1)
    Gaussian() : mu_(0.0), sigma_(1.0) {}

    /// Construct with mean μ and standard deviation σ
    Gaussian(double mu, double sigma) : mu_(mu), sigma_(sigma) {
        if (sigma <= 0.0)
            throw std::invalid_argument("Standard deviation must be positive.");
    }

    // ── Accessors ─────────────────────────────────────────────────────────────

    double mean()   const { return mu_; }
    double stddev() const { return sigma_; }
    double variance() const { return sigma_ * sigma_; }

    // ── Statistical functions ─────────────────────────────────────────────────

    /// Probability Density Function (PDF)
    double pdf(double x) const {
        static const double INV_SQRT_2PI = 1.0 / std::sqrt(2.0 * M_PI);
        double z = (x - mu_) / sigma_;
        return (INV_SQRT_2PI / sigma_) * std::exp(-0.5 * z * z);
    }

    /// Cumulative Distribution Function (CDF) using erfc
    double cdf(double x) const {
        return 0.5 * std::erfc(-(x - mu_) / (sigma_ * std::sqrt(2.0)));
    }

    /// Inverse CDF (quantile / probit function) via bisection
    double icdf(double p) const {
        if (p <= 0.0 || p >= 1.0)
            throw std::domain_error("Probability must be in (0, 1).");
        // Use a rational approximation (Beasley-Springer-Moro)
        return mu_ + sigma_ * rational_approx(p);
    }

    /// Log of the PDF (numerically stable)
    double log_pdf(double x) const {
        static const double LOG_SQRT_2PI = 0.5 * std::log(2.0 * M_PI);
        double z = (x - mu_) / sigma_;
        return -LOG_SQRT_2PI - std::log(sigma_) - 0.5 * z * z;
    }

    // ── Distribution algebra ──────────────────────────────────────────────────

    /// Sum of two independent Gaussians: N(μ₁+μ₂, σ₁²+σ₂²)
    Gaussian operator+(const Gaussian& other) const {
        return Gaussian(mu_ + other.mu_,
                        std::sqrt(variance() + other.variance()));
    }

    /// Scale a Gaussian by a scalar: N(a·μ, a²·σ²)
    Gaussian operator*(double a) const {
        return Gaussian(a * mu_, std::abs(a) * sigma_);
    }

    // ── Sampling ──────────────────────────────────────────────────────────────

    /// Draw a single sample
    double sample(std::mt19937& rng) const {
        std::normal_distribution<double> dist(mu_, sigma_);
        return dist(rng);
    }

    // ── Utility ───────────────────────────────────────────────────────────────

    /// Standardise a value to Z-score
    double z_score(double x) const { return (x - mu_) / sigma_; }

    /// KL divergence KL(this || other)
    double kl_divergence(const Gaussian& other) const {
        double ratio = sigma_ / other.sigma_;
        double delta = (other.mu_ - mu_) / other.sigma_;
        return std::log(other.sigma_ / sigma_)
             + (variance() + (mu_ - other.mu_) * (mu_ - other.mu_))
               / (2.0 * other.variance())
             - 0.5;
    }

    std::string to_string() const {
        return "N(" + std::to_string(mu_) + ", " + std::to_string(variance()) + ")";
    }

private:
    double mu_;     ///< Mean
    double sigma_;  ///< Standard deviation (σ > 0)

    /// Rational approximation of the inverse normal CDF (Beasley-Springer-Moro)
    static double rational_approx(double p) {
        static const double a[] = {-3.969683028665376e+01,  2.209460984245205e+02,
                                   -2.759285104469687e+02,  1.383577518672690e+02,
                                   -3.066479806614716e+01,  2.506628277459239e+00};
        static const double b[] = {-5.447609879822406e+01,  1.615858368580409e+02,
                                   -1.556989798598866e+02,  6.680131188771972e+01,
                                   -1.328068155288572e+01};
        static const double c[] = {-7.784894002430293e-03, -3.223964580411365e-01,
                                   -2.400758277161838e+00, -2.549732539343734e+00,
                                    4.374664141464968e+00,  2.938163982698783e+00};
        static const double d[] = { 7.784695709041462e-03,  3.224671290700398e-01,
                                    2.445134137142996e+00,  3.754408661907416e+00};

        double q, r;
        if (p < 0.02425) {
            q = std::sqrt(-2.0 * std::log(p));
            return (((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) /
                   ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1.0);
        } else if (p <= 0.97575) {
            q = p - 0.5;
            r = q * q;
            return (((((a[0]*r+a[1])*r+a[2])*r+a[3])*r+a[4])*r+a[5])*q /
                   (((((b[0]*r+b[1])*r+b[2])*r+b[3])*r+b[4])*r+1.0);
        } else {
            q = std::sqrt(-2.0 * std::log(1.0 - p));
            return -(((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) /
                    ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1.0);
        }
    }
};

// Scalar multiplication from the left: a * G
inline Gaussian operator*(double a, const Gaussian& g) { return g * a; }