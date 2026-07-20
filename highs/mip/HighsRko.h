#ifndef HIGHS_RKO_H_
#define HIGHS_RKO_H_

#include <vector>
#include <string>
#include <functional>
#include <random>
#include <chrono>
#include <queue>
#include <tuple>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>

// ============================================================================
// 1. ALL PARAMETERS IN ONE PLACE
// ============================================================================

struct RKOConfig {
    // ----- Problem -----
    int p = 7;
    double alpha = 0.6;
    std::string data_file = "cab25.txt";
    
    // ----- Control -----
    double max_time = 2.0;
    int num_algorithms = 8;
    int num_runs = 4;
    bool verbose = true;
    bool enable_ipr = true;
    
    // ----- Elite Pool -----
    int elite_pool_size = 10;
    
    // ----- SA -----
    double sa_T0 = 1000000.0;
    double sa_alpha = 0.99;
    int sa_max_iter = 1500;
    double sa_beta_min = 0.05;
    double sa_beta_max = 0.20;
    
    // ----- ILS -----
    double ils_beta_min = 0.05;
    double ils_beta_max = 0.20;
    int ils_local_search_iter = 50;
    
    // ----- VNS -----
    int vns_k_max = 10;
    double vns_beta_min = 0.005;
    int vns_local_search_iter = 50;
    
    // ----- LNS -----
    double lns_beta_min = 0.005;
    double lns_beta_max = 0.10;
    double lns_T0 = 1000000.0;
    double lns_alpha = 0.97;
    int lns_farey_level = 7;
    
    // ----- GRASP -----
    double grasp_hs = 0.125;
    double grasp_he = 0.00012;
    int grasp_local_search_iter = 50;
    
    // ----- GA -----
    int ga_population_size = 600;
    double ga_crossover_prob = 0.99;
    double ga_mutation_prob = 0.005;
    int ga_tournament_size = 3;
    
    // ----- PSO -----
    int pso_particle_size = 50;
    double pso_c1 = 2.05;
    double pso_c2 = 2.05;
    double pso_w = 0.73;
    int pso_local_search_iter = 50;
    
    // ----- BRKGA -----
    int brkga_population_size = 1597;
    double brkga_pe = 0.15;
    double brkga_pm = 0.20;
    double brkga_rhoe = 0.70;
    int brkga_local_search_iter = 50;
    
    // ----- IPR -----
    double ipr_block_ratio = 0.10;
    int ipr_direction = 1;
    
    // ----- Local Search -----
    int ls_simple_iterations = 30;
    int ls_strong_iterations = 200;
    double ls_step_size = 0.05;
    
    void print() const;
};

// ============================================================================
// 2. DATA STRUCTURES
// ============================================================================

struct RkoSolution {
    std::vector<double> rk;
    double ofv = 0.0;
    int generation = 0;
    
    RkoSolution() {}
    explicit RkoSolution(int n) : rk(n, 0.0) {}
    void resize(int n) { rk.resize(n, 0.0); }
    int size() const { return static_cast<int>(rk.size()); }
    void randomize();
};

struct THLPSolution {
    std::vector<int> hubs;
    std::vector<int> assignment;
    std::vector<std::pair<int, int> > hub_tree;
    double total_cost = 0.0;
    bool feasible = false;
    
    THLPSolution() {}
    explicit THLPSolution(int n) : assignment(n, -1) {}
    void clear();
};

struct THLPData {
    int n = 0;
    int p = 0;
    double alpha = 0.0;
    std::vector<std::vector<double> > W;
    std::vector<std::vector<double> > C;
    std::vector<double> O;
    std::vector<double> D;
    std::vector<std::string> node_names;
    
    void init(int size);
    void precompute();
    void printSummary() const;
};

// ============================================================================
// 3. THLP DECODER
// ============================================================================

class THLPDecoder {
public:
    THLPDecoder(const THLPData& data, const RKOConfig& config);
    
    THLPSolution decode(const std::vector<double>& rk, 
                        std::vector<double>& binary_solution);
    double evaluate(RkoSolution& sol, std::vector<double>& binary_solution);
    std::function<double(RkoSolution&)> createEvaluator(
        std::vector<double>& binary_solution);
    int getChromosomeSize() const { return chromosome_size_; }
    
private:
    std::vector<int> selectHubs(const std::vector<double>& scores, int p);
    std::vector<int> assignNodes(const std::vector<int>& hubs,
                                 const std::vector<double>& keys);
    std::vector<std::pair<int, int> > buildHubTree(
        const std::vector<int>& hubs,
        const std::vector<double>& edge_priorities,
        const std::vector<int>& assignment);
    double calculateCost(const THLPSolution& solution);
    
    const THLPData& data_;
    const RKOConfig& config_;
    int chromosome_size_;
};

// ============================================================================
// 4. ALGORITHM BASE
// ============================================================================

class RKOAlgorithm {
public:
    virtual ~RKOAlgorithm() {}
    virtual bool run(RkoSolution& sol,
                     std::function<double(RkoSolution&)> evaluator,
                     const RKOConfig& config,
                     bool verbose = true) = 0;
    virtual std::string name() const = 0;
};

class SAAlgorithm : public RKOAlgorithm {
public:
    std::string name() const { return "SA"; }
    bool run(RkoSolution& sol,
             std::function<double(RkoSolution&)> evaluator,
             const RKOConfig& config,
             bool verbose = true);
};

class ILSAlgorithm : public RKOAlgorithm {
public:
    std::string name() const { return "ILS"; }
    bool run(RkoSolution& sol,
             std::function<double(RkoSolution&)> evaluator,
             const RKOConfig& config,
             bool verbose = true);
};

class VNSAlgorithm : public RKOAlgorithm {
public:
    std::string name() const { return "VNS"; }
    bool run(RkoSolution& sol,
             std::function<double(RkoSolution&)> evaluator,
             const RKOConfig& config,
             bool verbose = true);
};

class LNSAlgorithm : public RKOAlgorithm {
public:
    std::string name() const { return "LNS"; }
    bool run(RkoSolution& sol,
             std::function<double(RkoSolution&)> evaluator,
             const RKOConfig& config,
             bool verbose = true);
};

class GRASPAlgorithm : public RKOAlgorithm {
public:
    std::string name() const { return "GRASP"; }
    bool run(RkoSolution& sol,
             std::function<double(RkoSolution&)> evaluator,
             const RKOConfig& config,
             bool verbose = true);
};

class GAAlgorithm : public RKOAlgorithm {
public:
    std::string name() const { return "GA"; }
    bool run(RkoSolution& sol,
             std::function<double(RkoSolution&)> evaluator,
             const RKOConfig& config,
             bool verbose = true);
};

class PSOAlgorithm : public RKOAlgorithm {
public:
    std::string name() const { return "PSO"; }
    bool run(RkoSolution& sol,
             std::function<double(RkoSolution&)> evaluator,
             const RKOConfig& config,
             bool verbose = true);
};

class BRKGAAlgorithm : public RKOAlgorithm {
public:
    std::string name() const { return "BRKGA"; }
    bool run(RkoSolution& sol,
             std::function<double(RkoSolution&)> evaluator,
             const RKOConfig& config,
             bool verbose = true);
};

// ============================================================================
// 5. RANDOM GENERATOR
// ============================================================================

class RandomGenerator {
public:
    static RandomGenerator& instance() {
        static RandomGenerator instance;
        return instance;
    }
    double uniform(double min = 0.0, double max = 1.0);
    int uniformInt(int min, int max);
    std::mt19937& getEngine() { return rng_; }
private:
    RandomGenerator() : rng_(std::chrono::steady_clock::now().time_since_epoch().count()) {}
    std::mt19937 rng_;
};

// ============================================================================
// 6. RKO OPTIMIZER
// ============================================================================

class RKOOptimizer {
public:
    RKOOptimizer();
    explicit RKOOptimizer(const RKOConfig& config);
    ~RKOOptimizer();
    
    RKOConfig& config() { return config_; }
    bool solveTHLP(const THLPData& data, std::vector<double>& solution);
    bool solveAllAlgorithms(const THLPData& data);
    THLPSolution getBestSolution() const { return best_solution_; }
    void printSummary() const;
    
    // Public methods for algorithms to use
    void shakeSolution(RkoSolution& sol, double beta_min, double beta_max);
    void rvnd(RkoSolution& sol, std::function<double(RkoSolution&)> evaluator);
    
private:
    // Core
    RkoSolution createInitialSolution(int chromosome_size);
    bool runSingleAlgorithm(RKOAlgorithm& algorithm,
                            RkoSolution& sol,
                            std::function<double(RkoSolution&)> evaluator);
    void updateBest(const RkoSolution& sol);
    
    // Elite pool
    void updateElitePool(const RkoSolution& sol);
    
    // IPR
    void runIPR(std::function<double(RkoSolution&)> evaluator);
    RkoSolution iprPath(const RkoSolution& atual, const RkoSolution& guia,
                        std::function<double(RkoSolution&)> evaluator);
    
    // Local search variants
    void simpleLocalSearch(RkoSolution& sol,
                           std::function<double(RkoSolution&)> evaluator,
                           int iterations = 30);
    void strongLocalSearch(RkoSolution& sol,
                           std::function<double(RkoSolution&)> evaluator,
                           int iterations = 200);
    void swapLocalSearch(RkoSolution& sol,
                         std::function<double(RkoSolution&)> evaluator);
    void fareyLocalSearch(RkoSolution& sol,
                          std::function<double(RkoSolution&)> evaluator);
    void mirrorLocalSearch(RkoSolution& sol,
                           std::function<double(RkoSolution&)> evaluator);
    void nelderMeadSearch(RkoSolution& sol,
                          std::function<double(RkoSolution&)> evaluator,
                          const std::vector<RkoSolution>& pool);
    RkoSolution blend(const RkoSolution& a, const RkoSolution& b, double rho);
    
    RKOConfig config_;
    THLPSolution best_solution_;
    RkoSolution best_rk_solution_;
    THLPDecoder* decoder_;
    std::vector<double> best_binary_solution_;
    std::vector<RkoSolution> elite_pool_;
    
    SAAlgorithm sa_;
    ILSAlgorithm ils_;
    VNSAlgorithm vns_;
    LNSAlgorithm lns_;
    GRASPAlgorithm grasp_;
    GAAlgorithm ga_;
    PSOAlgorithm pso_;
    BRKGAAlgorithm brkga_;
};

// ============================================================================
// 7. ENTRY POINTS
// ============================================================================

bool solveTHLP(const THLPData& data, std::vector<double>& solution);
bool testAllAlgorithmsOnTHLP();

#endif // HIGHS_RKO_H_